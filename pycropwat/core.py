"""
Core module for effective precipitation calculations.

This module provides the main :class:`EffectivePrecipitation` class for calculating
effective precipitation from climate datasets available on Google Earth Engine (GEE)
**or** from precipitation rasters/NetCDF files that already exist on disk.

Two precipitation sources are supported:

- **Google Earth Engine** (default) - pass ``asset_id`` and ``precip_band``.
- **Local files** - pass ``local_precip`` pointing at a directory of monthly rasters,
  a NetCDF file, or a glob string. See :mod:`pycropwat.local`.

With local precipitation, available water capacity (AWC), reference evapotranspiration
(ETo) and GEE ``FeatureCollection`` geometries are **still** downloaded from Earth
Engine, so the ``usda_scs``, ``ensemble`` and ``suet`` methods continue to require GEE
access. The precipitation-only methods (``cropwat``, ``fao_aglw``, ``fixed_percentage``,
``dependable_rainfall``, ``farmwest``) need **no** Earth Engine at all - Earth Engine is
never initialised for them, so those runs work completely offline. The ``pcml`` method
is a pre-computed GEE product and is therefore incompatible with local precipitation.

The module supports multiple effective precipitation methods:

- **Ensemble**: Mean of 6 methods (default)
- **CROPWAT**: FAO CROPWAT method
- **FAO/AGLW**: FAO Dependable Rainfall (80% exceedance)
- **Fixed Percentage**: Simple fixed percentage method
- **Dependable Rainfall**: FAO Dependable Rainfall method
- **FarmWest**: FarmWest method
- **USDA-SCS**: Soil moisture depletion method (requires AWC and ETo)
- **TAGEM-SuET**: Turkish irrigation method (requires ETo)
- **PCML**: Physics-Constrained ML for Western U.S. (pre-computed GEE asset)

Example
-------
```python
from pycropwat import EffectivePrecipitation
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2015,
    end_year=2020,
    precip_scale_factor=1000,
    method='ensemble'
)
results = ep.process(output_dir='./output', n_workers=4)
```

Local precipitation input (no Earth Engine required):

```python
from pycropwat import EffectivePrecipitation
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    method='cropwat',
    start_year=2005,
    end_year=2005
)
results = ep.process_sequential(output_dir='./output', months=[7, 8])
```

See Also
--------
pycropwat.local : Readers for local precipitation rasters and NetCDF files.
pycropwat.methods : Individual effective precipitation calculation functions.
pycropwat.analysis : Post-processing and analysis tools.
pycropwat.utils : Utility functions for GEE and file operations.
"""

import contextlib
import logging
import threading
import warnings
from pathlib import Path
from typing import Union, Optional, List, Tuple
import ee
import geopandas as gpd
import numpy as np
import xarray as xr
import rioxarray
from rasterio.crs import CRS
import dask
from dask import delayed, compute
from dask.diagnostics import ProgressBar

from .utils import (
    load_geometry,
    get_date_range,
    get_monthly_dates,
    initialize_gee,
    is_gee_asset
)
from .local import LocalPrecipitationSource
from .methods import (
    get_method_function,
    PeffMethod
)

logger = logging.getLogger(__name__)

# Maximum pixels per tile for GEE sampleRectangle (conservative limit)
MAX_PIXELS_PER_TILE = 65536  # 256 x 256

# Guards the lazy creation of an instance's NetCDF read lock. A threading.Lock
# cannot be pickled, so instances create theirs on first use rather than in
# __init__ (see EffectivePrecipitation._get_local_lock).
_LOCK_FACTORY = threading.Lock()

# Earth Engine's initialisation is process-local: a Dask worker process that
# unpickles an instance starts with no initialised client, and every ee call there
# raises until ee.Initialize() runs again in that process. The GEE projects already
# initialised in THIS process are recorded here, so re-initialising costs a set
# lookup rather than an ee.Initialize() round trip per task.
# See EffectivePrecipitation._ensure_gee.
_GEE_INIT_LOCK = threading.Lock()
_GEE_INITIALIZED_PROJECTS = set()

# google-auth is a dependency of earthengine-api, but keep the import defensive so
# a stripped installation still imports pycropwat.
try:
    from google.auth.exceptions import GoogleAuthError as _GoogleAuthError
    _GEE_AUTH_ERROR_TYPES: Tuple[type, ...] = (_GoogleAuthError,)
except Exception:  # pragma: no cover - only when google-auth is unavailable
    _GEE_AUTH_ERROR_TYPES = ()

# Message fragments (matched case-insensitively) that identify "Earth Engine could
# not be reached or authenticated" rather than "Earth Engine has no data here".
# The two must not share a code path: the second is a legitimate reason to fall
# back to a default, the first never is.
_GEE_UNAVAILABLE_MARKERS = (
    'not initialized',
    'not initialised',
    'ee.initialize',
    'ee.authenticate',
    'please authorize',
    'credentials',
    'unauthenticated',
    'invalid_grant',
)


class _NoTilesDownloadedError(ValueError):
    """
    Every Earth Engine tile request for a region came back empty.

    Subclasses :class:`ValueError` so existing callers that catch one keep
    working. Raised only when there is no target shape to fall back to, i.e.
    when a constant-filled array of the right size cannot be substituted.
    """


class _GeeUnavailableError(RuntimeError):
    """
    Earth Engine could not be initialised or authenticated in this process.

    Distinct from "the asset has no data for this region": nothing was ever asked
    of Earth Engine, so no answer - least of all a default constant - can be
    substituted. Raised instead of falling back to the AWC/ETo defaults, so a
    distributed run never writes effective precipitation derived from ancillary
    values that were invented rather than downloaded.
    """


def _is_gee_unavailable(exc: BaseException) -> bool:
    """
    Report whether an exception means Earth Engine itself is unusable here.

    Parameters
    ----------

    exc : BaseException
        Exception raised while talking to Earth Engine.

    Returns
    -------
    bool
        True for an initialisation or authentication failure, False for an
        ordinary data error such as an empty or missing asset.
    """
    if isinstance(exc, _GeeUnavailableError):
        return True
    if _GEE_AUTH_ERROR_TYPES and isinstance(exc, _GEE_AUTH_ERROR_TYPES):
        return True
    message = str(exc).lower()
    return any(marker in message for marker in _GEE_UNAVAILABLE_MARKERS)


def _mark_gee_initialized(project: Optional[str]) -> None:
    """
    Record that Earth Engine is initialised for ``project`` in this process.

    Parameters
    ----------

    project : str or None
        GEE project ID, or None for the default credentials.
    """
    with _GEE_INIT_LOCK:
        _GEE_INITIALIZED_PROJECTS.add(project or '')


# RuntimeWarning numpy raises from np.nanmean when every ensemble member is NaN,
# i.e. on a nodata pixel. Expected, so it is filtered where the effective
# precipitation function is called (see _process_single_month).
_EMPTY_SLICE_WARNING = 'Mean of empty slice'

# PCML (Physics-Constrained Machine Learning) default settings for Western U.S.
PCML_DEFAULT_ASSET = 'projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked'
PCML_DEFAULT_BAND = 'pcml'  # Special marker - actual bands are bYYYY_M format (e.g., b2015_9, b2016_10)
PCML_DEFAULT_SCALE = None  # Retrieved dynamically from asset using nominalScale()

# Western U.S. bounding box (17 states: AZ, CA, CO, ID, KS, MT, NE, NV, NM, ND, OK, OR, SD, TX, UT, WA, WY)
# The PCML image geometry is not bounded in the asset, so we use a predefined extent
# Note: Only Western U.S. vectors overlapping the 17-state extent can be used with PCML
PCML_WESTERN_US_BOUNDS = [
    [-125.0, 49.5],   # Northwest corner
    [-93.0, 49.5],    # Northeast corner
    [-93.0, 25.5],    # Southeast corner
    [-125.0, 25.5],   # Southwest corner
    [-125.0, 49.5]    # Close polygon
]

# PCML annual fraction asset (effective_precip / total_precip, available WY 2000-2024)
# Note: Only annual (water year, Oct-Sep) fractions are available for PCML, not monthly. Band format: bYYYY
PCML_FRACTION_ASSET = 'projects/ee-peff-westus-unmasked/assets/effective_precip_fraction_unmasked'


# Precipitation assets pyCropWat documents, as (marker found in the asset ID,
# asset ID, precipitation band). Used only to make the "you forgot the band" and
# "you forgot the asset" errors concrete; the authoritative band list for any
# asset is the Earth Engine Data Catalog.
_EXAMPLE_PRECIP_SOURCES = (
    ('ERA5', 'ECMWF/ERA5_LAND/MONTHLY_AGGR', 'total_precipitation_sum'),
    ('TERRACLIMATE', 'IDAHO_EPSCOR/TERRACLIMATE', 'pr'),
    ('GRIDMET', 'IDAHO_EPSCOR/GRIDMET', 'pr'),
    ('PRISM', 'OREGONSTATE/PRISM/AN81m', 'ppt'),
    ('CHIRPS', 'UCSB-CHG/CHIRPS/DAILY', 'precipitation'),
    ('IMERG', 'NASA/GPM_L3/IMERG_MONTHLY_V06', 'precipitation'),
)
# Fallback pair quoted when neither the asset nor the band is recognised.
_EXAMPLE_PRECIP_ASSET = 'IDAHO_EPSCOR/TERRACLIMATE'
_EXAMPLE_PRECIP_BAND = 'pr'


def _example_band_for_asset(asset_id: str) -> Optional[str]:
    """
    Get the documented precipitation band for a known asset ID.

    Parameters
    ----------

    asset_id : str
        GEE asset ID the user supplied.

    Returns
    -------

    str or None
        Band name if the asset is one pyCropWat documents, else None.
    """
    upper = str(asset_id).upper()
    for marker, _asset, band in _EXAMPLE_PRECIP_SOURCES:
        if marker in upper:
            return band
    return None


def _example_asset_for_band(precip_band: str) -> Optional[str]:
    """
    Get a documented asset ID that carries the given precipitation band.

    Parameters
    ----------

    precip_band : str
        Band name the user supplied.

    Returns
    -------

    str or None
        Asset ID of the first documented dataset using that band, else None.
    """
    band_name = str(precip_band).strip().lower()
    for _marker, asset, band in _EXAMPLE_PRECIP_SOURCES:
        if band.lower() == band_name:
            return asset
    return None


def _grid_crs(template_da: xr.DataArray) -> Union[str, CRS]:
    """
    Get the CRS attached to a DataArray, defaulting to EPSG:4326.
    
    Local precipitation grids may use any projection, whereas everything
    downloaded from Earth Engine is reprojected to EPSG:4326.
    
    Parameters
    ----------
    
    template_da : xr.DataArray
        DataArray whose CRS should be reused.
        
    Returns
    -------
    rasterio.crs.CRS or str
        The DataArray CRS, or ``'EPSG:4326'`` when it carries none.
    """
    try:
        crs = template_da.rio.crs
    except Exception:  # pragma: no cover - defensive
        crs = None
    return crs if crs is not None else "EPSG:4326"


def get_pcml_band_name(year: int, month: int) -> str:
    """Get PCML band name for a specific year and month.
    
    PCML bands are formatted as bYYYY_M where months 1-9 do not have a preceding zero.
    Examples: b2015_9, b2016_10
    
    Parameters
    ----------
    year : int
        Year (e.g., 2015)
    month : int
        Month (1-12)
        
    Returns
    -------
    str
        Band name in format bYYYY_M
    """
    return f"b{year}_{month}"


class EffectivePrecipitation:
    r"""
    Calculate effective precipitation from GEE or local climate data.
    
    Supports multiple effective precipitation calculation methods including
    CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall, FarmWest,
    and USDA-SCS (which requires AWC and ETo data).
    
    Precipitation is downloaded from Google Earth Engine by default. Pass
    ``local_precip`` to read monthly precipitation from GeoTIFFs or NetCDF files
    on disk instead. AWC, ETo and GEE ``FeatureCollection`` geometries are still
    read from Earth Engine in that case, so ``usda_scs``, ``ensemble`` and ``suet``
    still need GEE credentials. Precipitation-only methods (``cropwat``,
    ``fao_aglw``, ``fixed_percentage``, ``dependable_rainfall``, ``farmwest``)
    combined with local files need no Earth Engine at all - Earth Engine is never
    initialised and the run works offline.
    
    Parameters
    ----------
    
    asset_id : str, optional
        GEE ImageCollection asset ID for precipitation data. Required unless
        ``local_precip`` is given, or ``method='pcml'`` (which supplies its own
        asset). Common options: 
        
        * ``ECMWF/ERA5_LAND/MONTHLY_AGGR`` (ERA5-Land, global, ~11km),
        * ``IDAHO_EPSCOR/TERRACLIMATE`` (TerraClimate, global, ~4km),
        * ``IDAHO_EPSCOR/GRIDMET`` (GridMET, CONUS, ~4km),
        * ``OREGONSTATE/PRISM/AN81m`` (PRISM, CONUS, ~4km),
        * ``UCSB-CHG/CHIRPS/DAILY`` (CHIRPS, 50°S-50°N, ~5km),
        * ``NASA/GPM_L3/IMERG_MONTHLY_V06`` (GPM IMERG, global, ~11km).
    
    precip_band : str, optional
        Name of the precipitation band in the asset. Required unless
        ``local_precip`` is given, or ``method='pcml'``. Examples:
        
        * ERA5-Land: ``total_precipitation_sum``
        * TerraClimate: ``pr``
        * GridMET: ``pr``
        * PRISM: ``ppt``
        * CHIRPS: ``precipitation``
        * GPM IMERG: ``precipitation``
        
    geometry_path : str, Path, or None
        Path to shapefile or GeoJSON file defining the region of interest.
        Can also be a GEE FeatureCollection asset ID. Set to None if using
        gee_geometry_asset instead. With ``local_precip`` the geometry is
        optional (the extent of the local files is used when it is omitted) and
        a local vector file is additionally used to clip the local rasters, see
        ``clip_to_geometry``.

    start_year : int, optional
        Start year for processing (inclusive). Optional with ``local_precip``;
        inferred from the available files when omitted.
    
    end_year : int, optional
        End year for processing (inclusive). Optional with ``local_precip``;
        inferred from the available files when omitted.
    
    scale : float, optional
        Output resolution in meters. If None (default), uses native resolution
        of the dataset.
    
    precip_scale_factor : float, optional
        Factor to convert precipitation to mm. Default is 1.0.
        Common values: ERA5-Land (m to mm) = 1000, TerraClimate = 1.0, GridMET = 1.0.
        Also applied to local precipitation files.
    
    gee_project : str, optional
        GEE project ID for authentication. Required for cloud-based GEE access.
    
    gee_geometry_asset : str, optional
        GEE FeatureCollection asset ID for the region of interest.
        Takes precedence over geometry_path if both are provided.

    method : str, optional
        Effective precipitation calculation method. Default is 'ensemble'.
        Options:
        
        - ``'ensemble'`` - Mean of 6 methods (default, requires AWC and ETo)
        - ``'cropwat'`` - CROPWAT method (FAO standard)
        - ``'fao_aglw'`` - FAO Dependable Rainfall (80% exceedance)
        - ``'fixed_percentage'`` - Simple fixed percentage method
        - ``'dependable_rainfall'`` - FAO Dependable Rainfall method
        - ``'farmwest'`` - FarmWest method
        - ``'usda_scs'`` - USDA-SCS soil moisture depletion method
          (requires AWC and ETo data via method_params)
        - ``'suet'`` - TAGEM-SuET method (Turkish Irrigation Management System)
          (requires ETo data via method_params)
        - ``'pcml'`` - Physics-Constrained ML (Western U.S. only, Jan 2000 - Sep 2024)
          Uses default GEE asset: projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked
          
    method_params : dict, optional
        Additional parameters for the selected method:
        
        For ``'fixed_percentage'``:
            - ``percentage`` (float): Fraction 0-1. Default 0.7.
            
        For ``'dependable_rainfall'``:
            - ``probability`` (float): Probability level 0.5-0.9. Default 0.75.
            
        For ``'usda_scs'``:
            - ``awc_asset`` (str): GEE Image asset ID for AWC data. Required.
              U.S.: projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite
              Global: projects/sat-io/open-datasets/FAO/HWSD_V2_SMU
            - ``awc_band`` (str): Band name for AWC. Default 'AWC'.
            - ``eto_asset`` (str): GEE ImageCollection asset ID for ETo. Required.
              U.S.: projects/openet/assets/reference_et/conus/gridmet/monthly/v1
              Global: projects/climate-engine-pro/assets/ce-ag-era5-v2/daily
            - ``eto_band`` (str): Band name for ETo. Default 'eto'.
              U.S. (GridMET): 'eto', Global (AgERA5): 'ReferenceET_PenmanMonteith_FAO56'
            - ``eto_is_daily`` (bool): Whether ETo is daily. Default False.
              Set True for AgERA5 daily data.
            - ``awc_scale_factor`` (float): Scale factor for AWC. Default 1.0.
              FAO HWSD AWC is in mm/m; set to 0.001 to convert to volumetric
              fraction (required for correct USDA-SCS calculation).
              SSURGO AWC is already in inches/inch (volumetric fraction),
              so use 1.0 (default).
            - ``eto_scale_factor`` (float): Scale factor for ETo. Default 1.0.
            - ``rooting_depth`` (float): Rooting depth in meters. Default 1.0.
            - ``mad_factor`` (float): Management Allowed Depletion factor (0-1).
              Controls what fraction of soil water storage is available. Default 0.5.
    
    local_precip : str or Path, optional
        Local precipitation data to use **instead of** downloading precipitation
        from GEE. Accepts:
        
        * a directory of monthly rasters (globbed with ``local_precip_pattern``),
        * a single NetCDF file, or a directory of NetCDF files,
        * a glob string such as ``'./Precip/Precip_*.tif'``.
        
        Files must contain monthly totals (converted to mm with
        ``precip_scale_factor`` if needed) and must be datable, either from the
        file name (``Precip_2005_07.tif``) or from a NetCDF time coordinate.
        AWC, ETo and GEE ``FeatureCollection`` geometries are still read from
        Earth Engine; precipitation-only methods need no Earth Engine at all.
        Not compatible with ``method='pcml'``, which is a pre-computed GEE product.
    
    local_precip_pattern : str, optional
        Glob used when ``local_precip`` is a directory. Default ``'*.tif'``.
    
    local_precip_variable : str, optional
        NetCDF variable holding precipitation. Auto-detected when None (default).
    
    local_precip_nodata : float, optional
        Additional nodata sentinel, applied on top of the value stored in the file
        metadata (e.g. ``-9999``). Matching pixels become NaN and propagate as NaN
        to both output rasters.
    
    local_precip_crs : str, optional
        CRS to assign when the local files carry none, or to override the stored
        one, e.g. ``'EPSG:4326'``.
    
    local_precip_date_regex : str, optional
        Regular expression with named groups ``year`` and ``month``, used to date
        files the built-in parser cannot read, e.g.
        ``r'(?P<month>\d{2})_(?P<year>\d{4})'``.
    
    clip_to_geometry : bool, optional
        Clip local precipitation rasters to ``geometry_path`` when it points at a
        local vector file. Default True. Has no effect on GEE precipitation, which
        is always clipped to the geometry server-side.
        
    Attributes
    ----------
    geometry : ee.Geometry or None
        The loaded geometry for the region of interest. ``None`` when local
        precipitation is used and Earth Engine is not required.
    
    collection : ee.ImageCollection or None
        The filtered and scaled precipitation image collection. ``None`` for the
        PCML method and when local precipitation files are used.
    
    bounds : list
        Bounding box coordinates of the geometry (lon/lat ring).
        
    Examples
    --------
    Basic usage with Ensemble method (default):
    
    ```python
    from pycropwat import EffectivePrecipitation
    ep = EffectivePrecipitation(
        asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
        precip_band='total_precipitation_sum',
        geometry_path='roi.geojson',
        start_year=2015,
        end_year=2020,
        precip_scale_factor=1000
    )
    ep.process(output_dir='./output', n_workers=4)
    ```
    
    Using GEE FeatureCollection asset:
    
    ```python
    ep = EffectivePrecipitation(
        asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
        precip_band='total_precipitation_sum',
        gee_geometry_asset='projects/my-project/assets/study_area',
        start_year=2015,
        end_year=2020,
        precip_scale_factor=1000,
        gee_project='my-gee-project'
    )
    ```
    
    Using FAO/AGLW method:
    
    ```python
    ep = EffectivePrecipitation(
        asset_id='IDAHO_EPSCOR/TERRACLIMATE',
        precip_band='pr',
        geometry_path='study_area.geojson',
        start_year=2000,
        end_year=2020,
        method='fao_aglw'
    )
    ```
    
    Using fixed percentage method (80%):
    
    ```python
    ep = EffectivePrecipitation(
        asset_id='IDAHO_EPSCOR/GRIDMET',
        precip_band='pr',
        geometry_path='farm.geojson',
        start_year=2010,
        end_year=2020,
        method='fixed_percentage',
        method_params={'percentage': 0.8}
    )
    ```
    
    Using USDA-SCS method with AWC and ETo data:
    
    ```python
    ep = EffectivePrecipitation(
        asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
        precip_band='total_precipitation_sum',
        geometry_path='arizona.geojson',
        start_year=2015,
        end_year=2020,
        precip_scale_factor=1000,
        method='usda_scs',
        method_params={
            'awc_asset': 'projects/my-project/assets/soil_awc',
            'awc_band': 'AWC',
            'eto_asset': 'IDAHO_EPSCOR/GRIDMET',
            'eto_band': 'eto',
            'eto_is_daily': True,
            'rooting_depth': 1.0
        }
    )
    ```
    
    Local precipitation input (no Earth Engine required):
    
    ```python
    ep = EffectivePrecipitation(
        local_precip='../pyCropWat_Data/Precip',
        local_precip_pattern='Precip_*.tif',
        local_precip_nodata=-9999,
        method='cropwat'
    )
    # start_year/end_year inferred from the files when omitted
    ep.process_sequential(output_dir='./output', months=[7, 8])
    ```
    
    Local precipitation with AWC and ETo still coming from Earth Engine:
    
    ```python
    ep = EffectivePrecipitation(
        local_precip='../pyCropWat_Data/Precip',
        local_precip_pattern='Precip_*.tif',
        local_precip_nodata=-9999,
        geometry_path='basin.geojson',
        start_year=2005,
        end_year=2010,
        method='ensemble',
        gee_project='my-gee-project',
        method_params={
            'awc_asset': 'projects/sat-io/open-datasets/FAO/HWSD_V2_SMU',
            'awc_band': 'AWC',
            'awc_scale_factor': 0.001,
            'eto_asset': 'IDAHO_EPSCOR/TERRACLIMATE',
            'eto_band': 'pet',
            'eto_scale_factor': 0.1,
            'eto_is_daily': False,
            'rooting_depth': 2.0,
            'mad_factor': 1.0
        }
    )
    ep.process(output_dir='./output', n_workers=4)
    ```
    
    Local NetCDF input with an explicit variable and date pattern:
    
    ```python
    ep = EffectivePrecipitation(
        local_precip='./wrf_precip.nc',
        local_precip_variable='RAINNC',
        local_precip_crs='EPSG:4326',
        method='fao_aglw'
    )
    ```

    See Also
    --------
    pycropwat.local : Readers for local precipitation rasters and NetCDF files.
    pycropwat.methods : Individual effective precipitation calculation functions.
    pycropwat.analysis : Post-processing and analysis tools.
    """
    
    def __init__(
        self,
        asset_id: Optional[str] = None,
        precip_band: Optional[str] = None,
        geometry_path: Optional[Union[str, Path]] = None,
        start_year: int = None,
        end_year: int = None,
        scale: Optional[float] = None,
        precip_scale_factor: float = 1.0,
        gee_project: Optional[str] = None,
        gee_geometry_asset: Optional[str] = None,
        method: PeffMethod = 'ensemble',
        method_params: Optional[dict] = None,
        local_precip: Optional[Union[str, Path]] = None,
        local_precip_pattern: str = '*.tif',
        local_precip_variable: Optional[str] = None,
        local_precip_nodata: Optional[float] = None,
        local_precip_crs: Optional[str] = None,
        local_precip_date_regex: Optional[str] = None,
        clip_to_geometry: bool = True,
    ):
        self.asset_id = asset_id
        self.precip_band = precip_band
        self.geometry_path = geometry_path
        self.gee_geometry_asset = gee_geometry_asset
        self.start_year = start_year
        self.end_year = end_year
        self.scale = scale  # None means use native resolution
        self.precip_scale_factor = precip_scale_factor
        self.gee_project = gee_project
        self.method = method
        self.method_params = method_params or {}
        
        # Local precipitation configuration
        self.local_precip = local_precip
        self.local_precip_pattern = local_precip_pattern
        self.local_precip_variable = local_precip_variable
        self.local_precip_nodata = local_precip_nodata
        self.local_precip_crs = local_precip_crs
        self.local_precip_date_regex = local_precip_date_regex
        self.clip_to_geometry = clip_to_geometry
        self._is_local = local_precip is not None
        self._local_source = None
        self._local_scale_meters = None
        self._clip_gdf = None
        # get_month() is called from Dask worker threads in process(); NetCDF handles
        # are shared and are not thread-safe, so reads are serialised with a lock.
        # The same lock makes the coverage-warning check atomic (_warn_partial_coverage).
        # It is created on first use (_get_local_lock) because a threading.Lock
        # cannot be pickled and a GEE-only instance never needs one.
        self._local_lock = None
        # Ancillary fields already reported as only partly covering the local grid.
        # Process-local: it is pickled with the instance, so each Dask worker
        # process warns at most once per field (see _warn_partial_coverage).
        self._regrid_coverage_warned = set()
        # Whether written rasters declare NaN as their nodata value. Decided once per
        # run from the configuration instead of per month from the data, so a run
        # never emits a time series with mixed nodata metadata. Only local
        # precipitation carries nodata; GEE downloads use defaultValue and never do.
        self._declare_nodata = self._is_local
        
        # Get the effective precipitation function
        self._peff_function = get_method_function(method)
        
        # USDA-SCS specific: cache for AWC data (loaded once)
        self._awc_cache = None
        
        # Input directory for saving downloaded data (set during process())
        self._input_dir = None
        
        # Check if this is PCML method (uses single multi-band Image instead of ImageCollection)
        self._is_pcml = (method == 'pcml' or self.precip_band == PCML_DEFAULT_BAND)
        
        # PCML is a pre-computed GEE product; it cannot be derived from local rasters
        if self._is_local and self._is_pcml:
            raise ValueError(
                "local_precip cannot be combined with method='pcml'. PCML effective "
                "precipitation is a pre-computed Google Earth Engine product for the "
                "Western U.S. and is not calculated from precipitation. Use a "
                "precipitation-based method (e.g. 'cropwat', 'ensemble') with local files."
            )
        
        # A precipitation source is mandatory: either GEE or local files.
        # Each way of getting it wrong reports what is actually missing, so a
        # half-specified GEE source is not mistaken for no source at all.
        if not self._is_local and not self._is_pcml:
            if asset_id is None and precip_band is None:
                raise ValueError(
                    "No precipitation source provided. Either pass asset_id and "
                    "precip_band to download precipitation from Google Earth Engine, "
                    "or pass local_precip pointing at a directory of monthly rasters, "
                    "a NetCDF file, or a glob to read precipitation from disk."
                )
            if precip_band is None:
                band_example = _example_band_for_asset(asset_id)
                if band_example is not None:
                    hint = f", e.g. precip_band='{band_example}' for '{asset_id}'."
                else:
                    hint = (
                        f"; its band list is in the Earth Engine Data Catalog "
                        f"(e.g. precip_band='{_EXAMPLE_PRECIP_BAND}' for "
                        f"'{_EXAMPLE_PRECIP_ASSET}')."
                    )
                raise ValueError(
                    f"precip_band is required when asset_id is given: asset "
                    f"'{asset_id}' was provided but no precipitation band was named. "
                    f"Pass the band holding precipitation in that asset{hint} "
                    f"To read precipitation from disk instead, pass local_precip."
                )
            if asset_id is None:
                asset_example = _example_asset_for_band(precip_band) or _EXAMPLE_PRECIP_ASSET
                raise ValueError(
                    f"asset_id is required when precip_band is given: band "
                    f"'{precip_band}' was provided but no Google Earth Engine asset "
                    f"to read it from. Pass the ImageCollection holding that band, "
                    f"e.g. asset_id='{asset_example}'. To read precipitation from "
                    f"disk instead, pass local_precip."
                )
        
        # For PCML, use default asset if placeholder provided
        if self._is_pcml:
            if self.asset_id == 'PLACEHOLDER' or self.asset_id is None:
                self.asset_id = PCML_DEFAULT_ASSET
                logger.info(f"Using default PCML asset: {self.asset_id}")
            self.precip_band = PCML_DEFAULT_BAND
        
        # Validate that at least one geometry source is provided
        # (not required for PCML, and optional when reading local precipitation)
        if (geometry_path is None and gee_geometry_asset is None
                and not self._is_pcml and not self._is_local):
            raise ValueError("Either geometry_path or gee_geometry_asset must be provided")
        
        # Open the local precipitation source up front so a bad path fails immediately
        if self._is_local:
            self._local_source = self._build_local_source()
            self._local_scale_meters = self._local_native_scale()
        
        # Earth Engine is only needed for GEE precipitation, GEE geometries, or the
        # AWC/ETo driven methods. Precipitation-only methods on local files run offline.
        self._needs_gee = (
            (not self._is_local)
            or gee_geometry_asset is not None
            or (geometry_path is not None and is_gee_asset(str(geometry_path)))
            or method in ('usda_scs', 'ensemble', 'suet')
        )
        
        # Initialize GEE (only when something actually needs it)
        if self._needs_gee:
            initialize_gee(self.gee_project)
            # Remember it for this process, so _ensure_gee() is a no-op here and
            # only a worker process that unpickled the instance pays for a real
            # ee.Initialize().
            _mark_gee_initialized(self.gee_project)
        else:
            logger.info(
                f"Earth Engine not required: local precipitation with method "
                f"'{self.method}' needs no GEE assets. Skipping initialization."
            )
        
        # Load a local vector geometry for clipping the local rasters
        if self._is_local and geometry_path is not None and not is_gee_asset(str(geometry_path)):
            self._clip_gdf = self._read_clip_geometry(geometry_path)
            logger.info(f"Loaded clip geometry from {geometry_path}")
        
        # Resolve the region of interest
        if self._is_pcml and geometry_path is None and gee_geometry_asset is None:
            # Load PCML image first
            self._pcml_image = ee.Image(self.asset_id)
            # Use predefined Western U.S. bounding box since PCML image geometry is unbounded
            self.geometry = ee.Geometry.Polygon([PCML_WESTERN_US_BOUNDS])
            self.bounds = PCML_WESTERN_US_BOUNDS
            logger.info("Using predefined Western U.S. bounding box for PCML")
        elif self._is_local and geometry_path is None and gee_geometry_asset is None:
            # No geometry given: use the extent of the local precipitation files
            min_lon, min_lat, max_lon, max_lat = self._local_source.bounds
            self.bounds = [
                [min_lon, min_lat],
                [max_lon, min_lat],
                [max_lon, max_lat],
                [min_lon, max_lat],
                [min_lon, min_lat]
            ]
            self.geometry = (
                ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])
                if self._needs_gee else None
            )
            logger.info(
                f"No geometry provided; using local precipitation extent "
                f"({min_lon:.4f}, {min_lat:.4f}, {max_lon:.4f}, {max_lat:.4f})"
            )
        elif self._needs_gee:
            # Load geometry from GEE asset or local file
            self.geometry = load_geometry(geometry_path, gee_asset=gee_geometry_asset)
            self.bounds = self.geometry.bounds().getInfo()['coordinates'][0]
        else:
            # Local vector geometry, no Earth Engine needed
            self.geometry = None
            min_lon, min_lat, max_lon, max_lat = self._clip_gdf.to_crs('EPSG:4326').total_bounds
            self.bounds = [
                [min_lon, min_lat],
                [max_lon, min_lat],
                [max_lon, max_lat],
                [min_lon, max_lat],
                [min_lon, min_lat]
            ]
        
        # Infer / clamp the year range against the local files
        if self._is_local:
            self._resolve_local_years()
        
        # Get date range
        self.start_date, self.end_date = get_date_range(self.start_year, self.end_year)
        
        # Load and filter image collection (or load PCML image)
        if self._is_local:
            self.collection = None
            self._pcml_image = None
        else:
            self._load_collection()
    
    def _build_local_source(self) -> LocalPrecipitationSource:
        """
        Open the configured local precipitation files.
        
        Called from :meth:`__init__` and again by the :attr:`_local_source`
        property after unpickling, since open file handles cannot be sent to a
        Dask worker process.
        
        Returns
        -------
        LocalPrecipitationSource
            Source for the configured files.
        """
        return LocalPrecipitationSource(
            self.local_precip,
            pattern=self.local_precip_pattern,
            variable=self.local_precip_variable,
            scale_factor=self.precip_scale_factor,
            nodata=self.local_precip_nodata,
            crs=self.local_precip_crs,
            date_regex=self.local_precip_date_regex
        )
    
    @property
    def _local_source(self) -> Optional[LocalPrecipitationSource]:
        """
        The open local precipitation source, re-opened on demand.
        
        :meth:`__getstate__` drops the source because the xarray/NetCDF handles it
        caches cannot be pickled. The first access after unpickling rebuilds it
        from the constructor arguments, so a Dask worker process opens the files
        itself. ``None`` when precipitation comes from Earth Engine.
        """
        source = self.__dict__.get('_local_source_obj')
        if source is None and self.__dict__.get('_is_local'):
            source = self._build_local_source()
            self.__dict__['_local_source_obj'] = source
        return source
    
    @_local_source.setter
    def _local_source(self, source: Optional[LocalPrecipitationSource]) -> None:
        self.__dict__['_local_source_obj'] = source
    
    def _get_local_lock(self) -> threading.Lock:
        """
        Return this instance's local-mode lock, creating it on first use.
        
        Serialises the two things :meth:`process` does concurrently in local
        precipitation mode: reads from a shared NetCDF handle, which is not
        thread-safe, and the check-then-record in :meth:`_warn_partial_coverage`,
        which must be atomic for the warning to fire exactly once.
        
        A :class:`threading.Lock` cannot be pickled, so it is never stored in
        :meth:`__getstate__` and never allocated by a GEE-only run.
        
        Returns
        -------
        threading.Lock
            Lock guarding shared NetCDF handles and the coverage-warning set.
        """
        lock = self.__dict__.get('_local_lock')
        if lock is None:
            with _LOCK_FACTORY:
                lock = self.__dict__.get('_local_lock')
                if lock is None:
                    lock = threading.Lock()
                    self._local_lock = lock
        return lock
    
    def _ensure_gee(self) -> None:
        """
        Initialise Earth Engine in the current process, once.
        
        Earth Engine keeps its client state per process, so an instance that was
        pickled to a Dask worker process (``processes=True``, ``LocalCluster``,
        ``distributed``) arrives with ``_needs_gee`` set but no usable client:
        every ``ee`` call there raises "Earth Engine client library not
        initialized". This is called before each Earth Engine access and
        re-initialises the worker on first use.
        
        Cheap and idempotent: a project already initialised in this process is a
        set lookup, not an :func:`ee.Initialize` round trip, so per-task cost is
        negligible. A run that needs no Earth Engine at all - local precipitation
        with a precipitation-only method - has ``_needs_gee`` False and this
        initialises nothing.
        
        Raises
        ------
        _GeeUnavailableError
            Earth Engine could not be initialised. Raised rather than swallowed,
            so callers cannot mistake it for missing data and substitute a
            default.
        
        See Also
        --------
        __setstate__ : Re-initialises the worker as soon as the instance lands.
        """
        if not self.__dict__.get('_needs_gee'):
            return
        
        project = self.__dict__.get('gee_project')
        key = project or ''
        if key in _GEE_INITIALIZED_PROJECTS:
            return
        
        with _GEE_INIT_LOCK:
            # Re-check: another thread of this process may have just done it.
            if key in _GEE_INITIALIZED_PROJECTS:
                return
            where = f" for project {project}" if project else ""
            logger.info(f"Initializing Earth Engine in this process{where}")
            try:
                initialize_gee(project)
            except Exception as e:
                raise _GeeUnavailableError(
                    f"Earth Engine could not be initialized in this process: {e}. "
                    f"Ancillary data (AWC/ETo) will not be replaced by defaults."
                ) from e
            _GEE_INITIALIZED_PROJECTS.add(key)
    
    def __getstate__(self) -> dict:
        """
        Return picklable state, so instances survive Dask/``distributed``.
        
        The NetCDF read lock and the open local precipitation source are dropped:
        a lock and an open file handle cannot be pickled. Both are recreated on
        demand by :meth:`_get_local_lock` and :attr:`_local_source`.
        
        Returns
        -------
        dict
            Instance state without the unpicklable entries.
        """
        state = self.__dict__.copy()
        state.pop('_local_lock', None)
        state.pop('_local_source_obj', None)
        return state
    
    def __setstate__(self, state: dict) -> None:
        """
        Restore state, leaving the lock and local source to be rebuilt lazily.
        
        Earth Engine is re-initialised here when the restored instance needs it
        (:meth:`_ensure_gee`), because ``ee`` state is process-local and does not
        travel with the pickle: without this a Dask worker process would fail
        every ``ee`` call with "Earth Engine client library not initialized". A
        run that needs no Earth Engine initialises nothing, and a process that is
        already initialised pays only a set lookup.
        
        A failure here is logged rather than raised: unpickling must not fail, and
        the first real Earth Engine call retries and raises
        :class:`_GeeUnavailableError` so the failure is reported against the task
        that needed it.
        
        Parameters
        ----------
        
        state : dict
            State produced by :meth:`__getstate__`.
        """
        self.__dict__.update(state)
        self.__dict__['_local_lock'] = None
        if state.get('_needs_gee'):
            try:
                self._ensure_gee()
            except Exception as e:
                logger.error(
                    f"Earth Engine could not be initialized after unpickling: {e}"
                )
    
    @staticmethod
    def _read_clip_geometry(geometry_path: Union[str, Path]) -> gpd.GeoDataFrame:
        """
        Read the local vector geometry used to clip local precipitation.
        
        Applies the same existence and format checks as
        :func:`~pycropwat.utils.load_geometry_from_file` so a typo in the path
        raises a clear error instead of a driver-level one. Any format GeoPandas
        can open is accepted (``.shp``, ``.geojson``, ``.json``, ``.gpkg``, ...).
        
        Parameters
        ----------
        
        geometry_path : str or Path
            Path to the vector file.
            
        Returns
        -------
        geopandas.GeoDataFrame
            The loaded features.
            
        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        
        ValueError
            If the file cannot be read as a vector dataset, or holds no features.
        """
        path = Path(geometry_path)
        
        if not path.exists():
            raise FileNotFoundError(f"Geometry file not found: {path}")
        
        try:
            gdf = gpd.read_file(path)
        except Exception as e:
            raise ValueError(
                f"Could not read geometry file {path}: {e}. Use a vector format "
                f"GeoPandas can open, such as .shp, .geojson, .json or .gpkg."
            ) from e
        
        if gdf.empty:
            raise ValueError(f"Geometry file contains no features: {path}")
        
        return gdf
    
    def _local_native_scale(self) -> Optional[float]:
        """
        Approximate the local precipitation grid resolution in meters.
        
        Used as the download scale for GEE-sourced AWC and ETo so that those
        rasters are requested at roughly the resolution of the local grid.
        
        Returns
        -------
        float or None
            Approximate cell size in meters, or None if it cannot be determined.
        """
        try:
            res_x, res_y = self._local_source.resolution
            res = float(min(abs(res_x), abs(res_y)))
            crs = self._local_source.crs
            if crs is not None and crs.is_geographic:
                res = res * 111320.0
            logger.debug(f"Local grid native scale: {res:.2f} m")
            return res
        except Exception as e:  # pragma: no cover - defensive
            logger.warning(f"Could not determine local grid scale: {e}")
            return None
    
    def _resolve_local_years(self) -> None:
        """
        Infer or clamp ``start_year``/``end_year`` from the local precipitation files.
        
        When either bound is None it is taken from the source's year range. A
        user-supplied range is clamped to the years that actually exist on disk.
        
        Raises
        ------
        ValueError
            If the requested range does not overlap the available years.
        """
        src_start, src_end = self._local_source.year_range
        
        if self.start_year is None:
            self.start_year = src_start
            logger.info(f"start_year not given; inferred {self.start_year} from local files")
        if self.end_year is None:
            self.end_year = src_end
            logger.info(f"end_year not given; inferred {self.end_year} from local files")
        
        clamped_start = max(int(self.start_year), src_start)
        clamped_end = min(int(self.end_year), src_end)
        
        if clamped_start > clamped_end:
            raise ValueError(
                f"Requested years {self.start_year}-{self.end_year} do not overlap the "
                f"local precipitation data ({src_start}-{src_end})."
            )
        
        if (clamped_start, clamped_end) != (int(self.start_year), int(self.end_year)):
            logger.warning(
                f"Requested years {self.start_year}-{self.end_year} clamped to "
                f"{clamped_start}-{clamped_end} (local data covers {src_start}-{src_end})"
            )
        
        self.start_year = clamped_start
        self.end_year = clamped_end
    
    def _filter_local_dates(
        self,
        all_dates: List[Tuple[int, int]]
    ) -> List[Tuple[int, int]]:
        """
        Drop (and report) months that have no local precipitation file.
        
        Parameters
        ----------
        
        all_dates : list of tuple
            Candidate ``(year, month)`` pairs.
            
        Returns
        -------
        list of tuple
            The subset of ``all_dates`` that exists on disk. Returned unchanged
            when precipitation comes from Earth Engine.
        """
        if not self._is_local:
            return all_dates
        
        available = set(self._local_source.available_months())
        missing = [(y, m) for y, m in all_dates if (y, m) not in available]
        
        if missing:
            shown = ', '.join(f"{y}-{m:02d}" for y, m in missing[:12])
            if len(missing) > 12:
                shown = f"{shown}, ... (+{len(missing) - 12} more)"
            logger.warning(
                f"No local precipitation file for {len(missing)} requested month(s): {shown}"
            )
        
        return [(y, m) for y, m in all_dates if (y, m) in available]
    
    def _load_collection(self) -> None:
        """Load and prepare the GEE ImageCollection (or PCML Image)."""
        if self._is_pcml:
            # PCML is a single multi-band Image, not an ImageCollection
            # Bands are named bYYYY_M (e.g., b2015_9, b2016_10)
            # May already be loaded if geometry was derived from it
            if not hasattr(self, '_pcml_image') or self._pcml_image is None:
                self._pcml_image = ee.Image(self.asset_id)
            self.collection = None  # Not used for PCML
            
            # Get PCML native scale from asset using nominalScale()
            self._pcml_scale = self._pcml_image.projection().nominalScale().getInfo()
            logger.info(f"Loaded PCML image with dynamic band selection (bYYYY_M format)")
            logger.info(f"PCML native scale from asset: {self._pcml_scale:.2f}m")
        else:
            # Standard ImageCollection processing
            self.collection = (
                ee.ImageCollection(self.asset_id)
                .select(self.precip_band)
                .filterDate(self.start_date, self.end_date)
                .filterBounds(self.geometry)
            )
            
            # Apply scale factor and rename band to 'pr'
            def scale_and_rename(img):
                return (
                    img.multiply(self.precip_scale_factor)
                    .rename('pr')
                    .copyProperties(img, ['system:time_start', 'system:time_end'])
                )
            
            self.collection = self.collection.map(scale_and_rename)
            self._pcml_image = None
        
    @staticmethod
    def cropwat_effective_precip(pr: np.ndarray) -> np.ndarray:
        """
        Calculate CROPWAT effective precipitation.
        
        Parameters
        ----------
        
        pr : np.ndarray
            Precipitation in mm.
            
        Returns
        -------
        np.ndarray
            Effective precipitation in mm.
        """
        ep = np.where(
            pr <= 250,
            pr * (125 - 0.2 * pr) / 125,
            0.1 * pr + 125
        )
        return ep
    
    def _get_native_scale(self) -> float:
        """
        Get the native scale (resolution) of the image collection in meters.
        
        Returns
        -------
        float
            Native scale in meters.
        """
        try:
            # For local precipitation, use the resolution of the files on disk
            if self._is_local and self._local_scale_meters is not None:
                return self._local_scale_meters
            
            # For PCML, use the pre-computed scale from the asset
            if self._is_pcml and hasattr(self, '_pcml_scale') and self._pcml_scale is not None:
                return self._pcml_scale
            
            # Get the first image from the collection to determine native scale
            first_img = self.collection.first()
            projection = first_img.projection()
            native_scale = projection.nominalScale().getInfo()
            logger.info(f"Native scale: {native_scale} meters")
            return native_scale
        except Exception as e:
            logger.warning(f"Could not determine native scale, defaulting to 10000m: {e}")
            return 10000.0
    
    def _estimate_pixel_count(self, bounds_coords: List, scale_meters: float) -> int:
        """
        Estimate the number of pixels for a given bounds and scale.
        
        Parameters
        ----------
        
        bounds_coords : list
            Bounding box coordinates [[min_lon, min_lat], [max_lon, min_lat], ...]
        
        scale_meters : float
            Resolution in meters.
            
        Returns
        -------
        int
            Estimated number of pixels.
        """
        min_lon = min(c[0] for c in bounds_coords)
        max_lon = max(c[0] for c in bounds_coords)
        min_lat = min(c[1] for c in bounds_coords)
        max_lat = max(c[1] for c in bounds_coords)
        
        # Approximate width and height in meters (at mid-latitude)
        mid_lat = (min_lat + max_lat) / 2
        lat_meters_per_degree = 111320  # meters per degree latitude
        lon_meters_per_degree = 111320 * np.cos(np.radians(mid_lat))
        
        width_meters = (max_lon - min_lon) * lon_meters_per_degree
        height_meters = (max_lat - min_lat) * lat_meters_per_degree
        
        n_cols = int(np.ceil(width_meters / scale_meters))
        n_rows = int(np.ceil(height_meters / scale_meters))
        
        return n_cols * n_rows
    
    def _create_tile_grid(self, bounds_coords: List, scale_meters: float) -> List[ee.Geometry]:
        """
        Create a grid of tiles that cover the bounding box.
        
        Parameters
        ----------
        
        bounds_coords : list
            Bounding box coordinates.
        
        scale_meters : float
            Resolution in meters.
            
        Returns
        -------
        list
            List of ee.Geometry.Rectangle tiles.
        """
        min_lon = min(c[0] for c in bounds_coords)
        max_lon = max(c[0] for c in bounds_coords)
        min_lat = min(c[1] for c in bounds_coords)
        max_lat = max(c[1] for c in bounds_coords)
        
        # Calculate tile size in degrees based on MAX_PIXELS_PER_TILE
        tile_pixels = int(np.sqrt(MAX_PIXELS_PER_TILE))  # e.g., 256 pixels per side
        
        mid_lat = (min_lat + max_lat) / 2
        lat_meters_per_degree = 111320
        lon_meters_per_degree = 111320 * np.cos(np.radians(mid_lat))
        
        # Tile size in degrees
        tile_height_deg = (tile_pixels * scale_meters) / lat_meters_per_degree
        tile_width_deg = (tile_pixels * scale_meters) / lon_meters_per_degree
        
        tiles = []
        lat = min_lat
        while lat < max_lat:
            lon = min_lon
            while lon < max_lon:
                tile_max_lat = min(lat + tile_height_deg, max_lat)
                tile_max_lon = min(lon + tile_width_deg, max_lon)
                
                tile = ee.Geometry.Rectangle([lon, lat, tile_max_lon, tile_max_lat])
                tiles.append(tile)
                
                lon += tile_width_deg
            lat += tile_height_deg
        
        logger.info(f"Created {len(tiles)} tiles for download")
        return tiles
    
    def _download_tile(
        self,
        img: ee.Image,
        tile: ee.Geometry,
        band_name: str = 'pr',
        default_value: float = 0
    ) -> Optional[Tuple[np.ndarray, List]]:
        """
        Download a single tile from GEE.
        
        Parameters
        ----------
        
        img : ee.Image
            Image to download.
        
        tile : ee.Geometry
            Tile geometry.
        
        band_name : str
            Name of the band to extract from the sampled rectangle.
        
        default_value : float
            Default value for missing data.
            
        Returns
        -------
        tuple or None
            Tuple of (array, coordinates) or None if download fails.
        """
        try:
            arr = img.sampleRectangle(
                region=tile,
                defaultValue=default_value
            ).get(band_name).getInfo()
            
            if arr is None:
                return None
            
            arr = np.array(arr, dtype=np.float32)
            coords = tile.getInfo()['coordinates'][0]
            
            return arr, coords
            
        except Exception as e:
            logger.warning(f"Failed to download tile: {e}")
            return None
    
    @staticmethod
    def _describe_download(asset_id: Optional[str], bounds_coords: List) -> str:
        """
        Name the asset and region of a download, for error messages.
        
        Parameters
        ----------
        
        asset_id : str or None
            Earth Engine asset the image came from, if known.
        
        bounds_coords : list
            Bounding box coordinates ``[[lon, lat], ...]`` of the request.
            
        Returns
        -------
        str
            Phrase such as ``"asset 'X' over bounds [w, s, e, n]"``.
        """
        try:
            box = (f"[{min(c[0] for c in bounds_coords):.4f}, "
                   f"{min(c[1] for c in bounds_coords):.4f}, "
                   f"{max(c[0] for c in bounds_coords):.4f}, "
                   f"{max(c[1] for c in bounds_coords):.4f}]")
        except (TypeError, ValueError, IndexError):  # pragma: no cover - defensive
            box = str(bounds_coords)
        source = f"asset '{asset_id}'" if asset_id else "the requested image"
        return f"{source} over bounds {box}"
    
    def _download_image_chunked(
        self,
        img: ee.Image,
        bounds_coords: List,
        scale_meters: float,
        band_name: str,
        default_value: float,
        target_shape: Optional[tuple] = None,
        data_name: str = "data",
        asset_id: Optional[str] = None
    ) -> np.ndarray:
        """
        Generic chunked download for large regions.
        
        Downloads an image in tiles and mosaics them together. Works for
        precipitation, AWC, ETo, or any single-band image.
        
        Parameters
        ----------
        img : ee.Image
            Image to download (should have a single band with given name).
        bounds_coords : list
            Bounding box coordinates [[lon, lat], ...].
        scale_meters : float
            Resolution in meters.
        band_name : str
            Name of the band to extract (e.g., 'pr', 'awc', 'eto').
        default_value : float
            Default value for missing data and failed downloads.
        target_shape : tuple, optional
            Target shape (rows, cols) to resize output. If None, returns
            at native mosaic resolution.
        data_name : str
            Human-readable name for logging (e.g., "AWC", "ETo", "precipitation").
        asset_id : str, optional
            Earth Engine asset the image came from. Only used to name the source
            in error messages when no tile could be downloaded.
            
        Returns
        -------
        np.ndarray
            Mosaicked array, resized to target_shape if specified.
            
        Raises
        ------
        _NoTilesDownloadedError
            If every tile request came back empty and ``target_shape`` is None,
            so there is no shape to fall back to with ``default_value``.
        """
        try:
            # Create tiles
            tiles = self._create_tile_grid(bounds_coords, scale_meters)
            logger.info(f"Downloading {data_name} in {len(tiles)} tiles...")
            
            # Download each tile
            tile_arrays = []
            tile_coords = []
            
            for i, tile in enumerate(tiles):
                result = self._download_tile(img, tile, band_name, default_value)
                if result is not None:
                    arr, coords = result
                    tile_arrays.append(arr)
                    tile_coords.append(coords)
                else:
                    logger.debug(f"{data_name} tile {i+1}/{len(tiles)} failed")
            
            if not tile_arrays:
                if target_shape:
                    logger.warning(f"All {data_name} tiles failed, using default {default_value}")
                    return np.full(target_shape, default_value, dtype=np.float32)
                raise _NoTilesDownloadedError(
                    f"No {data_name} data could be downloaded from Earth Engine: all "
                    f"{len(tiles)} tile request(s) for "
                    f"{self._describe_download(asset_id, bounds_coords)} at "
                    f"{scale_meters:.0f} m came back empty. Check that the asset "
                    f"exists and is readable by this account, and that it has data "
                    f"over the geometry."
                )
            
            # Mosaic tiles together
            min_lon = min(c[0] for c in bounds_coords)
            max_lon = max(c[0] for c in bounds_coords)
            min_lat = min(c[1] for c in bounds_coords)
            max_lat = max(c[1] for c in bounds_coords)
            
            # Calculate output dimensions
            lat_range = max_lat - min_lat
            lon_range = max_lon - min_lon
            scale_deg = scale_meters / 111320  # approx meters per degree
            out_rows = int(np.ceil(lat_range / scale_deg))
            out_cols = int(np.ceil(lon_range / scale_deg))
            
            # Create output array
            output = np.full((out_rows, out_cols), np.nan, dtype=np.float32)
            
            # Place tiles in output
            for arr, coords in zip(tile_arrays, tile_coords):
                tile_min_lon = min(c[0] for c in coords)
                tile_max_lat = max(c[1] for c in coords)
                
                # Calculate pixel indices
                col_start = int((tile_min_lon - min_lon) / scale_deg)
                row_start = int((max_lat - tile_max_lat) / scale_deg)
                
                # Place data
                rows = min(arr.shape[0], out_rows - row_start)
                cols = min(arr.shape[1], out_cols - col_start)
                
                if rows > 0 and cols > 0:
                    output[row_start:row_start+rows, col_start:col_start+cols] = arr[:rows, :cols]
            
            # Fill NaN with default
            output = np.nan_to_num(output, nan=default_value)
            
            # Resize to match target shape if specified
            if target_shape and output.shape != target_shape:
                from scipy.ndimage import zoom
                zoom_factors = (target_shape[0] / output.shape[0],
                               target_shape[1] / output.shape[1])
                output = zoom(output, zoom_factors, order=1)
            
            logger.info(f"Successfully mosaicked {len(tile_arrays)} {data_name} tiles")
            return output
            
        except _NoTilesDownloadedError:
            raise
        except Exception as e:
            if target_shape:
                logger.warning(
                    f"{data_name} chunked download failed: {e}. Using default {default_value}"
                )
                return np.full(target_shape, default_value, dtype=np.float32)
            logger.warning(
                f"{data_name} chunked download failed for "
                f"{self._describe_download(asset_id, bounds_coords)}: {e}"
            )
            raise
    
    def _get_monthly_image(self, year: int, month: int) -> ee.Image:
        """
        Get a single monthly image from the collection.
        
        Parameters
        ----------
        
        year : int
            Year.
        
        month : int
            Month (1-12).
            
        Returns
        -------
        ee.Image
            Monthly precipitation image (sum of all images in that month).
            For PCML, returns the specific band for that year/month.
        """
        if self._is_pcml:
            # PCML: select band by name bYYYY_M (e.g., b2015_9, b2016_10)
            band_name = get_pcml_band_name(year, month)
            monthly_img = self._pcml_image.select([band_name]).rename('pr')
            logger.debug(f"PCML: Selected band {band_name}")
        else:
            # Standard ImageCollection: filter and sum
            monthly_img = (
                self.collection
                .filter(ee.Filter.calendarRange(year, year, 'year'))
                .filter(ee.Filter.calendarRange(month, month, 'month'))
                .sum()  # Sum all images to get monthly total precipitation
            )
        return monthly_img.clip(self.geometry)
    
    def _download_pcml_annual_fraction(self, year: int, template_da: xr.DataArray) -> np.ndarray:
        """
        Download PCML annual effective precipitation fraction from GEE asset.
        
        The PCML annual fraction asset contains pre-computed effective_precip / total_precip
        ratios for each water year (Oct-Sep, 2000-2024).
        
        Parameters
        ----------
        year : int
            Year (2000-2024).
        template_da : xr.DataArray
            Template DataArray to match spatial extent.
            
        Returns
        -------
        np.ndarray
            Annual effective precipitation fraction at PCML scale.
        """
        logger.info(f"Loading PCML annual fraction for {year}")
        
        try:
            # Load PCML annual fraction image
            pcml_fraction_image = ee.Image(PCML_FRACTION_ASSET)
            
            # Select the band for this year (format: bYYYY)
            band_name = f"b{year}"
            fraction_img = pcml_fraction_image.select([band_name]).rename('fraction').clip(self.geometry)
            logger.debug(f"PCML Fraction: Selected band {band_name}")
            
            # Use PCML native scale
            scale_meters = self._pcml_scale
            
            # Reproject to PCML scale
            fraction_img = fraction_img.reproject(
                crs='EPSG:4326',
                scale=scale_meters
            )
            
            # Download - use chunked download for large regions
            region = self.geometry.bounds()
            bounds_coords = region.getInfo()['coordinates'][0]
            
            # Estimate pixel count
            estimated_pixels = self._estimate_pixel_count(bounds_coords, scale_meters)
            logger.debug(f"PCML fraction download: estimated pixels={estimated_pixels}, max={MAX_PIXELS_PER_TILE}")
            
            if estimated_pixels <= MAX_PIXELS_PER_TILE:
                # Direct download (small region)
                arr = fraction_img.sampleRectangle(
                    region=region,
                    defaultValue=0
                ).get('fraction').getInfo()
                
                if arr is None:
                    logger.warning(f"No PCML fraction data for {year}")
                    return np.full(template_da.shape, 0, dtype=np.float32)
                fraction_arr = np.array(arr, dtype=np.float32)
            else:
                # Chunked download for large regions
                logger.info(f"Large region for PCML fraction ({estimated_pixels} pixels), using chunked download...")
                fraction_arr = self._download_image_chunked(
                    fraction_img, bounds_coords, scale_meters,
                    band_name='fraction', default_value=0.0,
                    target_shape=template_da.shape, data_name="PCML_fraction"
                )
            
            return fraction_arr
            
        except Exception as e:
            logger.warning(f"Error loading PCML annual fraction data: {e}. Returning zeros.")
            return np.full(template_da.shape, 0, dtype=np.float32)

    def _download_monthly_precip(self, year: int, month: int, temp_dir: Optional[Path] = None) -> Optional[xr.DataArray]:
        """
        Get monthly precipitation, from local files or from GEE.
        
        When ``local_precip`` was supplied this delegates to
        :meth:`_load_local_monthly_precip` and never touches Earth Engine.
        Otherwise the month is downloaded from GEE, using a chunked download
        for large regions.
        
        Parameters
        ----------
        
        year : int
            Year.
        
        month : int
            Month (1-12).
        
        temp_dir : Path, optional
            Temporary directory for tile files. If None, uses system temp.
            
        Returns
        -------
        xr.DataArray or None
            Precipitation data array, or None if the month is unavailable.
        """
        if self._is_local:
            return self._load_local_monthly_precip(year, month)
        
        # A Dask worker process unpickled this instance with no Earth Engine
        # client; initialise it before the first ee call (no-op in the process
        # that built the instance).
        self._ensure_gee()
        
        try:
            img = self._get_monthly_image(year, month)
            region = self.geometry.bounds()
            bounds_coords = region.getInfo()['coordinates'][0]
            
            # Determine the scale to use (native or specified)
            if self.scale is not None:
                scale_meters = self.scale
            else:
                # Use native scale from the dataset
                scale_meters = self._get_native_scale()
            
            # Always reproject to ensure consistent resolution
            img = img.reproject(
                crs='EPSG:4326',
                scale=scale_meters
            )
            
            # Estimate pixel count
            estimated_pixels = self._estimate_pixel_count(bounds_coords, scale_meters)
            logger.debug(f"Estimated pixels: {estimated_pixels}, max allowed: {MAX_PIXELS_PER_TILE}")
            
            # Check if we need chunked download
            if estimated_pixels <= MAX_PIXELS_PER_TILE:
                # Direct download (small region)
                return self._download_single_tile(img, region, year, month)
            else:
                # Chunked download (large region)
                logger.info(f"Large region detected ({estimated_pixels} pixels), using chunked download...")
                return self._download_chunked(img, bounds_coords, scale_meters, year, month, temp_dir)
            
        except Exception as e:
            logger.error(f"Error downloading data for {year}-{month:02d}: {e}")
            return None
    
    def _load_local_monthly_precip(self, year: int, month: int) -> Optional[xr.DataArray]:
        """
        Read one month of precipitation from the local source.
        
        Values are already in mm (``precip_scale_factor`` is applied by
        :class:`~pycropwat.local.LocalPrecipitationSource`) and nodata pixels are
        NaN. The array is optionally clipped to ``geometry_path`` and saved to the
        input directory when ``save_inputs`` is enabled.
        
        Parameters
        ----------
        
        year : int
            Year.
        
        month : int
            Month (1-12).
            
        Returns
        -------
        xr.DataArray or None
            Precipitation in mm with dims ``('y', 'x')``, or None when the month
            has no file.
            
        Notes
        -----
        :meth:`process` fans months out over Dask worker threads. NetCDF datasets
        are held open and shared between those threads, so NetCDF reads are
        serialised with an instance-level :class:`threading.Lock`, created on first
        use by :meth:`_get_local_lock`. Raster reads open and close their own handle
        per month and run concurrently.
        """
        source = self._local_source
        lock = self._get_local_lock() if source.kind == 'netcdf' else contextlib.nullcontext()
        
        with lock:
            pr_da = source.get_month(year, month)
        
        if pr_da is None:
            logger.warning(f"No local precipitation data for {year}-{month:02d}, skipping")
            return None
        
        if self.clip_to_geometry and self._clip_gdf is not None:
            pr_da = self._clip_local(pr_da, year, month)
        
        if pr_da.rio.crs is None:
            pr_da = pr_da.rio.write_crs(source.crs)
        
        logger.info(
            f"Loaded local precipitation for {year}-{month:02d}: shape={pr_da.shape}, "
            f"source={Path(pr_da.attrs.get('source', str(source.path))).name}"
        )
        
        # Save precipitation input if input_dir is set
        if self._input_dir is not None:
            pr_path = self._input_dir / f"precip_{year}_{month:02d}.tif"
            if not pr_path.exists():
                out_da = pr_da
                if self._declare_nodata:
                    out_da = out_da.rio.write_nodata(np.nan, encoded=False)
                out_da.rio.to_raster(pr_path)
                logger.info(f"Saved input precipitation: {pr_path.name}")
        
        return pr_da
    
    def _clip_local(self, pr_da: xr.DataArray, year: int, month: int) -> xr.DataArray:
        """
        Clip a local precipitation array to the vector geometry.
        
        Parameters
        ----------
        
        pr_da : xr.DataArray
            Precipitation array to clip.
        
        year : int
            Year (for logging).
        
        month : int
            Month (for logging).
            
        Returns
        -------
        xr.DataArray
            Clipped array, or the original array if clipping fails.
        """
        try:
            gdf = self._clip_gdf
            target_crs = pr_da.rio.crs
            if target_crs is not None and gdf.crs is not None:
                gdf = gdf.to_crs(target_crs)
            clipped = pr_da.rio.clip(gdf.geometry, drop=True, all_touched=True)
            logger.debug(
                f"Clipped {year}-{month:02d} to geometry: {pr_da.shape} -> {clipped.shape}"
            )
            return clipped
        except Exception as e:
            logger.warning(
                f"Could not clip local precipitation for {year}-{month:02d} to geometry: "
                f"{e}. Using the unclipped array."
            )
            return pr_da
    
    def _regrid_to_template(
        self,
        arr: np.ndarray,
        src_bounds: List,
        template_da: xr.DataArray,
        name: str = "data"
    ) -> np.ndarray:
        """
        Land a GEE-downloaded array exactly on the local precipitation grid.
        
        AWC and ETo are downloaded as bare arrays covering the EPSG:4326 bounding
        box of ``self.geometry``. When precipitation comes from local files the
        target grid can use any CRS, extent and resolution, so the array is wrapped
        in a georeferenced DataArray and reprojected onto the precipitation grid
        with :meth:`rioxarray.raster_array.RasterArray.reproject_match`.
        
        Parameters
        ----------
        
        arr : np.ndarray
            Downloaded values, north-up, covering ``src_bounds``.
        
        src_bounds : list
            Bounding box coordinates ``[[lon, lat], ...]`` of the download region.
        
        template_da : xr.DataArray
            Precipitation grid to match.
        
        name : str, optional
            Human-readable name for logging. Default ``'data'``.
            
        Returns
        -------
        np.ndarray
            Values on the template grid. Cells of the precipitation grid that the
            downloaded region does not cover are NaN.
            
        Notes
        -----
        Uncovered cells are **not** backfilled. Fabricating a value there - the
        mean of the covered cells, say - produces plausible-looking effective
        precipitation over ground no ancillary data was ever downloaded for, so
        those cells stay NaN and propagate to NaN in the outputs. The uncovered
        fraction is reported once per field per process by
        :meth:`_warn_partial_coverage`.
        
        Only used in local precipitation mode. The GEE-to-GEE path keeps its
        original :func:`scipy.ndimage.zoom` behaviour, which is also the fallback
        here if ``reproject_match`` fails.
        """
        arr = np.asarray(arr, dtype=np.float32)
        target_shape = (len(template_da.coords['y']), len(template_da.coords['x']))
        
        try:
            min_lon = min(c[0] for c in src_bounds)
            max_lon = max(c[0] for c in src_bounds)
            min_lat = min(c[1] for c in src_bounds)
            max_lat = max(c[1] for c in src_bounds)
            
            src_da = xr.DataArray(
                arr,
                dims=['y', 'x'],
                coords={
                    'y': np.linspace(max_lat, min_lat, arr.shape[0]),
                    'x': np.linspace(min_lon, max_lon, arr.shape[1])
                }
            )
            src_da = src_da.rio.write_crs("EPSG:4326")
            src_da = src_da.rio.write_nodata(np.nan, encoded=False)
            
            matched = src_da.rio.reproject_match(template_da)
            out = np.asarray(matched.values, dtype=np.float32)
            
            self._warn_partial_coverage(out, name)
            
            logger.debug(f"Regridded {name} {arr.shape} -> {out.shape} via reproject_match")
            return out
            
        except Exception as e:
            logger.warning(
                f"reproject_match failed for {name} ({e}); falling back to scipy zoom. "
                f"Values may be misaligned if the local grid does not share bounds "
                f"with the downloaded region."
            )
            if arr.shape != target_shape:
                from scipy.ndimage import zoom
                zoom_factors = (target_shape[0] / arr.shape[0],
                                target_shape[1] / arr.shape[1])
                arr = zoom(arr, zoom_factors, order=1)
            return arr
    
    def _warn_partial_coverage(self, arr: np.ndarray, name: str) -> None:
        """
        Warn once per field per process when a regridded field is not fully covered.
        
        AWC and ETo are downloaded for the bounding box of ``self.geometry``. When
        the local precipitation rasters are not clipped to that geometry - either
        ``clip_to_geometry=False`` or a GEE ``FeatureCollection`` geometry, which
        never clips local files - most of the precipitation grid can fall outside
        the download. Those cells are NaN, and a single warning per field says how
        much of the grid they are, so a mostly-NaN result is never silent.
        
        "Once" means once per field per **process**. The set of already-warned
        fields is instance state that travels with the pickle, so months that run
        in different Dask worker processes each warn once; within a process the
        check-and-record is atomic under :meth:`_get_local_lock`, so concurrent
        worker threads cannot both slip through and warn twice.
        
        Parameters
        ----------
        
        arr : np.ndarray
            Regridded field, NaN where the download did not reach.
        
        name : str
            Human-readable field name, e.g. ``'AWC'`` or ``'ETo'``.
        """
        n_nan = int(np.isnan(arr).sum())
        if not n_nan or arr.size == 0:
            return
        
        # ETo is regridded every month, and process() regrids months concurrently:
        # test-and-record must be one atomic step or two threads both warn. Only
        # reached in local precipitation mode, where the lock already exists.
        with self._get_local_lock():
            if name in self._regrid_coverage_warned:
                return
            self._regrid_coverage_warned.add(name)
        
        covered_pct = 100.0 * (arr.size - n_nan) / arr.size
        logger.warning(
            f"{name} covers only {covered_pct:.1f}% of the precipitation grid; "
            f"{n_nan} pixel(s) will be NaN in the output. Clip the precipitation "
            f"to the geometry (clip_to_geometry=True) or widen the geometry."
        )
    
    @staticmethod
    def _extend_nodata_mask(mask: np.ndarray, *fields: Optional[np.ndarray]) -> np.ndarray:
        """
        Add the NaN cells of ancillary fields to the output nodata mask.
        
        Keeps uncovered AWC/ETo cells (see :meth:`_regrid_to_template`) NaN in the
        outputs even for methods that would otherwise average them away, such as
        ``ensemble``, whose :func:`numpy.nanmean` ignores a NaN member. Fields
        whose shape does not match the mask are skipped, so the GEE-to-GEE path -
        where a directly sampled array need not match the precipitation grid
        exactly - is never affected.
        
        Parameters
        ----------
        
        mask : np.ndarray
            Boolean nodata mask of the precipitation grid.
        
        *fields : np.ndarray or None
            Ancillary arrays to fold in.
            
        Returns
        -------
        np.ndarray
            The mask, with the ancillary NaN cells added.
        """
        for field in fields:
            if field is None:
                continue
            field = np.asarray(field)
            if field.shape != mask.shape or not np.issubdtype(field.dtype, np.floating):
                continue
            mask = mask | np.isnan(field)
        return mask
    
    def _download_single_tile(
        self,
        img: ee.Image,
        region: ee.Geometry,
        year: int,
        month: int
    ) -> Optional[xr.DataArray]:
        """
        Download a single tile directly (for small regions).
        """
        arr = img.sampleRectangle(
            region=region,
            defaultValue=0
        ).get('pr').getInfo()
        
        if arr is None:
            logger.warning(f"No data for {year}-{month:02d}")
            return None
        
        arr = np.array(arr, dtype=np.float32)
        
        # Get coordinates
        coords = region.getInfo()['coordinates'][0]
        min_lon = min(c[0] for c in coords)
        max_lon = max(c[0] for c in coords)
        min_lat = min(c[1] for c in coords)
        max_lat = max(c[1] for c in coords)
        
        # Create coordinate arrays
        lats = np.linspace(max_lat, min_lat, arr.shape[0])
        lons = np.linspace(min_lon, max_lon, arr.shape[1])
        
        # Create xarray DataArray
        da = xr.DataArray(
            arr,
            dims=['y', 'x'],
            coords={
                'y': lats,
                'x': lons
            },
            attrs={
                'units': 'mm',
                'long_name': 'precipitation',
                'year': year,
                'month': month
            }
        )
        da = da.rio.write_crs("EPSG:4326")
        
        # Save precipitation input if input_dir is set
        if self._input_dir is not None:
            pr_path = self._input_dir / f"precip_{year}_{month:02d}.tif"
            if not pr_path.exists():
                da.rio.to_raster(pr_path)
                logger.info(f"Saved input precipitation: {pr_path.name}")
        
        return da
    
    def _load_awc_data(self, template_da: xr.DataArray) -> np.ndarray:
        """
        Load Available Water Capacity (AWC) data for USDA-SCS method.
        
        Downloads AWC data from a GEE Image asset and resamples it to match
        the template DataArray's spatial extent and resolution. Results are
        cached for efficiency across multiple months.
        
        Parameters
        ----------
        
        template_da : xr.DataArray
            Template DataArray to match spatial extent and resolution.
            Typically a precipitation DataArray from the same month.
            
        Returns
        -------
        np.ndarray
            AWC values (fraction, 0-1) resampled to match template grid.
            
        Notes
        -----
        AWC data is loaded from the GEE asset specified in ``method_params['awc_asset']``.
        A scale factor (``method_params['awc_scale_factor']``) is applied after loading
        to convert to volumetric fraction. For FAO HWSD data (mm/m), use 0.001.
        
        If the asset has no data for the region, a default value of 0.15 (15% AWC)
        is used and the substitution is logged at WARNING. An Earth Engine
        initialisation or authentication failure is **not** treated that way: it
        raises :class:`_GeeUnavailableError`, because a client that was never able
        to ask Earth Engine anything has no grounds to answer with a constant.
        
        Earth Engine is always used for AWC, including when precipitation comes from
        local files. In that case the downloaded array is reprojected onto the local
        precipitation grid with :meth:`_regrid_to_template` rather than resampled by
        pixel count, since the local grid may use a different CRS and extent.
        
        The AWC data is cached after first load to avoid repeated downloads
        for subsequent months.
        
        If ``save_inputs`` is True during processing, the AWC data is saved
        as ``awc.tif`` in the input directory.
        
        See Also
        --------
        _load_monthly_eto : Load reference evapotranspiration data.
        """
        if self._awc_cache is not None:
            return self._awc_cache
        
        awc_asset = self.method_params.get('awc_asset')
        awc_band = self.method_params.get('awc_band')  # None for single-band SSURGO, 'AWC' for HWSD
        
        # AWC always comes from Earth Engine, including in local precipitation
        # mode, so make sure this process has an initialised client first.
        self._ensure_gee()
        
        logger.info(f"Loading AWC data from {awc_asset}")
        
        try:
            # Load AWC image (static, not a time series)
            awc_img = ee.Image(awc_asset)
            
            # Select band if specified
            if awc_band:
                awc_img = awc_img.select(awc_band)
            
            # Get the scale to use
            if self.scale is not None:
                scale_meters = self.scale
            else:
                scale_meters = self._get_native_scale()
            
            # Reproject to match template
            awc_img = awc_img.reproject(
                crs='EPSG:4326',
                scale=scale_meters
            ).rename('awc')
            
            # Download AWC data - use chunked download for large regions
            region = self.geometry.bounds()
            bounds_coords = region.getInfo()['coordinates'][0]
            
            # Estimate pixel count
            estimated_pixels = self._estimate_pixel_count(bounds_coords, scale_meters)
            logger.debug(f"AWC download: estimated pixels={estimated_pixels}, max={MAX_PIXELS_PER_TILE}")
            
            if estimated_pixels <= MAX_PIXELS_PER_TILE:
                # Direct download (small region)
                arr = awc_img.sampleRectangle(
                    region=region,
                    defaultValue=0.15  # Default AWC if missing
                ).get('awc').getInfo()
                
                if arr is None:
                    logger.warning("No AWC data available, using default value of 0.15")
                    awc_arr = np.full(template_da.shape, 0.15, dtype=np.float32)
                else:
                    awc_arr = np.array(arr, dtype=np.float32)
            else:
                # Chunked download for large regions
                logger.info(f"Large region for AWC ({estimated_pixels} pixels), using chunked download...")
                awc_arr = self._download_image_chunked(
                    awc_img, bounds_coords, scale_meters,
                    band_name='awc', default_value=0.15,
                    target_shape=None if self._is_local else template_da.shape,
                    data_name="AWC", asset_id=awc_asset
                )
            
            # Apply scale factor (e.g., 0.001 for FAO HWSD mm/m -> volumetric fraction)
            awc_scale_factor = self.method_params.get('awc_scale_factor', 1.0)
            if awc_scale_factor != 1.0:
                awc_arr = awc_arr * awc_scale_factor
                logger.info(f"Applied AWC scale factor: {awc_scale_factor}")
            
            # Local precipitation grids can use any CRS/extent, so reproject the
            # downloaded array onto the precipitation grid instead of a naive zoom.
            if self._is_local:
                awc_arr = self._regrid_to_template(awc_arr, bounds_coords, template_da, "AWC")
            
            # Cache the result
            self._awc_cache = awc_arr
            
            # Save AWC input if input_dir is set
            if self._input_dir is not None:
                awc_path = self._input_dir / "awc.tif"
                if not awc_path.exists():
                    # Create DataArray for saving - use template coords but ensure shape matches
                    # Resize array if needed to exactly match template coordinates
                    if awc_arr.shape != (len(template_da.coords['y']), len(template_da.coords['x'])):
                        from scipy.ndimage import zoom
                        target_shape = (len(template_da.coords['y']), len(template_da.coords['x']))
                        zoom_factors = (target_shape[0] / awc_arr.shape[0],
                                       target_shape[1] / awc_arr.shape[1])
                        awc_arr = zoom(awc_arr, zoom_factors, order=1)
                    
                    awc_da = xr.DataArray(
                        awc_arr,
                        dims=template_da.dims,
                        coords=template_da.coords,
                        attrs={
                            'units': 'fraction',
                            'long_name': 'available_water_capacity',
                            'source': awc_asset
                        }
                    )
                    awc_da = awc_da.rio.write_crs(_grid_crs(template_da))
                    if self._declare_nodata:
                        awc_da = awc_da.rio.write_nodata(np.nan, encoded=False)
                    awc_da.rio.to_raster(awc_path)
                    logger.info(f"Saved input AWC: {awc_path.name}")
            
            return awc_arr
            
        except _GeeUnavailableError:
            # Earth Engine itself is unusable here - nothing was asked of the
            # asset, so there is no missing-data case to fall back for.
            raise
        except Exception as e:
            if _is_gee_unavailable(e):
                raise _GeeUnavailableError(
                    f"Earth Engine is not usable in this process, so AWC could not be "
                    f"downloaded from {awc_asset}: {e}. Refusing to substitute the 0.15 "
                    f"default, which would fabricate ancillary data."
                ) from e
            logger.warning(f"Error loading AWC data: {e}. Using default value of 0.15")
            awc_arr = np.full(template_da.shape, 0.15, dtype=np.float32)
            self._awc_cache = awc_arr
            return awc_arr

    def _load_monthly_eto(self, year: int, month: int, template_da: xr.DataArray) -> np.ndarray:
        """
        Load Reference Evapotranspiration (ETo) data for USDA-SCS method.
        
        Downloads ETo data from a GEE ImageCollection and aggregates to monthly
        totals. The data is resampled to match the template DataArray's spatial
        extent and resolution.
        
        Parameters
        ----------
        
        year : int
            Year to load.
        
        month : int
            Month to load (1-12).
        
        template_da : xr.DataArray
            Template DataArray to match spatial extent and resolution.
            Typically a precipitation DataArray from the same month.
            
        Returns
        -------
        np.ndarray
            Monthly ETo values in mm, resampled to match template grid.
            
        Notes
        -----
        ETo data is loaded from the GEE asset specified in ``method_params['eto_asset']``.
        
        If ``method_params['eto_is_daily']`` is True, daily values are summed
        to monthly totals. Otherwise, the first image in the month is used.
        
        A scale factor can be applied via ``method_params['eto_scale_factor']``
        for unit conversion (e.g., 0.1 if ETo is stored in 0.1 mm units).
        
        Earth Engine is always used for ETo, including when precipitation comes from
        local files. In that case the downloaded array is reprojected onto the local
        precipitation grid with :meth:`_regrid_to_template`.
        
        If the asset has no data for the month, a default value of 100 mm/month is
        used and the substitution is logged at WARNING. An Earth Engine
        initialisation or authentication failure is **not** treated that way: it
        raises :class:`_GeeUnavailableError` rather than fabricating ETo.
        
        If ``save_inputs`` is True during processing, the ETo data is saved
        as ``eto_YYYY_MM.tif`` in the input directory.
        
        See Also
        --------
        _load_awc_data : Load available water capacity data.
        """
        eto_asset = self.method_params.get('eto_asset')
        eto_band = self.method_params.get('eto_band', 'eto')
        eto_is_daily = self.method_params.get('eto_is_daily', False)
        eto_scale_factor = self.method_params.get('eto_scale_factor', 1.0)
        
        # ETo always comes from Earth Engine, including in local precipitation
        # mode, so make sure this process has an initialised client first.
        self._ensure_gee()
        
        logger.info(f"Loading ETo data from {eto_asset} for {year}-{month:02d}")
        
        try:
            # Get date range for this month
            import calendar
            _, days_in_month = calendar.monthrange(year, month)
            start_date = f"{year}-{month:02d}-01"
            end_date = f"{year}-{month:02d}-{days_in_month}"
            
            # Load ETo collection
            eto_coll = (
                ee.ImageCollection(eto_asset)
                .select(eto_band)
                .filterDate(start_date, end_date)
                .filterBounds(self.geometry)
            )
            
            # Sum to monthly (for daily data) or just get the image (for monthly)
            if eto_is_daily:
                eto_img = eto_coll.sum().rename('eto')
                logger.debug(f"Aggregated {days_in_month} daily ETo images to monthly")
            else:
                eto_img = eto_coll.first().rename('eto')
            
            # Get the scale to use
            if self.scale is not None:
                scale_meters = self.scale
            else:
                scale_meters = self._get_native_scale()
            
            # Reproject to match template
            eto_img = eto_img.reproject(
                crs='EPSG:4326',
                scale=scale_meters
            )
            
            # Download ETo data - use chunked download for large regions
            region = self.geometry.bounds()
            bounds_coords = region.getInfo()['coordinates'][0]
            
            # Estimate pixel count
            estimated_pixels = self._estimate_pixel_count(bounds_coords, scale_meters)
            logger.debug(f"ETo download: estimated pixels={estimated_pixels}, max={MAX_PIXELS_PER_TILE}")
            
            if estimated_pixels <= MAX_PIXELS_PER_TILE:
                # Direct download (small region)
                arr = eto_img.sampleRectangle(
                    region=region,
                    defaultValue=0
                ).get('eto').getInfo()
                
                if arr is None:
                    logger.warning(f"No ETo data for {year}-{month:02d}, using default 100 mm")
                    eto_arr = np.full(template_da.shape, 100.0, dtype=np.float32)
                else:
                    eto_arr = np.array(arr, dtype=np.float32)
            else:
                # Chunked download for large regions
                logger.info(f"Large region for ETo ({estimated_pixels} pixels), using chunked download...")
                eto_arr = self._download_image_chunked(
                    eto_img, bounds_coords, scale_meters,
                    band_name='eto', default_value=100.0,
                    target_shape=None if self._is_local else template_da.shape,
                    data_name="ETo", asset_id=eto_asset
                )
            
            # Apply scale factor if needed (for both direct and chunked downloads)
            if eto_scale_factor != 1.0:
                eto_arr = eto_arr * eto_scale_factor
                logger.debug(f"Applied ETo scale factor: {eto_scale_factor}")
            
            # Local precipitation grids can use any CRS/extent, so reproject the
            # downloaded array onto the precipitation grid instead of a naive zoom.
            if self._is_local:
                eto_arr = self._regrid_to_template(eto_arr, bounds_coords, template_da, "ETo")
            
            # Save ETo input if input_dir is set
            if self._input_dir is not None:
                eto_path = self._input_dir / f"eto_{year}_{month:02d}.tif"
                if not eto_path.exists():
                    # Create DataArray for saving - use template coords but ensure shape matches
                    # Resize array if needed to exactly match template coordinates
                    eto_arr_save = eto_arr
                    if eto_arr_save.shape != (len(template_da.coords['y']), len(template_da.coords['x'])):
                        from scipy.ndimage import zoom
                        target_shape = (len(template_da.coords['y']), len(template_da.coords['x']))
                        zoom_factors = (target_shape[0] / eto_arr_save.shape[0],
                                       target_shape[1] / eto_arr_save.shape[1])
                        eto_arr_save = zoom(eto_arr_save, zoom_factors, order=1)
                    
                    eto_da = xr.DataArray(
                        eto_arr_save,
                        dims=template_da.dims,
                        coords=template_da.coords,
                        attrs={
                            'units': 'mm',
                            'long_name': 'reference_evapotranspiration',
                            'year': year,
                            'month': month,
                            'source': eto_asset
                        }
                    )
                    eto_da = eto_da.rio.write_crs(_grid_crs(template_da))
                    if self._declare_nodata:
                        eto_da = eto_da.rio.write_nodata(np.nan, encoded=False)
                    eto_da.rio.to_raster(eto_path)
                    logger.info(f"Saved input ETo: {eto_path.name}")
            
            return eto_arr
            
        except _GeeUnavailableError:
            # Earth Engine itself is unusable here - nothing was asked of the
            # asset, so there is no missing-data case to fall back for.
            raise
        except Exception as e:
            if _is_gee_unavailable(e):
                raise _GeeUnavailableError(
                    f"Earth Engine is not usable in this process, so ETo for "
                    f"{year}-{month:02d} could not be downloaded from {eto_asset}: {e}. "
                    f"Refusing to substitute the 100 mm default, which would fabricate "
                    f"ancillary data."
                ) from e
            logger.warning(f"Error loading ETo data: {e}. Using default value of 100 mm")
            return np.full(template_da.shape, 100.0, dtype=np.float32)

    def _download_chunked(
        self,
        img: ee.Image,
        bounds_coords: List,
        scale_meters: float,
        year: int,
        month: int,
        temp_dir: Optional[Path] = None
    ) -> Optional[xr.DataArray]:
        """
        Download precipitation image in chunks and return as DataArray.
        
        Wrapper around _download_image_chunked that returns an xr.DataArray
        with proper coordinates for precipitation data.
        
        Parameters
        ----------
        img : ee.Image
            Precipitation image to download.
        bounds_coords : list
            Bounding box coordinates.
        scale_meters : float
            Resolution in meters.
        year : int
            Year (for metadata).
        month : int
            Month (for metadata).
        temp_dir : Path, optional
            Unused, kept for API compatibility.
            
        Returns
        -------
        xr.DataArray or None
            Precipitation data array with coordinates, or None if download fails.
        """
        try:
            # Use the generic chunked download
            arr = self._download_image_chunked(
                img, bounds_coords, scale_meters,
                band_name='pr', default_value=0,
                target_shape=None, data_name="precipitation",
                asset_id=self.asset_id
            )
            
            # Get bounds for coordinate creation
            min_lon = min(c[0] for c in bounds_coords)
            max_lon = max(c[0] for c in bounds_coords)
            min_lat = min(c[1] for c in bounds_coords)
            max_lat = max(c[1] for c in bounds_coords)
            
            # Create coordinate arrays
            lats = np.linspace(max_lat, min_lat, arr.shape[0])
            lons = np.linspace(min_lon, max_lon, arr.shape[1])
            
            # Create xarray DataArray
            da = xr.DataArray(
                arr,
                dims=['y', 'x'],
                coords={'y': lats, 'x': lons},
                attrs={
                    'units': 'mm',
                    'long_name': 'precipitation',
                    'year': year,
                    'month': month
                }
            )
            da = da.rio.write_crs("EPSG:4326")
            
            # Save precipitation input if input_dir is set
            if self._input_dir is not None:
                pr_path = self._input_dir / f"precip_{year}_{month:02d}.tif"
                if not pr_path.exists():
                    da.rio.to_raster(pr_path)
                    logger.info(f"Saved input precipitation: {pr_path.name}")
            
            return da
            
        except Exception as e:
            logger.warning(f"Chunked download failed for {year}-{month:02d}: {e}")
            return None
    
    def _process_single_month(
        self,
        year: int,
        month: int,
        output_dir: Path
    ) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Process a single month: download, calculate effective precipitation, save.
        
        Parameters
        ----------
        
        year : int
            Year.
        
        month : int
            Month (1-12).
        
        output_dir : Path
            Output directory.
            
        Returns
        -------
        tuple
            Tuple of (ep_path, epf_path) or (None, None) if processing fails.
            
        Notes
        -----
        Nodata pixels (NaN) in the precipitation input propagate to NaN in both the
        effective precipitation and the fraction output, and the written rasters
        declare NaN as their nodata value. GEE downloads use ``defaultValue=0`` and
        never produce NaN, so this only affects local precipitation files.
        """
        logger.info(f"Processing {year}-{month:02d}")
        
        # Download precipitation data (or read it from local files)
        pr_da = self._download_monthly_precip(year, month)
        
        if pr_da is None:
            return None, None
        
        pr_values = pr_da.values
        # Local rasters carry real nodata; keep numpy quiet about NaN comparisons
        # and 0/0 divisions, then mask those pixels out of both outputs below.
        # AWC/ETo cells the Earth Engine download did not cover are NaN too
        # (see _regrid_to_template) and are folded into the mask as they load.
        nan_mask = np.isnan(pr_values)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            # Calculate effective precipitation using the configured method
            if self.method == 'pcml':
                # PCML: pr_da already contains effective precipitation (PCML Peff)
                # Download annual fraction directly from the GEE asset (once per year)
                ep_arr = pr_values  # PCML Peff is directly downloaded
                
                # For PCML, check if annual fraction file already exists
                pcml_frac_path = output_dir / f"effective_precip_fraction_{year}.tif"
                if pcml_frac_path.exists():
                    # Load existing annual fraction
                    epf_da_existing = rioxarray.open_rasterio(
                        pcml_frac_path).squeeze('band', drop=True)
                    epf_arr = epf_da_existing.values
                    logger.info(f"PCML annual fraction for {year} already exists, reusing")
                else:
                    epf_arr = self._download_pcml_annual_fraction(year, pr_da)
                    logger.info(f"PCML annual fraction loaded from GEE asset for {year}")
            elif self.method == 'usda_scs':
                # USDA-SCS method requires AWC and ETo data
                awc_arr = self._load_awc_data(pr_da)
                eto_arr = self._load_monthly_eto(year, month, pr_da)
                nan_mask = self._extend_nodata_mask(nan_mask, awc_arr, eto_arr)
                rooting_depth = self.method_params.get('rooting_depth', 1.0)
                mad_factor = self.method_params.get('mad_factor', 0.5)
                
                ep_arr = self._peff_function(
                    pr_values,
                    eto_arr,
                    awc_arr,
                    rooting_depth,
                    mad_factor
                )
            elif self.method == 'ensemble':
                # Ensemble method requires AWC and ETo data (same as USDA-SCS)
                awc_arr = self._load_awc_data(pr_da)
                eto_arr = self._load_monthly_eto(year, month, pr_da)
                nan_mask = self._extend_nodata_mask(nan_mask, awc_arr, eto_arr)
                rooting_depth = self.method_params.get('rooting_depth', 1.0)
                fixed_percentage = self.method_params.get('percentage', 0.7)
                dependable_probability = self.method_params.get('probability', 0.75)
                
                # np.nanmean over the 6 members warns on nodata pixels, where every
                # member is NaN. np.errstate does not cover it - it is a warnings
                # module RuntimeWarning, not a floating point error - so filter that
                # one message here. Any other warning still reaches the user.
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        'ignore',
                        message=_EMPTY_SLICE_WARNING,
                        category=RuntimeWarning
                    )
                    ep_arr = self._peff_function(
                        pr_values,
                        eto_arr,
                        awc_arr,
                        rooting_depth,
                        fixed_percentage,
                        dependable_probability
                    )
            elif self.method == 'suet':
                # SuET method requires ETo data
                eto_arr = self._load_monthly_eto(year, month, pr_da)
                nan_mask = self._extend_nodata_mask(nan_mask, eto_arr)
                ep_arr = self._peff_function(pr_values, eto_arr)
            elif self.method_params:
                # Other methods with parameters (fixed_percentage, dependable_rainfall)
                # Filter to only pass valid parameters for the method
                valid_params = {}
                if self.method == 'fixed_percentage':
                    valid_params['percentage'] = self.method_params.get('percentage', 0.7)
                elif self.method == 'dependable_rainfall':
                    valid_params['probability'] = self.method_params.get('probability', 0.75)
                ep_arr = self._peff_function(pr_values, **valid_params)
            else:
                ep_arr = self._peff_function(pr_values)
            
            # Calculate effective precipitation fraction (skip for PCML - already calculated above)
            if self.method != 'pcml':
                epf_arr = np.where(pr_values > 0, ep_arr / pr_values, 0)
            
            # Propagate nodata: NaN precipitation, or an AWC/ETo cell the download
            # never covered -> NaN effective precipitation and NaN fraction (the
            # np.where above would otherwise report a fraction of 0).
            if nan_mask.any():
                ep_arr = np.where(nan_mask, np.nan, ep_arr)
                epf_arr = np.where(nan_mask, np.nan, epf_arr)
        
        # Create effective precipitation DataArray
        ep_da = xr.DataArray(
            ep_arr,
            dims=pr_da.dims,
            coords=pr_da.coords,
            attrs={
                'units': 'mm',
                'long_name': 'effective_precipitation',
                'year': year,
                'month': month,
                'method': self.method.upper()
            }
        )
        ep_da = ep_da.rio.write_crs(_grid_crs(pr_da))
        
        # Create effective precipitation fraction DataArray
        # For PCML: annual fraction loaded directly from GEE asset
        # For others: fraction = peff / precip
        epf_da = xr.DataArray(
            epf_arr.astype(np.float32),
            dims=pr_da.dims,
            coords=pr_da.coords,
            attrs={
                'units': 'fraction',
                'long_name': 'effective_precipitation_fraction',
                'year': year,
                'month': month,
                'method': self.method.upper(),
                'note': 'PCML annual fraction from GEE asset' if self.method == 'pcml' else 'peff / precip'
            }
        )
        epf_da = epf_da.rio.write_crs(_grid_crs(pr_da))
        
        # Declare NaN as nodata so downstream readers mask it out. The decision is
        # self._declare_nodata, taken once per run from the configuration, so every
        # raster of a time series carries the same metadata whether or not this
        # particular month happens to contain a NaN.
        if self._declare_nodata:
            ep_da = ep_da.rio.write_nodata(np.nan, encoded=False)
            epf_da = epf_da.rio.write_nodata(np.nan, encoded=False)
        
        # Save to GeoTIFF
        ep_path = output_dir / f"effective_precip_{year}_{month:02d}.tif"
        
        # For PCML, save annual fraction only once per year (without month suffix)
        if self.method == 'pcml':
            epf_path = output_dir / f"effective_precip_fraction_{year}.tif"
        else:
            epf_path = output_dir / f"effective_precip_fraction_{year}_{month:02d}.tif"
        
        ep_da.rio.to_raster(ep_path)
        
        # Only save fraction file if it doesn't exist (for PCML) or always (for others)
        if self.method != 'pcml' or not epf_path.exists():
            epf_da.rio.to_raster(epf_path)
            logger.info(f"Saved: {ep_path.name}, {epf_path.name}")
        else:
            logger.info(f"Saved: {ep_path.name} (fraction already exists)")
        
        return ep_path, epf_path
    
    def process(
        self,
        output_dir: Union[str, Path],
        n_workers: int = 4,
        months: Optional[List[int]] = None,
        input_dir: Optional[Union[str, Path]] = None,
        save_inputs: bool = False
    ) -> List[Tuple[Optional[Path], Optional[Path]]]:
        """
        Process all months and save effective precipitation rasters.
        
        Downloads precipitation data from Google Earth Engine - or reads it from
        the local files given by ``local_precip`` - calculates effective
        precipitation using the configured method, and saves results as GeoTIFF
        files. Uses Dask for parallel processing of multiple months.
        
        Parameters
        ----------
        
        output_dir : str or Path
            Directory to save output rasters. Will be created if it
            doesn't exist.
        
        n_workers : int, optional
            Number of parallel workers for Dask. Default is 4.
            Set to 1 for sequential processing.
        
        months : list of int, optional
            List of months to process (1-12). If None, processes all months
            in the date range. Useful for seasonal analyses.
        
        input_dir : str or Path, optional
            Directory to save downloaded input data (precipitation, AWC, ETo).
            If None and save_inputs is True, uses ``output_dir/../analysis_inputs``.
        
        save_inputs : bool, optional
            Whether to save downloaded input data as GeoTIFF files.
            Default is False. Useful for debugging or further analysis.
            
        Returns
        -------
        list of tuple
            List of tuples containing paths to saved files:
            ``(effective_precip_path, effective_precip_fraction_path)``.
            Returns ``(None, None)`` for months that failed to process.
            
        Notes
        -----
        Output files are named:
        
        - ``effective_precip_YYYY_MM.tif`` - Effective precipitation in mm
        - ``effective_precip_fraction_YYYY_MM.tif`` - Effective/total ratio (non-PCML methods)
        - ``effective_precip_fraction_YYYY.tif`` - Annual (water year) fraction (PCML method only)
        
        For the USDA-SCS method, AWC and ETo data are automatically downloaded
        and cached for efficiency. This is true for local precipitation too:
        AWC and ETo always come from Earth Engine.
        
        With local precipitation, months that have no file are skipped with a
        warning instead of being reported as failures, and nodata pixels are
        written as NaN in both outputs.
        
        Examples
        --------
        Process all months in parallel:
        
        ```python
        ep = EffectivePrecipitation(...)
        results = ep.process(output_dir='./output', n_workers=8)
        ```
        
        Process only summer months:
        
        ```python
        results = ep.process(
            output_dir='./output',
            months=[6, 7, 8]  # June, July, August
        )
        ```
        
        Save input data for debugging:
        
        ```python
        results = ep.process(
            output_dir='./output',
            save_inputs=True,
            input_dir='./inputs'
        )
        ```

        See Also
        --------
            process_sequential: Sequential processing for debugging.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up input directory for saving downloaded data
        if save_inputs:
            if input_dir is not None:
                self._input_dir = Path(input_dir)
            else:
                # Default: parallel to output_dir in analysis_inputs
                self._input_dir = output_dir.parent / 'analysis_inputs' / output_dir.name
            self._input_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Input data will be saved to: {self._input_dir}")
        else:
            self._input_dir = None
        
        # Generate list of (year, month) to process
        all_dates = get_monthly_dates(self.start_year, self.end_year)
        
        if months is not None:
            all_dates = [(y, m) for y, m in all_dates if m in months]
        
        # Drop months that have no local precipitation file
        all_dates = self._filter_local_dates(all_dates)
        
        logger.info(f"Processing {len(all_dates)} months with {n_workers} workers")
        
        # Create delayed tasks
        tasks = [
            delayed(self._process_single_month)(year, month, output_dir)
            for year, month in all_dates
        ]
        
        # Execute in parallel with progress bar
        with ProgressBar():
            results = compute(*tasks, num_workers=n_workers)
        
        return list(results)
    
    def process_sequential(
        self,
        output_dir: Union[str, Path],
        months: Optional[List[int]] = None,
        input_dir: Optional[Union[str, Path]] = None,
        save_inputs: bool = False
    ) -> List[Tuple[Optional[Path], Optional[Path]]]:
        """
        Process all months sequentially (useful for debugging).
        
        Same as :meth:`process` but without parallel processing. Useful for
        debugging issues, testing on small datasets, or when GEE rate limits
        are a concern. Also the simplest way to run a local precipitation
        workflow with no Earth Engine involvement.
        
        Parameters
        ----------
        
        output_dir : str or Path
            Directory to save output rasters. Will be created if it
            doesn't exist.
        
        months : list of int, optional
            List of months to process (1-12). If None, processes all months
            in the date range.
        
        input_dir : str or Path, optional
            Directory to save downloaded input data (precipitation, AWC, ETo).
            If None and save_inputs is True, uses ``output_dir/../analysis_inputs``.
        
        save_inputs : bool, optional
            Whether to save downloaded input data. Default is False.
            
        Returns
        -------
        list of tuple
            List of tuples containing paths to saved files:
            ``(effective_precip_path, effective_precip_fraction_path)``.
            Returns ``(None, None)`` for months that failed to process.
            
        Examples
        --------
        Debug a single month:
        
        ```python
        ep = EffectivePrecipitation(...)
        results = ep.process_sequential(
            output_dir='./output',
            months=[1]  # Process only January
        )
        ```

        See Also
        --------
            process: Parallel processing method (recommended for production).
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up input directory for saving downloaded data
        if save_inputs:
            if input_dir is not None:
                self._input_dir = Path(input_dir)
            else:
                # Default: parallel to output_dir in analysis_inputs
                self._input_dir = output_dir.parent / 'analysis_inputs' / output_dir.name
            self._input_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Input data will be saved to: {self._input_dir}")
        else:
            self._input_dir = None
        
        all_dates = get_monthly_dates(self.start_year, self.end_year)
        
        if months is not None:
            all_dates = [(y, m) for y, m in all_dates if m in months]
        
        # Drop months that have no local precipitation file
        all_dates = self._filter_local_dates(all_dates)
        
        logger.info(f"Processing {len(all_dates)} months sequentially")
        
        results = []
        for year, month in all_dates:
            result = self._process_single_month(year, month, output_dir)
            results.append(result)
        
        return results

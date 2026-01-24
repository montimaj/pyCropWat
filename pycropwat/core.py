"""
Core module for effective precipitation calculations using Google Earth Engine.

This module provides the main :class:`EffectivePrecipitation` class for calculating
effective precipitation from various climate datasets available on Google Earth Engine.

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

See Also
--------
pycropwat.methods : Individual effective precipitation calculation functions.
pycropwat.analysis : Post-processing and analysis tools.
pycropwat.utils : Utility functions for GEE and file operations.
"""

import logging
from pathlib import Path
from typing import Union, Optional, List, Tuple
import ee
import numpy as np
import xarray as xr
import rioxarray
import dask
from dask import delayed, compute
from dask.diagnostics import ProgressBar

from .utils import load_geometry, get_date_range, get_monthly_dates, initialize_gee
from .methods import (
    get_method_function,
    PeffMethod
)

logger = logging.getLogger(__name__)

# Maximum pixels per tile for GEE sampleRectangle (conservative limit)
MAX_PIXELS_PER_TILE = 65536  # 256 x 256

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
    """
    Calculate effective precipitation from GEE climate data.
    
    Supports multiple effective precipitation calculation methods including
    CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall, FarmWest,
    and USDA-SCS (which requires AWC and ETo data).
    
    Parameters
    ----------
    
    asset_id : str
        GEE ImageCollection asset ID for precipitation data.
        Common options: 
        
        * ``ECMWF/ERA5_LAND/MONTHLY_AGGR`` (ERA5-Land, global, ~11km),
        * ``IDAHO_EPSCOR/TERRACLIMATE`` (TerraClimate, global, ~4km),
        * ``IDAHO_EPSCOR/GRIDMET`` (GridMET, CONUS, ~4km),
        * ``OREGONSTATE/PRISM/AN81m`` (PRISM, CONUS, ~4km),
        * ``UCSB-CHG/CHIRPS/DAILY`` (CHIRPS, 50°S-50°N, ~5km),
        * ``NASA/GPM_L3/IMERG_MONTHLY_V06`` (GPM IMERG, global, ~11km).
    
    precip_band : str
        Name of the precipitation band in the asset. Examples:
        
        * ERA5-Land: ``total_precipitation_sum``
        * TerraClimate: ``pr``
        * GridMET: ``pr``
        * PRISM: ``ppt``
        * CHIRPS: ``precipitation``
        * GPM IMERG: ``precipitation``
        
    geometry_path : str, Path, or None
        Path to shapefile or GeoJSON file defining the region of interest.
        Can also be a GEE FeatureCollection asset ID. Set to None if using
        gee_geometry_asset instead.

    start_year : int
        Start year for processing (inclusive).
    
    end_year : int
        End year for processing (inclusive).
    
    scale : float, optional
        Output resolution in meters. If None (default), uses native resolution
        of the dataset.
    
    precip_scale_factor : float, optional
        Factor to convert precipitation to mm. Default is 1.0.
        Common values: ERA5-Land (m to mm) = 1000, TerraClimate = 1.0, GridMET = 1.0.
    
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
            - ``eto_scale_factor`` (float): Scale factor for ETo. Default 1.0.
            - ``rooting_depth`` (float): Rooting depth in meters. Default 1.0.
            - ``mad_factor`` (float): Management Allowed Depletion factor (0-1).
              Controls what fraction of soil water storage is available. Default 0.5.
        
    Attributes
    ----------
    geometry : ee.Geometry
        The loaded geometry for the region of interest.
    
    collection : ee.ImageCollection
        The filtered and scaled precipitation image collection.
    
    bounds : list
        Bounding box coordinates of the geometry.
        
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

    See Also
    --------
    pycropwat.methods : Individual effective precipitation calculation functions.
    pycropwat.analysis : Post-processing and analysis tools.
    """
    
    def __init__(
        self,
        asset_id: str,
        precip_band: str,
        geometry_path: Optional[Union[str, Path]] = None,
        start_year: int = None,
        end_year: int = None,
        scale: Optional[float] = None,
        precip_scale_factor: float = 1.0,
        gee_project: Optional[str] = None,
        gee_geometry_asset: Optional[str] = None,
        method: PeffMethod = 'ensemble',
        method_params: Optional[dict] = None,
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
        
        # Get the effective precipitation function
        self._peff_function = get_method_function(method)
        
        # USDA-SCS specific: cache for AWC data (loaded once)
        self._awc_cache = None
        
        # Input directory for saving downloaded data (set during process())
        self._input_dir = None
        
        # Check if this is PCML method (uses single multi-band Image instead of ImageCollection)
        self._is_pcml = (method == 'pcml' or self.precip_band == PCML_DEFAULT_BAND)
        
        # For PCML, use default asset if placeholder provided
        if self._is_pcml:
            if self.asset_id == 'PLACEHOLDER' or self.asset_id is None:
                self.asset_id = PCML_DEFAULT_ASSET
                logger.info(f"Using default PCML asset: {self.asset_id}")
            self.precip_band = PCML_DEFAULT_BAND
        
        # Validate that at least one geometry source is provided (not required for PCML)
        if geometry_path is None and gee_geometry_asset is None and not self._is_pcml:
            raise ValueError("Either geometry_path or gee_geometry_asset must be provided")
        
        # Initialize GEE
        initialize_gee(self.gee_project)
        
        # For PCML, use the asset's own geometry if no geometry provided
        if self._is_pcml and geometry_path is None and gee_geometry_asset is None:
            # Load PCML image first
            self._pcml_image = ee.Image(self.asset_id)
            # Use predefined Western U.S. bounding box since PCML image geometry is unbounded
            self.geometry = ee.Geometry.Polygon([PCML_WESTERN_US_BOUNDS])
            self.bounds = PCML_WESTERN_US_BOUNDS
            logger.info("Using predefined Western U.S. bounding box for PCML")
        else:
            # Load geometry from GEE asset or local file
            self.geometry = load_geometry(geometry_path, gee_asset=gee_geometry_asset)
            self.bounds = self.geometry.bounds().getInfo()['coordinates'][0]
        
        # Get date range
        self.start_date, self.end_date = get_date_range(start_year, end_year)
        
        # Load and filter image collection (or load PCML image)
        self._load_collection()
        
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
    
    def _download_image_chunked(
        self,
        img: ee.Image,
        bounds_coords: List,
        scale_meters: float,
        band_name: str,
        default_value: float,
        target_shape: Optional[tuple] = None,
        data_name: str = "data"
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
            
        Returns
        -------
        np.ndarray
            Mosaicked array, resized to target_shape if specified.
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
                logger.warning(f"All {data_name} tiles failed, using default {default_value}")
                if target_shape:
                    return np.full(target_shape, default_value, dtype=np.float32)
                else:
                    raise ValueError(f"No tiles downloaded and no target_shape specified")
            
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
            
        except Exception as e:
            logger.warning(f"{data_name} chunked download failed: {e}. Using default {default_value}")
            if target_shape:
                return np.full(target_shape, default_value, dtype=np.float32)
            else:
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
        Download monthly precipitation data from GEE, using chunked download if needed.
        
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
            Precipitation data array, or None if download fails.
        """
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
        If loading fails, a default value of 0.15 (15% AWC) is used.
        
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
                    target_shape=template_da.shape, data_name="AWC"
                )
            
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
                    awc_da = awc_da.rio.write_crs("EPSG:4326")
                    awc_da.rio.to_raster(awc_path)
                    logger.info(f"Saved input AWC: {awc_path.name}")
            
            return awc_arr
            
        except Exception as e:
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
        
        If loading fails, a default value of 100 mm/month is used.
        
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
                    target_shape=template_da.shape, data_name="ETo"
                )
            
            # Apply scale factor if needed (for both direct and chunked downloads)
            if eto_scale_factor != 1.0:
                eto_arr = eto_arr * eto_scale_factor
                logger.debug(f"Applied ETo scale factor: {eto_scale_factor}")
            
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
                    eto_da = eto_da.rio.write_crs("EPSG:4326")
                    eto_da.rio.to_raster(eto_path)
                    logger.info(f"Saved input ETo: {eto_path.name}")
            
            return eto_arr
            
        except Exception as e:
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
                target_shape=None, data_name="precipitation"
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
        """
        logger.info(f"Processing {year}-{month:02d}")
        
        # Download precipitation data
        pr_da = self._download_monthly_precip(year, month)
        
        if pr_da is None:
            return None, None
        
        # Calculate effective precipitation using the configured method
        if self.method == 'pcml':
            # PCML: pr_da already contains effective precipitation (PCML Peff)
            # Download annual fraction directly from the GEE asset (once per year)
            ep_arr = pr_da.values  # PCML Peff is directly downloaded
            
            # For PCML, check if annual fraction file already exists
            pcml_frac_path = output_dir / f"effective_precip_fraction_{year}.tif"
            if pcml_frac_path.exists():
                # Load existing annual fraction
                epf_da_existing = rioxarray.open_rasterio(pcml_frac_path).squeeze('band', drop=True)
                epf_arr = epf_da_existing.values
                logger.info(f"PCML annual fraction for {year} already exists, reusing")
            else:
                epf_arr = self._download_pcml_annual_fraction(year, pr_da)
                logger.info(f"PCML annual fraction loaded from GEE asset for {year}")
        elif self.method == 'usda_scs':
            # USDA-SCS method requires AWC and ETo data
            awc_arr = self._load_awc_data(pr_da)
            eto_arr = self._load_monthly_eto(year, month, pr_da)
            rooting_depth = self.method_params.get('rooting_depth', 1.0)
            mad_factor = self.method_params.get('mad_factor', 0.5)
            
            ep_arr = self._peff_function(
                pr_da.values,
                eto_arr,
                awc_arr,
                rooting_depth,
                mad_factor
            )
        elif self.method == 'ensemble':
            # Ensemble method requires AWC and ETo data (same as USDA-SCS)
            awc_arr = self._load_awc_data(pr_da)
            eto_arr = self._load_monthly_eto(year, month, pr_da)
            rooting_depth = self.method_params.get('rooting_depth', 1.0)
            fixed_percentage = self.method_params.get('percentage', 0.7)
            dependable_probability = self.method_params.get('probability', 0.75)
            
            ep_arr = self._peff_function(
                pr_da.values,
                eto_arr,
                awc_arr,
                rooting_depth,
                fixed_percentage,
                dependable_probability
            )
        elif self.method == 'suet':
            # SuET method requires ETo data
            eto_arr = self._load_monthly_eto(year, month, pr_da)
            ep_arr = self._peff_function(pr_da.values, eto_arr)
        elif self.method_params:
            # Other methods with parameters (fixed_percentage, dependable_rainfall)
            # Filter to only pass valid parameters for the method
            valid_params = {}
            if self.method == 'fixed_percentage':
                valid_params['percentage'] = self.method_params.get('percentage', 0.7)
            elif self.method == 'dependable_rainfall':
                valid_params['probability'] = self.method_params.get('probability', 0.75)
            ep_arr = self._peff_function(pr_da.values, **valid_params)
        else:
            ep_arr = self._peff_function(pr_da.values)
        
        # Calculate effective precipitation fraction (skip for PCML - already calculated above)
        if self.method != 'pcml':
            with np.errstate(divide='ignore', invalid='ignore'):
                epf_arr = np.where(pr_da.values > 0, ep_arr / pr_da.values, 0)
        
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
        ep_da = ep_da.rio.write_crs("EPSG:4326")
        
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
        epf_da = epf_da.rio.write_crs("EPSG:4326")
        
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
        
        Downloads precipitation data from Google Earth Engine, calculates
        effective precipitation using the configured method, and saves
        results as GeoTIFF files. Uses Dask for parallel processing of
        multiple months.
        
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
        and cached for efficiency.
        
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
        are a concern.
        
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
        
        results = []
        for year, month in all_dates:
            result = self._process_single_month(year, month, output_dir)
            results.append(result)
        
        return results

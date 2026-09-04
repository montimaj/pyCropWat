"""
pyCropWat - Calculate effective precipitation from GEE or local climate data.

pyCropWat is a Python package for calculating effective precipitation using
various methods from climate data available on Google Earth Engine (GEE), or
from precipitation rasters and NetCDF files you already have on disk. It
supports multiple precipitation datasets and effective precipitation
calculation methods.

Main Features
-------------
- Calculate effective precipitation from GEE climate datasets
- Calculate effective precipitation from local GeoTIFF/NetCDF precipitation
- Support for multiple methods (CROPWAT, FAO/AGLW, USDA-SCS, etc.)
- Temporal aggregation (seasonal, annual, growing season)
- Statistical analysis (trends, anomalies, climatology)
- Publication-quality visualization
- Export to NetCDF and Cloud-Optimized GeoTIFF formats

Quick Start
-----------
>>> from pycropwat import EffectivePrecipitation
>>> 
>>> # Calculate effective precipitation using ERA5-Land data
>>> ep = EffectivePrecipitation(
...     asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
...     precip_band='total_precipitation_sum',
...     geometry_path='study_area.geojson',
...     start_year=2015,
...     end_year=2020,
...     precip_scale_factor=1000,  # Convert m to mm
...     method='ensemble'
... )
>>> results = ep.process(output_dir='./output', n_workers=4)

Local Precipitation Input
-------------------------
Pass ``local_precip`` to read monthly precipitation from disk instead of
downloading it from Earth Engine. It accepts a directory of monthly rasters, a
single NetCDF file, or a glob string:

>>> ep = EffectivePrecipitation(
...     local_precip='../pyCropWat_Data/Precip',
...     local_precip_pattern='Precip_*.tif',
...     local_precip_nodata=-9999,
...     method='cropwat'
... )
>>> results = ep.process_sequential(output_dir='./output', months=[7, 8])

``start_year``/``end_year`` are inferred from the files when omitted, and nodata
pixels are propagated as NaN to both outputs. Available water capacity (AWC),
reference evapotranspiration (ETo) and GEE ``FeatureCollection`` geometries are
still downloaded from Earth Engine, so ``usda_scs``, ``ensemble`` and ``suet``
still require GEE access. The precipitation-only methods (``cropwat``,
``fao_aglw``, ``fixed_percentage``, ``dependable_rainfall``, ``farmwest``) need
no Earth Engine at all - it is never initialised and the run works offline. The
``pcml`` method is a pre-computed GEE product and cannot use local precipitation.

Supported Precipitation Datasets
--------------------------------
- ERA5-Land (global, ~11km): ``'ECMWF/ERA5_LAND/MONTHLY_AGGR'``
- TerraClimate (global, ~4km): ``'IDAHO_EPSCOR/TERRACLIMATE'``
- GridMET (CONUS, ~4km): ``'IDAHO_EPSCOR/GRIDMET'``
- PRISM (CONUS, ~4km): ``'OREGONSTATE/PRISM/AN81m'``
- CHIRPS (50°S-50°N, ~5km): ``'UCSB-CHG/CHIRPS/DAILY'``
- GPM IMERG (global, ~11km): ``'NASA/GPM_L3/IMERG_MONTHLY_V06'``

Effective Precipitation Methods
-------------------------------
- ``'ensemble'`` - Ensemble mean of 6 methods (default, requires AWC and ETo)
- ``'cropwat'`` - CROPWAT method (FAO standard)
- ``'fao_aglw'`` - FAO Dependable Rainfall (80% exceedance)
- ``'fixed_percentage'`` - Simple fixed percentage method
- ``'dependable_rainfall'`` - FAO Dependable Rainfall method
- ``'farmwest'`` - FarmWest method
- ``'usda_scs'`` - USDA-SCS method (requires AWC and ETo assets)
- ``'suet'`` - TAGEM-SuET method (requires ETo asset)
- ``'pcml'`` - Physics-Constrained ML (Western U.S. only, pre-computed GEE asset)

Modules
-------
core
    Main :class:`EffectivePrecipitation` class for calculations.
local
    :class:`LocalPrecipitationSource` and helpers for reading local
    precipitation rasters and NetCDF files.
methods
    Individual effective precipitation calculation functions.
analysis
    :class:`TemporalAggregator`, :class:`StatisticalAnalyzer`, 
    :class:`Visualizer`, and export functions.
utils
    Utility functions for GEE and file operations.

See Also
--------
Documentation: https://pycropwat.readthedocs.io/
GitHub: https://github.com/username/pycropwat
"""

from .core import EffectivePrecipitation
from .local import (
    LocalPrecipitationSource,
    open_local_precipitation,
    parse_year_month,
)
from .utils import load_geometry, load_geometry_from_gee_asset, get_date_range, is_gee_asset
from .methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip,
    farmwest_effective_precip,
    usda_scs_effective_precip,
    suet_effective_precip,
    pcml_effective_precip,
    ensemble_effective_precip,
    get_method_function,
    list_available_methods,
)
from .analysis import (
    TemporalAggregator,
    StatisticalAnalyzer,
    Visualizer,
    export_to_netcdf,
    export_to_cog,
)

__version__ = "1.3.0"
__all__ = [
    # Core
    "EffectivePrecipitation",
    # Local precipitation input
    "LocalPrecipitationSource",
    "open_local_precipitation",
    "parse_year_month",
    # Utils
    "load_geometry",
    "load_geometry_from_gee_asset",
    "get_date_range",
    "is_gee_asset",
    # Methods
    "cropwat_effective_precip",
    "fao_aglw_effective_precip",
    "fixed_percentage_effective_precip",
    "dependable_rainfall_effective_precip",
    "farmwest_effective_precip",
    "usda_scs_effective_precip",
    "suet_effective_precip",
    "pcml_effective_precip",
    "ensemble_effective_precip",
    "get_method_function",
    "list_available_methods",
    # Analysis
    "TemporalAggregator",
    "StatisticalAnalyzer",
    "Visualizer",
    "export_to_netcdf",
    "export_to_cog",
]

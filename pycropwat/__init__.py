"""
pyCropWat - Calculate effective precipitation from Google Earth Engine climate data.

pyCropWat is a Python package for calculating effective precipitation using
various methods from climate data available on Google Earth Engine (GEE).
It supports multiple precipitation datasets and effective precipitation
calculation methods.

Main Features
-------------
- Calculate effective precipitation from GEE climate datasets
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
- ``'usda_scs'`` - USDA-SCS soil moisture depletion method
- ``'suet'`` - TAGEM-SuET method (Turkish Irrigation Management System)

Modules
-------
core
    Main :class:`EffectivePrecipitation` class for calculations.
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
from .utils import load_geometry, load_geometry_from_gee_asset, get_date_range, is_gee_asset
from .methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip,
    farmwest_effective_precip,
    usda_scs_effective_precip,
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

__version__ = "1.2"
__all__ = [
    # Core
    "EffectivePrecipitation",
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
    "get_method_function",
    "list_available_methods",
    # Analysis
    "TemporalAggregator",
    "StatisticalAnalyzer",
    "Visualizer",
    "export_to_netcdf",
    "export_to_cog",
]

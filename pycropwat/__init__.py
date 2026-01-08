"""
pyCropWat - A Python package for calculating CROPWAT effective precipitation from GEE climate data.
"""

from .core import EffectivePrecipitation
from .utils import load_geometry, load_geometry_from_gee_asset, get_date_range, is_gee_asset
from .methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip,
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

__version__ = "1.0.0"
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
    "get_method_function",
    "list_available_methods",
    # Analysis
    "TemporalAggregator",
    "StatisticalAnalyzer",
    "Visualizer",
    "export_to_netcdf",
    "export_to_cog",
]

"""
pyCropWat - A Python package for calculating CROPWAT effective precipitation from GEE climate data.
"""

from .core import EffectivePrecipitation
from .utils import load_geometry, load_geometry_from_gee_asset, get_date_range, is_gee_asset

__version__ = "0.1.0"
__all__ = [
    "EffectivePrecipitation",
    "load_geometry",
    "load_geometry_from_gee_asset",
    "get_date_range",
    "is_gee_asset"
]

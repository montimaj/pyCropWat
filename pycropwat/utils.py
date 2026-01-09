"""
Utility functions for pyCropWat package.

This module provides utility functions for working with Google Earth Engine,
loading geometries from various sources, and managing date ranges.

Functions
---------
initialize_gee
    Initialize Google Earth Engine with optional project ID.
    
load_geometry
    Load geometry from local file (shapefile/GeoJSON) or GEE asset.
    
load_geometry_from_file
    Load geometry from local shapefile or GeoJSON file.
    
load_geometry_from_gee_asset
    Load geometry from GEE FeatureCollection asset.
    
is_gee_asset
    Check if a path string represents a GEE asset.
    
get_date_range
    Generate date range strings for GEE filtering.
    
get_monthly_dates
    Generate list of (year, month) tuples for iteration.

Example
-------
```python
from pycropwat.utils import initialize_gee, load_geometry

# Initialize GEE with project
initialize_gee(project='my-gee-project')

# Load geometry from local file
geom = load_geometry('study_area.geojson')

# Load geometry from GEE asset
geom = load_geometry(gee_asset='projects/my-project/assets/boundary')
```
"""

from pathlib import Path
from typing import Union, Tuple, List, Optional
import ee
import geopandas as gpd
from shapely.geometry import mapping


def is_gee_asset(path: str) -> bool:
    """
    Check if a path looks like a GEE asset path.
    
    Parameters
    ----------
    
    path : str
        Path to check.
        
    Returns
    -------
    bool
        True if the path appears to be a GEE asset path.
    """
    # GEE assets typically start with 'projects/', 'users/', or contain '/'
    # without a file extension
    path_str = str(path)
    if path_str.startswith(('projects/', 'users/', 'TIGER/', 'USDOS/', 'FAO/', 'WCMC/', 'RESOLVE/')):
        return True
    # Check if it has no file extension and contains forward slashes
    if '/' in path_str and not Path(path_str).suffix:
        return True
    return False


def load_geometry_from_gee_asset(asset_id: str) -> ee.Geometry:
    """
    Load a geometry from a GEE FeatureCollection asset.
    
    Parameters
    ----------
    
    asset_id : str
        GEE FeatureCollection asset ID (e.g., 'projects/my-project/assets/my_boundary').
        
    Returns
    -------
    ee.Geometry
        Earth Engine Geometry object (dissolved/union of all features).
    """
    fc = ee.FeatureCollection(asset_id)
    # Dissolve all features into a single geometry
    geometry = fc.geometry().dissolve()
    return geometry


def load_geometry_from_file(geometry_path: Union[str, Path]) -> ee.Geometry:
    """
    Load a geometry from a local shapefile or GeoJSON file.
    
    Parameters
    ----------
    
    geometry_path : str or Path
        Path to the shapefile (.shp) or GeoJSON (.geojson, .json) file.
        
    Returns
    -------
    ee.Geometry
        Earth Engine Geometry object.
    """
    geometry_path = Path(geometry_path)
    
    if not geometry_path.exists():
        raise FileNotFoundError(f"Geometry file not found: {geometry_path}")
    
    suffix = geometry_path.suffix.lower()
    
    if suffix == ".shp":
        gdf = gpd.read_file(geometry_path)
    elif suffix in [".geojson", ".json"]:
        gdf = gpd.read_file(geometry_path)
    else:
        raise ValueError(f"Unsupported file format: {suffix}. Use .shp, .geojson, or .json")
    
    # Ensure CRS is WGS84
    if gdf.crs is not None and gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)
    
    # Dissolve all geometries into a single geometry
    dissolved = gdf.dissolve()
    geom = dissolved.geometry.values[0]
    
    # Convert to GeoJSON and then to ee.Geometry
    geojson = mapping(geom)
    return ee.Geometry(geojson)


def load_geometry(
    geometry_source: Union[str, Path],
    gee_asset: Optional[str] = None
) -> ee.Geometry:
    """
    Load a geometry from a local file or GEE FeatureCollection asset.
    
    Parameters
    ----------
    
    geometry_source : str or Path
        Path to a local shapefile (.shp) or GeoJSON (.geojson, .json) file,
        OR a GEE FeatureCollection asset ID.
    
    gee_asset : str, optional
        Explicit GEE FeatureCollection asset ID. If provided, this takes
        precedence over geometry_source.
        
    Returns
    -------
    ee.Geometry
        Earth Engine Geometry object.
        
    Examples
    --------
    Load from local file:
    
    ```python
    geom = load_geometry('boundary.geojson')
    ```
    
    Load from GEE asset:
    
    ```python
    geom = load_geometry(gee_asset='projects/my-project/assets/boundary')
    ```
    
    Auto-detect GEE asset from path:
    
    ```python
    geom = load_geometry('users/username/my_boundary')
    ```
    """
    # If explicit GEE asset provided, use it
    if gee_asset is not None:
        return load_geometry_from_gee_asset(gee_asset)
    
    # Check if geometry_source looks like a GEE asset
    if is_gee_asset(str(geometry_source)):
        return load_geometry_from_gee_asset(str(geometry_source))
    
    # Otherwise, treat as local file
    return load_geometry_from_file(geometry_source)


def get_date_range(start_year: int, end_year: int) -> Tuple[str, str]:
    """
    Generate date range strings for filtering.
    
    Parameters
    ----------
    
    start_year : int
        Start year (inclusive).
    
    end_year : int
        End year (inclusive).
        
    Returns
    -------
    tuple
        Tuple of (start_date, end_date) strings in 'YYYY-MM-DD' format.
    """
    start_date = f"{start_year}-01-01"
    end_date = f"{end_year + 1}-01-01"
    return start_date, end_date


def get_monthly_dates(start_year: int, end_year: int) -> List[Tuple[int, int]]:
    """
    Generate list of (year, month) tuples for the given date range.
    
    Parameters
    ----------
    
    start_year : int
        Start year (inclusive).
    
    end_year : int
        End year (inclusive).
        
    Returns
    -------
    list
        List of (year, month) tuples.
    """
    dates = []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            dates.append((year, month))
    return dates


def initialize_gee(project: str = None) -> None:
    """
    Initialize Google Earth Engine.
    
    Parameters
    ----------
    
    project : str, optional
        GEE project ID. If None, uses default authentication.
    """
    try:
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
    except Exception:
        ee.Authenticate()
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()

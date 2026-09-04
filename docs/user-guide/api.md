# Python API

This guide covers the Python API for pyCropWat.

## EffectivePrecipitation Class

The main interface for calculating effective precipitation.

### Basic Usage

```python
from pycropwat import EffectivePrecipitation

# Initialize with geometry file
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2020,
    end_year=2023,
    precip_scale_factor=1000
)

# Process all months in parallel
results = ep.process(output_dir='./outputs', n_workers=4)
```

### Common Precipitation Datasets

```python
# ERA5-Land Monthly (~11 km global)
ep_era5 = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2020, end_year=2023,
    precip_scale_factor=1000  # meters to mm
)

# TerraClimate (~4 km global)
ep_terra = EffectivePrecipitation(
    asset_id='IDAHO_EPSCOR/TERRACLIMATE',
    precip_band='pr',
    geometry_path='study_area.geojson',
    start_year=2020, end_year=2023,
    precip_scale_factor=1  # already in mm
)

# GridMET (~4 km U.S. only)
ep_gridmet = EffectivePrecipitation(
    asset_id='IDAHO_EPSCOR/GRIDMET',
    precip_band='pr',
    geometry_path='study_area.geojson',
    start_year=2020, end_year=2023,
    precip_scale_factor=1
)

# PRISM (~4 km U.S. only)
ep_prism = EffectivePrecipitation(
    asset_id='OREGONSTATE/PRISM/AN81m',
    precip_band='ppt',
    geometry_path='study_area.geojson',
    start_year=2020, end_year=2023,
    precip_scale_factor=1
)

# CHIRPS (~5.5 km global, 50°S-50°N)
ep_chirps = EffectivePrecipitation(
    asset_id='UCSB-CHG/CHIRPS/PENTAD',
    precip_band='precipitation',
    geometry_path='study_area.geojson',
    start_year=2020, end_year=2023,
    precip_scale_factor=1
)
```

### Using GEE Vector Asset

```python
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    gee_geometry_asset='projects/my-project/assets/study_area',
    start_year=2020,
    end_year=2023,
    precip_scale_factor=1000
)
```

### Local Precipitation Input

Pass `local_precip` to read monthly precipitation from disk instead of downloading it from Earth
Engine. `asset_id` and `precip_band` are then unused (they are `Optional[str] = None` since v1.3.0),
and `start_year`/`end_year` are inferred from the files when omitted:

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',   # directory, NetCDF file, or glob
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    method='cropwat'
)

# No Earth Engine at all for precipitation-only methods
results = ep.process(output_dir='./outputs', n_workers=4)
```

#### New Constructor Arguments

All local-precipitation arguments are appended after the existing ones, so existing positional and
keyword calls are unchanged.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `local_precip` | `str` or `Path` | `None` | Local precipitation to use **instead of** GEE: a directory of monthly rasters, a single NetCDF file (or a directory of them), or a glob string such as `'./Precip/Precip_*.tif'`. Not compatible with `method='pcml'` |
| `local_precip_pattern` | `str` | `'*.tif'` | Glob used when `local_precip` is a directory |
| `local_precip_variable` | `str` | `None` | NetCDF variable holding precipitation. Auto-detected when `None` |
| `local_precip_nodata` | `float` | `None` | Additional nodata sentinel, applied on top of the value stored in the file metadata (e.g. `-9999`). Matching pixels become NaN and propagate as NaN to both output rasters |
| `local_precip_crs` | `str` | `None` | CRS **override**, e.g. `'EPSG:4326'`. Replaces whatever CRS the files declare (logged at INFO) and supplies one when they declare none. It relabels the grid in place - values and transform are untouched, nothing is reprojected |
| `local_precip_date_regex` | `str` | `None` | Regular expression with named groups `year` and `month`, used to date files the built-in parser cannot read, e.g. `r'(?P<month>\d{2})_(?P<year>\d{4})'` |
| `clip_to_geometry` | `bool` | `True` | Clip local precipitation rasters to `geometry_path` when it points at a local vector file. No effect on GEE precipitation, which is always clipped server-side. Setting it to `False` keeps the full grid, but pixels outside the AWC/ETo download region then stay NaN in `usda_scs`/`ensemble`/`suet` outputs |

`precip_scale_factor` applies to local files too, so use it to convert your units to millimetres per
month (e.g. `1000` for metres). It is a single multiplier, so it cannot turn a rate such as
kg m⁻² s⁻¹ into a monthly total - accumulate those to monthly totals before pyCropWat reads them.

Related changes to existing arguments:

| Parameter | Change |
|-----------|--------|
| `asset_id` | Now `Optional[str] = None`. Required only when `local_precip` is `None` and `method != 'pcml'` |
| `precip_band` | Now `Optional[str] = None`. Same rule as `asset_id` |
| `geometry_path` / `gee_geometry_asset` | Optional with `local_precip`. Without either, the extent of the local files is used |

Omitting a precipitation source altogether raises `ValueError("No precipitation source provided...")`.

#### What Still Comes From Earth Engine

| Method | Earth Engine required with `local_precip`? |
|--------|--------------------------------------------|
| `cropwat`, `fao_aglw`, `fixed_percentage`, `dependable_rainfall`, `farmwest` | **No** - `ee.Initialize()` is never called; runs offline |
| `usda_scs`, `ensemble`, `suet` | **Yes** - AWC and/or ETo are still downloaded from GEE |
| `pcml` | Not supported - raises `ValueError` (pre-computed GEE product) |

A GEE `FeatureCollection` geometry (`gee_geometry_asset`, or a `geometry_path` that is a GEE asset)
also requires Earth Engine, regardless of method.

Local precipitation with AWC and ETo still coming from Earth Engine:

```python
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    geometry_path='basin.geojson',      # clips the local rasters
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
ep.process(output_dir='./outputs', n_workers=4)
```

Local NetCDF input with an explicit variable and CRS:

```python
ep = EffectivePrecipitation(
    local_precip='./wrf_precip.nc',
    local_precip_variable='RAINNC',
    local_precip_crs='EPSG:4326',   # override: relabels the grid, never reprojects
    method='fao_aglw'
)
```

!!! warning "AWC/ETo coverage"
    `usda_scs`, `ensemble` and `suet` still download AWC and ETo from Earth Engine for the
    geometry you give them. Any part of the local precipitation grid the download did not cover
    stays **NaN** in the outputs - pyCropWat does not backfill it - and the run logs a WARNING
    naming the uncovered pixel count, once per field per process (twice for `usda_scs` and
    `ensemble`: once for AWC, once for ETo):

    ```
    AWC covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
    ```

    Clipping to the geometry (`clip_to_geometry=True`, the default) or widening the geometry gives
    full coverage.

See [Local Precipitation Data](local-data.md) for the complete guide.

### Custom Resolution

```python
ep = EffectivePrecipitation(
    asset_id='IDAHO_EPSCOR/TERRACLIMATE',
    precip_band='pr',
    geometry_path='study_area.geojson',
    start_year=2020,
    end_year=2023,
    scale=5000  # 5 km resolution in meters
)
```

### Processing Options

```python
# Process all months with parallel workers
results = ep.process(output_dir='./outputs', n_workers=4)

# Process specific months only
results = ep.process(
    output_dir='./outputs',
    n_workers=4,
    months=[6, 7, 8, 9]  # June through September
)

# Sequential processing (useful for debugging)
results = ep.process_sequential(output_dir='./outputs')
```

## LocalPrecipitationSource Class

`pycropwat.local.LocalPrecipitationSource` is the reader behind `local_precip`. Use it directly when
you want to inspect a dataset, or to pull monthly grids without running the full workflow.

```python
from pycropwat import LocalPrecipitationSource, open_local_precipitation, parse_year_month
```

### Constructor

```python
LocalPrecipitationSource(
    path: Union[str, Path],                                   # directory of rasters, a raster/NetCDF file, or a glob string
    pattern: str = '*.tif',                                   # glob applied when `path` is a directory
    variable: Optional[str] = None,                           # NetCDF variable; auto-detected when None
    scale_factor: float = 1.0,                                # multiplier applied to every returned array (source units -> mm)
    nodata: Optional[Union[float, Sequence[float]]] = None,   # extra nodata sentinel(s), in addition to the file metadata
    crs: Optional[Union[str, CRS]] = None,                    # CRS override: relabels the grid, replacing any declared CRS
    date_regex: Optional[str] = None,                         # custom regex with named groups 'year' and 'month'
    time_dim: str = 'time',                                   # NetCDF time coordinate name
    x_dim: Optional[str] = None,                              # NetCDF X dimension; auto-detected when None
    y_dim: Optional[str] = None                               # NetCDF Y dimension; auto-detected when None
)
```

Files are indexed up front (cheap); pixel data is read only when `get_month()` is called.

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `kind` | `str` | `'raster'` for one-file-per-month rasters, `'netcdf'` otherwise |
| `files` | `List[Path]` | The data files backing this source, sorted |
| `crs` | `CRS` | Coordinate reference system of the source data |
| `bounds` | `Tuple[float, float, float, float]` | `(minx, miny, maxx, maxy)` in EPSG:4326 lon/lat degrees |
| `native_bounds` | `Tuple[float, float, float, float]` | `(minx, miny, maxx, maxy)` in the source CRS |
| `year_range` | `Tuple[int, int]` | `(first_year, last_year)` present in the source |
| `resolution` | `Tuple[float, float]` | `(x_resolution, y_resolution)` as positive values in CRS units |
| `shape` | `Tuple[int, int]` | `(n_rows, n_cols)` of every monthly grid |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `available_months()` | `List[Tuple[int, int]]` | Sorted `(year, month)` pairs present in the source |
| `has_month(year, month)` | `bool` | Whether that month can be read |
| `get_month(year, month)` | `Optional[xr.DataArray]` | The monthly grid, or `None` when the month is missing |
| `close()` | `None` | Release open NetCDF handles (a no-op in raster mode); safe to call twice |

`len(src)` gives the number of months available, and `repr(src)` summarises kind, files, months,
year range, CRS and shape.

Every array returned by `get_month()` has dims exactly `('y', 'x')`, dtype `float32`, values in
millimetres (`scale_factor` applied), nodata as `NaN`, a descending (north-up) `y` coordinate, the
source CRS attached via rioxarray, and attrs `{'units', 'long_name', 'year', 'month', 'source'}`.

### Usage

```python
from pycropwat import LocalPrecipitationSource

src = LocalPrecipitationSource(
    '../pyCropWat_Data/Precip',
    pattern='Precip_*.tif',
    nodata=-9999
)
print(src.kind, len(src), src.year_range)   # raster 264 (2000, 2021)
print(src.crs, src.shape, src.resolution)   # EPSG:4326 (689, 799) (0.0397..., 0.0397...)
print(src.bounds)                           # (-71.3018, -40.3510, -39.5434, -12.9648)

da = src.get_month(2005, 7)                 # xr.DataArray, dims ('y', 'x'), mm
print(da.dims, da.shape, float(da.mean()))
src.close()
```

The source is also a context manager, which is the safest way to use NetCDF inputs:

```python
from pycropwat import open_local_precipitation

# open_local_precipitation() is a thin factory; kwargs are forwarded unchanged
with open_local_precipitation('precip_2000_2020.nc',
                              variable='tp',
                              scale_factor=1000.0) as src:   # metres -> mm
    for year, month in src.available_months():
        da = src.get_month(year, month)
```

### parse_year_month()

```python
parse_year_month(name: Union[str, Path], date_regex: Optional[str] = None) -> Optional[Tuple[int, int]]
```

Extracts `(year, month)` from a file name, stem, path or arbitrary string. Built-in layouts are
`YYYY_MM`, `YYYY-MM`, `YYYYMM` and `YYYY.MM`; the last plausible match in the string wins, so
prefixes such as `Precip_` never confuse it. Returns `None` when nothing matched.

```python
from pycropwat import parse_year_month

parse_year_month('Precip_2005_07.tif')       # (2005, 7)
parse_year_month('pr200507')                 # (2005, 7)
parse_year_month('/data/2019/x-2005-07.nc')  # (2005, 7)
parse_year_month('no_date_here.tif')         # None

# Custom layout: MM_YYYY
parse_year_month('rain_07_2005.tif',
                 date_regex=r'(?P<month>\d{2})_(?P<year>\d{4})')   # (2005, 7)
```

`date_regex` must define both named groups, otherwise a `ValueError` is raised.

Full generated reference: [Local Data API](../api/local.md). Task-oriented guide:
[Local Precipitation Data](local-data.md).

## Effective Precipitation Methods

pyCropWat supports multiple methods for calculating effective precipitation:

```python
from pycropwat.methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip,
    farmwest_effective_precip,
    usda_scs_effective_precip,
    list_available_methods,
    get_method_function
)
import numpy as np

# List all available methods with descriptions
for method, description in list_available_methods().items():
    print(f"{method}: {description}")

# Get a method function by name
peff_func = get_method_function('cropwat')

precip = np.array([50, 100, 200, 300, 400])

# CROPWAT method
peff_cropwat = cropwat_effective_precip(precip)
# [46.  72.  136.  155.  165.]

# FAO/AGLW method
peff_fao = fao_aglw_effective_precip(precip)

# Fixed percentage (e.g., 80%)
peff_fixed = fixed_percentage_effective_precip(precip, percentage=0.8)
# [40.  80.  160.  240.  320.]

# Dependable rainfall at 75% probability
peff_depend = dependable_rainfall_effective_precip(precip, probability=0.75)

# FarmWest method
peff_farmwest = farmwest_effective_precip(precip)
# [33.75  71.25  146.25  221.25  296.25]

# USDA-SCS method (requires AWC and ETo arrays)
# AWC is in volumetric fraction (cm/cm or inch/inch)
# ETo is reference evapotranspiration in mm
eto = np.array([80, 120, 180, 220, 260])  # mm
awc = np.array([0.15, 0.15, 0.15, 0.15, 0.15])  # volumetric fraction
peff_usda = usda_scs_effective_precip(precip, eto, awc, rooting_depth=1.0, mad_factor=0.5)

# SuET method (requires ETo array)
from pycropwat.methods import suet_effective_precip
peff_suet = suet_effective_precip(precip, eto)
```

### Using USDA-SCS Method with GEE Assets

The USDA-SCS method can be used with GEE assets for AWC and ETo:

```python
# For U.S. regions (SSURGO AWC + GridMET ETo)
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='arizona.geojson',
    start_year=2015,
    end_year=2020,
    precip_scale_factor=1000,
    method='usda_scs',
    method_params={
        'awc_asset': 'projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite',
        'awc_band': 'AWC',
        'eto_asset': 'projects/openet/assets/reference_et/conus/gridmet/monthly/v1',
        'eto_band': 'eto',
        'rooting_depth': 1.0,
        'mad_factor': 0.5
    }
)

# For global regions (FAO HWSD AWC + AgERA5 ETo)
ep_global = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2015,
    end_year=2020,
    precip_scale_factor=1000,
    method='usda_scs',
    method_params={
        'awc_asset': 'projects/sat-io/open-datasets/FAO/HWSD_V2_SMU',
        'awc_band': 'AWC',
        'eto_asset': 'projects/climate-engine-pro/assets/ce-ag-era5-v2/daily',
        'eto_band': 'ReferenceET_PenmanMonteith_FAO56',
        'eto_is_daily': True,  # Will aggregate daily to monthly
        'rooting_depth': 1.0,
        'mad_factor': 0.5
    }
)
```

## Temporal Aggregation

The `TemporalAggregator` class provides functions for aggregating monthly data:

```python
from pycropwat.analysis import TemporalAggregator

agg = TemporalAggregator('./outputs')

# Annual totals
annual = agg.annual_aggregate(
    year=2020,
    method='sum',
    output_path='./annual_2020.tif'
)

# Seasonal aggregation (DJF, MAM, JJA, SON)
summer = agg.seasonal_aggregate(
    year=2020,
    season='JJA',
    method='sum',
    output_path='./summer_2020.tif'
)

# Growing season (customizable months)
# Northern Hemisphere (April-October, same year)
growing_nh = agg.growing_season_aggregate(
    year=2020,
    start_month=4,
    end_month=10,
    method='sum',
    output_path='./growing_2020.tif'
)

# Southern Hemisphere (October-March, cross-year)
# When start_month > end_month, automatically loads from two years
# This aggregates Oct 2020 - Mar 2021
growing_sh = agg.growing_season_aggregate(
    year=2020,
    start_month=10,
    end_month=3,
    method='sum',
    output_path='./growing_2020_2021.tif'
)

# Custom month range
custom = agg.custom_aggregate(
    year=2020,
    months=[5, 6, 7, 8],
    method='mean',
    output_path='./may_aug_mean_2020.tif'
)

# Multi-year climatology (per-month long-term mean)
climatology = agg.multi_year_climatology(
    start_year=1990,
    end_year=2020,
    months=None,  # All months; or specify list [6, 7, 8]
    output_dir='./climatology/'
)
```

## Statistical Analysis

The `StatisticalAnalyzer` class provides anomaly detection, trend analysis, and zonal statistics:

```python
from pycropwat.analysis import StatisticalAnalyzer

stats = StatisticalAnalyzer('./outputs')

# Calculate anomaly relative to climatology
anomaly = stats.calculate_anomaly(
    year=2023,
    month=6,
    clim_start=1990,
    clim_end=2020,
    anomaly_type='percent',  # 'absolute', 'percent', or 'standardized'
    output_path='./anomaly_2023_06.tif'
)

# Trend analysis with linear regression
trend_linear = stats.calculate_trend(
    start_year=2000,
    end_year=2020,
    month=None,  # Annual trends; specify 1-12 for specific month
    method='linear',
    output_dir='./trend_linear/'
)
# Returns a (slope, pvalue) tuple of DataArrays; output_dir also writes
# trend_slope_*.tif and trend_pvalue_*.tif

# Trend analysis with Theil-Sen slope and Mann-Kendall test
trend_sen = stats.calculate_trend(
    start_year=2000,
    end_year=2020,
    method='sen',
    output_dir='./trend_sen/'
)
# Returns the same (slope, pvalue) tuple as 'linear'

# Zonal statistics by polygon
zonal_df = stats.zonal_statistics(
    geometry_path='./regions.shp',
    start_year=2000,
    end_year=2020,
    months=None,  # All months; or specify list [6, 7, 8]
    stats=['mean', 'sum', 'min', 'max', 'std'],
    output_path='./zonal_stats.csv'
)
```

## Export Functions

Export data to different formats:

```python
from pycropwat.analysis import export_to_netcdf, export_to_cog

# Export time series to NetCDF with time dimension
export_to_netcdf(
    input_dir='./outputs',
    output_path='./data.nc',
    variable_name='effective_precip',
    pattern='effective_precip_[0-9]*.tif',  # Excludes fraction files
    compression=True
)

# Convert to Cloud-Optimized GeoTIFF
export_to_cog(
    input_path='./effective_precip_2020_06.tif',
    output_path='./cog_2020_06.tif'
)

# Batch convert to COGs (excludes fraction files)
from pathlib import Path
for tif in Path('./outputs').glob('effective_precip_[0-9]*.tif'):
    export_to_cog(str(tif), f'./cogs/{tif.name}')
```

## Visualization

The `Visualizer` class provides plotting functions:

```python
from pycropwat.analysis import Visualizer

viz = Visualizer('./outputs')

# Time series plot
fig = viz.plot_time_series(
    start_year=2000,
    end_year=2020,
    stat='mean',  # Spatial statistic: 'mean', 'sum', 'min', 'max'
    title='Effective Precipitation Time Series',
    figsize=(12, 6),
    output_path='./timeseries.png'
)

# Monthly climatology bar chart
fig = viz.plot_monthly_climatology(
    start_year=2000,
    end_year=2020,
    stat='mean',
    title='Monthly Climatology',
    output_path='./climatology.png'
)

# Single raster map
fig = viz.plot_raster(
    year=2020,
    month=6,
    cmap='YlGnBu',
    vmin=0,
    vmax=150,
    title='Effective Precipitation - June 2020',
    output_path='./map_2020_06.png'
)

# Interactive HTML map (requires leafmap or folium)
viz.plot_interactive_map(
    year=2020,
    month=6,
    cmap='YlGnBu',
    opacity=0.7,
    zoom_start=6,
    output_path='./map.html'
)

# Compare two datasets side-by-side with difference map
fig = viz.plot_comparison(
    other_dir='./terraclimate_outputs',
    year=2020,
    month=6,
    labels=('ERA5', 'TerraClimate'),
    cmap='YlGnBu',
    diff_cmap='RdBu',
    output_path='./comparison.png'
)

# Scatter plot comparison with statistics (R², RMSE, bias)
fig = viz.plot_scatter_comparison(
    other_dir='./terraclimate_outputs',
    start_year=2000,
    end_year=2020,
    months=None,  # All months
    sample_size=10000,  # Max points to plot
    labels=('ERA5 (mm)', 'TerraClimate (mm)'),
    output_path='./scatter.png'
)

# Annual totals comparison bar chart
fig = viz.plot_annual_comparison(
    other_dir='./terraclimate_outputs',
    start_year=2000,
    end_year=2020,
    labels=('ERA5', 'TerraClimate'),
    output_path='./annual_compare.png'
)
```

## Working with Results

```python
import rioxarray

# Process and get file paths
results = ep.process(output_dir='./outputs', n_workers=4)

# Results is a list of tuples: [(ep_path, epf_path), ...]
for ep_path, epf_path in results:
    if ep_path is not None:
        # Read effective precipitation
        da = rioxarray.open_rasterio(ep_path)
        print(f"Mean effective precip: {da.mean().values:.2f} mm")
```

## Advanced: Accessing Internal Methods

```python
# Get native scale of the dataset
native_scale = ep._get_native_scale()
print(f"Native scale: {native_scale} meters")

# Get a monthly image directly
monthly_img = ep._get_monthly_image(year=2020, month=6)

# Download monthly precipitation as xarray
pr_data = ep._download_monthly_precip(year=2020, month=6)
print(pr_data)
```

## Error Handling

```python
from pycropwat import EffectivePrecipitation

try:
    ep = EffectivePrecipitation(
        asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
        precip_band='total_precipitation_sum',
        geometry_path='nonexistent.geojson',
        start_year=2020,
        end_year=2023
    )
except FileNotFoundError:
    print("Geometry file not found!")
except ValueError as e:
    print(f"Invalid parameters: {e}")
```

## Logging

pyCropWat uses Python's logging module:

```python
import logging

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)

# Or just for pycropwat
logger = logging.getLogger('pycropwat')
logger.setLevel(logging.DEBUG)
```

## Integration with Dask

pyCropWat uses Dask for parallel processing. You can customize the scheduler:

```python
from dask.distributed import Client

# Create a Dask client with custom settings
client = Client(n_workers=8, threads_per_worker=2)

# Process (will use the active client)
results = ep.process(output_dir='./outputs', n_workers=8)

# Close client when done
client.close()
```

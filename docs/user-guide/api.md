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

## Effective Precipitation Methods

pyCropWat supports multiple methods for calculating effective precipitation:

```python
from pycropwat.methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip
)
import numpy as np

precip = np.array([50, 100, 200, 300, 400])

# CROPWAT (USDA SCS method - default)
peff_cropwat = cropwat_effective_precip(precip)
# [46.  72.  136.  155.  165.]

# FAO/AGLW method
peff_fao = fao_aglw_effective_precip(precip)

# Fixed percentage (e.g., 80%)
peff_fixed = fixed_percentage_effective_precip(precip, percentage=0.8)
# [40.  80.  160.  240.  320.]

# Dependable rainfall at 75% probability
peff_depend = dependable_rainfall_effective_precip(precip, probability=0.75)
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
growing = agg.growing_season_aggregate(
    year=2020,
    start_month=4,
    end_month=10,
    method='sum',
    output_path='./growing_2020.tif'
)

# Custom month range
custom = agg.custom_aggregate(
    year=2020,
    months=[5, 6, 7, 8],
    method='mean',
    output_path='./may_aug_mean_2020.tif'
)

# Multi-year climatology (long-term average)
climatology = agg.climatology(
    start_year=1990,
    end_year=2020,
    method='mean',
    output_dir='./climatology/'
)
```

## Statistical Analysis

The `StatisticalAnalyzer` class provides anomaly detection, trend analysis, and zonal statistics:

```python
from pycropwat.analysis import StatisticalAnalyzer, TemporalAggregator

agg = TemporalAggregator('./outputs')
stats = StatisticalAnalyzer(agg)

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
# Returns dict with 'slope', 'intercept', 'r_squared', 'p_value' rasters

# Trend analysis with Theil-Sen slope and Mann-Kendall test
trend_sen = stats.calculate_trend(
    start_year=2000,
    end_year=2020,
    method='sen',
    output_dir='./trend_sen/'
)
# Returns dict with 'slope', 'p_value' rasters

# Zonal statistics by polygon
zonal_df = stats.zonal_statistics(
    zones_path='./regions.shp',
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
    variable='effective_precip',
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
    basemap='OpenStreetMap',
    output_path='./map.html'
)

# Compare two datasets side-by-side with difference map
fig = viz.plot_comparison(
    other_input_dir='./terraclimate_outputs',
    year=2020,
    month=6,
    label1='ERA5',
    label2='TerraClimate',
    cmap='YlGnBu',
    diff_cmap='RdBu',
    output_path='./comparison.png'
)

# Scatter plot comparison with statistics (R², RMSE, bias)
fig = viz.plot_scatter_comparison(
    other_input_dir='./terraclimate_outputs',
    start_year=2000,
    end_year=2020,
    months=None,  # All months
    sample_size=10000,  # Max points to plot
    label1='ERA5 (mm)',
    label2='TerraClimate (mm)',
    output_path='./scatter.png'
)

# Annual totals comparison bar chart
fig = viz.plot_annual_comparison(
    other_input_dir='./terraclimate_outputs',
    start_year=2000,
    end_year=2020,
    label1='ERA5',
    label2='TerraClimate',
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

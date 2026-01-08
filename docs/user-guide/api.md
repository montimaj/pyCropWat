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

## Static Method: CROPWAT Formula

The CROPWAT effective precipitation formula is available as a static method:

```python
import numpy as np
from pycropwat import EffectivePrecipitation

# Calculate effective precipitation for an array
precip = np.array([50, 100, 200, 300, 400])
eff_precip = EffectivePrecipitation.cropwat_effective_precip(precip)

print(eff_precip)
# [46.  72.  136.  155.  165.]
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

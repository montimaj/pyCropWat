# Examples

## Example 1: ERA5-Land for South America

Calculate effective precipitation for a large river basin using ERA5-Land monthly data.

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    gee_geometry_asset='projects/my-project/assets/amazon_basin',
    start_year=2010,
    end_year=2023,
    precip_scale_factor=1000  # ERA5 precipitation is in meters
)

results = ep.process(
    output_dir='./amazon_effective_precip',
    n_workers=8
)
```

## Example 2: TerraClimate Long-term Analysis

Process 40+ years of TerraClimate data for drought analysis.

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    asset_id='IDAHO_EPSCOR/TERRACLIMATE',
    precip_band='pr',
    geometry_path='western_us.geojson',
    start_year=1980,
    end_year=2023
    # No scale_factor needed - TerraClimate is already in mm
)

# Process only growing season months
results = ep.process(
    output_dir='./terraclimate_growing_season',
    n_workers=4,
    months=[4, 5, 6, 7, 8, 9, 10]  # April-October
)
```

## Example 3: High-Resolution Custom Output

Generate 1km resolution outputs from ERA5-Land.

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='small_watershed.geojson',
    start_year=2020,
    end_year=2023,
    scale=1000,  # 1 km resolution
    precip_scale_factor=1000
)

results = ep.process(output_dir='./high_res_outputs', n_workers=4)
```

## Example 4: CHIRPS for Africa

Use CHIRPS daily data for African agriculture analysis.

```bash
pycropwat \
    --asset UCSB-CHG/CHIRPS/DAILY \
    --band precipitation \
    --gee-geometry projects/my-project/assets/africa_croplands \
    --start-year 2015 \
    --end-year 2023 \
    --months 3 4 5 6 7 8 9 10 \
    --workers 8 \
    --output ./africa_chirps_outputs
```

!!! note
    Daily CHIRPS data is automatically summed to monthly totals.

## Example 5: Analyzing Results with xarray

```python
import xarray as xr
import rioxarray
import numpy as np
from pathlib import Path

# Load all effective precipitation files
output_dir = Path('./outputs')
files = sorted(output_dir.glob('effective_precip_*.tif'))

# Stack into a time series
datasets = []
for f in files:
    da = rioxarray.open_rasterio(f).squeeze('band', drop=True)
    # Extract year and month from filename
    year = int(f.stem.split('_')[2])
    month = int(f.stem.split('_')[3])
    da = da.expand_dims({'time': [f'{year}-{month:02d}-01']})
    datasets.append(da)

# Combine
ds = xr.concat(datasets, dim='time')
ds['time'] = pd.to_datetime(ds['time'])

# Calculate annual totals
annual = ds.resample(time='Y').sum()
print(f"Mean annual effective precipitation: {annual.mean().values:.1f} mm")

# Calculate long-term mean
ltm = ds.groupby('time.month').mean()
print("Long-term monthly means:")
print(ltm.mean(dim=['x', 'y']).values)
```

## Example 6: Comparing Datasets

```python
from pycropwat import EffectivePrecipitation

# ERA5-Land
ep_era5 = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2015,
    end_year=2020,
    precip_scale_factor=1000
)
ep_era5.process(output_dir='./era5_outputs', n_workers=4)

# TerraClimate
ep_tc = EffectivePrecipitation(
    asset_id='IDAHO_EPSCOR/TERRACLIMATE',
    precip_band='pr',
    geometry_path='study_area.geojson',
    start_year=2015,
    end_year=2020
)
ep_tc.process(output_dir='./terraclimate_outputs', n_workers=4)

# Now compare the outputs...
```

## Example 7: Batch Processing Multiple Regions

```python
from pycropwat import EffectivePrecipitation
from pathlib import Path

regions = [
    'projects/my-project/assets/region_north',
    'projects/my-project/assets/region_south',
    'projects/my-project/assets/region_east',
    'projects/my-project/assets/region_west',
]

for region in regions:
    region_name = region.split('/')[-1]
    print(f"Processing {region_name}...")
    
    ep = EffectivePrecipitation(
        asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
        precip_band='total_precipitation_sum',
        gee_geometry_asset=region,
        start_year=2020,
        end_year=2023,
        precip_scale_factor=1000
    )
    
    ep.process(
        output_dir=f'./outputs/{region_name}',
        n_workers=4
    )
```

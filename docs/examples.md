# Examples

## Example 1: ERA5-Land for South America

Calculate effective precipitation for a large river basin using ERA5-Land monthly data.

**Python:**

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

**CLI:**

```bash
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --gee-geometry projects/my-project/assets/amazon_basin \
    --start-year 2010 \
    --end-year 2023 \
    --scale-factor 1000 \
    --workers 8 \
    --output ./amazon_effective_precip
```

## Example 2: TerraClimate Long-term Analysis

Process 40+ years of TerraClimate data for drought analysis.

**Python:**

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

**CLI:**

```bash
pycropwat process \
    --asset IDAHO_EPSCOR/TERRACLIMATE \
    --band pr \
    --geometry western_us.geojson \
    --start-year 1980 \
    --end-year 2023 \
    --months 4 5 6 7 8 9 10 \
    --workers 4 \
    --output ./terraclimate_growing_season
```

## Example 3: High-Resolution Custom Output

Generate 1km resolution outputs from ERA5-Land.

**Python:**

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

**CLI:**

```bash
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry small_watershed.geojson \
    --start-year 2020 \
    --end-year 2023 \
    --scale 1000 \
    --scale-factor 1000 \
    --output ./high_res_outputs
```

## Example 4: CHIRPS for Africa

Use CHIRPS daily data for African agriculture analysis.

**Python:**

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    asset_id='UCSB-CHG/CHIRPS/DAILY',
    precip_band='precipitation',
    gee_geometry_asset='projects/my-project/assets/africa_croplands',
    start_year=2015,
    end_year=2023
)

results = ep.process(
    output_dir='./africa_chirps_outputs',
    n_workers=8,
    months=[3, 4, 5, 6, 7, 8, 9, 10]  # Growing season
)
```

**CLI:**

```bash
pycropwat process \
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

## Example 5: Using Different Peff Methods

Compare different effective precipitation methods.

**Python:**

```python
from pycropwat import EffectivePrecipitation
from pycropwat.methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip
)
import numpy as np

# Direct calculation with numpy arrays
precip = np.array([50, 100, 150, 200, 250, 300])

print("Monthly Precip:", precip)
print("CROPWAT:", cropwat_effective_precip(precip))
print("FAO/AGLW:", fao_aglw_effective_precip(precip))
print("Fixed 80%:", fixed_percentage_effective_precip(precip, 0.8))
print("Depend 75%:", dependable_rainfall_effective_precip(precip, 0.75))
```

**CLI:**

```bash
# CROPWAT method (default)
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2015 --end-year 2020 \
    --method cropwat --output ./output_cropwat

# FAO/AGLW method
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2015 --end-year 2020 \
    --method fao_aglw --output ./output_fao

# Fixed percentage (80%)
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2015 --end-year 2020 \
    --method fixed_percentage --percentage 0.8 --output ./output_fixed

# Dependable rainfall (80% probability)
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2015 --end-year 2020 \
    --method dependable_rainfall --probability 0.8 --output ./output_depend
```

## Example 6: Temporal Aggregation

Aggregate monthly data into seasonal, annual, and climatology products.

**Python:**

```python
from pycropwat.analysis import TemporalAggregator

agg = TemporalAggregator('./outputs')

# Annual totals for each year
for year in range(2010, 2024):
    agg.annual_aggregate(
        year=year,
        method='sum',
        output_path=f'./annual/annual_{year}.tif'
    )

# Seasonal aggregations
for year in range(2010, 2024):
    for season in ['DJF', 'MAM', 'JJA', 'SON']:
        agg.seasonal_aggregate(
            year=year,
            season=season,
            method='sum',
            output_path=f'./seasonal/{season}_{year}.tif'
        )

# Growing season (April-October)
for year in range(2010, 2024):
    agg.growing_season_aggregate(
        year=year,
        start_month=4,
        end_month=10,
        method='sum',
        output_path=f'./growing/growing_{year}.tif'
    )

# 30-year climatology (long-term monthly means)
climatology = agg.climatology(
    start_year=1990,
    end_year=2020,
    method='mean',
    output_dir='./climatology/'
)
```

**CLI:**

```bash
# Annual total
pycropwat aggregate --input ./outputs --type annual \
    --year 2020 --output ./annual_2020.tif

# Summer (JJA) aggregate
pycropwat aggregate --input ./outputs --type seasonal \
    --year 2020 --season JJA --output ./summer_2020.tif

# Growing season (April-October)
pycropwat aggregate --input ./outputs --type growing-season \
    --year 2020 --start-month 4 --end-month 10 \
    --output ./growing_2020.tif

# 30-year climatology
pycropwat aggregate --input ./outputs --type climatology \
    --start-year 1990 --end-year 2020 --output ./climatology/
```

## Example 7: Statistical Analysis

Perform anomaly detection, trend analysis, and zonal statistics.

**Python:**

```python
from pycropwat.analysis import TemporalAggregator, StatisticalAnalyzer

agg = TemporalAggregator('./outputs')
stats = StatisticalAnalyzer(agg)

# Percent anomaly for June 2023 relative to 1990-2020 climatology
anomaly = stats.calculate_anomaly(
    year=2023,
    month=6,
    clim_start=1990,
    clim_end=2020,
    anomaly_type='percent',
    output_path='./anomaly_2023_06.tif'
)

# Standardized anomaly (z-score)
anomaly_std = stats.calculate_anomaly(
    year=2023,
    month=6,
    clim_start=1990,
    clim_end=2020,
    anomaly_type='standardized',
    output_path='./anomaly_std_2023_06.tif'
)

# Linear trend analysis (2000-2020)
trend = stats.calculate_trend(
    start_year=2000,
    end_year=2020,
    method='linear',
    output_dir='./trend_linear/'
)
print(f"Trend outputs: slope, intercept, r_squared, p_value")

# Theil-Sen trend with Mann-Kendall test
trend_sen = stats.calculate_trend(
    start_year=2000,
    end_year=2020,
    method='sen',
    output_dir='./trend_sen/'
)

# Zonal statistics by administrative region
zonal_df = stats.zonal_statistics(
    zones_path='./admin_regions.shp',
    start_year=2010,
    end_year=2020,
    stats=['mean', 'sum', 'min', 'max', 'std'],
    output_path='./zonal_stats.csv'
)
print(zonal_df.head())
```

**CLI:**

```bash
# Percent anomaly
pycropwat analyze anomaly --input ./outputs \
    --year 2023 --month 6 \
    --clim-start 1990 --clim-end 2020 \
    --anomaly-type percent \
    --output ./anomaly_2023_06.tif

# Theil-Sen trend analysis
pycropwat analyze trend --input ./outputs \
    --start-year 2000 --end-year 2020 \
    --trend-method sen \
    --output ./trend/

# Zonal statistics
pycropwat analyze zonal --input ./outputs \
    --zones ./admin_regions.shp \
    --start-year 2010 --end-year 2020 \
    --output ./zonal_stats.csv
```

## Example 8: Export to Different Formats

Export data to NetCDF and Cloud-Optimized GeoTIFF.

**Python:**

```python
from pycropwat.analysis import export_to_netcdf, export_to_cog
from pathlib import Path

# Export all monthly data to a single NetCDF file
export_to_netcdf(
    input_dir='./outputs',
    output_path='./effective_precip.nc',
    variable='effective_precip',
    pattern='effective_precip_[0-9]*.tif',  # Excludes fraction files
    compression=True
)

# Batch convert to Cloud-Optimized GeoTIFFs (excludes fraction files)
output_cog_dir = Path('./cogs')
output_cog_dir.mkdir(exist_ok=True)

for tif in Path('./outputs').glob('effective_precip_[0-9]*.tif'):
    export_to_cog(
        input_path=str(tif),
        output_path=str(output_cog_dir / tif.name)
    )
```

**CLI:**

```bash
# Export to NetCDF
pycropwat export netcdf --input ./outputs --output ./effective_precip.nc

# Convert single file to COG
pycropwat export cog \
    --input ./outputs/effective_precip_2020_06.tif \
    --output ./cogs/effective_precip_2020_06.tif

# Batch convert with shell loop (excludes fraction files)
mkdir -p ./cogs
for f in ./outputs/effective_precip_[0-9]*.tif; do
    pycropwat export cog --input "$f" --output "./cogs/$(basename $f)"
done
```

## Example 9: Visualization

Create static plots and interactive maps.

**Python:**

```python
from pycropwat.analysis import Visualizer

viz = Visualizer('./outputs')

# Time series plot
fig = viz.plot_time_series(
    start_year=2000,
    end_year=2023,
    stat='mean',
    title='Mean Effective Precipitation (2000-2023)',
    figsize=(14, 6),
    output_path='./timeseries.png'
)

# Monthly climatology bar chart
fig = viz.plot_monthly_climatology(
    start_year=2000,
    end_year=2023,
    stat='mean',
    title='Monthly Climatology (2000-2023)',
    output_path='./climatology.png'
)

# Single month raster map
fig = viz.plot_raster(
    year=2020,
    month=6,
    cmap='YlGnBu',
    vmin=0,
    vmax=150,
    title='Effective Precipitation - June 2020',
    output_path='./map_2020_06.png'
)

# Interactive map (requires: pip install leafmap)
viz.plot_interactive_map(
    year=2020,
    month=6,
    cmap='YlGnBu',
    opacity=0.7,
    basemap='OpenStreetMap',
    output_path='./interactive_map.html'
)

# Anomaly map from GeoTIFF (diverging colormap, centered at 100%)
fig = viz.plot_anomaly_map(
    anomaly_path='./analysis_outputs/anomalies/ERA5Land/anomaly_2023_08.tif',
    title='Peff Anomaly - August 2023 (% of Climatology)',
    output_path='./anomaly_map_2023_08.png'
)

# Climatology map from GeoTIFF
fig = viz.plot_climatology_map(
    climatology_path='./analysis_outputs/climatology/ERA5Land/climatology_06.tif',
    title='Peff Climatology - June (2000-2020)',
    output_path='./climatology_map_06.png'
)

# Trend map with significance stippling
fig = viz.plot_trend_map(
    slope_path='./analysis_outputs/trend/ERA5Land/slope.tif',
    pvalue_path='./analysis_outputs/trend/ERA5Land/pvalue.tif',
    title='Peff Trend (2000-2023)',
    show_significance=True,  # Adds stippling for p < 0.05
    output_path='./trend_map.png'
)

# Trend panel (slope and p-value side-by-side)
fig = viz.plot_trend_panel(
    slope_path='./analysis_outputs/trend/ERA5Land/slope.tif',
    pvalue_path='./analysis_outputs/trend/ERA5Land/pvalue.tif',
    title='Peff Trend Analysis (2000-2023)',
    output_path='./trend_panel.png'
)
```

**CLI:**

```bash
# Time series plot
pycropwat plot timeseries --input ./outputs \
    --start-year 2000 --end-year 2023 \
    --title "Mean Effective Precipitation (2000-2023)" \
    --output ./timeseries.png

# Monthly climatology
pycropwat plot climatology --input ./outputs \
    --start-year 2000 --end-year 2023 \
    --output ./climatology.png

# Raster map
pycropwat plot map --input ./outputs \
    --year 2020 --month 6 \
    --cmap YlGnBu --vmin 0 --vmax 150 \
    --output ./map_2020_06.png

# Interactive map
pycropwat plot interactive --input ./outputs \
    --year 2020 --month 6 \
    --output ./interactive_map.html
```

## Example 10: Comparing Datasets

Compare ERA5-Land vs TerraClimate effective precipitation.

**Python:**

```python
from pycropwat import EffectivePrecipitation
from pycropwat.analysis import Visualizer

# Process ERA5-Land
ep_era5 = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2000,
    end_year=2020,
    precip_scale_factor=1000
)
ep_era5.process(output_dir='./era5_outputs', n_workers=4)

# Process TerraClimate
ep_tc = EffectivePrecipitation(
    asset_id='IDAHO_EPSCOR/TERRACLIMATE',
    precip_band='pr',
    geometry_path='study_area.geojson',
    start_year=2000,
    end_year=2020
)
ep_tc.process(output_dir='./terraclimate_outputs', n_workers=4)

# Compare the outputs
viz = Visualizer('./era5_outputs')

# Side-by-side comparison with difference map
fig = viz.plot_comparison(
    other_input_dir='./terraclimate_outputs',
    year=2020,
    month=6,
    label1='ERA5-Land',
    label2='TerraClimate',
    cmap='YlGnBu',
    diff_cmap='RdBu',
    output_path='./comparison_2020_06.png'
)

# Scatter plot with statistics (R², RMSE, bias)
fig = viz.plot_scatter_comparison(
    other_input_dir='./terraclimate_outputs',
    start_year=2000,
    end_year=2020,
    sample_size=10000,
    label1='ERA5-Land (mm)',
    label2='TerraClimate (mm)',
    output_path='./scatter_comparison.png'
)

# Annual totals comparison
fig = viz.plot_annual_comparison(
    other_input_dir='./terraclimate_outputs',
    start_year=2000,
    end_year=2020,
    label1='ERA5-Land',
    label2='TerraClimate',
    output_path='./annual_comparison.png'
)
```

**CLI:**

```bash
# Process both datasets
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry study_area.geojson \
    --start-year 2000 --end-year 2020 \
    --scale-factor 1000 --output ./era5_outputs

pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE \
    --band pr \
    --geometry study_area.geojson \
    --start-year 2000 --end-year 2020 \
    --output ./terraclimate_outputs

# Side-by-side comparison
pycropwat plot compare \
    --input ./era5_outputs \
    --other-input ./terraclimate_outputs \
    --year 2020 --month 6 \
    --label1 ERA5-Land --label2 TerraClimate \
    --output ./comparison_2020_06.png

# Scatter plot
pycropwat plot scatter \
    --input ./era5_outputs \
    --other-input ./terraclimate_outputs \
    --start-year 2000 --end-year 2020 \
    --output ./scatter_comparison.png

# Annual comparison
pycropwat plot annual-compare \
    --input ./era5_outputs \
    --other-input ./terraclimate_outputs \
    --start-year 2000 --end-year 2020 \
    --label1 ERA5-Land --label2 TerraClimate \
    --output ./annual_comparison.png
```

## Example 11: Batch Processing Multiple Regions

**Python:**

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

**CLI:**

```bash
# Using a shell loop
for region in region_north region_south region_east region_west; do
    echo "Processing $region..."
    pycropwat process \
        --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
        --band total_precipitation_sum \
        --gee-geometry "projects/my-project/assets/$region" \
        --start-year 2020 \
        --end-year 2023 \
        --scale-factor 1000 \
        --output "./outputs/$region"
done
```

## Example 12: Complete Workflow

A complete analysis workflow from data processing to visualization.

!!! tip "Ready-to-Run Scripts"
    Comprehensive workflow scripts are available in the `Examples/` directory:
    
    **South America (Rio de la Plata):** [`Examples/south_america_cropwat_example.py`](https://github.com/montimaj/pyCropWat/blob/main/Examples/south_america_cropwat_example.py)
    
    ```bash
    # Run with existing data (analysis only)
    python Examples/south_america_cropwat_example.py --analysis-only
    
    # Run full workflow with GEE processing
    python Examples/south_america_cropwat_example.py --gee-project your-project-id --workers 8
    ```
    
    **Arizona (USDA-SCS Method):** [`Examples/arizona_usda_scs_example.py`](https://github.com/montimaj/pyCropWat/blob/main/Examples/arizona_usda_scs_example.py)
    
    ```bash
    # Run with existing data (analysis only)
    python Examples/arizona_usda_scs_example.py --analysis-only
    
    # Run full workflow with GEE processing
    python Examples/arizona_usda_scs_example.py --gee-project your-project-id --workers 8
    ```

```python
from pycropwat import EffectivePrecipitation
from pycropwat.analysis import (
    TemporalAggregator,
    StatisticalAnalyzer,
    Visualizer,
    export_to_netcdf
)
from pathlib import Path

# 1. Process effective precipitation
print("Step 1: Processing effective precipitation...")
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=1990,
    end_year=2023,
    precip_scale_factor=1000
)
ep.process(output_dir='./outputs', n_workers=8)

# 2. Temporal aggregation
print("Step 2: Creating temporal aggregations...")
agg = TemporalAggregator('./outputs')

# Annual totals
Path('./annual').mkdir(exist_ok=True)
for year in range(1990, 2024):
    agg.annual_aggregate(year, output_path=f'./annual/annual_{year}.tif')

# Climatology
agg.climatology(1990, 2020, output_dir='./climatology/')

# 3. Statistical analysis
print("Step 3: Running statistical analysis...")
stats = StatisticalAnalyzer(agg)

# Anomaly for recent years
Path('./anomalies').mkdir(exist_ok=True)
for year in range(2021, 2024):
    for month in range(1, 13):
        try:
            stats.calculate_anomaly(
                year=year, month=month,
                clim_start=1990, clim_end=2020,
                anomaly_type='percent',
                output_path=f'./anomalies/anomaly_{year}_{month:02d}.tif'
            )
        except:
            pass  # Skip missing months

# Trend analysis
stats.calculate_trend(1990, 2023, method='sen', output_dir='./trend/')

# Zonal statistics
zonal_df = stats.zonal_statistics(
    zones_path='./admin_regions.shp',
    start_year=1990,
    end_year=2023,
    output_path='./zonal_stats.csv'
)

# 4. Visualization
print("Step 4: Creating visualizations...")
viz = Visualizer('./outputs')

viz.plot_time_series(1990, 2023, output_path='./figures/timeseries.png')
viz.plot_monthly_climatology(1990, 2020, output_path='./figures/climatology.png')
viz.plot_raster(2023, 6, output_path='./figures/map_2023_06.png')
viz.plot_interactive_map(2023, 6, output_path='./figures/map_2023_06.html')

# 5. Export to NetCDF
print("Step 5: Exporting to NetCDF...")
export_to_netcdf('./outputs', './effective_precip_1990_2023.nc')

print("Workflow complete!")
```

### Example Outputs

The complete workflow example generates the following visualizations using real Rio de la Plata basin data:

#### Time Series & Climatology

<p align="center">
  <img src="../assets/examples/figures/ERA5Land/time_series.png" width="48%" alt="Time Series">
  <img src="../assets/examples/figures/ERA5Land/monthly_climatology.png" width="48%" alt="Monthly Climatology">
</p>

*Left: Monthly effective precipitation time series (2000-2025). Right: Monthly climatology showing seasonal patterns.*

#### Spatial Maps

<p align="center">
  <img src="../assets/examples/figures/ERA5Land/map_2023_06.png" width="32%" alt="Winter Map">
  <img src="../assets/examples/figures/ERA5Land/map_2023_01.png" width="32%" alt="Summer Map">
  <img src="../assets/examples/figures/ERA5Land/map_notable_2015_12.png" width="32%" alt="El Niño Event">
</p>

*Left: Winter dry season (June 2023). Center: Summer wet season (January 2023). Right: El Niño event (December 2015).*

#### Dataset Comparison (ERA5-Land vs TerraClimate)

<p align="center">
  <img src="../assets/examples/comparisons/comparison_2023_06.png" width="100%" alt="Side-by-side Comparison">
</p>

*Side-by-side comparison of ERA5-Land and TerraClimate effective precipitation with difference map.*

<p align="center">
  <img src="../assets/examples/comparisons/scatter_comparison.png" width="48%" alt="Scatter Plot">
  <img src="../assets/examples/comparisons/annual_comparison.png" width="48%" alt="Annual Comparison">
</p>

*Left: Scatter plot comparison with R², RMSE, and bias statistics. Right: Annual totals comparison.*

<p align="center">
  <img src="../assets/examples/comparisons/zonal_comparison.png" width="70%" alt="Zonal Comparison">
</p>

*Zonal statistics comparison between ERA5-Land and TerraClimate for Eastern and Western Rio de la Plata regions.*

#### Anomaly, Climatology & Trend Maps

<p align="center">
  <img src="../assets/examples/figures/ERA5Land/anomaly_map_2023_01.png" width="32%" alt="Anomaly Map">
  <img src="../assets/examples/figures/ERA5Land/climatology_map_01.png" width="32%" alt="Climatology Map">
  <img src="../assets/examples/figures/ERA5Land/trend_map_with_significance.png" width="32%" alt="Trend Map">
</p>

*Left: Percent anomaly for January 2023 (diverging colormap centered at 100%). Center: January climatology (2000-2020). Right: Long-term trend with significance stippling (p < 0.05).*

<p align="center">
  <img src="../assets/examples/figures/ERA5Land/trend_maps.png" width="70%" alt="Trend Panel">
</p>

*Trend analysis panel showing slope (mm/year) and p-value maps side by side.*

#### Method Comparison

<p align="center">
  <img src="../assets/examples/method_comparison/ERA5Land_method_maps_2020_01.png" width="100%" alt="Method Comparison Maps">
</p>

*Comparison of effective precipitation methods: CROPWAT, FAO/AGLW, Fixed Percentage (70%), Dependable Rainfall (75%), and FarmWest.*

<p align="center">
  <img src="../assets/examples/method_comparison/ERA5Land_method_curves.png" width="60%" alt="Method Curves">
</p>

*Theoretical response curves for different effective precipitation methods.*

---

## Arizona USDA-SCS Example

This example demonstrates the **USDA-SCS method** for calculating effective precipitation using U.S.-specific datasets for Arizona, with comparisons to global datasets and multiple Peff methods.

!!! tip "Ready-to-Run Script"
    A comprehensive Arizona workflow is available at [`Examples/arizona_usda_scs_example.py`](https://github.com/montimaj/pyCropWat/blob/main/Examples/arizona_usda_scs_example.py):
    
    ```bash
    # Run with existing data
    python Examples/arizona_usda_scs_example.py --analysis-only
    
    # Run full workflow with GEE processing
    python Examples/arizona_usda_scs_example.py --gee-project your-project-id
    ```

### Dataset Configuration

The Arizona example processes **4 precipitation datasets** (2 U.S., 2 Global):

| Dataset | Type | GEE Asset ID | Band | Resolution |
|---------|------|-------------|------|------------|
| **GridMET** | U.S. | `IDAHO_EPSCOR/GRIDMET` | `pr` | ~4 km |
| **PRISM** | U.S. | `projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY` | `ppt` | ~800 m |
| **ERA5-Land** | Global | `ECMWF/ERA5_LAND/MONTHLY_AGGR` | `total_precipitation_sum` | ~11 km |
| **TerraClimate** | Global | `IDAHO_EPSCOR/TERRACLIMATE` | `pr` | ~4 km |

For **USDA-SCS method**, dataset-specific AWC and ETo are used:

**U.S. Datasets (High-Resolution):**

| Dataset | GEE Asset ID | Description |
|---------|-------------|-------------|
| **AWC (SSURGO)** | `projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite` | USDA SSURGO soil water holding capacity |
| **ETo (GridMET)** | `projects/openet/assets/reference_et/conus/gridmet/monthly/v1` | OpenET GridMET reference ET (band: `eto`) |

**Global Datasets:**

| Dataset | GEE Asset ID | Description |
|---------|-------------|-------------|
| **AWC (FAO HWSD)** | `projects/sat-io/open-datasets/FAO/HWSD_V2_SMU` | FAO Harmonized World Soil Database v2 (band: `AWC`) |
| **ETo (AgERA5)** | `projects/climate-engine-pro/assets/ce-ag-era5-v2/daily` | AgERA5 FAO-56 Penman-Monteith ETo (band: `ReferenceET_PenmanMonteith_FAO56`) |

### Python Example

```python
from pycropwat import EffectivePrecipitation

# Arizona with GridMET precipitation and USDA-SCS method
ep = EffectivePrecipitation(
    asset_id='IDAHO_EPSCOR/GRIDMET',
    precip_band='pr',
    gee_geometry_asset='users/montimajumdar/AZ',
    start_year=2000,
    end_year=2024,
    scale=4000,  # 4km resolution
    # USDA-SCS specific parameters
    method='usda_scs',
    awc_asset='projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite',
    eto_asset='projects/openet/assets/reference_et/conus/gridmet/monthly/v1',
    eto_band='eto',
    rooting_depth=1.0  # 1 meter
)

ep.process(output_dir='./Arizona/AZ_GridMET_USDA_SCS', n_workers=8)
```

### CLI Example

```bash
# GridMET precipitation with USDA-SCS method for Arizona
pycropwat process --asset IDAHO_EPSCOR/GRIDMET --band pr \
    --gee-geometry users/montimajumdar/AZ \
    --start-year 2000 --end-year 2024 --scale 4000 \
    --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto --rooting-depth 1.0 \
    --workers 8 --output ./Arizona/AZ_GridMET_USDA_SCS

# PRISM precipitation (higher resolution ~800m)
pycropwat process --asset projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY --band ppt \
    --gee-geometry users/montimajumdar/AZ \
    --start-year 2000 --end-year 2024 --scale 4000 \
    --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto --rooting-depth 1.0 \
    --workers 8 --output ./Arizona/AZ_PRISM_USDA_SCS
```

### Key Features

- **USDA-SCS Method**: Accounts for soil water holding capacity (AWC) and evaporative demand (ETo)
- **U.S. High-Resolution Data**: GridMET (~4km) and PRISM (~800m) precipitation with SSURGO AWC and GridMET ETo
- **Global Coverage Data**: ERA5-Land and TerraClimate with FAO HWSD AWC and AgERA5 ETo
- **ETo Scale Factor**: Support for different ETo units via `eto_scale_factor` parameter
- **U.S. vs Global Comparison**: Compare U.S. datasets (GridMET, PRISM) against global datasets (ERA5-Land, TerraClimate)
- **Method Comparison**: Compare all 6 Peff methods (CROPWAT, FAO/AGLW, Fixed %, Dependable Rain, FarmWest, USDA-SCS)
- **Arizona-Specific Analysis**: Monsoon season (Jul-Sep), winter season aggregations
- **Regional Zones**: Central AZ (Phoenix), Southern AZ (Tucson), Northern AZ (Flagstaff)

### Generated Outputs

The workflow generates comprehensive analysis outputs:

```
Examples/Arizona/
├── AZ_GridMET_USDA_SCS/    # U.S. dataset outputs
├── AZ_PRISM_USDA_SCS/
├── AZ_ERA5Land_USDA_SCS/   # Global dataset outputs
├── AZ_TerraClimate_USDA_SCS/
└── analysis_outputs/
    ├── annual/             # Temporal aggregations
    ├── climatology/
    ├── anomalies/          # Statistical analysis & CSV exports
    ├── trend/
    ├── figures/            # Raster maps & time series plots
    ├── comparisons/
    ├── zonal_stats/
    ├── us_vs_global/       # U.S. vs Global dataset comparison
    │   ├── us_global_spatial_*.png
    │   ├── multi_dataset_summary.png
    │   └── dataset_statistics.csv
    ├── method_comparison/  # All 6 methods comparison
    │   ├── method_spatial_*.png
    │   ├── method_curves_arizona.png
    │   ├── method_comparison_summary.png
    │   └── method_statistics.csv
    ├── az_zones.geojson
    └── *_usda_scs_peff_2000_2024.nc
```

### Arizona Example Outputs

The Arizona workflow generates the following visualizations:

#### Time Series & Climatology

<p align="center">
  <img src="../assets/examples/arizona/figures/GridMET/time_series.png" width="48%" alt="Arizona Time Series">
  <img src="../assets/examples/arizona/figures/GridMET/monthly_climatology.png" width="48%" alt="Arizona Climatology">
</p>

*GridMET USDA-SCS effective precipitation for Arizona: time series (left) and monthly climatology (right).*

#### Anomaly, Climatology & Trend Maps

<p align="center">
  <img src="../assets/examples/arizona/figures/GridMET/anomaly_map_2023_08.png" width="32%" alt="Arizona Anomaly">
  <img src="../assets/examples/arizona/figures/GridMET/climatology_map_08.png" width="32%" alt="Arizona Climatology Map">
  <img src="../assets/examples/arizona/figures/GridMET/trend_map_with_significance.png" width="32%" alt="Arizona Trend">
</p>

*Left: August 2023 anomaly (monsoon season). Center: August climatology (monsoon peak). Right: Long-term trend with significance stippling (p < 0.05).*

<p align="center">
  <img src="../assets/examples/arizona/figures/GridMET/trend_maps.png" width="70%" alt="Arizona Trend Panel">
</p>

*Trend analysis panel for Arizona showing slope (mm/year) and p-value maps side by side.*

#### U.S. vs Global Dataset Comparison

<p align="center">
  <img src="../assets/examples/arizona/us_vs_global/multi_dataset_summary.png" width="70%" alt="U.S. vs Global">
</p>

*Comparison of U.S. datasets (GridMET, PRISM) vs Global datasets (ERA5-Land, TerraClimate) for Arizona.*

#### Method Comparison

<p align="center">
  <img src="../assets/examples/arizona/method_comparison/method_curves_arizona.png" width="60%" alt="Arizona Method Curves">
</p>

*Theoretical Peff response curves for Arizona conditions (AWC=150 mm/m, ETo=180 mm/month).*


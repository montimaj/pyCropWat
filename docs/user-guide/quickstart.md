# Quick Start

This guide will help you get started with pyCropWat in minutes.

pyCropWat works with any GEE ImageCollection containing precipitation data from the [GEE Data Catalog](https://developers.google.com/earth-engine/datasets) or [Community Catalog](https://gee-community-catalog.org/).

## Basic Usage

### Using the Command Line

The simplest way to use pyCropWat is through the CLI:

```bash
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry study_area.geojson \
    --start-year 2020 \
    --end-year 2023 \
    --scale-factor 1000 \
    --output ./outputs
```

### Using the Python API

```python
from pycropwat import EffectivePrecipitation

# Initialize with a local geometry file
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2020,
    end_year=2023,
    precip_scale_factor=1000  # Convert ERA5 meters to mm
)

# Process all months
results = ep.process(output_dir='./outputs', n_workers=4)
```

## Common Precipitation Datasets

pyCropWat works with any monthly precipitation dataset in Google Earth Engine:

| Dataset | Asset ID | Band | Scale Factor | Coverage |
|---------|----------|------|--------------|----------|
| ERA5-Land | `ECMWF/ERA5_LAND/MONTHLY_AGGR` | `total_precipitation_sum` | 1000 | Global (~11 km) |
| TerraClimate | `IDAHO_EPSCOR/TERRACLIMATE` | `pr` | 1 | Global (~4 km) |
| GridMET | `IDAHO_EPSCOR/GRIDMET` | `pr` | 1 | U.S. (~4 km) |
| PRISM | `OREGONSTATE/PRISM/AN81m` | `ppt` | 1 | U.S. (~4 km) |
| CHIRPS | `UCSB-CHG/CHIRPS/PENTAD` | `precipitation` | 1 | 50°S-50°N (~5.5 km) |
| GPM IMERG | `NASA/GPM_L3/IMERG_MONTHLY_V06` | `precipitation` | 1 | Global (~10 km) |

## Using GEE Vector Assets

You can use a GEE FeatureCollection asset instead of a local file:

**CLI:**

```bash
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --gee-geometry projects/my-project/assets/study_area \
    --start-year 2020 \
    --end-year 2023 \
    --scale-factor 1000 \
    --output ./outputs
```

**Python:**

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

## Effective Precipitation Methods

pyCropWat supports eight methods for calculating effective precipitation:

| Method | Formula | Use Case |
|--------|---------|----------|
| `cropwat` | P ≤ 250: Peff = P × (125 - 0.2P) / 125 <br> P > 250: Peff = 0.1P + 125 | Default FAO CROPWAT method |
| `fao_aglw` | P ≤ 250: Peff = 0.6P - 10 <br> P > 250: Peff = 0.8P - 25 | FAO Land and Water Division |
| `fixed_percentage` | Peff = P × percentage | Simple, adjustable (default 70%) |
| `dependable_rainfall` | Probability-scaled FAO method | Conservative estimates |
| `farmwest` | Peff = (P - 5) × 0.75 | Pacific Northwest irrigation |
| `usda_scs` | Uses AWC, ETo, and soil storage | Site-specific with soil data |
| `suet` | P ≤ ETo: 0 <br> P - ETo < 75: P - ETo <br> else: 75 + f(P - ETo - 75) | TAGEM-SuET method (P - ETo with 75mm cap) |
| `ensemble` | Mean of 6 methods (excludes TAGEM-SuET) | Robust multi-method estimate |

**CLI:**

```bash
# CROPWAT (default)
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2020 --end-year 2023 \
    --method cropwat --output ./output

# FAO/AGLW method
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2020 --end-year 2023 \
    --method fao_aglw --output ./output

# Fixed percentage (80%)
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2020 --end-year 2023 \
    --method fixed_percentage --percentage 0.8 --output ./output

# Dependable rainfall (80% probability)
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2020 --end-year 2023 \
    --method dependable_rainfall --probability 0.8 --output ./output

# FarmWest method
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson --start-year 2020 --end-year 2023 \
    --method farmwest --output ./output

# USDA-SCS method (requires AWC and ETo assets)
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR --band total_precipitation_sum \
    --geometry roi.geojson --start-year 2020 --end-year 2023 \
    --scale-factor 1000 --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --rooting-depth 1.0 --output ./output

# SuET method (requires ETo asset)
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR --band total_precipitation_sum \
    --geometry roi.geojson --start-year 2020 --end-year 2023 \
    --scale-factor 1000 --method suet \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --output ./output
```

**Python:**

```python
from pycropwat import EffectivePrecipitation
from pycropwat.methods import (
    cropwat_effective_precip, fao_aglw_effective_precip,
    fixed_percentage_effective_precip, dependable_rainfall_effective_precip,
    farmwest_effective_precip, usda_scs_effective_precip, suet_effective_precip
)

# Using different methods with numpy arrays
import numpy as np
precip = np.array([50, 100, 200, 300])

peff_cropwat = cropwat_effective_precip(precip)
peff_fao = fao_aglw_effective_precip(precip)
peff_fixed = fixed_percentage_effective_precip(precip, percentage=0.8)
peff_depend = dependable_rainfall_effective_precip(precip, probability=0.8)
peff_farmwest = farmwest_effective_precip(precip)

# USDA-SCS method requires ETo and AWC arrays
eto = np.array([80, 120, 180, 220])  # mm
awc = np.array([0.15, 0.15, 0.15, 0.15])  # volumetric fraction
peff_usda = usda_scs_effective_precip(precip, eto, awc, rooting_depth=1.0)

# SuET method requires ETo array
peff_suet = suet_effective_precip(precip, eto)
```

## Processing Specific Months

To process only certain months (e.g., growing season):

**CLI:**

```bash
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry study_area.geojson \
    --start-year 2020 \
    --end-year 2023 \
    --months 4 5 6 7 8 9 \
    --scale-factor 1000 \
    --output ./outputs
```

**Python:**

```python
results = ep.process(
    output_dir='./outputs',
    n_workers=4,
    months=[4, 5, 6, 7, 8, 9]  # April through September
)
```

## Temporal Aggregation

Aggregate monthly rasters into annual, seasonal, or custom periods:

**CLI:**

```bash
# Annual total
pycropwat aggregate --input ./outputs --type annual --year 2020 \
    --output ./annual_2020.tif

# Summer (JJA) aggregate
pycropwat aggregate --input ./outputs --type seasonal --year 2020 \
    --season JJA --output ./summer_2020.tif

# Growing season (April-October)
pycropwat aggregate --input ./outputs --type growing-season --year 2020 \
    --start-month 4 --end-month 10 --output ./growing_2020.tif

# Multi-year climatology (long-term mean)
pycropwat aggregate --input ./outputs --type climatology \
    --start-year 2000 --end-year 2020 --output ./climatology/
```

**Python:**

```python
from pycropwat.analysis import TemporalAggregator

agg = TemporalAggregator('./outputs')

# Annual totals
annual = agg.annual_aggregate(2020, method='sum', output_path='./annual_2020.tif')

# Seasonal aggregation
summer = agg.seasonal_aggregate(2020, 'JJA', method='sum')

# Growing season - Northern Hemisphere (April-October)
growing_nh = agg.growing_season_aggregate(2020, start_month=4, end_month=10)

# Growing season - Southern Hemisphere (October-March, cross-year)
# When start_month > end_month, automatically aggregates across two years
# This aggregates Oct 2020 - Mar 2021
growing_sh = agg.growing_season_aggregate(2020, start_month=10, end_month=3)

# Multi-year climatology
climatology = agg.climatology(start_year=2000, end_year=2020)
```

## Statistical Analysis

Perform anomaly detection, trend analysis, and zonal statistics:

**CLI:**

```bash
# Calculate anomaly relative to climatology
pycropwat analyze anomaly --input ./outputs --year 2023 --month 6 \
    --clim-start 1990 --clim-end 2020 \
    --anomaly-type percent --output ./anomaly_2023_06.tif

# Trend analysis with Theil-Sen slope
pycropwat analyze trend --input ./outputs \
    --start-year 2000 --end-year 2020 \
    --trend-method sen --output ./trend/

# Zonal statistics by region
pycropwat analyze zonal --input ./outputs --zones ./regions.shp \
    --start-year 2000 --end-year 2020 --output ./zonal_stats.csv
```

**Python:**

```python
from pycropwat.analysis import StatisticalAnalyzer, TemporalAggregator

agg = TemporalAggregator('./outputs')
stats = StatisticalAnalyzer(agg)

# Anomaly calculation
anomaly = stats.calculate_anomaly(
    year=2023, month=6,
    clim_start=1990, clim_end=2020,
    anomaly_type='percent'
)

# Trend analysis
trend = stats.calculate_trend(
    start_year=2000, end_year=2020,
    method='sen'
)

# Zonal statistics
zonal = stats.zonal_statistics(
    zones_path='./regions.shp',
    start_year=2000, end_year=2020
)
```

## Export Options

Export data to NetCDF or Cloud-Optimized GeoTIFF:

**CLI:**

```bash
# Export time series to NetCDF
pycropwat export netcdf --input ./outputs --output ./data.nc

# Convert to Cloud-Optimized GeoTIFF
pycropwat export cog --input ./effective_precip_2020_06.tif \
    --output ./cog_2020_06.tif
```

**Python:**

```python
from pycropwat.analysis import export_to_netcdf, export_to_cog

# Export to NetCDF with time dimension
export_to_netcdf('./outputs', './data.nc', variable='effective_precip')

# Convert to Cloud-Optimized GeoTIFF
export_to_cog('./effective_precip_2020_06.tif', './cog_2020_06.tif')
```

## Visualization

Create static and interactive visualizations:

**CLI:**

```bash
# Time series plot
pycropwat plot timeseries --input ./outputs \
    --start-year 2000 --end-year 2020 --output ./timeseries.png

# Monthly climatology bar chart
pycropwat plot climatology --input ./outputs \
    --start-year 2000 --end-year 2020 --output ./climatology.png

# Single month map
pycropwat plot map --input ./outputs --year 2020 --month 6 \
    --output ./map_2020_06.png

# Interactive HTML map (requires: pip install leafmap)
pycropwat plot interactive --input ./outputs --year 2020 --month 6 \
    --output ./map.html

# Compare two datasets side-by-side
pycropwat plot compare --input ./era5_outputs --other-input ./terraclimate_outputs \
    --year 2020 --month 6 --label1 ERA5 --label2 TerraClimate \
    --output ./comparison.png

# Scatter plot comparison with statistics
pycropwat plot scatter --input ./era5_outputs --other-input ./terraclimate_outputs \
    --start-year 2000 --end-year 2020 --output ./scatter.png

# Annual comparison bar chart
pycropwat plot annual-compare --input ./era5_outputs --other-input ./terraclimate_outputs \
    --start-year 2000 --end-year 2020 --output ./annual_compare.png
```

**Python:**

```python
from pycropwat.analysis import Visualizer

viz = Visualizer('./outputs')

# Time series plot
viz.plot_time_series(2000, 2020, output_path='./timeseries.png')

# Monthly climatology
viz.plot_monthly_climatology(2000, 2020, output_path='./climatology.png')

# Raster map
viz.plot_raster(2020, 6, output_path='./map_2020_06.png')

# Interactive map (requires leafmap or folium)
viz.plot_interactive_map(2020, 6, output_path='./map.html')

# Compare two datasets
viz.plot_comparison(
    other_input_dir='./terraclimate_outputs',
    year=2020, month=6,
    label1='ERA5', label2='TerraClimate',
    output_path='./comparison.png'
)

# Scatter plot with R², RMSE, bias
viz.plot_scatter_comparison(
    other_input_dir='./terraclimate_outputs',
    start_year=2000, end_year=2020,
    output_path='./scatter.png'
)

# Annual comparison
viz.plot_annual_comparison(
    other_input_dir='./terraclimate_outputs',
    start_year=2000, end_year=2020,
    output_path='./annual_compare.png'
)
```

## Custom Resolution

By default, pyCropWat uses the native resolution of the dataset. You can specify a custom resolution:

```bash
pycropwat process \
    --asset IDAHO_EPSCOR/TERRACLIMATE \
    --band pr \
    --geometry study_area.geojson \
    --start-year 2020 \
    --end-year 2023 \
    --scale 1000 \
    --output ./outputs
```

The `--scale` parameter is in meters.

## Output Files

pyCropWat generates two GeoTIFF files per month:

| File | Description | Units |
|------|-------------|-------|
| `effective_precip_YYYY_MM.tif` | Effective precipitation | mm |
| `effective_precip_fraction_YYYY_MM.tif` | Effective/Total precipitation ratio | fraction (0-1) |

Example output directory:

```
outputs/
├── effective_precip_2020_01.tif
├── effective_precip_2020_02.tif
├── ...
├── effective_precip_fraction_2020_01.tif
├── effective_precip_fraction_2020_02.tif
└── ...
```

## Next Steps

- **Try the Complete Workflow Example**: Run `python Examples/south_america_example.py --analysis-only` to see all features in action using real Rio de la Plata basin data
- See [CLI Reference](cli.md) for complete command documentation
- See [Python API](api.md) for advanced programmatic usage
- See [Examples](../examples.md) for real-world use cases including the [Complete Workflow](../examples.md#example-12-complete-workflow)

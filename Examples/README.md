# pyCropWat Complete Workflow Example

This directory contains a comprehensive example script demonstrating the full capabilities of the pyCropWat package for effective precipitation analysis.

## Overview

The `complete_workflow_example.py` script processes effective precipitation data from two different climate datasets (ERA5-Land and TerraClimate) for the Rio de la Plata Basin in South America, then performs extensive temporal aggregation, statistical analysis, visualization, and data comparison.

## Study Area

**Rio de la Plata Basin** - A major drainage basin in South America covering parts of:
- Argentina
- Brazil  
- Paraguay
- Uruguay
- Bolivia

GEE Asset: `projects/ssebop-471916/assets/Riodelaplata`

## Data Sources

| Dataset | GEE Asset | Precipitation Band | Scale Factor |
|---------|-----------|-------------------|--------------|
| ERA5-Land | `ECMWF/ERA5_LAND/MONTHLY_AGGR` | `total_precipitation_sum` | 1000 (m → mm) |
| TerraClimate | `IDAHO_EPSCOR/TERRACLIMATE` | `pr` | 1.0 (already mm) |

## Workflow Steps

### Step 1: Process Effective Precipitation from GEE

**Function:** `process_effective_precipitation()`

Calculates monthly effective precipitation using the CROPWAT method:
- Downloads monthly precipitation data from Google Earth Engine
- Applies the USDA SCS/CROPWAT formula to convert total precipitation to effective precipitation
- Exports monthly GeoTIFF rasters for both effective precipitation (mm) and effective precipitation fraction

**Outputs:**
- `RDP_ERA5Land/effective_precip_YYYY_MM.tif`
- `RDP_ERA5Land/effective_precip_fraction_YYYY_MM.tif`
- `RDP_TerraClimate/effective_precip_YYYY_MM.tif`
- `RDP_TerraClimate/effective_precip_fraction_YYYY_MM.tif`

---

### Step 2: Temporal Aggregation

**Function:** `run_temporal_aggregation()`

Aggregates monthly effective precipitation data into different temporal scales:

| Aggregation Type | Description | Output |
|-----------------|-------------|--------|
| **Annual Totals** | Sum of monthly values per year | `analysis_outputs/annual/{dataset}/` |
| **Seasonal** | Apr-Sep (growing season) totals | `analysis_outputs/growing_season/{dataset}/` |
| **Monthly Climatology** | Multi-year monthly means (2000-2020) | `analysis_outputs/climatology/{dataset}/` |

**Key Features:**
- Uses `TemporalAggregator` class from pyCropWat
- Pattern matching excludes fraction files (only processes `effective_precip_[0-9]*.tif`)
- Supports custom date ranges

---

### Step 3: Statistical Analysis

**Function:** `run_statistical_analysis()`

Performs advanced statistical analysis on the effective precipitation data:

#### 3a. Anomaly Calculation
- Computes percent anomalies for recent years (2021-2023)
- Reference period: 2000-2020 climatology
- Output: `analysis_outputs/anomalies/{dataset}/`

#### 3b. Trend Analysis
- Calculates pixel-wise trends using Mann-Kendall test and Sen's slope
- Produces trend magnitude and p-value rasters
- Output: `analysis_outputs/trend/{dataset}/`

#### 3c. Zonal Statistics
- Creates sample zones (Eastern RDP and Western RDP)
- Computes zonal statistics (mean, min, max, std) for each zone
- Generates summary plots and CSV exports
- Output: `analysis_outputs/zonal/{dataset}/`

---

### Step 4: Visualization

**Function:** `run_visualization()`

Creates comprehensive visualizations:

#### Time Series Plot
- Annual mean effective precipitation over the study period
- Output: `analysis_outputs/figures/{dataset}/time_series.png`

#### Monthly Climatology Chart
- Bar chart showing average monthly effective precipitation
- Output: `analysis_outputs/figures/{dataset}/monthly_climatology.png`

#### Spatial Maps
Static maps for different seasons/months:
- `map_2023_01.png` - January (Southern Hemisphere summer/wet season)
- `map_2023_06.png` - June (winter/dry season)
- `map_2023_10.png` - October (spring)

#### Notable Climate Event Maps
Maps highlighting significant climate events:
- `map_notable_2015_12.png` - El Niño event (Dec 2015)
- `map_notable_2010_08.png` - La Niña event (Aug 2010)

#### Interactive Maps
- HTML-based interactive maps using Leafmap/Folium
- Output: `analysis_outputs/figures/{dataset}/interactive_map_*.html`

---

### Step 5: Dataset Comparison

**Function:** `compare_datasets()` and `compare_zonal_statistics()`

Compares ERA5-Land and TerraClimate datasets:

#### Side-by-Side Spatial Comparison
- Dual-panel maps comparing the same month from both datasets
- Output: `analysis_outputs/comparisons/comparison_YYYY_MM.png`

#### Scatter Plot Comparison
- Pixel-by-pixel scatter plot between datasets
- Includes 1:1 line and correlation statistics
- Output: `analysis_outputs/comparisons/scatter_comparison.png`

#### Annual Bar Chart Comparison
- Side-by-side annual totals for both datasets
- Output: `analysis_outputs/comparisons/annual_comparison.png`

#### Zonal Statistics Comparison
- Compares zonal statistics between datasets for Eastern and Western RDP
- Annual time series and monthly climatology by zone
- Output: `analysis_outputs/comparisons/zonal_comparison.png`

---

### Step 6: Effective Precipitation Method Comparison

**Function:** `compare_peff_methods()`

Compares different effective precipitation calculation methods to help users understand how each method transforms total precipitation into effective precipitation.

#### Methods Compared
| Method | Formula | Description |
|--------|---------|-------------|
| **CROPWAT** | Peff = 0.6×P - 10 (P ≤ 70mm) or 0.8×P - 24 (P > 70mm) | USDA SCS method (default) |
| **FAO/AGLW** | Peff = 0.8×P - 25 (P > 75mm) or 0.6×P - 10 (P ≤ 75mm) | FAO Land and Water Division formula |
| **Fixed 70%** | Peff = 0.70 × P | Simple fixed percentage method |
| **Dependable Rain** | Peff = 0.75 × (P - 5) where P > 5mm | FAO Dependable Rainfall at 75% probability |

#### Method Comparison Maps
- 2x2 spatial maps showing Peff calculated by each method for the same month
- Common color scale across all methods for direct comparison
- Shared colorbar positioned outside the figure
- Output: `analysis_outputs/method_comparison/{dataset}_method_maps_YYYY_MM.png`

#### Theoretical Method Curves
- Line plots showing how each method transforms precipitation to effective precipitation (0-400 mm range)
- Efficiency curves showing Peff/P ratio across precipitation values
- Highlights breakpoints where methods change behavior
- Output: `analysis_outputs/method_comparison/{dataset}_method_curves.png`

#### Cross-Dataset Summary
- Box plots of Peff distribution by method
- Bar charts comparing methods across datasets
- Rainfall effectiveness comparison
- Output: `analysis_outputs/method_comparison/method_comparison_summary.png`

#### Statistics Export
- CSV file with mean, std, min, max for each dataset-method combination
- Output: `analysis_outputs/method_comparison/method_statistics.csv`

---

### Step 7: Export to NetCDF

**Function:** `export_netcdf()`

Exports the complete time series to CF-compliant NetCDF format:
- Includes metadata and coordinate reference system
- Single file per dataset containing all months
- Output: `analysis_outputs/{dataset}_effective_precip_2000_2023.nc`

---

## Output Directory Structure

```
Examples/
├── complete_workflow_example.py
├── README.md                       # This file
├── RDP_ERA5Land/                   # ERA5-Land monthly outputs
│   ├── effective_precip_YYYY_MM.tif
│   └── effective_precip_fraction_YYYY_MM.tif
├── RDP_TerraClimate/               # TerraClimate monthly outputs
│   ├── effective_precip_YYYY_MM.tif
│   └── effective_precip_fraction_YYYY_MM.tif
└── analysis_outputs/
    ├── annual/
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── climatology/
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── growing_season/
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── anomalies/
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── trend/
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── zonal/
    │   ├── sample_zones.geojson
    │   ├── ERA5Land/
    │   │   ├── zonal_stats.csv
    │   │   └── zonal_summary.png
    │   └── TerraClimate/
    │       ├── zonal_stats.csv
    │       └── zonal_summary.png
    ├── figures/
    │   ├── ERA5Land/
    │   │   ├── time_series.png
    │   │   ├── monthly_climatology.png
    │   │   ├── map_2023_01.png
    │   │   ├── map_2023_06.png
    │   │   ├── map_2023_10.png
    │   │   ├── map_notable_2015_12.png
    │   │   ├── map_notable_2010_08.png
    │   │   └── interactive_map_*.html
    │   └── TerraClimate/
    ├── comparisons/
    │   ├── comparison_2023_06.png
    │   ├── scatter_comparison.png
    │   ├── annual_comparison.png
    │   └── zonal_comparison.png
    ├── method_comparison/          # Peff method comparison
    │   ├── ERA5Land_method_maps_2020_01.png
    │   ├── ERA5Land_method_curves.png
    │   ├── TerraClimate_method_maps_2020_01.png
    │   ├── TerraClimate_method_curves.png
    │   ├── method_comparison_summary.png
    │   └── method_statistics.csv
    ├── ERA5Land_effective_precip_2000_2023.nc
    └── TerraClimate_effective_precip_2000_2023.nc
```

## Usage

### Full Workflow (includes GEE processing)

```bash
cd Examples
python complete_workflow_example.py
```

### Analysis Only (skip GEE processing)

If you already have the monthly effective precipitation rasters:

```bash
python complete_workflow_example.py --analysis-only
```

### Force Reprocessing

To reprocess data even if outputs exist:

```bash
python complete_workflow_example.py --force-reprocess
```

### Specify GEE Project

```bash
python complete_workflow_example.py --gee-project your-project-id
```

## Requirements

- Python 3.9+
- pyCropWat package installed
- Google Earth Engine account
- Access to the study area GEE asset

### Python Dependencies

Install pyCropWat with all dependencies:

```bash
# Basic installation
pip install -e .

# With interactive map support (leafmap, localtileserver)
pip install -e ".[interactive]"
```

Core dependencies (installed automatically):
- `earthengine-api`
- `rasterio`
- `xarray`
- `dask`
- `matplotlib`
- `pandas`
- `scipy`
- `folium`
- `rasterstats`
- `netCDF4`

Optional dependencies (for interactive maps):
- `leafmap`
- `localtileserver`

## Configuration

Edit the following variables in the script to customize:

```python
# GEE Project ID
GEE_PROJECT = "your-gee-project-id"

# Time period
START_YEAR = 2000
END_YEAR = 2023

# Climatology reference period
CLIM_START = 2000
CLIM_END = 2020

# Study area asset
STUDY_AREA_ASSET = "projects/ssebop-471916/assets/Riodelaplata"
```

## Sample Zones

The script creates two sample zones for demonstrating zonal statistics:

| Zone | Region | Coordinates |
|------|--------|-------------|
| Eastern RDP | Uruguay, SE Brazil | Lon: -54° to -50°, Lat: -35° to -28° |
| Western RDP | N Argentina, Paraguay | Lon: -62° to -54°, Lat: -35° to -25° |

## Notes

- The script uses pattern matching (`effective_precip_[0-9]*.tif`) to exclude fraction files from analysis
- Zero values are masked as nodata in visualizations
- Maps include degree symbols (°) in latitude/longitude labels
- Interactive maps require a web browser to view

## See Also

- [pyCropWat Documentation](https://montimaj.github.io/pyCropWat/)
- [API Reference](https://montimaj.github.io/pyCropWat/api/core/)
- [Installation Guide](https://montimaj.github.io/pyCropWat/installation/)

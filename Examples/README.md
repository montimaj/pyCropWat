# pyCropWat Workflow Examples

This directory contains comprehensive example scripts demonstrating the full capabilities of the pyCropWat package for effective precipitation analysis.

## Available Examples

| Script | Region | Method | Precipitation Data | Comparisons |
|--------|--------|--------|-------------------|-------------|
| `south_america_example.py` | Rio de la Plata (South America) | 8 methods | ERA5-Land, TerraClimate | Method comparison |
| `arizona_example.py` | Arizona (U.S.) | 8 methods | GridMET, PRISM, ERA5-Land, TerraClimate | U.S. vs Global, Method comparison |
| `new_mexico_example.py` | New Mexico (U.S.) | 8 methods | PRISM | Method comparison |
| `western_us_pcml_example.py` | Western U.S. (17 states) | PCML | PCML GEE Asset | Water year totals, Mean fractions |
| `ucrb_example.py` | Upper Colorado River Basin | 8 methods | GeoPackage (pre-existing) | Field-scale aggregation |

---

## 1. Western U.S. PCML Example

### Overview

The `western_us_pcml_example.py` script demonstrates the Physics-Constrained Machine Learning (PCML) method for effective precipitation estimation in the Western United States. This example downloads pre-computed PCML Peff values and calculates water year totals.

### Study Area

**Western United States** - 17 states: AZ, CA, CO, ID, KS, MT, NE, NV, NM, ND, OK, OR, SD, TX, UT, WA, WY

### PCML Method

| Attribute | Value |
|-----------|-------|
| **GEE Asset** | `projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked` |
| **Temporal** | January 2000 - September 2024 (monthly) |
| **Resolution** | ~2 km (native scale retrieved dynamically from GEE asset) |
| **Band Format** | `bYYYY_M` (e.g., `b2015_9`, `b2016_10`) |
| **Fraction Calculation** | Annual (water year, Oct-Sep) fractions only from GEE asset (WY 2000-2024, band format: `bYYYY`) |

### PCML Geometry Options

- **No geometry provided**: Downloads the entire PCML asset (full Western U.S. - 17 states)
- **User provides geometry**: PCML data is clipped/subsetted to that geometry. **Note:** Only Western U.S. vectors that overlap with the 17-state extent can be used (e.g., AZ.geojson, pacific_northwest.geojson)

### Water Year Definition

- Western U.S. water year runs from **October to September**
- Water Year 2001 = October 2000 - September 2001
- Water Year 2024 = October 2023 - September 2024

### Workflow Steps

1. **Download PCML Data** - Downloads monthly PCML Peff and annual fractions from GEE assets
2. **Calculate Water Year Totals** - Aggregates Oct-Sep totals for each water year
3. **Calculate Mean Fractions** - Computes long-term mean Peff fractions (WY 2001-2024)
4. **Visualization** - Creates time series, maps, and summary figures

**Note:** For PCML, `save_inputs=False` is recommended since the PCML effective precipitation comes directly from the GEE asset (not calculated from raw precipitation). Setting `save_inputs=True` would save redundant copies of the data.

### Outputs

| Output Type | Location | Description |
|-------------|----------|-------------|
| Monthly Peff | `WesternUS_PCML/effective_precip_YYYY_MM.tif` | Monthly PCML effective precipitation |
| Annual Fractions | `WesternUS_PCML/effective_precip_fraction_YYYY.tif` | Annual (water year, Oct-Sep) Peff/Precip fractions |
| Water Year Totals | `WesternUS_PCML/water_year/pcml_peff_wy*.tif` | Oct-Sep total Peff |
| Mean Fractions | `WesternUS_PCML/fractions/pcml_mean_fraction_*.tif` | Long-term mean fractions |
| Figures | `WesternUS_PCML/figures/` | Time series, maps, summary plots |

### Usage

```bash
# Basic usage (4 workers)
python western_us_pcml_example.py

# Use 8 workers
python western_us_pcml_example.py -w 8

# Skip download, run analysis only
python western_us_pcml_example.py --analysis-only

# Specify GEE project
python western_us_pcml_example.py --gee-project your-project-id
```

### Reference

Hasan, M. F., Smith, R. G., Majumdar, S., Huntington, J. L., Alves Meira Neto, A., & Minor, B. A. (2025). Satellite data and physics-constrained machine learning for estimating effective precipitation in the Western United States and application for monitoring groundwater irrigation. *Agricultural Water Management*, 319, 109821. https://doi.org/10.1016/j.agwat.2025.109821

---

## 2. Upper Colorado River Basin (UCRB) Field-Scale Example

### Overview

The `ucrb_example.py` script demonstrates **field-scale effective precipitation aggregation** using existing precipitation volumes stored in a GeoPackage file. This example shows how to apply pyCropWat methods to pre-computed field-level data without requiring new GEE downloads.

### Key Features

- **Field-scale analysis**: Works with individual agricultural fields/parcels
- **GeoPackage input**: Reads existing precipitation volumes from vector data
- **No GEE precipitation download**: Uses pre-existing monthly precipitation volumes
- **Zonal AWC statistics**: Downloads SSURGO AWC data and calculates mean AWC per field
- **Unit conversion**: Handles acre-feet ↔ mm conversions for method application
- **Parallel processing**: Uses Dask for efficient multi-method calculation

### Study Area

**Upper Colorado River Basin (UCRB)** - Agricultural fields in the upper Colorado River watershed

### Data Requirements

| Data | Source | Description |
|------|--------|-------------|
| **Precipitation** | GeoPackage layer `PPT_VOLUME_ET_Demands` | Monthly volumes in acre-feet (columns: `PPT_VOLUME_mm_yy_acft`) |
| **ETo** | GeoPackage layer `ETDa_VOLUME_ET_Demands` | Monthly ET volumes in acre-feet (columns: `ETDa_VOLUME_mm_yy_acft`) |
| **AWC** | GEE SSURGO asset | Downloaded and cached as zonal means per field |
| **Field geometry** | GeoPackage | Field polygons with `DRI_ID`, `ACRES_FTR_GEOM` |

### Effective Precipitation Methods (8 methods)

| Method | Description | Requires AWC/ETo |
|--------|-------------|------------------|
| CROPWAT | FAO CROPWAT method | No |
| FAO/AGLW | FAO Dependable Rainfall (80% exceedance) | No |
| Fixed 70% | Simple empirical method | No |
| Dependable Rainfall | 75% probability level | No |
| FarmWest | WSU irrigation scheduling formula | No |
| USDA-SCS (MAD=0.5) | Site-specific method, default MAD factor | Yes |
| USDA-SCS (MAD=1.0) | Site-specific method, full soil storage | Yes |
| USDA-SCS (MAD=1.0, RD=2m) | Full soil storage, 2m rooting depth | Yes |

### Unit Conversion Workflow

```
Input: Precipitation volume (acre-feet)
    ↓ Divide by field area (acres) → Depth (feet)
    ↓ Convert feet to mm (× 304.8)
    ↓ Apply pyCropWat Peff method
    ↓ Convert mm back to feet (÷ 304.8)
    ↓ Multiply by field area (acres)
Output: Effective precipitation volume (acre-feet)
```

### Workflow Steps

1. **Load GeoPackage** - Read precipitation and ETo volumes from existing layers
2. **Download AWC** - Get SSURGO AWC data from GEE, calculate zonal mean per field
3. **Apply Methods** - Calculate Peff for all 8 methods across all time periods
4. **Save Results** - Add new layers to GeoPackage with Peff volumes
5. **Generate Statistics** - Create CSV summaries by year and method

### Outputs

| Output Type | Location | Description |
|-------------|----------|-------------|
| Peff Layers | Input GeoPackage | New layers: `Peff_CROPWAT`, `Peff_FAO_AGLW`, etc. |
| AWC Cache | `UCRB/awc_data/field_awc_values.csv` | Mean AWC per field (DRI_ID) |
| Annual Totals | `UCRB/summaries/annual_totals_by_method.csv` | Annual Peff totals by method |
| Method Stats | `UCRB/summaries/method_comparison_statistics.csv` | Method comparison statistics |

### Usage

```bash
# Basic usage (8 workers)
python ucrb_example.py

# Specify GEE project
python ucrb_example.py --gee-project your-project-id

# Skip AWC download if cached
python ucrb_example.py --skip-awc-download

# Use 16 workers
python ucrb_example.py -w 16

# Verbose logging
python ucrb_example.py -v
```

### Python API for Field-Scale Aggregation

For custom field-scale workflows, you can use pyCropWat methods directly:

```python
import numpy as np
import geopandas as gpd
from pycropwat.methods import (
    cropwat_effective_precip,
    usda_scs_effective_precip
)

# Load field data
gdf = gpd.read_file('fields.gpkg', layer='precipitation')

# Convert volume to depth
FT_TO_MM = 304.8
ppt_volume_acft = gdf['PPT_VOLUME_01_20_acft'].values
area_acres = gdf['ACRES'].values
ppt_mm = (ppt_volume_acft / area_acres) * FT_TO_MM

# Apply pyCropWat method
peff_mm = cropwat_effective_precip(ppt_mm)

# Convert back to volume
peff_acft = (peff_mm / FT_TO_MM) * area_acres

# For USDA-SCS with AWC and ETo
awc_values = gdf['mean_awc'].values  # mm/m
eto_mm = (gdf['ETO_VOLUME_01_20_acft'].values / area_acres) * FT_TO_MM
peff_scs_mm = usda_scs_effective_precip(
    ppt_mm, eto_mm, awc_values, 
    rooting_depth=1.0, mad_factor=0.5
)
```

---

## 3. Rio de la Plata Example (Global)

### Overview

The `south_america_example.py` script processes effective precipitation data from two different climate datasets (ERA5-Land and TerraClimate) for the Rio de la Plata Basin in South America, then performs extensive temporal aggregation, statistical analysis, visualization, and data comparison.

### Study Area

**Rio de la Plata Basin** - A major drainage basin in South America covering parts of:
- Argentina
- Brazil  
- Paraguay
- Uruguay
- Bolivia

GEE Asset: `projects/ssebop-471916/assets/Riodelaplata`

### Data Sources

| Dataset | GEE Asset | Precipitation Band | Scale Factor |
|---------|-----------|-------------------|--------------|
| ERA5-Land | `ECMWF/ERA5_LAND/MONTHLY_AGGR` | `total_precipitation_sum` | 1000 (m → mm) |
| TerraClimate | `IDAHO_EPSCOR/TERRACLIMATE` | `pr` | 1.0 (already mm) |

### Workflow Steps

### Step 1: Process Effective Precipitation from GEE

**Function:** `process_effective_precipitation()`

Calculates monthly effective precipitation using the CROPWAT method:
- Downloads monthly precipitation data from Google Earth Engine
- Optionally saves downloaded input data to `analysis_inputs/` folder
- Applies the CROPWAT formula to convert total precipitation to effective precipitation
- Exports monthly GeoTIFF rasters for both effective precipitation (mm) and effective precipitation fraction

**Input Data (saved to `analysis_inputs/`):**
- `RioDelaPlata/analysis_inputs/RDP_ERA5Land/precip_YYYY_MM.tif`
- `RioDelaPlata/analysis_inputs/RDP_TerraClimate/precip_YYYY_MM.tif`

**Outputs:**
- `RioDelaPlata/RDP_ERA5Land/effective_precip_YYYY_MM.tif`
- `RioDelaPlata/RDP_ERA5Land/effective_precip_fraction_YYYY_MM.tif`
- `RioDelaPlata/RDP_TerraClimate/effective_precip_YYYY_MM.tif`
- `RioDelaPlata/RDP_TerraClimate/effective_precip_fraction_YYYY_MM.tif`

---

### Step 2: Temporal Aggregation

**Function:** `run_temporal_aggregation()`

Aggregates monthly effective precipitation data into different temporal scales:

| Aggregation Type | Description | Output |
|-----------------|-------------|--------|
| **Annual Totals** | Sum of monthly values per year | `RioDelaPlata/analysis_outputs/annual/{dataset}/` |
| **Seasonal** | Apr-Sep (growing season) totals | `RioDelaPlata/analysis_outputs/growing_season/{dataset}/` |
| **Monthly Climatology** | Multi-year monthly means (2000-2020) | `RioDelaPlata/analysis_outputs/climatology/{dataset}/` |

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
- Output: `RioDelaPlata/analysis_outputs/anomalies/{dataset}/`

#### 3b. Trend Analysis
- Calculates pixel-wise trends using Mann-Kendall test and Sen's slope
- Produces trend magnitude and p-value rasters
- Output: `RioDelaPlata/analysis_outputs/trend/{dataset}/`

#### 3c. Zonal Statistics
- Creates sample zones (Eastern RDP and Western RDP)
- Computes zonal statistics (mean, min, max, std) for each zone
- Generates summary plots and CSV exports
- Output: `RioDelaPlata/analysis_outputs/zonal/{dataset}/`

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

### Step 4b: Anomaly, Climatology, and Trend Maps

**Functions:** `plot_anomaly_maps()`, `plot_climatology_maps()`, `plot_trend_maps()`

Creates specialized maps from statistical analysis outputs using the `Visualizer` class methods:

#### Anomaly Maps
- Percent anomalies using diverging RdBu colormap (red=dry, blue=wet)
- Centered at 100% (normal conditions)
- Plotted for recent years (2021-2023) during wet/monsoon season
- Output: `analysis_outputs/figures/{dataset}/anomaly_map_YYYY_MM.png`

#### Climatology Maps
- Monthly climatology spatial distribution
- Shows typical monthly effective precipitation patterns
- Key months: January (wet), June (dry), October (spring), December (late spring)
- Output: `analysis_outputs/figures/{dataset}/climatology_map_MM.png`

#### Trend Maps
- **Trend Panel**: Two-panel figure with slope and p-value maps side by side
  - Slope map using RdBu colormap (red=decreasing, blue=increasing)
  - P-value map using RdYlGn_r colormap (green=significant)
  - Output: `analysis_outputs/figures/{dataset}/trend_maps.png`
- **Trend with Significance**: Single map with stippling overlay for significant trends (p < 0.05)
  - Output: `analysis_outputs/figures/{dataset}/trend_map_with_significance.png`

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
| **Ensemble** | Mean of 6 methods | Ensemble mean (default, excludes TAGEM-SuET and PCML) |
| **CROPWAT** | Peff = P×(125-0.2P)/125 (P ≤ 250mm) or 0.1P+125 (P > 250mm) | CROPWAT method |
| **FAO/AGLW** | Peff = 0.6×P - 10 (P ≤ 70mm) or 0.8×P - 24 (P > 70mm) | FAO Dependable Rainfall (80% exceedance) |
| **Fixed 70%** | Peff = 0.70 × P | Simple fixed percentage method |
| **Dependable Rain** | Peff = 0.6×P - 10 (P ≤ 70mm) or 0.8×P - 24 (P > 70mm) | Same as FAO/AGLW (80% probability) |
| **FarmWest** | Peff = (P - 5) × 0.75 | Pacific Northwest method |
| **USDA-SCS** | Complex formula with AWC and ETo | Site-specific method |
| **TAGEM-SuET** | P ≤ ETo: 0; P - ETo < 75: P - ETo; else: 75 + f(excess) | Turkish Irrigation Management System |

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
├── south_america_example.py
├── arizona_example.py
├── new_mexico_example.py
├── AZ.geojson                      # Arizona boundary GeoJSON
├── NM.geojson                      # New Mexico boundary GeoJSON
├── README.md                       # This file
├── RioDelaPlata/                   # Rio de la Plata region outputs
│   ├── RDP_ERA5Land/               # ERA5-Land monthly outputs
│   │   ├── effective_precip_YYYY_MM.tif
│   │   └── effective_precip_fraction_YYYY_MM.tif
│   ├── RDP_TerraClimate/           # TerraClimate monthly outputs
│   │   ├── effective_precip_YYYY_MM.tif
│   │   └── effective_precip_fraction_YYYY_MM.tif
│   ├── analysis_inputs/            # Downloaded input data (optional)
│   └── analysis_outputs/
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
    │   │   ├── interactive_map_*.html
    │   │   ├── anomaly_map_2023_01.png     # Anomaly maps
    │   │   ├── anomaly_map_2022_06.png
    │   │   ├── climatology_map_01.png      # Climatology maps
    │   │   ├── climatology_map_06.png
    │   │   ├── trend_maps.png              # Trend panel (slope + p-value)
    │   │   └── trend_map_with_significance.png
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
python south_america_example.py
```

### Analysis Only (skip GEE processing)

If you already have the monthly effective precipitation rasters:

```bash
python south_america_example.py --analysis-only
```

### Force Reprocessing

To reprocess data even if outputs exist:

```bash
python south_america_example.py --force-reprocess
```

### Specify GEE Project

```bash
python south_america_example.py --gee-project your-project-id
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

1. Processing large datasets requires significant time and computational resources
2. GEE quotas may limit concurrent requests
3. Analysis-only mode is much faster if you already have processed data
4. Interactive maps require `leafmap` installation

---

## 4. Arizona USDA-SCS Example (U.S.)

### Overview

The `arizona_example.py` script demonstrates the **USDA-SCS method** for effective precipitation using U.S.-specific high-resolution datasets for Arizona, with **comparisons to global datasets** and **8 Peff methods** (excludes PCML).

### Study Area

**Arizona** - A southwestern U.S. state with diverse precipitation patterns:
- Summer monsoon (July-September) 
- Winter Pacific storms (December-February)
- Arid climate with spatial variability

GEE Asset: `users/montimajumdar/AZ`

### Data Sources

#### Precipitation Datasets (U.S. + Global)

| Dataset | Type | GEE Asset | Band | Resolution |
|---------|------|-----------|------|------------|
| GridMET | U.S. | `IDAHO_EPSCOR/GRIDMET` | `pr` | ~4 km |
| PRISM | U.S. | `projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY` | `ppt` | ~800 m |
| ERA5-Land | Global | `ECMWF/ERA5_LAND/MONTHLY_AGGR` | `total_precipitation_sum` | ~11 km |
| TerraClimate | Global | `IDAHO_EPSCOR/TERRACLIMATE` | `pr` | ~4 km |

#### USDA-SCS Required Data (Dataset-Specific)

**U.S. Datasets (High-Resolution):**

| Dataset | GEE Asset | Band | Description |
|---------|-----------|------|-------------|
| AWC (SSURGO) | `projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite` | (single band) | USDA SSURGO soil water holding capacity |
| ETo (GridMET) | `projects/openet/assets/reference_et/conus/gridmet/monthly/v1` | `eto` | OpenET GridMET monthly reference ET |

**Global Datasets:**

| Dataset | GEE Asset | Band | Description |
|---------|-----------|------|-------------|
| AWC (FAO HWSD) | `projects/sat-io/open-datasets/FAO/HWSD_V2_SMU` | `AWC` | FAO Harmonized World Soil Database v2 |
| ETo (AgERA5) | `projects/climate-engine-pro/assets/ce-ag-era5-v2/daily` | `ReferenceET_PenmanMonteith_FAO56` | AgERA5 daily FAO-56 Penman-Monteith reference ET |

### USDA-SCS Method

The USDA-SCS method calculates effective precipitation using:
1. **AWC (Available Water Capacity)**: U.S. uses SSURGO, Global uses FAO HWSD
2. **ETo (Reference Evapotranspiration)**: U.S. uses GridMET, Global uses AgERA5 daily (auto-aggregated to monthly)
3. **Rooting Depth**: User-specified crop root zone depth (default 1.0 m)
4. **ETo Scale Factor**: Optional scaling for unit conversion via `eto_scale_factor`

This method is more site-specific than CROPWAT because it accounts for local soil conditions and evaporative demand.

### Workflow Steps

1. **Process Effective Precipitation** - Uses USDA-SCS method with SSURGO AWC and GridMET ETo
   - Downloads precipitation, AWC, and ETo data from GEE
   - Optionally saves input data to `analysis_inputs/` folder (enabled by default)
2. **Process All Precipitation Sources** - U.S. (GridMET, PRISM) and Global (ERA5-Land, TerraClimate)
3. **U.S. Dataset Comparison** - GridMET (~4km) vs PRISM (~800m) scatter plots and maps
4. **U.S. vs Global Comparison** - Compares U.S. datasets against global datasets
5. **Method Comparison** - Compares all 8 Peff methods across all datasets
6. **Arizona-Specific Aggregation** - Monsoon season (Jul-Sep), winter season (Jan-Feb)
7. **Zonal Statistics** - Central AZ, Southern AZ, Northern AZ regions
8. **Anomaly, Climatology, and Trend Maps** - Spatial visualization of statistical outputs
9. **Export to NetCDF** - CF-compliant NetCDF for each dataset

### Visualization Outputs

The workflow produces comprehensive visualizations using the `Visualizer` class:

#### Standard Visualizations
- Time series plots
- Monthly climatology bar charts
- Static raster maps for key months
- Interactive HTML maps

#### Statistical Maps (New)
- **Anomaly Maps**: Percent anomalies using RdBu colormap (red=dry, blue=wet), centered at 100%
  - `figures/{dataset}/anomaly_map_2023_08.png`
- **Climatology Maps**: Monthly climatology spatial patterns
  - `figures/{dataset}/climatology_map_08.png`
- **Trend Maps**: Slope and significance maps
  - `figures/{dataset}/trend_maps.png` (slope + p-value panel)
  - `figures/{dataset}/trend_map_with_significance.png` (stippling for p < 0.05)

### Input Data Organization

When `save_inputs=True` (default), downloaded input data is saved to region-specific folders:

```
Arizona/analysis_inputs/
├── AZ_GridMET_USDA_SCS/
│   ├── precip_YYYY_MM.tif      # Monthly precipitation
│   ├── awc.tif                  # Available Water Capacity (static)
│   └── eto_YYYY_MM.tif         # Monthly Reference ET
├── AZ_PRISM_USDA_SCS/
├── AZ_ERA5Land_USDA_SCS/
└── AZ_TerraClimate_USDA_SCS/
```

### U.S. vs Global Dataset Comparison

The workflow includes a comprehensive comparison between:
- **U.S. Datasets**: GridMET (~4km), PRISM (~800m)
- **Global Datasets**: ERA5-Land (~11km), TerraClimate (~4km)

Outputs include:
- Side-by-side spatial maps for each dataset
- Box plots comparing distributions
- Monthly climatology by dataset region
- Annual time series by U.S. vs Global

### Effective Precipitation Method Comparison

Compares **all 8 methods** available in pyCropWat (excludes PCML):

| Method | Description |
|--------|-------------|
| Ensemble | Mean of 6 methods (default, excludes TAGEM-SuET and PCML) |
| CROPWAT | CROPWAT formula |
| FAO/AGLW | FAO Dependable Rainfall (80% exceedance) |
| Fixed 70% | Simple fixed percentage |
| Dependable Rain | Same as FAO/AGLW (80% probability) |
| FarmWest | Pacific Northwest method |
| USDA-SCS | Full USDA-SCS with AWC and ETo |
| TAGEM-SuET | Turkish Irrigation Management System (P - ETo) |
| PCML | Physics-Constrained ML (Western U.S. only, 2000-2024) |

Outputs include:
- 2×3 spatial comparison maps for each dataset
- Theoretical method curves (using Arizona-typical AWC=150mm/m, ETo=180mm/month)
- Cross-dataset method statistics

### Running the Script

#### Full Workflow (with GEE processing)
```bash
python arizona_example.py --gee-project your-project-id --workers 8
```

#### Analysis Only
```bash
python arizona_example.py --analysis-only
```

#### Force Reprocessing
```bash
python arizona_example.py --force-reprocess --gee-project your-project-id
```

### CLI Equivalent

```bash
# GridMET with USDA-SCS
pycropwat process --asset IDAHO_EPSCOR/GRIDMET --band pr \
    --gee-geometry users/montimajumdar/AZ \
    --start-year 2000 --end-year 2024 \
    --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto --rooting-depth 1.0 --mad-factor 0.5 \
    --workers 8 --output ./Arizona/AZ_GridMET_USDA_SCS

# PRISM with USDA-SCS  
pycropwat process --asset projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY --band ppt \
    --gee-geometry users/montimajumdar/AZ \
    --start-year 2000 --end-year 2024 \
    --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto --rooting-depth 1.0 --mad-factor 0.5 \
    --workers 8 --output ./Arizona/AZ_PRISM_USDA_SCS
```

### Sample Zones

| Zone | Region | Description |
|------|--------|-------------|
| Central AZ | Phoenix Metropolitan | Salt River Valley agricultural area |
| Southern AZ | Tucson/Santa Cruz Valley | Semi-arid with summer monsoon |
| Northern AZ | Flagstaff/Colorado Plateau | Higher elevation with winter precip |

### Output Structure

```
Examples/
├── south_america_example.py   # Rio de la Plata workflow
├── arizona_example.py         # Arizona workflow
├── new_mexico_example.py      # New Mexico workflow
├── AZ.geojson                         # Arizona boundary geometry
├── NM.geojson                         # New Mexico boundary geometry
├── RioDelaPlata/                      # Rio de la Plata region outputs
│   ├── RDP_ERA5Land/
│   │   ├── effective_precip_YYYY_MM.tif
│   │   └── effective_precip_fraction_YYYY_MM.tif
│   ├── RDP_TerraClimate/
│   │   ├── effective_precip_YYYY_MM.tif
│   │   └── effective_precip_fraction_YYYY_MM.tif
│   ├── analysis_inputs/               # Downloaded input data
│   │   ├── RDP_ERA5Land/
│   │   │   └── precip_YYYY_MM.tif
│   │   └── RDP_TerraClimate/
│   │       └── precip_YYYY_MM.tif
│   └── analysis_outputs/
│       ├── annual/
│       ├── climatology/
│       ├── anomalies/
│       ├── trend/
│       ├── figures/
│       │   └── {dataset}/
│       │       ├── time_series.png
│       │       ├── monthly_distribution.png
│       │       ├── spatial_comparison.png
│       │       ├── anomaly_map_YYYY_MM.png       # NEW
│       │       ├── climatology_map_MM.png        # NEW
│       │       ├── trend_maps.png                # NEW
│       │       └── trend_map_with_significance.png  # NEW
│       ├── comparisons/
│       ├── zonal_stats/
│       ├── method_comparison/
│       └── rdp_zones.geojson
└── Arizona/                           # Arizona region outputs
    ├── AZ_GridMET_USDA_SCS/
    │   ├── effective_precip_YYYY_MM.tif
    │   └── effective_precip_fraction_YYYY_MM.tif
    ├── AZ_PRISM_USDA_SCS/
    │   ├── effective_precip_YYYY_MM.tif
    │   └── effective_precip_fraction_YYYY_MM.tif
    ├── AZ_ERA5Land_USDA_SCS/
    │   ├── effective_precip_YYYY_MM.tif
    │   └── effective_precip_fraction_YYYY_MM.tif
    ├── AZ_TerraClimate_USDA_SCS/
    │   ├── effective_precip_YYYY_MM.tif
    │   └── effective_precip_fraction_YYYY_MM.tif
    ├── analysis_inputs/               # Downloaded input data (precip, AWC, ETo)
    │   ├── AZ_GridMET_USDA_SCS/
    │   │   ├── precip_YYYY_MM.tif
    │   │   ├── awc.tif
    │   │   └── eto_YYYY_MM.tif
    │   ├── AZ_PRISM_USDA_SCS/
    │   ├── AZ_ERA5Land_USDA_SCS/
    │   └── AZ_TerraClimate_USDA_SCS/
    └── analysis_outputs/
        ├── annual/
        ├── climatology/
        ├── anomalies/
        ├── trend/
        ├── figures/
        │   └── {dataset}/
        │       ├── time_series.png
        │       ├── monthly_distribution.png
        │       ├── spatial_comparison.png
        │       ├── anomaly_map_YYYY_MM.png       # NEW
        │       ├── climatology_map_MM.png        # NEW
        │       ├── trend_maps.png                # NEW
        │       └── trend_map_with_significance.png  # NEW
        ├── comparisons/
        ├── zonal_stats/
        ├── us_vs_global/
        ├── method_comparison/
        ├── az_zones.geojson
        └── *_usda_scs_peff_2000_2024.nc
```

---

## Additional Notes

- The script uses pattern matching (`effective_precip_[0-9]*.tif`) to exclude fraction files from analysis
- Zero values are masked as nodata in visualizations
- Maps include degree symbols (°) in latitude/longitude labels
- Interactive maps require a web browser to view

## See Also

- [pyCropWat Documentation](https://montimaj.github.io/pyCropWat/)
- [API Reference](https://montimaj.github.io/pyCropWat/api/core/)
- [Installation Guide](https://montimaj.github.io/pyCropWat/installation/)

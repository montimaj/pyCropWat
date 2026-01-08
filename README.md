# pyCropWat

[![Release](https://img.shields.io/badge/release-v1.0.0-green.svg)](https://github.com/montimaj/pyCropWat/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://montimaj.github.io/pyCropWat)
[![GEE](https://img.shields.io/badge/Google%20Earth%20Engine-4285F4?logo=google-earth&logoColor=white)](https://earthengine.google.com/)

A Python package for calculating CROPWAT effective precipitation from Google Earth Engine climate data.

## Project Structure

```
pyCropWat/
├── pycropwat/               # Main package
│   ├── __init__.py          # Package exports
│   ├── core.py              # EffectivePrecipitation class
│   ├── methods.py           # Effective precipitation methods (CROPWAT, FAO/AGLW, etc.)
│   ├── analysis.py          # Temporal aggregation, statistics, visualization
│   ├── utils.py             # Utility functions (geometry loading, GEE init)
│   └── cli.py               # Command-line interface
├── tests/                   # Unit tests
│   ├── __init__.py
│   └── test_core.py
├── docs/                    # MkDocs documentation
│   ├── index.md             # Documentation home
│   ├── installation.md      # Installation guide
│   ├── examples.md          # Usage examples
│   ├── contributing.md      # Contribution guidelines
│   ├── api/                 # API reference
│   │   ├── cli.md
│   │   ├── core.md
│   │   └── utils.md
│   └── user-guide/          # User guide
│       ├── api.md
│       ├── cli.md
│       └── quickstart.md
├── Examples/                # Example outputs and scripts
│   ├── README.md               # Detailed workflow documentation
│   ├── complete_workflow_example.py  # Complete workflow script
│   ├── RDP_ERA5Land.zip           # ERA5-Land outputs for Rio de la Plata
│   ├── RDP_TerraClimate.zip       # TerraClimate outputs for Rio de la Plata
│   └── analysis_outputs/       # Analysis outputs from running the workflow script
├── .github/                 # GitHub configuration
│   └── workflows/
│       └── docs.yml         # GitHub Pages deployment workflow
├── CHANGELOG.md             # Release notes
├── mkdocs.yml               # MkDocs configuration
├── environment.yml          # Conda environment file
├── pyproject.toml           # Package configuration
├── requirements.txt         # pip dependencies
├── LICENSE
└── README.md
```

**Note:** The `Examples/` folder contains:
- **`README.md`**: Detailed step-by-step documentation of the complete workflow example.
- **`complete_workflow_example.py`**: A comprehensive Python script demonstrating the complete pyCropWat workflow including data processing, temporal aggregation, statistical analysis, visualization, and dataset comparison using real Rio de la Plata data.
- **Sample output rasters**: Monthly TerraClimate and ERA5-Land-based effective precipitation and effective precipitation fraction TIF files (2000-2025) for the Rio de la Plata (RDP) region of South America.

See the [Complete Workflow Example](#complete-workflow-example) section below for details.

**Changelog:** See [CHANGELOG.md](CHANGELOG.md) for release notes and version history.

## Overview

pyCropWat converts precipitation data from any GEE climate dataset into effective precipitation and effective precipitation fraction rasters. It supports:

- Any GEE ImageCollection with precipitation data
- Shapefile, GeoJSON, or GEE FeatureCollection asset for region of interest
- **Multiple effective precipitation methods**: CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall
- Parallel processing using Dask
- Monthly output rasters in GeoTIFF format
- **Temporal aggregation**: Seasonal, annual, growing season, custom date ranges
- **Statistical analysis**: Climatology, anomalies, trend analysis
- **Enhanced exports**: NetCDF, Cloud-Optimized GeoTIFF (COG), zonal statistics CSV
- **Visualization**: Time series plots, maps, climatology charts

### Effective Precipitation Methods

pyCropWat supports multiple methods for calculating effective precipitation:

| Method | Description |
|--------|-------------|
| `cropwat` | USDA SCS method as implemented in FAO CROPWAT (default) |
| `fao_aglw` | FAO/AGLW formula from FAO Irrigation Paper No. 33 |
| `fixed_percentage` | Simple fixed percentage method (configurable, default 70%) |
| `dependable_rainfall` | FAO Dependable Rainfall at specified probability level |

### CROPWAT Formula (Default)

The effective precipitation is calculated using the USDA SCS method as implemented in FAO CROPWAT (Smith, 1992; Muratoglu et al., 2023):

- If precipitation ≤ 250 mm: `Peff = P × (125 - 0.2 × P) / 125`
- If precipitation > 250 mm: `Peff = 0.1 × P + 125`

### FAO/AGLW Formula

The FAO Land and Water Division (AGLW) formula from FAO Irrigation and Drainage Paper No. 33:

- If precipitation ≤ 250 mm: `Peff = max(0.6 × P - 10, 0)`
- If precipitation > 250 mm: `Peff = 0.8 × P - 25`

### Fixed Percentage Method

A simple method assuming a constant fraction of precipitation is effective:

- `Peff = P × f` where `f` is the effectiveness fraction (default: 0.7 or 70%)

### Dependable Rainfall Method

The FAO Dependable Rainfall method estimates rainfall at a given probability level (default 75%):

- If precipitation < 100 mm: `Peff = max(0.6 × P - 10, 0)`
- If precipitation ≥ 100 mm: `Peff = 0.8 × P - 25`

A probability scaling factor is applied:
- 50% probability: ~1.2× base estimate (less conservative)
- 75% probability: 1.0× base estimate (default)
- 90% probability: ~0.8× base estimate (more conservative)

### Method Comparison

| Method | Use Case | Characteristics |
|--------|----------|-----------------|
| **CROPWAT** | General irrigation planning | Balanced, widely validated |
| **FAO/AGLW** | Yield response studies | Similar to CROPWAT, slightly different curve |
| **Fixed Percentage** | Quick estimates, calibration | Simple, requires local calibration |
| **Dependable Rainfall** | Risk-averse planning | Conservative, probability-based |

## Installation

### Using Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/pyCropWat.git
cd pyCropWat

# Create conda environment from environment.yml
conda env create -f environment.yml

# Activate the environment
conda activate pycropwat

# Install the package (registers the 'pycropwat' CLI command)
pip install -e .

# Or with interactive map support (leafmap, localtileserver)
pip install -e ".[interactive]"

# Verify installation
pycropwat --help
```

### Using pip

```bash
# Clone the repository
git clone https://github.com/yourusername/pyCropWat.git
cd pyCropWat

# Create and activate a virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install the package (registers the 'pycropwat' CLI command)
pip install -e .

# Or with interactive map support (leafmap, localtileserver)
pip install -e ".[interactive]"

# Verify installation
pycropwat --help
```

**Note:** After running `pip install -e .`, the `pycropwat` command will be available globally in your environment. Do not use `./pycropwat` - just use `pycropwat` directly.

**Optional Dependencies:**

- `pip install -e ".[interactive]"` - Adds leafmap and localtileserver for interactive HTML maps
- `pip install -e ".[dev]"` - Adds development tools (pytest, black, ruff)
- `pip install -e ".[docs]"` - Adds documentation tools (mkdocs)
```

## Requirements

- Python >= 3.9
- Google Earth Engine account and authentication
- Dependencies: earthengine-api, numpy, xarray, rioxarray, geopandas, shapely, dask

### Conda Environment

The `environment.yml` file provides a complete conda environment with all dependencies:

```bash
# Create environment
conda env create -f environment.yml

# Activate environment
conda activate pycropwat

# Update existing environment
conda env update -f environment.yml --prune

# Remove environment
conda env remove -n pycropwat
```

## Usage

### Python API

```python
from pycropwat import EffectivePrecipitation

# Initialize the processor with a local file
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='path/to/region.geojson',
    start_year=2015,
    end_year=2020,
    precip_scale_factor=1000,  # ERA5 precipitation is in meters, convert to mm
    gee_project='your-gee-project'  # Optional
)

# Or use a GEE FeatureCollection asset for the study area
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    gee_geometry_asset='projects/my-project/assets/study_boundary',
    start_year=2015,
    end_year=2020,
    precip_scale_factor=1000
)

# Use alternative effective precipitation methods
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    gee_geometry_asset='projects/my-project/assets/study_boundary',
    start_year=2015,
    end_year=2020,
    precip_scale_factor=1000,
    method='fao_aglw'  # Options: 'cropwat', 'fao_aglw', 'fixed_percentage', 'dependable_rainfall'
)

# Process with parallel execution (using dask)
results = ep.process(
    output_dir='./output',
    n_workers=4,
    months=[6, 7, 8]  # Optional: process only specific months
)

# Or process sequentially (useful for debugging)
results = ep.process_sequential(output_dir='./output')
```

### Temporal Aggregation & Analysis

```python
from pycropwat import TemporalAggregator, StatisticalAnalyzer, Visualizer

# Temporal aggregation
agg = TemporalAggregator('./output')

# Annual total
annual = agg.annual_aggregate(2020, method='sum', output_path='./annual_2020.tif')

# Seasonal aggregate (JJA = June-July-August)
summer = agg.seasonal_aggregate(2020, 'JJA', method='sum')

# Growing season (April-October)
growing = agg.growing_season_aggregate(2020, start_month=4, end_month=10)

# Multi-year climatology
climatology = agg.multi_year_climatology(2000, 2020, output_dir='./climatology')

# Statistical analysis
stats = StatisticalAnalyzer('./output')

# Calculate anomaly
anomaly = stats.calculate_anomaly(2020, 6, clim_start=1990, clim_end=2020, 
                                   anomaly_type='percent')

# Trend analysis (returns slope in mm/year and p-value)
slope, pvalue = stats.calculate_trend(start_year=2000, end_year=2020, month=6)

# Zonal statistics
zonal_df = stats.zonal_statistics('./zones.shp', 2000, 2020, output_path='./zonal_stats.csv')

# Visualization
viz = Visualizer('./output')
viz.plot_time_series(2000, 2020, output_path='./timeseries.png')
viz.plot_monthly_climatology(2000, 2020, output_path='./climatology.png')
viz.plot_raster(2020, 6, output_path='./map_2020_06.png')

# Interactive map (requires leafmap or folium: pip install leafmap)
viz.plot_interactive_map(2020, 6, output_path='./interactive_map.html')

# Dataset comparison
viz.plot_comparison(2020, 6, other_dir='./terraclimate_output', 
                    labels=('ERA5', 'TerraClimate'), output_path='./comparison.png')
viz.plot_scatter_comparison(2000, 2020, other_dir='./terraclimate_output',
                            labels=('ERA5', 'TerraClimate'), output_path='./scatter.png')
viz.plot_annual_comparison(2000, 2020, other_dir='./terraclimate_output',
                           labels=('ERA5', 'TerraClimate'), output_path='./annual_comparison.png')
```

### Export Options

```python
from pycropwat import export_to_netcdf, export_to_cog

# Export to NetCDF (single file with time dimension)
export_to_netcdf('./output', './effective_precip.nc')

# Convert to Cloud-Optimized GeoTIFF
export_to_cog('./output/effective_precip_2020_06.tif', './cog_2020_06.tif')
```

### Command Line Interface

pyCropWat provides a subcommand-based CLI for all functionality:

```bash
pycropwat <command> [OPTIONS]
```

**Available Commands:**

| Command | Description |
|---------|-------------|
| `process` | Calculate effective precipitation from GEE climate data |
| `aggregate` | Temporal aggregation (annual, seasonal, growing season) |
| `analyze` | Statistical analysis (anomaly, trend, zonal statistics) |
| `export` | Export to NetCDF or Cloud-Optimized GeoTIFF |
| `plot` | Create visualizations (time series, climatology, maps) |

#### Process Command Examples

```bash
# Process ERA5-Land data (actual working example)
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
                  --band total_precipitation_sum \
                  --gee-geometry projects/ssebop-471916/assets/Riodelaplata \
                  --start-year 2000 --end-year 2025 \
                  --scale-factor 1000 --scale 4000 \
                  --workers 32 --output ./Examples/RDP_ERA5Land

# Use alternative effective precipitation method
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
                  --band total_precipitation_sum \
                  --gee-geometry projects/my-project/assets/study_area \
                  --start-year 2020 --end-year 2023 \
                  --scale-factor 1000 \
                  --method fao_aglw --output ./outputs

# List available methods
pycropwat --list-methods

# Process TerraClimate data (actual working example)
pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE \
                  --band pr \
                  --gee-geometry projects/ssebop-471916/assets/Riodelaplata \
                  --start-year 2000 --end-year 2025 \
                  --workers 32 --output ./Examples/RDP_TerraClimate

# Process with local shapefile
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
                  --band total_precipitation_sum \
                  --geometry roi.geojson \
                  --start-year 2015 --end-year 2020 \
                  --scale-factor 1000 --output ./output
```

#### Aggregate Command Examples

```bash
# Annual total
pycropwat aggregate --input ./output --type annual --year 2020 --output ./annual_2020.tif

# Seasonal (summer)
pycropwat aggregate --input ./output --type seasonal --year 2020 --season JJA --output ./summer_2020.tif

# Growing season (April-October)
pycropwat aggregate --input ./output --type growing-season --year 2020 \
                    --start-month 4 --end-month 10 --output ./growing_2020.tif

# Multi-year climatology
pycropwat aggregate --input ./output --type climatology \
                    --start-year 2000 --end-year 2020 --output ./climatology/
```

#### Analyze Command Examples

```bash
# Calculate anomaly
pycropwat analyze anomaly --input ./output --year 2020 --month 6 \
                          --clim-start 1990 --clim-end 2020 --output ./anomaly_2020_06.tif

# Calculate trend
pycropwat analyze trend --input ./output --start-year 2000 --end-year 2020 \
                        --trend-method sen --output ./trend/

# Zonal statistics
pycropwat analyze zonal --input ./output --zones ./regions.shp \
                        --start-year 2000 --end-year 2020 --output ./zonal_stats.csv
```

#### Export Command Examples

```bash
# Export to NetCDF
pycropwat export netcdf --input ./output --output ./data.nc

# Convert to Cloud-Optimized GeoTIFF
pycropwat export cog --input ./effective_precip_2020_06.tif --output ./cog_2020_06.tif
```

#### Plot Command Examples

```bash
# Time series plot
pycropwat plot timeseries --input ./output --start-year 2000 --end-year 2020 --output ./timeseries.png

# Monthly climatology bar chart
pycropwat plot climatology --input ./output --start-year 2000 --end-year 2020 --output ./climatology.png

# Single month map
pycropwat plot map --input ./output --year 2020 --month 6 --output ./map_2020_06.png

# Interactive map (requires leafmap: pip install leafmap)
pycropwat plot interactive --input ./output --year 2020 --month 6 --output ./map.html

# Compare two datasets (e.g., ERA5 vs TerraClimate)
pycropwat plot compare --input ./era5_output --other-input ./terraclimate_output \
                       --year 2020 --month 6 --label1 ERA5 --label2 TerraClimate \
                       --output ./comparison.png

# Scatter plot for validation
pycropwat plot scatter --input ./era5_output --other-input ./terraclimate_output \
                       --start-year 2000 --end-year 2020 --output ./scatter.png

# Annual comparison bar chart
pycropwat plot annual-compare --input ./era5_output --other-input ./terraclimate_output \
                              --start-year 2000 --end-year 2020 --output ./annual.png
```

### CLI Arguments

#### Global Options

| Argument | Description |
|----------|-------------|
| `--help` | Show help message |
| `--version` | Show version number |
| `--list-methods` | List available effective precipitation methods |

#### Process Command Arguments

| Argument | Short | Required | Default | Description |
|----------|-------|----------|---------|-------------|
| `--asset` | `-a` | Yes | - | GEE ImageCollection asset ID |
| `--band` | `-b` | Yes | - | Precipitation band name |
| `--geometry` | `-g` | No* | - | Path to shapefile or GeoJSON |
| `--gee-geometry` | `-G` | No* | - | GEE FeatureCollection asset ID |
| `--start-year` | `-s` | Yes | - | Start year (inclusive) |
| `--end-year` | `-e` | Yes | - | End year (inclusive) |
| `--output` | `-o` | Yes | - | Output directory |
| `--scale-factor` | `-f` | No | 1.0 | Conversion factor to mm |
| `--scale` | `-r` | No | Native | Output resolution in meters |
| `--workers` | `-w` | No | 4 | Number of parallel workers |
| `--months` | `-m` | No | All | Specific months to process |
| `--project` | `-p` | No | None | GEE project ID |
| `--method` | - | No | cropwat | Peff method: cropwat, fao_aglw, fixed_percentage, dependable_rainfall |
| `--percentage` | - | No | 0.7 | Percentage for fixed_percentage method |
| `--probability` | - | No | 0.75 | Probability for dependable_rainfall method |
| `--sequential` | - | No | False | Process sequentially |
| `--verbose` | `-v` | No | False | Verbose output |

\* Either `--geometry` or `--gee-geometry` must be provided.

For full CLI documentation, run `pycropwat <command> --help` or see the [CLI Reference](https://montimaj.github.io/pyCropWat/user-guide/cli/).

## Output Files

The package generates two GeoTIFF files per month:

1. `effective_precip_YYYY_MM.tif` - Effective precipitation (mm)
2. `effective_precip_fraction_YYYY_MM.tif` - Effective precipitation fraction (0-1)

### Output Resolution

- **Default (no `--scale`):** Uses the native resolution of the input dataset
  - ERA5-Land: ~11 km (0.1°)
  - TerraClimate: ~4 km (1/24°)
  - CHIRPS: ~5.5 km (0.05°)
- **With `--scale`:** Reprojects to the specified resolution in meters (e.g., `--scale 1000` for 1 km)

### Large Region Handling

For large study areas or high-resolution outputs that exceed GEE's pixel limits, pyCropWat automatically:
1. Splits the region into smaller tiles (max 256×256 pixels per tile)
2. Downloads each tile separately
3. Mosaics the tiles back together using `rioxarray`

This is handled automatically with no additional dependencies required.

### Important: Units

The CROPWAT formula is calibrated for precipitation in **millimeters (mm)**. The output effective precipitation is always in mm, provided you use the correct `--scale-factor` to convert input precipitation to mm first.

The formula constants (125, 250, 0.2, 0.1) are specifically designed for mm units:
- If P ≤ 250mm: `Peff = P × (125 - 0.2 × P) / 125`
- If P > 250mm: `Peff = 0.1 × P + 125`

**Warning:** If you pass precipitation in wrong units (e.g., ERA5 in meters without `--scale-factor 1000`), the results will be incorrect because the 250mm threshold won't match properly.

### Temporal Aggregation

pyCropWat automatically **sums** all images within each month to compute monthly total precipitation, regardless of the input data's temporal resolution:

- **Monthly data (ERA5, TerraClimate):** Uses the single monthly image directly
- **Daily data (CHIRPS/DAILY):** Sums all ~30 daily images → monthly total
- **Sub-daily data (GPM IMERG):** Sums all timesteps → monthly total

This ensures the CROPWAT formula always receives the correct monthly precipitation totals.

## Common GEE Climate Assets

| Asset ID | Precipitation Band | Scale Factor | Spatial Resolution | Temporal Resolution |
|----------|-------------------|--------------|-------------------|---------------------|
| `ECMWF/ERA5_LAND/MONTHLY_AGGR` | `total_precipitation_sum` | 1000 | ~11 km (0.1°) | Monthly |
| `ECMWF/ERA5/MONTHLY` | `total_precipitation` | 1000 | ~27 km (0.25°) | Monthly |
| `IDAHO_EPSCOR/TERRACLIMATE` | `pr` | 1 | ~4 km (1/24°) | Monthly |
| `UCSB-CHG/CHIRPS/DAILY` | `precipitation` | 1 | ~5.5 km (0.05°) | Daily |
| `UCSB-CHG/CHIRPS/PENTAD` | `precipitation` | 1 | ~5.5 km (0.05°) | 5-day (Pentad) |
| `NASA/GPM_L3/IMERG_V06` | `precipitation` | 1 | ~11 km (0.1°) | Half-hourly |
| `projects/climate-engine-pro/assets/ce-ag-era5-v2/daily` | `Precipitation_Flux` | 1 | ~9 km (0.1°) | Daily |

## Complete Workflow Example

The `Examples/complete_workflow_example.py` script demonstrates the complete pyCropWat workflow using real data from the Rio de la Plata basin. This is an excellent starting point to understand all package capabilities.

📖 **For detailed step-by-step documentation, see [Examples/README.md](Examples/README.md)**

### What the Example Does

The script performs a comprehensive 6-step workflow:

1. **Process Effective Precipitation** - Downloads and calculates effective precipitation from ERA5-Land and TerraClimate via GEE
2. **Temporal Aggregation** - Creates annual totals, growing season aggregations (Apr-Sep), and monthly climatology
3. **Statistical Analysis** - Computes percent anomalies, trends (Sen's slope), and zonal statistics
4. **Visualization** - Generates time series plots, climatology charts, static maps, and interactive HTML maps
5. **Dataset Comparison** - Creates side-by-side comparison plots, scatter plots, annual charts, and zonal comparisons
6. **NetCDF Export** - Exports data to CF-compliant NetCDF format

### Running the Example

```bash
# Navigate to the Examples directory
cd Examples/

# Run analysis only (using existing pre-processed data)
python complete_workflow_example.py --analysis-only

# Run full workflow with GEE processing (requires authentication)
python complete_workflow_example.py --gee-project your-project-id

# Force reprocess all data from GEE
python complete_workflow_example.py --force-reprocess --gee-project your-project-id
```

### Configuration

The script uses these real-world settings:

| Parameter | Value |
|-----------|-------|
| Study Area | Rio de la Plata Basin (GEE Asset: `projects/ssebop-471916/assets/Riodelaplata`) |
| Time Period | 2000-2023 |
| Climatology Period | 2000-2020 |
| Datasets | ERA5-Land, TerraClimate |
| Sample Zones | Eastern RDP (Uruguay, SE Brazil), Western RDP (N Argentina, Paraguay) |

### Output Structure

Running the example creates organized outputs:

```
Examples/
├── complete_workflow_example.py    # The workflow script
├── README.md                       # Detailed documentation
├── RDP_ERA5Land/                   # Monthly effective precipitation (ERA5-Land)
│   ├── effective_precip_YYYY_MM.tif
│   └── effective_precip_fraction_YYYY_MM.tif
├── RDP_TerraClimate/               # Monthly effective precipitation (TerraClimate)
│   ├── effective_precip_YYYY_MM.tif
│   └── effective_precip_fraction_YYYY_MM.tif
└── analysis_outputs/
    ├── annual/                     # Annual totals by dataset
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── climatology/                # Monthly climatology (2000-2020)
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── growing_season/             # Apr-Sep aggregations
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── anomalies/                  # Percent anomalies (2021-2023)
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── trend/                      # Trend analysis (Sen's slope)
    │   ├── ERA5Land/
    │   └── TerraClimate/
    ├── zonal/                      # Zonal statistics
    │   ├── sample_zones.geojson    # Eastern & Western RDP zones
    │   ├── ERA5Land/
    │   │   ├── zonal_stats.csv
    │   │   └── zonal_summary.png
    │   └── TerraClimate/
    │       ├── zonal_stats.csv
    │       └── zonal_summary.png
    ├── figures/                    # Visualizations by dataset
    │   ├── ERA5Land/
    │   │   ├── time_series.png
    │   │   ├── monthly_climatology.png
    │   │   ├── map_2023_01.png         # Summer/Wet season
    │   │   ├── map_2023_06.png         # Winter/Dry season
    │   │   ├── map_2023_10.png         # Spring
    │   │   ├── interactive_map_*.html  # Interactive maps
    │   │   └── map_notable_*.png       # El Niño/La Niña events
    │   └── TerraClimate/
    ├── comparisons/                # ERA5-Land vs TerraClimate
    │   ├── comparison_2023_06.png      # Side-by-side maps
    │   ├── scatter_comparison.png      # Scatter plot
    │   ├── annual_comparison.png       # Annual bar chart
    │   └── zonal_comparison.png        # Zonal statistics comparison
    ├── ERA5Land_effective_precip_2000_2023.nc
    └── TerraClimate_effective_precip_2000_2023.nc
```

For more examples and detailed API usage, see the [Examples documentation](https://montimaj.github.io/pyCropWat/examples/).

## Google Earth Engine Authentication

Before using pyCropWat, authenticate with Google Earth Engine:

```bash
earthengine authenticate
```

Or in Python:

```python
import ee
ee.Authenticate()
ee.Initialize(project='your-project-id')
```

## Citation

If you use pyCropWat in your research, please cite:

```bibtex
@software{pycropwat,
  author = {Majumdar, Sayantan and ReVelle, Peter and Nozari, Soheil and Huntington, Justin and Smith, Ryan},
  title = {pyCropWat: Python implementation of FAO CROPWAT effective precipitation using Google Earth Engine},
  year = {2026},
  url = {https://github.com/montimaj/pyCropWat}
}
```

### Effective Precipitation Method References

- Muratoglu, A., Bilgen, G. K., Angin, I., & Kodal, S. (2023). Performance analyses of effective rainfall estimation methods for accurate quantification of agricultural water footprint. *Water Research*, *238*, 120011. https://doi.org/10.1016/j.watres.2023.120011

- Smith, M. (1992). *CROPWAT: A computer program for irrigation planning and management* (FAO Irrigation and Drainage Paper No. 46). Food and Agriculture Organization of the United Nations. https://www.fao.org/sustainable-development-goals-helpdesk/champion/article-detail/cropwat/en


## Funding

This work was supported by a U.S. Army Corps of Engineers grant (W912HZ25C0016) for the project *"Improved Characterization of Groundwater Resources in Transboundary Watersheds using Satellite Data and Integrated Models."*

**Principal Investigator:** Dr. Ryan Smith (Colorado State University)

**Co-Principal Investigators:**
- Dr. Ryan Bailey (Colorado State University)
- Dr. Justin Huntington (Desert Research Institute)
- Dr. Sayantan Majumdar (Desert Research Institute)
- Mr. Peter ReVelle (Desert Research Institute)

**Research Scientist:**
- Dr. Soheil Nozari (Colorado State University)

## Roadmap

The following features have been implemented or are under consideration for future releases:

### ✅ Implemented in v1.0.0

#### 📊 Temporal Aggregation
- ✅ Seasonal summaries (DJF, MAM, JJA, SON)
- ✅ Annual totals with statistics
- ✅ Growing season aggregations based on crop calendars
- ✅ Custom date range aggregations

#### 📈 Statistical Analysis
- ✅ Long-term climatology (e.g., 30-year normals)
- ✅ Anomaly detection (absolute, percent, standardized)
- ✅ Trend analysis using linear regression or Theil-Sen/Mann-Kendall
- ✅ Zonal statistics for polygons

#### 🔄 Additional Effective Precipitation Methods
- ✅ FAO/AGLW method
- ✅ Fixed percentage method (configurable percentage)
- ✅ Dependable rainfall (FAO method at different probability levels)

#### 📤 Enhanced Export Options
- ✅ NetCDF output for time-series analysis
- ✅ Cloud-Optimized GeoTIFFs (COGs)
- ✅ Zonal statistics CSV export by polygon

#### 📉 Visualization
- ✅ Built-in plotting functions for time series
- ✅ Monthly climatology bar charts
- ✅ Raster map visualization
- ✅ Interactive maps (folium/leafmap integration)
- ✅ Comparison plots between datasets (side-by-side, scatter, annual)

### 🚧 Planned for Future Releases

#### 🌾 Crop Water Requirements
- Crop coefficient (Kc) integration for different crop types/stages
- Net irrigation requirement (ETc - Peff) calculations
- Crop calendar support for region-specific growing seasons

#### 📈 Advanced Analysis
- Drought indices (SPI, SPEI) integration
- Direct cloud storage export (GCS, S3)

#### ✅ Validation Tools
- Station data comparison module
- Cross-dataset validation (e.g., ERA5 vs TerraClimate)
- Uncertainty quantification

#### 💧 Water Balance Extension
- Evapotranspiration integration (from MODIS, SSEBop, OpenET)
- Simple water balance (P - ET - Runoff)

**Have a feature request?** Please submit your ideas via [GitHub Issues](https://github.com/montimaj/pyCropWat/issues). We welcome community contributions and feedback!

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

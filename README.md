# pyCropWat

A Python package for calculating CROPWAT effective precipitation from Google Earth Engine climate data.

## Overview

pyCropWat converts precipitation data from any GEE climate dataset into effective precipitation and effective precipitation fraction rasters using the CROPWAT methodology. It supports:

- Any GEE ImageCollection with precipitation data
- Shapefile, GeoJSON, or GEE FeatureCollection asset for region of interest
- Parallel processing using Dask
- Monthly output rasters in GeoTIFF format

### CROPWAT Effective Precipitation Formula

- If precipitation ≤ 250mm: `Peff = P × (125 - 0.2 × P) / 125`
- If precipitation > 250mm: `Peff = 0.1 × P + 125`

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

# Verify installation
pycropwat --help
```

**Note:** After running `pip install -e .`, the `pycropwat` command will be available globally in your environment. Do not use `./pycropwat` - just use `pycropwat` directly.

# Or install dependencies directly
pip install -r requirements.txt
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

# Process with parallel execution (using dask)
results = ep.process(
    output_dir='./output',
    n_workers=4,
    months=[6, 7, 8]  # Optional: process only specific months
)

# Or process sequentially (useful for debugging)
results = ep.process_sequential(output_dir='./output')
```

### Command Line Interface

```bash
# Process ERA5-Land data with a GEE FeatureCollection asset (actual working example)
pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
          --band total_precipitation_sum \
          --gee-geometry projects/ssebop-471916/assets/Riodelaplata \
          --start-year 2025 \
          --end-year 2025 \
          --months 6 7 8 \
          --scale-factor 1000 \
          --output ./Outputs

# Process TerraClimate data with a GEE FeatureCollection asset (actual working example)
# Note: TerraClimate precipitation is already in mm, so no scale-factor needed
pycropwat --asset IDAHO_EPSCOR/TERRACLIMATE \
          --band pr \
          --gee-geometry projects/ssebop-471916/assets/Riodelaplata \
          --start-year 2000 \
          --end-year 2025 \
          --output ./Outputs_TerraClimate

# Process ERA5-Land data with a local shapefile
pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
          --band total_precipitation_sum \
          --geometry roi.geojson \
          --start-year 2015 \
          --end-year 2020 \
          --output ./output \
          --scale-factor 1000 \
          --workers 4

# Process specific months only
pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
          --band total_precipitation_sum \
          --geometry roi.geojson \
          --start-year 2020 \
          --end-year 2020 \
          --output ./output \
          --scale-factor 1000 \
          --months 6 7 8
```

### CLI Arguments

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
| `--scale` | `-r` | No | 11132 | Output resolution (meters) |
| `--workers` | `-w` | No | 4 | Number of parallel workers |
| `--months` | `-m` | No | All | Specific months to process |
| `--project` | `-p` | No | None | GEE project ID |
| `--sequential` | - | No | False | Process sequentially |
| `--verbose` | `-v` | No | False | Verbose output |

\* Either `--geometry` or `--gee-geometry` must be provided.

## Output Files

The package generates two GeoTIFF files per month:

1. `effective_precip_YYYY_MM.tif` - Effective precipitation (mm)
2. `effective_precip_fraction_YYYY_MM.tif` - Effective precipitation fraction (0-1)

### Important: Units

The CROPWAT formula is calibrated for precipitation in **millimeters (mm)**. The output effective precipitation is always in mm, provided you use the correct `--scale-factor` to convert input precipitation to mm first.

The formula constants (125, 250, 0.2, 0.1) are specifically designed for mm units:
- If P ≤ 250mm: `Peff = P × (125 - 0.2 × P) / 125`
- If P > 250mm: `Peff = 0.1 × P + 125`

**Warning:** If you pass precipitation in wrong units (e.g., ERA5 in meters without `--scale-factor 1000`), the results will be incorrect because the 250mm threshold won't match properly.

## Common GEE Climate Assets

| Asset ID | Precipitation Band | Scale Factor | Notes |
|----------|-------------------|--------------|-------|
| `ECMWF/ERA5_LAND/MONTHLY_AGGR` | `total_precipitation_sum` | 1000 | Monthly, meters to mm |
| `ECMWF/ERA5/MONTHLY` | `total_precipitation` | 1000 | Monthly, meters to mm |
| `IDAHO_EPSCOR/TERRACLIMATE` | `pr` | 1 | Monthly, already in mm |
| `UCSB-CHG/CHIRPS/DAILY` | `precipitation` | 1 | Daily, already in mm |
| `UCSB-CHG/CHIRPS/PENTAD` | `precipitation` | 1 | 5-day, already in mm |
| `NASA/GPM_L3/IMERG_V06` | `precipitation` | 1 | Half-hourly, mm/hr |

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
Python and Google Earth Engine Python API-based implementation of FAO CropWat

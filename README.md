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
│   ├── user-guide/          # User guide
│   │   ├── api.md
│   │   ├── cli.md
│   │   └── quickstart.md
├── Examples/                # Example outputs (see note below)
│   ├── RDP_ERA5Land.zip     # ERA5-Land outputs for Rio de la Plata 
│   └── RDP_TerraClimate.zip # TerraClimate outputs for Rio de la Plata
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

**Note:** The `Examples/` folder contains sample output rasters generated from actual runs. These can be used as reference to verify your outputs or for testing visualization workflows. Currently, this folder has monthly TerraClimate and ERA5-Land-based effective precipitation and effective precipitation fraction TIF files (as zips; 2000-2025) for the Rio de la Plata (RDP) region of South America.

**Changelog:** See [CHANGELOG.md](CHANGELOG.md) for release notes and version history.

## Overview

pyCropWat converts precipitation data from any GEE climate dataset into effective precipitation and effective precipitation fraction rasters using the CROPWAT methodology. It supports:

- Any GEE ImageCollection with precipitation data
- Shapefile, GeoJSON, or GEE FeatureCollection asset for region of interest
- Parallel processing using Dask
- Monthly output rasters in GeoTIFF format

### CROPWAT Effective Precipitation Formula

The effective precipitation is calculated using the USDA SCS method as implemented in FAO CROPWAT (Smith, 1992; Muratoglu et al., 2023):

- If precipitation ≤ 250 mm: `Peff = P × (125 - 0.2 × P) / 125`
- If precipitation > 250 mm: `Peff = 0.1 × P + 125`

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
# Process ERA5-Land data with a GEE FeatureCollection asset and output in 4 km resolution instead of the native 11 km resolution (actual working example)
pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
          --band total_precipitation_sum \
          --gee-geometry projects/ssebop-471916/assets/Riodelaplata \
          --start-year 2025 \
          --end-year 2025 \
          --months 6 7 8 \
          --scale-factor 1000 \
          --scale 4000 \
          --output ./Outputs

# Process TerraClimate data with a GEE FeatureCollection asset and output in the native resolution (actual working example)
# Note: TerraClimate precipitation is already in mm, so no scale-factor needed
pycropwat --asset IDAHO_EPSCOR/TERRACLIMATE \
          --band pr \
          --gee-geometry projects/ssebop-471916/assets/Riodelaplata \
          --start-year 2000 \
          --end-year 2025 \
          --workers 32 \
          --output ./Examples/RDP_TerraClimate

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
| `--scale` | `-r` | No | Native | Output resolution in meters |
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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

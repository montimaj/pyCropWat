# Installation

## Requirements

- Python 3.9 or higher
- Google Earth Engine account with authenticated access

## Install with pip

```bash
pip install pycropwat
```

## Install from source

```bash
git clone https://github.com/sayantan-majumdar/pyCropWat.git
cd pyCropWat
pip install -e .
```

## Install with conda

```bash
# Create environment from file
conda env create -f environment.yml
conda activate pycropwat

# Install package
pip install -e .
```

## Dependencies

The following packages are automatically installed:

| Package | Version | Description |
|---------|---------|-------------|
| earthengine-api | ≥0.1.370 | Google Earth Engine Python API |
| numpy | ≥1.21.0 | Numerical computing |
| xarray | ≥2022.3.0 | N-D labeled arrays |
| rioxarray | ≥0.14.0 | Rasterio xarray extension |
| geopandas | ≥0.12.0 | Geospatial pandas |
| shapely | ≥2.0.0 | Geometric operations |
| dask | ≥2022.1.0 | Parallel computing |
| distributed | ≥2022.1.0 | Dask distributed scheduler |
| rasterio | ≥1.3.0 | Raster I/O |

## Google Earth Engine Authentication

Before using pyCropWat, you must authenticate with Google Earth Engine:

```bash
earthengine authenticate
```

This will open a browser window to complete the OAuth flow. For server environments, use:

```bash
earthengine authenticate --auth_mode=notebook
```

### Using a GEE Project

If you have a specific GEE project, you can specify it:

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2020,
    end_year=2023,
    gee_project='my-gee-project'  # Your GEE project ID
)
```

Or via CLI:

```bash
export GOOGLE_CLOUD_PROJECT=my-gee-project
pycropwat --asset ... --gee-geometry ...
```

## Verifying Installation

```python
import pycropwat
print(pycropwat.__version__)

# Test GEE connection
import ee
ee.Initialize()
print("GEE connection successful!")
```

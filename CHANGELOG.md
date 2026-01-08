# Changelog

All notable changes to pyCropWat will be documented in this file.

## [1.0.0] - 2026-01-08

### 🎉 Initial Release

pyCropWat is a Python package for calculating CROPWAT effective precipitation from Google Earth Engine (GEE) climate data using the USDA SCS method as implemented in FAO CROPWAT.

### ✨ Features

- **Multiple GEE Dataset Support**: Works with ERA5-Land, TerraClimate, CHIRPS, GPM IMERG, AgERA5, and any GEE dataset with precipitation bands
- **Flexible Geometry Inputs**: Accepts Shapefiles, GeoJSON, or GEE FeatureCollection assets
- **Automatic Chunked Downloads**: Handles large regions exceeding GEE pixel limits with intelligent tiling (256×256 pixels)
- **Parallel Processing**: Leverages Dask for efficient multi-month/multi-year processing
- **Dual Interface**: Full-featured CLI and Python API for seamless integration
- **Output Products**:
  - Monthly effective precipitation rasters (GeoTIFF)
  - Effective precipitation fraction (ratio of effective to total)

### 📦 Installation

```bash
pip install pycropwat
```

Or with conda:

```bash
conda env create -f environment.yml
conda activate pycropwat
```

### 🚀 Quick Start

```bash
pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
          --band total_precipitation_sum \
          --gee-geometry projects/my-project/assets/study_area \
          --start-year 2020 --end-year 2023 \
          --scale-factor 1000 --output ./outputs
```

### 📚 Documentation

Full documentation available at: https://montimaj.github.io/pyCropWat

### 👥 Contributors

- Sayantan Majumdar (Desert Research Institute)
- Peter ReVelle (Desert Research Institute)
- Soheil Nozari (Colorado State University)
- Justin Huntington (Desert Research Institute)
- Ryan Smith (Colorado State University)

### 🙏 Acknowledgments

This work was supported by the U.S. Army Corps of Engineers (Grant W912HZ25C0016) for the project *"Improved Characterization of Groundwater Resources in Transboundary Watersheds using Satellite Data and Integrated Models."*

### 📄 References

- Smith, M. (1992). *CROPWAT: A computer program for irrigation planning and management* (FAO Irrigation and Drainage Paper No. 46). Food and Agriculture Organization of the United Nations. https://www.fao.org/3/t7202e/t7202e00.htm
- Muratoglu, A., Bilgen, G. K., Angin, I., & Kodal, S. (2023). Performance analyses of effective rainfall estimation methods for accurate quantification of agricultural water footprint. *Water Research*, *238*, 120011. https://doi.org/10.1016/j.watres.2023.120011

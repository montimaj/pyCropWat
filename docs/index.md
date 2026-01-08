# pyCropWat

**Calculate CROPWAT effective precipitation from Google Earth Engine climate data**

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://montimaj.github.io/pyCropWat)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

## Overview

pyCropWat is a Python package that calculates effective precipitation using the CROPWAT methodology from any Google Earth Engine (GEE) climate dataset. It supports:

- **Multiple GEE datasets**: ERA5-Land, TerraClimate, CHIRPS, GPM IMERG, and more
- **Flexible geometry inputs**: Shapefiles, GeoJSON, or GEE FeatureCollection assets
- **Automatic chunked downloads**: Handles large regions that exceed GEE pixel limits
- **Parallel processing**: Uses Dask for efficient multi-month processing
- **CLI and Python API**: Use from command line or integrate into your workflows

## CROPWAT Formula

The effective precipitation is calculated using the USDA SCS method as implemented in FAO CROPWAT ([Smith, 1992](https://www.fao.org/sustainable-development-goals-helpdesk/champion/article-detail/cropwat/en); [Muratoglu et al., 2023](https://doi.org/10.1016/j.watres.2023.120011)):

$$
P_{eff} = \begin{cases}
P \times \frac{125 - 0.2P}{125} & \text{if } P \leq 250 \text{ mm} \\
0.1P + 125 & \text{if } P > 250 \text{ mm}
\end{cases}
$$

## Quick Start

```bash
# Install
pip install pycropwat

# Run from CLI
pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
          --band total_precipitation_sum \
          --gee-geometry projects/my-project/assets/study_area \
          --start-year 2020 \
          --end-year 2023 \
          --scale-factor 1000 \
          --output ./outputs
```

## Features

- 📊 **Monthly effective precipitation rasters** (GeoTIFF)
- 📈 **Effective precipitation fraction** (ratio of effective to total)
- 🌍 **Any GEE climate dataset** with precipitation band
- 🗺️ **Flexible resolution control** (native or custom scale)
- ⚡ **Parallel processing** with Dask
- 🧩 **Automatic tiling** for large regions

## Documentation

For full documentation, visit [https://montimaj.github.io/pyCropWat](https://montimaj.github.io/pyCropWat)

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

MIT License - see [LICENSE](LICENSE) for details.

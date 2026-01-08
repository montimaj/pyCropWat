# pyCropWat

**Calculate CROPWAT effective precipitation from Google Earth Engine climate data**

[![Release](https://img.shields.io/badge/release-v1.0.0-green.svg)](https://github.com/montimaj/pyCropWat/releases)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://montimaj.github.io/pyCropWat)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

## Overview

pyCropWat is a Python package that calculates effective precipitation using multiple methodologies from any Google Earth Engine (GEE) climate dataset. It supports:

- **Multiple GEE datasets**: ERA5-Land, TerraClimate, CHIRPS, GPM IMERG, and more
- **Multiple Peff methods**: CROPWAT (USDA SCS), FAO/AGLW, Fixed Percentage, Dependable Rainfall
- **Flexible geometry inputs**: Shapefiles, GeoJSON, or GEE FeatureCollection assets
- **Automatic chunked downloads**: Handles large regions that exceed GEE pixel limits
- **Temporal aggregation**: Annual, seasonal, growing season, climatology
- **Statistical analysis**: Anomaly detection, trend analysis, zonal statistics
- **Visualization**: Time series, maps, interactive HTML maps, dataset comparison
- **Export options**: NetCDF, Cloud-Optimized GeoTIFFs
- **Parallel processing**: Uses Dask for efficient multi-month processing
- **CLI and Python API**: Use from command line or integrate into your workflows

## CROPWAT Formula

The default effective precipitation is calculated using the USDA SCS method as implemented in FAO CROPWAT ([Smith, 1992](https://www.fao.org/sustainable-development-goals-helpdesk/champion/article-detail/cropwat/en); [Muratoglu et al., 2023](https://doi.org/10.1016/j.watres.2023.120011)):

$$
P_{eff} = \begin{cases}
P \times \frac{125 - 0.2P}{125} & \text{if } P \leq 250 \text{ mm} \\
0.1P + 125 & \text{if } P > 250 \text{ mm}
\end{cases}
$$

## Effective Precipitation Methods

pyCropWat supports four different methods for calculating effective precipitation:

### 1. CROPWAT (USDA SCS) - Default

The USDA Soil Conservation Service method as implemented in FAO CROPWAT software. This is the most widely used method for irrigation planning.

| Condition | Formula |
|-----------|---------|
| P ≤ 250 mm | $P_{eff} = P \times \frac{125 - 0.2P}{125}$ |
| P > 250 mm | $P_{eff} = 0.1P + 125$ |

**Usage:**
```python
ep = EffectivePrecipitation(..., method='cropwat')
```

**Reference:** Smith, M. (1992). *CROPWAT: A computer program for irrigation planning and management*. FAO Irrigation and Drainage Paper No. 46.

---

### 2. FAO/AGLW

The FAO Land and Water Division (AGLW) formula from FAO Irrigation and Drainage Paper No. 33.

| Condition | Formula |
|-----------|---------|
| P ≤ 250 mm | $P_{eff} = \max(0.6P - 10, 0)$ |
| P > 250 mm | $P_{eff} = 0.8P - 25$ |

**Usage:**
```python
ep = EffectivePrecipitation(..., method='fao_aglw')
```

**Reference:** FAO. (1986). *Yield response to water*. FAO Irrigation and Drainage Paper No. 33.

---

### 3. Fixed Percentage

A simple method that assumes a constant fraction of precipitation is effective. Common values range from 70-80%.

$$P_{eff} = P \times f$$

Where $f$ is the effectiveness fraction (default: 0.7 or 70%).

**Usage:**
```python
ep = EffectivePrecipitation(..., method='fixed_percentage', method_params={'percentage': 0.7})
```

---

### 4. Dependable Rainfall

The FAO Dependable Rainfall method estimates the amount of rainfall that can be depended upon at a given probability level. More conservative than other methods.

| Condition | Formula (at 75% probability) |
|-----------|------------------------------|
| P < 100 mm | $P_{eff} = \max(0.6P - 10, 0)$ |
| P ≥ 100 mm | $P_{eff} = 0.8P - 25$ |

A probability scaling factor is applied for other probability levels:
- 50% probability: ~1.2× base estimate
- 75% probability: 1.0× base estimate (default)
- 90% probability: ~0.8× base estimate

**Usage:**
```python
ep = EffectivePrecipitation(..., method='dependable_rainfall', method_params={'probability': 0.75})
```

**Reference:** FAO. (1992). *CROPWAT - A computer program for irrigation planning and management*. FAO Irrigation and Drainage Paper No. 46.

---

### Method Comparison

| Method | Use Case | Characteristics |
|--------|----------|-----------------|
| **CROPWAT** | General irrigation planning | Balanced, widely validated |
| **FAO/AGLW** | Yield response studies | Similar to CROPWAT, slightly different curve |
| **Fixed Percentage** | Quick estimates, calibration | Simple, requires local calibration |
| **Dependable Rainfall** | Risk-averse planning | Conservative, probability-based |

!!! tip "Choosing a Method"
    - Use **CROPWAT** (default) for most irrigation planning applications
    - Use **FAO/AGLW** when following FAO Irrigation Paper No. 33 guidelines
    - Use **Fixed Percentage** when you have locally calibrated effectiveness values
    - Use **Dependable Rainfall** for drought-sensitive crops or risk-averse planning

## Quick Start

### CLI

```bash
# Install (basic)
pip install pycropwat

# Or with interactive map support
pip install pycropwat[interactive]

# Process effective precipitation
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --gee-geometry projects/my-project/assets/study_area \
    --start-year 2020 --end-year 2023 \
    --scale-factor 1000 \
    --output ./outputs

# Create annual aggregation
pycropwat aggregate --input ./outputs --type annual --year 2020 \
    --output ./annual_2020.tif

# Generate time series plot
pycropwat plot timeseries --input ./outputs \
    --start-year 2020 --end-year 2023 --output ./timeseries.png
```

### Python

```python
from pycropwat import EffectivePrecipitation
from pycropwat.analysis import TemporalAggregator, Visualizer

# Process effective precipitation
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    gee_geometry_asset='projects/my-project/assets/study_area',
    start_year=2020,
    end_year=2023,
    precip_scale_factor=1000
)
ep.process(output_dir='./outputs', n_workers=4)

# Create annual aggregation
agg = TemporalAggregator('./outputs')
agg.annual_aggregate(2020, output_path='./annual_2020.tif')

# Generate time series plot
viz = Visualizer('./outputs')
viz.plot_time_series(2020, 2023, output_path='./timeseries.png')
```

## Features

### ✅ Implemented

- 📊 **Multiple Peff methods**: CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall
- 🗓️ **Temporal aggregation**: Annual, seasonal, growing season, climatology
- 📈 **Statistical analysis**: Anomaly detection, trend analysis (linear, Theil-Sen), zonal statistics
- 📉 **Visualization**: Time series, climatology charts, maps, interactive HTML maps, dataset comparison (side-by-side, scatter, annual)
- 📤 **Export options**: NetCDF with time dimension, Cloud-Optimized GeoTIFFs
- 🌍 **Any GEE climate dataset** with precipitation band
- 🗺️ **Flexible resolution control** (native or custom scale)
- ⚡ **Parallel processing** with Dask
- 🧩 **Automatic tiling** for large regions

### 🚧 Planned

- 🌾 **Crop water requirements**: Kc integration, net irrigation requirement
- 📈 **Advanced analysis**: Drought indices (SPI, SPEI), direct cloud export
- ✅ **Validation tools**: Station comparison, uncertainty quantification
- 💧 **Water balance**: ET integration, simple water balance

## Example Outputs

The following visualizations are generated by the [complete workflow example](examples.md#example-12-complete-workflow) using real Rio de la Plata basin data:

### Time Series & Climatology

<p align="center">
  <img src="assets/examples/figures/ERA5Land/time_series.png" width="48%" alt="Time Series">
  <img src="assets/examples/figures/ERA5Land/monthly_climatology.png" width="48%" alt="Monthly Climatology">
</p>

*Left: Monthly effective precipitation time series (2000-2025). Right: Monthly climatology showing seasonal patterns.*

### Spatial Maps

<p align="center">
  <img src="assets/examples/figures/ERA5Land/map_2023_06.png" width="32%" alt="Winter Map">
  <img src="assets/examples/figures/ERA5Land/map_2023_01.png" width="32%" alt="Summer Map">
  <img src="assets/examples/figures/ERA5Land/map_notable_2015_12.png" width="32%" alt="El Niño Event">
</p>

*Left: Winter dry season (June 2023). Center: Summer wet season (January 2023). Right: El Niño event (December 2015).*

### Dataset Comparison (ERA5-Land vs TerraClimate)

<p align="center">
  <img src="assets/examples/comparisons/comparison_2023_06.png" width="100%" alt="Side-by-side Comparison">
</p>

*Side-by-side comparison of ERA5-Land and TerraClimate effective precipitation with difference map.*

<p align="center">
  <img src="assets/examples/comparisons/scatter_comparison.png" width="48%" alt="Scatter Plot">
  <img src="assets/examples/comparisons/annual_comparison.png" width="48%" alt="Annual Comparison">
</p>

*Left: Scatter plot comparison with R², RMSE, and bias statistics. Right: Annual totals comparison.*

### Method Comparison

<p align="center">
  <img src="assets/examples/method_comparison/ERA5Land_method_maps_2020_01.png" width="100%" alt="Method Comparison Maps">
</p>

*Comparison of effective precipitation methods: CROPWAT, FAO/AGLW, Fixed Percentage (70%), and Dependable Rainfall (75%).*

<p align="center">
  <img src="assets/examples/method_comparison/ERA5Land_method_curves.png" width="60%" alt="Method Curves">
</p>

*Theoretical response curves for different effective precipitation methods.*

## Documentation

For full documentation, visit [https://montimaj.github.io/pyCropWat](https://montimaj.github.io/pyCropWat)

- [Quick Start Guide](user-guide/quickstart.md)
- [CLI Reference](user-guide/cli.md)
- [Python API](user-guide/api.md)
- [Examples](examples.md)
- [Complete Workflow Example](examples.md#example-12-complete-workflow) - A comprehensive script demonstrating all features

!!! tip "Try the Complete Workflow Example"
    The `Examples/complete_workflow_example.py` script provides a ready-to-run demonstration of all pyCropWat features using real Rio de la Plata basin data:
    
    ```bash
    # Run with existing sample data
    python Examples/complete_workflow_example.py --analysis-only
    ```

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

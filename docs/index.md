# A Python Package for Computing Effective Precipitation Using Google Earth Engine Climate Data
<p align="left">
  <img src="assets/pyCropWat_logo.png" alt="pyCropWat Logo" width="200"><br>
</p>

[![Release](https://img.shields.io/badge/release-v1.3.0-green.svg)](https://github.com/montimaj/pyCropWat/releases)
[![PyPI](https://img.shields.io/pypi/v/pycropwat.svg)](https://pypi.org/project/pycropwat/)
[![Downloads](https://static.pepy.tech/badge/pycropwat/month)](https://pepy.tech/project/pycropwat)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18201619.svg)](https://doi.org/10.5281/zenodo.18201619)
[![GitHub stars](https://img.shields.io/github/stars/montimaj/pyCropWat)](https://github.com/montimaj/pyCropWat/stargazers)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://montimaj.github.io/pyCropWat)
[![GEE](https://img.shields.io/badge/Google%20Earth%20Engine-4285F4?logo=google-earth&logoColor=white)](https://earthengine.google.com/)

<p align="center">
  <img src="assets/pyCropWat.gif" alt="pyCropWat Demo">
</p>

## Quick Install

```bash
pip install pycropwat
```

## Overview

pyCropWat is a Python package that calculates effective precipitation using multiple methodologies from any Google Earth Engine (GEE) climate dataset. It supports:

- **Multiple GEE datasets**: Any ImageCollection from the [GEE Data Catalog](https://developers.google.com/earth-engine/datasets) or [Community Catalog](https://gee-community-catalog.org/) (e.g., ERA5-Land, TerraClimate, CHIRPS, GPM IMERG)
- **Local precipitation input**: Bring your own monthly GeoTIFFs or NetCDF files - AWC, ETo and geometry still come from Earth Engine, and the precipitation-only methods run with **no Earth Engine at all** (see [Local Precipitation Input](#local-precipitation-input))
- **Multiple Peff methods**: CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall, FarmWest, USDA-SCS, TAGEM-SuET, PCML, Ensemble
- **Flexible geometry inputs**: Shapefiles, GeoJSON, or GEE FeatureCollection assets
- **Automatic chunked downloads**: Handles large regions that exceed GEE pixel limits
- **Temporal aggregation**: Annual, seasonal, growing season (with cross-year support for Southern Hemisphere), climatology
- **Statistical analysis**: Anomaly detection, trend analysis, zonal statistics
- **Visualization**: Time series, maps, interactive HTML maps, dataset comparison
- **Export options**: NetCDF, Cloud-Optimized GeoTIFFs
- **Parallel processing**: Uses Dask for efficient multi-month processing
- **CLI and Python API**: Use from command line or integrate into your workflows

## CROPWAT Formula

The CROPWAT method calculates effective precipitation ([Smith, 1992](https://www.fao.org/sustainable-development-goals-helpdesk/champion/article-detail/cropwat/en); [Muratoglu et al., 2023](https://doi.org/10.1016/j.watres.2023.120011)):

$$
P_{eff} = \begin{cases}
P \times \frac{125 - 0.2P}{125} & \text{if } P \leq 250 \text{ mm} \\
0.1P + 125 & \text{if } P > 250 \text{ mm}
\end{cases}
$$

## Effective Precipitation Methods

pyCropWat supports nine different methods for calculating effective precipitation:

### 1. CROPWAT

The method used in FAO CROPWAT software. This is the most widely used method for irrigation planning.

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

### 2. FAO/AGLW (Dependable Rainfall)

The FAO Land and Water Division (AGLW) Dependable Rainfall formula, based on 80% probability exceedance, from FAO Irrigation and Drainage Paper No. 33.

| Condition | Formula |
|-----------|---------|
| P ≤ 70 mm | $P_{eff} = \max(0.6P - 10, 0)$ |
| P > 70 mm | $P_{eff} = 0.8P - 24$ |

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

The FAO Dependable Rainfall method is the same as the FAO/AGLW method, based on 80% probability exceedance. It estimates the amount of rainfall that can be depended upon at a given probability level.

| Condition | Formula (at 80% probability) |
|-----------|------------------------------|
| P ≤ 70 mm | $P_{eff} = \max(0.6P - 10, 0)$ |
| P > 70 mm | $P_{eff} = 0.8P - 24$ |

A probability scaling factor is applied for other probability levels:
- 50% probability: 1.3× base estimate
- 80% probability: 1.0× base estimate (default)
- 90% probability: 0.9× base estimate

**Usage:**
```python
ep = EffectivePrecipitation(..., method='dependable_rainfall', method_params={'probability': 0.80})
```

**Reference:** FAO. (1992). *CROPWAT - A computer program for irrigation planning and management*. FAO Irrigation and Drainage Paper No. 46.

---

### 5. FarmWest

A simple empirical formula used by [FarmWest program](https://farmwest.com/climate/calculator-information/et/effective-precipitation/) for irrigation scheduling in the Pacific Northwest.

$$P_{eff} = \max((P - 5) \times 0.75, 0)$$

The method assumes the first 5 mm is lost to interception/evaporation, and 75% of the remaining precipitation is effective.

**Usage:**
```python
ep = EffectivePrecipitation(..., method='farmwest')
```

**Reference:** FarmWest. [Effective Precipitation](https://farmwest.com/climate/calculator-information/et/effective-precipitation/).

---

### 6. USDA-SCS (with AWC and ETo)

The USDA Soil Conservation Service method that accounts for soil water holding capacity (AWC) and reference evapotranspiration (ETo). This method is more site-specific than other methods. Replace ETo with crop evapotranspiration (ETc) when available for more accurate results.

**Formula:**

1. Soil storage depth: $d = AWC \times MAD \times D_{root}$ (in inches), where $MAD$ is the Maximum Allowable Depletion factor (default: 0.5)
2. Storage factor: $sf = 0.531747 + 0.295164 \cdot d - 0.057697 \cdot d^2 + 0.003804 \cdot d^3$
3. Effective precipitation: $P_{eff} = sf \times (P^{0.82416} \times 0.70917 - 0.11556) \times 10^{ET_o \times 0.02426}$
4. Clamped: $P_{eff} = \min(P_{eff}, P, ET_o)$, $P_{eff} \geq 0$

!!! tip "Using ETc Instead of ETo"
    For more accurate crop-specific estimates, replace ETo (grass reference ET) with ETc (crop evapotranspiration) when available. ETc accounts for crop coefficients (Kc) and provides better estimates for specific crops. OpenET and similar platforms provide ETc data for agricultural fields.

**Required GEE Assets:**

| Region | Precipitation Asset | AWC Asset | ETo Asset |
|--------|---------------------|-----------|-----------|
| **U.S.** | `IDAHO_EPSCOR/GRIDMET` (band: `pr`, daily) | `projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite` (volumetric fraction) | `projects/openet/assets/reference_et/conus/gridmet/monthly/v1` (band: `eto`) |
| **Global** | `ECMWF/ERA5_LAND/MONTHLY_AGGR` (band: `total_precipitation_sum`) | `projects/sat-io/open-datasets/FAO/HWSD_V2_SMU` (band: `AWC`, mm/m) | `projects/climate-engine-pro/assets/ce-ag-era5-v2/daily` (band: `ReferenceET_PenmanMonteith_FAO56`) |

!!! warning "AWC Unit Conversion"
    FAO HWSD AWC is in **mm/m** and must be converted to volumetric fraction for the USDA-SCS formula. Use `--awc-scale-factor 0.001` (CLI) or `'awc_scale_factor': 0.001` (Python API). SSURGO AWC is already in volumetric fraction (inches/inch) and needs no conversion (default scale factor = 1.0).

**Usage (CLI - U.S.):**
```bash
pycropwat process --asset IDAHO_EPSCOR/GRIDMET --band pr \
    --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto --rooting-depth 1.0 ...
```

**Usage (CLI - Global):**
```bash
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR --band total_precipitation_sum \
    --method usda_scs \
    --awc-asset projects/sat-io/open-datasets/FAO/HWSD_V2_SMU --awc-band AWC \
    --awc-scale-factor 0.001 \
    --eto-asset projects/climate-engine-pro/assets/ce-ag-era5-v2/daily \
    --eto-band ReferenceET_PenmanMonteith_FAO56 --eto-is-daily \
    --rooting-depth 1.0 ...
```

**Reference:** USDA SCS. (1993). [Chapter 2 Irrigation Water Requirements](https://www.wcc.nrcs.usda.gov/ftpref/wntsc/waterMgt/irrigation/NEH15/ch2.pdf). In Part 623 National Engineering Handbook. USDA Soil Conservation Service.

---

### 7. TAGEM-SuET (with ETo)

The TAGEM-SuET (Türkiye'de Sulanan Bitkilerin Bitki Su Tüketimleri - Turkish Irrigation Management and Plant Water Consumption System) method calculates effective precipitation based on the difference between precipitation and reference evapotranspiration. When precipitation exceeds ETo, the excess becomes effective precipitation.

**Formula:**

$$
P_{eff} = \begin{cases}
0 & \text{if } P \leq ET_o \\
P - ET_o & \text{if } P > ET_o \text{ and } (P - ET_o) < 75 \\
75 + 0.0011(P - ET_o - 75)^2 + 0.44(P - ET_o - 75) & \text{otherwise}
\end{cases}
$$

!!! warning "Performance Note"
    Studies have shown that the TAGEM-SuET method tends to underperform compared to other methods, particularly in arid and semi-arid climates where ETo often exceeds precipitation. In our method comparison analyses, TAGEM-SuET consistently produced the lowest effective precipitation estimates. Users should consider this limitation when selecting a method for their application.

!!! tip "Using ETc Instead of ETo"
    For more accurate crop-specific estimates, replace ETo (grass reference ET) with ETc (crop evapotranspiration) when available. ETc = ETo × Kc, where Kc is the crop coefficient.

**Reference:** [Muratoglu, A., Bilgen, G. K., Angin, I., & Kodal, S. (2023). Performance analyses of effective rainfall estimation methods for accurate quantification of agricultural water footprint. *Water Research*, 238, 120011.](https://doi.org/10.1016/j.watres.2023.120011)

**Usage:**
```python
ep = EffectivePrecipitation(..., method='suet', method_params={
    'eto_asset': 'projects/openet/assets/reference_et/conus/gridmet/monthly/v1',
    'eto_band': 'eto'
})
```

**Usage (CLI):**
```bash
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR --band total_precipitation_sum \
    --method suet \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto ...
```

---

### 8. PCML (Physics-Constrained Machine Learning)

The PCML method uses pre-computed effective precipitation from a physics-constrained machine learning model trained specifically for the Western United States. Unlike other methods, PCML Peff is retrieved directly from a GEE asset rather than calculated from precipitation.

**Coverage:**

| Attribute | Value |
|-----------|-------|
| **Region** | Western U.S. (17 states: AZ, CA, CO, ID, KS, MT, NE, NV, NM, ND, OK, OR, SD, TX, UT, WA, WY) |
| **Temporal** | January 2000 - September 2024 (monthly) |
| **Resolution** | ~2 km (native scale retrieved dynamically from GEE asset) |
| **GEE Asset** | `projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked` |
| **Band Format** | `bYYYY_M` (e.g., `b2015_9`, `b2016_10`) - bands auto-selected by year/month |

!!! note "Pre-computed Data with Annual Fractions"
    PCML is different from other methods - it provides pre-computed Peff values from a trained ML model. When using `--method pcml`, the default PCML asset is automatically used and bands are dynamically selected based on the year/month being processed. The native scale (~2km) is retrieved from the asset using GEE's `nominalScale()` function. **Only annual (water year, Oct-Sep)** effective precipitation fractions are available for PCML (not monthly), loaded directly from a separate GEE asset (`projects/ee-peff-westus-unmasked/assets/effective_precip_fraction_unmasked`, WY 2000-2024, band format: `bYYYY`).

!!! tip "PCML Geometry Options"
    - **No geometry provided**: Downloads the entire PCML asset (full Western U.S. - 17 states)
    - **User provides geometry**: PCML data is clipped/subsetted to that geometry (e.g., just Pacific Northwest states), reducing download size and processing time

**Usage (CLI - full Western U.S.):**
```bash
pycropwat process \
    --method pcml \
    --start-year 2000 --end-year 2024 \
    --output ./WesternUS_PCML
```

**Usage (CLI - subset to specific region):**
```bash
pycropwat process \
    --method pcml \
    --geometry pacific_northwest.geojson \
    --start-year 2000 --end-year 2024 \
    --output ./PacificNW_PCML
```

**Reference:** [Hasan, M. F., Smith, R. G., Majumdar, S., Huntington, J. L., Alves Meira Neto, A., & Minor, B. A. (2025). Satellite data and physics-constrained machine learning for estimating effective precipitation in the Western United States and application for monitoring groundwater irrigation. *Agricultural Water Management*, 319, 109821.](https://doi.org/10.1016/j.agwat.2025.109821)

---

### 9. Ensemble - Default (Mean of Methods)

The ensemble method provides a robust estimate by calculating the mean of six methods, excluding TAGEM-SuET which has been shown to underperform in comparative analyses. This multi-method average reduces bias from any single method.

**Included Methods:**

1. CROPWAT - FAO standard method
2. FAO/AGLW - Dependable Rainfall (80% exceedance)
3. Fixed Percentage - 70% of precipitation
4. Dependable Rainfall - 75% probability level
5. FarmWest - Pacific Northwest method
6. USDA-SCS - Soil-based method

**Formula:**

$$P_{eff}^{ensemble} = \frac{P_{eff}^{cropwat} + P_{eff}^{fao\_aglw} + P_{eff}^{fixed} + P_{eff}^{dependable} + P_{eff}^{farmwest} + P_{eff}^{usda\_scs}}{6}$$

!!! info "Recommended Method"
    The ensemble method is recommended when users want a robust, multi-method average that reduces bias from any single method. It requires AWC and ETo assets (same as USDA-SCS) since it internally calculates all component methods.

**Required GEE Assets:** Same as USDA-SCS method (AWC and ETo).

**Usage:**
```python
ep = EffectivePrecipitation(..., method='ensemble', method_params={
    'awc_asset': 'projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite',
    'awc_band': 'AWC',
    'awc_scale_factor': 1.0,  # SSURGO is already volumetric fraction; use 0.001 for FAO HWSD (mm/m)
    'eto_asset': 'projects/openet/assets/reference_et/conus/gridmet/monthly/v1',
    'eto_band': 'eto',
    'rooting_depth': 1.0
})
```

**Usage (CLI):**
```bash
pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR --band total_precipitation_sum \
    --method ensemble \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite --awc-band AWC \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto --rooting-depth 1.0 ...
```

---

### Method Comparison

| Method | Use Case | Characteristics |
|--------|----------|-----------------|
| **CROPWAT** | General irrigation planning | Balanced, widely validated |
| **FAO/AGLW** | Yield response studies | FAO Dependable Rainfall (80% exceedance) |
| **Fixed Percentage** | Quick estimates, calibration | Simple, requires local calibration |
| **Dependable Rainfall** | Risk-averse planning | Same as FAO/AGLW, with probability scaling |
| **FarmWest** | Pacific Northwest irrigation | Simple, accounts for interception loss |
| **USDA-SCS** | Site-specific irrigation | Accounts for soil AWC and ETo |
| **TAGEM-SuET** | ET-based irrigation planning | Based on P - ETo difference |
| **PCML** | Western U.S. applications | ML-based, pre-computed (2000-2024) |
| **Ensemble** | Robust multi-method estimate | Mean of 6 methods (excludes TAGEM-SuET and PCML) |

!!! tip "Choosing a Method"
    - Use **Ensemble** (default) for robust multi-method estimates when soil and ETo data are available
    - Use **PCML** for Western U.S. applications (2000-2024) - trained specifically for that region
    - Use **CROPWAT** for most irrigation planning applications when soil data is unavailable
    - Use **FAO/AGLW** or **Dependable Rainfall** when following FAO Irrigation Paper No. 33 guidelines (they use the same base formula)
    - Use **Fixed Percentage** when you have locally calibrated effectiveness values
    - Use **FarmWest** for Pacific Northwest regions or when accounting for interception loss
    - Use **USDA-SCS** when soil AWC and ETo data are available for site-specific estimates
    - Use **TAGEM-SuET** with caution - this method tends to produce the lowest Peff estimates and may underperform in arid/semi-arid regions (see [Muratoglu et al., 2023](https://doi.org/10.1016/j.watres.2023.120011))

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

### ✅ Implemented (v1)

- 📊 **Multiple Peff methods**: CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall, FarmWest, USDA-SCS, TAGEM-SuET, PCML (Western U.S.), Ensemble
- 🗓️ **Temporal aggregation**: Annual, seasonal, growing season (with cross-year support), climatology
- 📈 **Statistical analysis**: Anomaly detection, trend analysis (linear, Theil-Sen), zonal statistics
- 📉 **Visualization**: Time series, climatology charts, maps, interactive HTML maps, dataset comparison (side-by-side, scatter, annual)
- 📤 **Export options**: NetCDF with time dimension, Cloud-Optimized GeoTIFFs
- 🌍 **Any GEE climate dataset** with precipitation band
- 💾 **Local precipitation input**: your own monthly GeoTIFFs, globs or NetCDF files; precipitation-only methods run offline with no Earth Engine
- 🗺️ **Flexible resolution control** (native or custom scale)
- ⚡ **Parallel processing** with Dask
- 🧩 **Automatic tiling** for large regions (in-memory mosaicking)
- 📁 **Example workflows**: South America, Arizona, New Mexico, Western U.S. PCML, UCRB field-scale, local precipitation (no GEE)

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

*Comparison of all 9 effective precipitation methods: CROPWAT, FAO/AGLW, Fixed Percentage (70%), Dependable Rainfall (80%), FarmWest, USDA-SCS, TAGEM-SuET, PCML (Western U.S. only), and Ensemble.*

<p align="center">
  <img src="assets/examples/method_comparison/ERA5Land_method_curves.png" width="60%" alt="Method Curves">
</p>

*Theoretical response curves for different effective precipitation methods.*

### Anomaly, Climatology & Trend Maps

<p align="center">
  <img src="assets/examples/figures/ERA5Land/anomaly_map_2023_01.png" width="32%" alt="Anomaly Map">
  <img src="assets/examples/figures/ERA5Land/climatology_map_01.png" width="32%" alt="Climatology Map">
  <img src="assets/examples/figures/ERA5Land/trend_map_with_significance.png" width="32%" alt="Trend Map">
</p>

*Left: Percent anomaly (January 2023). Center: January climatology (2000-2020). Right: Long-term trend with significance stippling (p < 0.05).*

---

### Arizona Example (USDA-SCS Method)

The [Arizona workflow](examples.md#arizona-usda-scs-example) demonstrates U.S.-specific datasets with the USDA-SCS method:

<p align="center">
  <img src="assets/examples/arizona/figures/GridMET/time_series.png" width="48%" alt="Arizona Time Series">
  <img src="assets/examples/arizona/figures/GridMET/monthly_climatology.png" width="48%" alt="Arizona Climatology">
</p>

*GridMET USDA-SCS effective precipitation for Arizona: time series (left) and monthly climatology (right).*

<p align="center">
  <img src="assets/examples/arizona/figures/GridMET/anomaly_map_2023_08.png" width="32%" alt="Arizona Anomaly">
  <img src="assets/examples/arizona/figures/GridMET/climatology_map_08.png" width="32%" alt="Arizona Climatology Map">
  <img src="assets/examples/arizona/figures/GridMET/trend_map_with_significance.png" width="32%" alt="Arizona Trend">
</p>

*Left: August 2023 anomaly (monsoon). Center: August climatology. Right: Long-term trend with significance stippling.*

### New Mexico Example (8-Method Comparison)

The [New Mexico workflow](examples.md#new-mexico-method-comparison-example) compares 8 effective precipitation methods using PRISM data (excludes PCML):

<p align="center">
  <img src="assets/examples/new_mexico/method_comparison/mean_annual_comparison_1986_2025.png" width="100%" alt="New Mexico Mean Annual Comparison">
</p>

*Mean annual effective precipitation (1986-2025) for 8 methods: CROPWAT, FAO/AGLW, Fixed 70%, Dependable Rainfall, FarmWest, USDA-SCS, TAGEM-SuET, and Ensemble.*

<p align="center">
  <img src="assets/examples/new_mexico/method_comparison/method_curves_new_mexico.png" width="60%" alt="New Mexico Method Curves">
</p>

*Theoretical response curves for 8 effective precipitation methods using New Mexico typical values (AWC=120mm/m, ETo=200mm/month).*

## Documentation

For full documentation, visit [https://montimaj.github.io/pyCropWat](https://montimaj.github.io/pyCropWat)

- [Quick Start Guide](user-guide/quickstart.md)
- [CLI Reference](user-guide/cli.md)
- [Python API](user-guide/api.md)
- [Local Precipitation Data](user-guide/local-data.md) - use your own GeoTIFFs/NetCDF instead of a GEE download
- [Examples](examples.md)
- [Complete Workflow Example](examples.md#example-12-complete-workflow) - A comprehensive script demonstrating all features

!!! tip "Try the Complete Workflow Examples"
    The `Examples/` directory provides ready-to-run demonstrations:
    
    **Rio de la Plata Example (Global datasets):**
    ```bash
    # Full workflow with GEE processing (generates ~5 GB of output data)
    python Examples/south_america_example.py --gee-project your-project-id --workers 8
    
    # Analysis only (requires previously generated data)
    python Examples/south_america_example.py --analysis-only
    ```
    
    **Arizona USDA-SCS Example (U.S. datasets):**
    ```bash
    # Full workflow with GEE processing (generates ~12 GB of output data)
    python Examples/arizona_example.py --gee-project your-project-id --workers 8
    
    # Analysis only (requires previously generated data)
    python Examples/arizona_example.py --analysis-only
    ```
    
    **New Mexico 8-Method Comparison:**
    ```bash
    # Full workflow with GEE processing (generates ~8 GB of output data)
    python Examples/new_mexico_example.py --gee-project your-project-id --workers 8
    ```
    
    **Western U.S. PCML Example (17 states):**
    ```bash
    # PCML effective precipitation with water year aggregation (generates ~3 GB of output data)
    python Examples/western_us_pcml_example.py --gee-project your-project-id --workers 8
    
    # Analysis only (requires previously generated data)
    python Examples/western_us_pcml_example.py --analysis-only
    ```
    
    **UCRB Field-Scale Example (GeoPackage):**
    ```bash
    # Field-scale Peff from existing precipitation volumes (~3 MB of output data)
    # Note: The GeoPackage file (~7 GB) is not included; contact authors for access
    python Examples/ucrb_example.py
    ```
    
    **Local Precipitation Example (no Earth Engine):**
    ```bash
    # Effective precipitation straight from monthly GeoTIFFs on disk - runs offline
    python Examples/LocalPrecip/local_precip_quickstart.py
    
    # A different method, a subset of years, and specific months
    python Examples/LocalPrecip/local_precip_quickstart.py --method farmwest \
        --start-year 2005 --end-year 2010 --months 6 7 8
    
    # NetCDF round trip (build a .nc from the GeoTIFFs, then process it) - offline
    python Examples/LocalPrecip/local_netcdf_example.py
    
    # Hybrid: local precipitation + GEE AWC/ETo (needs Earth Engine authentication)
    python Examples/LocalPrecip/wrf_south_america_local_example.py --gee-project your-project-id
    ```

## U.S.-Specific Datasets

For U.S.-based applications, pyCropWat supports high-resolution precipitation and the USDA-SCS method:

### Precipitation Datasets

| Dataset | Asset ID | Band | Resolution |
|---------|----------|------|------------|
| **GridMET** | `IDAHO_EPSCOR/GRIDMET` | `pr` | ~4 km |
| **PRISM** | `projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY` | `ppt` | ~800 m |

### USDA-SCS Required Assets

| Dataset | Asset ID | Band | Description |
|---------|----------|------|-------------|
| **AWC (SSURGO)** | `projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite` | (single) | USDA SSURGO soil data |
| **ETo (GridMET)** | `projects/openet/assets/reference_et/conus/gridmet/monthly/v1` | `eto` | OpenET GridMET monthly ETo |

### Example: Arizona with USDA-SCS

```bash
pycropwat process --asset IDAHO_EPSCOR/GRIDMET --band pr \
    --gee-geometry users/montimajumdar/AZ \
    --start-year 2000 --end-year 2024 --scale 4000 \
    --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --eto-band eto --rooting-depth 1.0 \
    --output ./AZ_GridMET_USDA_SCS
```

See [Examples](examples.md#arizona-usda-scs-example) for the complete Arizona workflow script.

## Local Precipitation Input

*New in v1.3.0.* Instead of downloading precipitation from Earth Engine, point pyCropWat at
precipitation you already have on disk - a directory of monthly GeoTIFFs, a glob, or one or more
NetCDF files. Everything else is unchanged: AWC, ETo and GEE `FeatureCollection` geometries are still
read from Earth Engine, all the same methods are available, and the output file names are identical
(`effective_precip_YYYY_MM.tif`, `effective_precip_fraction_YYYY_MM.tif`).

The precipitation-only methods need **no Earth Engine at all** - `ee.Initialize()` is never called,
so those runs work offline and without GEE credentials.

### Supported Layouts

| Layout | `local_precip` / `--local-precip` value | Notes |
|--------|------------------------------------------|-------|
| Directory of monthly rasters | `../pyCropWat_Data/Precip` | One file per month, dated from the file name. Select the files with `--local-precip-pattern` (e.g. `'Precip_*.tif'`; default `'*.tif'`) |
| Glob of monthly rasters | `'../pyCropWat_Data/Precip/Precip_*.tif'` | Same as above, without a separate pattern |
| Single monthly raster | `./Precip_2005_07.tif` | Processes just that month |
| Single NetCDF stack | `./wrf_precip.nc` | Months read from the `time` coordinate; set `--local-precip-variable` when auto-detection is ambiguous |
| Several NetCDF files | `'./nc/*.nc'` | Combined by coordinates when they carry a `time` axis, otherwise dated one-month-per-file from the file name |

- **Raster suffixes:** `.tif`, `.tiff`, `.vrt`, `.img`, `.bil`, `.bsq`, `.bip`, `.asc`, `.jp2`, `.grd`, `.dat`, `.hgt`
- **NetCDF/HDF suffixes:** `.nc`, `.nc4`, `.cdf`, `.netcdf`, `.h5`, `.hdf5`
- Anything else in the directory (`.aux.xml` sidecars, `.json` state files, `.csv` metadata) is ignored, so stray companion files never break indexing.
- File names are dated as `YYYY_MM`, `YYYY-MM`, `YYYY.MM` or `YYYYMM`. Any other convention works via `--local-precip-date-regex` / `local_precip_date_regex` (a regex with named groups `year` and `month`).
- One multi-band raster per year (12 bands = 12 months) is **not** supported - split it into monthly files.
- `--scale-factor` / `precip_scale_factor` converts the local units to millimetres exactly as it does for GEE data, and nodata (the value in the file metadata plus any `--local-precip-nodata` sentinel) is propagated as `NaN` to **both** output rasters.
- Missing months are skipped with a warning, and `--start-year`/`--end-year` are inferred from the files when omitted (an explicit range is clamped to what is available).

### Python Example

The example dataset used by the `Examples/LocalPrecip/` scripts is 264 monthly WRF regional climate
model GeoTIFFs (`Precip_YYYY_MM.tif`, 2000-01 to 2021-12, South America / Rio de la Plata domain,
EPSG:4326, 689 x 799 pixels at ~0.0397° (~4.4 km), float32, nodata `-9999`, already in mm/month):

```python
from pycropwat import EffectivePrecipitation

# Precipitation-only method: no Earth Engine, no authentication, works offline
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',   # directory, glob, or NetCDF file
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,                 # on top of the file's own nodata
    method='cropwat',
    start_year=2005,                           # optional - inferred from the files when omitted
    end_year=2005
)
results = ep.process(output_dir='./output', n_workers=4, months=[7])
```

Keep AWC and ETo coming from Earth Engine by choosing a method that needs them:

```python
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    geometry_path='roi.geojson',      # local rasters are clipped to this (clip_to_geometry=True)
    start_year=2005,
    end_year=2010,
    method='ensemble',                # needs AWC + ETo, so Earth Engine IS initialized
    gee_project='your-gee-project',
    method_params={
        'awc_asset': 'projects/sat-io/open-datasets/FAO/HWSD_V2_SMU',
        'awc_band': 'AWC',
        'awc_scale_factor': 0.001,    # mm/m -> volumetric fraction
        'eto_asset': 'IDAHO_EPSCOR/TERRACLIMATE',
        'eto_band': 'pet',
        'eto_scale_factor': 0.1,
        'eto_is_daily': False,
        'rooting_depth': 2.0,
        'mad_factor': 1.0
    }
)
results = ep.process(output_dir='./output', n_workers=4)
```

Read the monthly arrays directly, without running the workflow:

```python
from pycropwat import open_local_precipitation

with open_local_precipitation('../pyCropWat_Data/Precip',
                              pattern='Precip_*.tif',
                              nodata=-9999) as src:
    print(src.kind, len(src), src.year_range, src.crs, src.shape)
    # raster 264 (2000, 2021) EPSG:4326 (689, 799)
    da = src.get_month(2005, 7)     # xarray.DataArray, dims ('y', 'x'), mm, NaN nodata
```

### CLI Example

```bash
# No Earth Engine at all - years inferred from the file names
pycropwat process --local-precip ../pyCropWat_Data/Precip \
                  --local-precip-pattern 'Precip_*.tif' \
                  --local-precip-nodata -9999 \
                  --method cropwat --output ./output

# Local precipitation, with AWC and ETo still from Earth Engine
pycropwat process --local-precip ../pyCropWat_Data/Precip \
                  --local-precip-pattern 'Precip_*.tif' \
                  --local-precip-nodata -9999 \
                  --geometry roi.geojson \
                  --start-year 2005 --end-year 2010 \
                  --method ensemble \
                  --awc-asset projects/sat-io/open-datasets/FAO/HWSD_V2_SMU \
                  --awc-band AWC --awc-scale-factor 0.001 \
                  --eto-asset projects/climate-engine-pro/assets/ce-ag-era5-v2/daily \
                  --eto-band ReferenceET_PenmanMonteith_FAO56 --eto-is-daily \
                  --rooting-depth 2.0 --mad-factor 1.0 \
                  --output ./output

# NetCDF input with an explicit variable and CRS override
# (--local-precip-crs relabels the grid; it never reprojects it)
pycropwat process --local-precip ./wrf_precip.nc \
                  --local-precip-variable RAINNC \
                  --local-precip-crs EPSG:4326 \
                  --method fao_aglw --output ./output
```

New `process` flags: `--local-precip`, `--local-precip-pattern`, `--local-precip-variable`,
`--local-precip-nodata`, `--local-precip-crs`, `--local-precip-date-regex`. With `--local-precip`,
`--asset`/`--band` are not needed at all and `--start-year`/`--end-year` become optional.
`--local-precip-crs` is an **override**: it replaces whatever CRS the files declare (and supplies
one when they declare none), relabelling the grid in place without reprojecting it.

### Earth Engine Requirements

Which methods still touch Earth Engine when precipitation comes from local files:

| Method | Extra inputs | Earth Engine needed with local precipitation? |
|--------|--------------|-----------------------------------------------|
| `cropwat` | Precipitation only | ❌ No - fully offline |
| `fao_aglw` | Precipitation only | ❌ No - fully offline |
| `fixed_percentage` | Precipitation only | ❌ No - fully offline |
| `dependable_rainfall` | Precipitation only | ❌ No - fully offline |
| `farmwest` | Precipitation only | ❌ No - fully offline |
| `suet` | ETo | ✅ Yes - ETo is downloaded from GEE |
| `usda_scs` | AWC + ETo | ✅ Yes - AWC and ETo are downloaded from GEE |
| `ensemble` (default) | AWC + ETo | ✅ Yes - AWC and ETo are downloaded from GEE |
| `pcml` | Pre-computed GEE product | ⛔ Not supported - `local_precip` with `method='pcml'` is rejected |

!!! note "Geometry also decides whether Earth Engine starts"
    Earth Engine is initialized when the study area is a GEE `FeatureCollection` (`--gee-geometry`,
    or `--geometry` pointing at a GEE asset), whatever the method. A local shapefile/GeoJSON keeps an
    offline run offline and is used to clip the local rasters (`clip_to_geometry=True` by default);
    with no geometry at all, the full extent of the local files is used.

!!! warning "AWC/ETo coverage"
    For `usda_scs`, `ensemble` and `suet`, AWC and ETo are downloaded only for the geometry you
    supply. Any part of the local precipitation grid the download did not cover stays **NaN** in
    the outputs - pyCropWat does not backfill it with a stand-in value - and the run logs a
    WARNING naming the uncovered pixel count, once per field per process (so `usda_scs` and
    `ensemble` log it twice: once for AWC, once for ETo):

    ```
    AWC covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
    ```

    This bites when the grid is wider than the download region: a GEE `FeatureCollection` geometry
    never clips the local rasters, and `clip_to_geometry=False` deliberately keeps them whole.
    Clipping to the geometry (the default for a local vector file) or widening the geometry gives
    full coverage. Precipitation-only methods download nothing and are unaffected.

See the [Local Precipitation Data guide](user-guide/local-data.md) for the full walkthrough, the
[`pycropwat.local` API reference](api/local.md) for `LocalPrecipitationSource`,
`open_local_precipitation` and `parse_year_month`, and
[Examples/LocalPrecip/](https://github.com/montimaj/pyCropWat/blob/main/Examples/LocalPrecip/README.md)
for runnable scripts.

## Citation

If you use pyCropWat in your research, please cite:

```bibtex
@software{pycropwat,
  author = {Majumdar, Sayantan and ReVelle, Peter and Pearson, Christopher and Nozari, Soheil and Minor, Blake A. and Hasan, M. F. and Huntington, Justin L. and Smith, Ryan G.},
  title = {pyCropWat: A Python Package for Computing Effective Precipitation Using Google Earth Engine Climate Data},
  year = {2026},
  url = {https://github.com/montimaj/pyCropWat},
  doi = {10.5281/zenodo.18201619}
}
```

### Effective Precipitation Method References

- Muratoglu, A., Bilgen, G. K., Angin, I., & Kodal, S. (2023). Performance analyses of effective rainfall estimation methods for accurate quantification of agricultural water footprint. *Water Research*, *238*, 120011. [https://doi.org/10.1016/j.watres.2023.120011](https://doi.org/10.1016/j.watres.2023.120011)

- Smith, M. (1992). *CROPWAT: A computer program for irrigation planning and management* (FAO Irrigation and Drainage Paper No. 46). Food and Agriculture Organization of the United Nations. [https://www.fao.org/sustainable-development-goals-helpdesk/champion/article-detail/cropwat/en](https://www.fao.org/sustainable-development-goals-helpdesk/champion/article-detail/cropwat/en)

- USDA SCS. (1993). Chapter 2 Irrigation Water Requirements. In Part 623 National Engineering Handbook. USDA Soil Conservation Service. [https://www.wcc.nrcs.usda.gov/ftpref/wntsc/waterMgt/irrigation/NEH15/ch2.pdf](https://www.wcc.nrcs.usda.gov/ftpref/wntsc/waterMgt/irrigation/NEH15/ch2.pdf)

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

MIT License - see [LICENSE](https://github.com/montimaj/pyCropWat/blob/main/LICENSE) for details.
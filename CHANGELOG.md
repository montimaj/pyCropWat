# Changelog

All notable changes to pyCropWat will be documented in this file.

## [1.2] - 2026-01-23

### ✨ New Features
- **PCML Method**: Added Physics-Constrained Machine Learning (PCML) effective precipitation method for Western U.S.
  - Pre-computed Peff from GEE asset: `projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked`
  - Coverage: 17 Western U.S. states (AZ, CA, CO, ID, KS, MT, NE, NV, NM, ND, OK, OR, SD, TX, UT, WA, WY)
  - Temporal: January 2000 - September 2024 (monthly)
  - Resolution: ~2 km (native scale retrieved dynamically from GEE asset)
  - Band format: `bYYYY_M` (e.g., `b2015_9` for September 2015)
  - Annual (water year, Oct-Sep) fractions from separate GEE asset
  - CLI: `pycropwat process --method pcml --start-year 2000 --end-year 2024 --output ./output`
  - Reference: [Hasan et al. (2025)](https://doi.org/10.1016/j.agwat.2025.109821)
- **Simplified PCML CLI**: PCML method no longer requires `--asset`, `--band`, or `--geometry` arguments
  - Default asset and band are automatically set when `--method pcml` is used
  - Geometry defaults to full Western U.S. extent (17 states)
  - User can optionally provide geometry to subset the region
- **UCRB Field-Scale Example**: Added new Upper Colorado River Basin example for field-scale Peff calculations
  - Uses existing precipitation volumes from GeoPackage
  - Demonstrates AWC lookup from CSV and all 8 Peff methods

### 📁 New Files
- `Examples/western_us_pcml_example.py` - Western U.S. PCML workflow with water year aggregation
- `Examples/ucrb_example.py` - UCRB field-scale Peff calculation example

### 📚 Documentation
- Updated all PCML documentation to clarify that only Western U.S. vectors overlapping the 17-state extent can be used
- Added PCML CLI examples to help text and docstrings
- Synchronized badges between README.md and docs/index.md

---

## [1.1.1.post3] - 2026-01-12

### 📚 Documentation
- **FarmWest Reference Fix**: Removed incorrect Washington State University attribution from FarmWest method descriptions

---

## [1.1.1.post2] - 2026-01-11

### 🔧 Improvements
- **Default Method Changed**: Changed default method from `cropwat` to `ensemble` for more robust multi-method estimates
  - Ensemble requires AWC and ETo assets but provides superior results by averaging 6 methods
  - CROPWAT remains available for users without AWC/ETo data

### 🐛 Fixes
- **AWC Band Selection Fix**: Fixed issue where `--awc-band` defaulted to `'AWC'` causing errors with single-band SSURGO
  - Now defaults to `None` - SSURGO works without specifying band, HWSD uses `--awc-band AWC`

### 📚 Documentation
- Updated all documentation to reflect ensemble as the default method

---

## [1.1.0] - 2026-01-11

### ✨ New Features
- **Ensemble Method**: Added new robust effective precipitation method that calculates the mean of 6 methods (excludes TAGEM-SuET and PCML)
  - Formula: Peff_ensemble = (CROPWAT + FAO/AGLW + Fixed 70% + Dependable 75% + FarmWest + USDA-SCS) / 6
  - Requires AWC and ETo data (same as USDA-SCS) via `--awc-asset` and `--eto-asset` CLI flags
  - Recommended for robust multi-method estimates that reduce bias from any single method
- **TAGEM-SuET Method**: Added new effective precipitation method based on P - ETo difference (Turkish Irrigation Management and Plant Water Consumption System)
  - Formula: Peff = 0 if P ≤ ETo; Peff = P - ETo if (P - ETo) < 75; else Peff = 75 + 0.0011(P-ETo-75)² + 0.44(P-ETo-75)
  - Requires ETo data via `--eto-asset` CLI flag or `method_params`
  - ⚠️ Note: Studies show TAGEM-SuET tends to underperform in arid/semi-arid climates (see [Muratoglu et al., 2023](https://doi.org/10.1016/j.watres.2023.120011))
- **Cross-Year Season Aggregation**: Added support for temporal aggregations spanning two calendar years (e.g., Southern Hemisphere growing season Oct-Mar)
  - `custom_aggregate()` now accepts `cross_year=True` parameter for cross-year aggregations
  - `growing_season_aggregate()` auto-detects cross-year seasons when `start_month > end_month`
  - Example: `agg.growing_season_aggregate(2020, start_month=10, end_month=3)` aggregates Oct 2020 - Mar 2021
- **Chunked Download for Large Regions**: Automatic chunked download for AWC and ETo data when region exceeds GEE pixel limits
  - New generic `_download_image_chunked()` method consolidates all tiled downloads
  - In-memory mosaicking for faster performance (no temp files)

### 🔧 Improvements
- **Code Refactoring**: Consolidated three chunked download methods into one generic `_download_image_chunked()` method
  - Removed redundant `_download_awc_chunked()` and `_download_eto_chunked()` methods
  - `_download_chunked()` now wraps the generic method for precipitation DataArray output
  - Removed unused `_mosaic_tiles()` method and `tempfile`/`merge_arrays` imports
  - ~150 lines of code reduction with cleaner architecture

### 🐛 Fixes
- **FAO/AGLW Formula Correction**: Fixed threshold from 250mm to 70mm, formula from 0.8P-25 to 0.8P-24
- **Dependable Rainfall Formula Correction**: Fixed threshold from 100mm to 70mm, default probability from 75% to 80%
- **Method Descriptions**: Clarified that CROPWAT and USDA-SCS are different methods (removed "USDA SCS" from CROPWAT descriptions)
- **Large Region Downloads**: Fixed `_create_tiles` → `_create_tile_grid` method name in chunked downloads

### 📚 Documentation
- Updated all method formulas and descriptions across README, docs, and examples
- Added TAGEM-SuET to all method comparison tables and feature lists
- Updated CLI help text for ETo asset to mention both usda_scs and suet methods
- Added New Mexico method comparison example with efficient single-download workflow
- Fixed South America example to use correct Southern Hemisphere growing season (Oct-Mar)

---

## [1.0.5.post1] - 2026-01-10

### 🎨 UI/Branding
- Updated logo styling with green "Crop" and blue "Wat" text colors in MkDocs docs
- Added logo to Overview section in documentation and README
- Switched PyPI downloads badge to pepy.tech for reliability

---

## [1.0.4] - 2026-01-09

### 📦 Package & Distribution
- Added animated logo with PNG fallback for PyPI compatibility
- Added PyPI downloads and GitHub stars badges

---

## [1.0.3] - 2026-01-09

### 📦 Package & Distribution
- Added Zenodo DOI badge to README and documentation

---

## [1.0.2] - 2026-01-09

### 📦 Package & Distribution
- Excluded documentation assets and example figures from PyPI package to reduce size
- Converted all relative links and image paths to absolute GitHub URLs for PyPI README rendering

---

## [1.0.1] - 2026-01-09

### 📦 Package & Distribution
- Added PyPI publishing workflow via GitHub Actions with trusted publishing
- Added `MANIFEST.in` to exclude large example files from PyPI package
- Updated package description to "A Python Package for Computing Effective Precipitation Using Google Earth Engine Climate Data"
- Added Zenodo DOI (`10.5281/zenodo.18201619`) to citation

### 🔧 Fixes
- Fixed Git clone URLs to use correct repository path (`montimaj/pyCropWat`)

---

## [1.0.0] - 2026-01-08

### ✨ New Features

#### Multiple Effective Precipitation Methods
- **Ensemble (default)**: Ensemble mean of 6 methods (requires AWC and ETo)
- **CROPWAT**: CROPWAT method from FAO
- **FAO/AGLW**: FAO Dependable Rainfall (80% exceedance)
- **Fixed Percentage**: Configurable percentage method (default 70%)
- **Dependable Rainfall**: FAO method at specified probability levels (50-90%)
- **FarmWest**: Washington State University's simple empirical formula: `Peff = (P - 5) × 0.75`
- **USDA-SCS with AWC**: Site-specific method using Available Water Capacity and Reference ET from GEE assets
- **TAGEM-SuET**: Turkish Irrigation Management System method based on P - ETo difference

#### USDA-SCS Method with AWC and ETo
- Accounts for soil water holding capacity (AWC) and evaporative demand (ETo)
- Supports U.S. datasets:
  - AWC: `projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite`
  - ETo: `projects/openet/assets/reference_et/conus/gridmet/monthly/v1`
- Supports Global datasets:
  - AWC: `projects/sat-io/open-datasets/FAO/HWSD_V2_SMU`
  - ETo: `projects/climate-engine-pro/assets/ce-ag-era5-v2/daily`
- Configurable crop rooting depth (default: 1 meter)
- Daily ETo aggregation to monthly supported via `--eto-is-daily` flag

#### Temporal Aggregation (`pycropwat.analysis.TemporalAggregator`)
- Seasonal aggregation (DJF, MAM, JJA, SON)
- Annual totals with configurable statistics (sum, mean, min, max, std)
- Growing season aggregation with customizable start/end months
- Custom date range aggregation
- Multi-year climatology calculation

#### Statistical Analysis (`pycropwat.analysis.StatisticalAnalyzer`)
- Anomaly calculation (absolute, percent, standardized)
- Trend analysis with linear regression
- Theil-Sen slope with Mann-Kendall significance test
- Zonal statistics for polygon features (CSV export)

#### Visualization (`pycropwat.analysis.Visualizer`)
- Time series plots
- Monthly climatology bar charts
- Single raster map visualization
- Anomaly maps with diverging colormaps (`plot_anomaly_map()`)
- Climatology maps (`plot_climatology_map()`)
- Trend maps with significance stippling (`plot_trend_map()`)
- Trend panel with slope and p-value side by side (`plot_trend_panel()`)
- Interactive maps using leafmap or folium (`plot_interactive_map()`)
- Side-by-side dataset comparison with difference map (`plot_comparison()`)
- Scatter plot comparison with R², RMSE, bias statistics (`plot_scatter_comparison()`)
- Annual totals comparison bar chart (`plot_annual_comparison()`)

#### Enhanced Export Options
- NetCDF export with time dimension (`export_to_netcdf()`)
- Cloud-Optimized GeoTIFF conversion (`export_to_cog()`)
- Zonal statistics CSV export

#### CLI Enhancements
- **New subcommand structure**: `pycropwat <command> [OPTIONS]`
- `process` subcommand: Main effective precipitation calculation
- `aggregate` subcommand: Temporal aggregation (annual, seasonal, growing season, climatology)
- `analyze` subcommand: Statistical analysis (anomaly, trend, zonal)
- `export` subcommand: Export to NetCDF or Cloud-Optimized GeoTIFF
- `plot` subcommand: Visualization (timeseries, climatology, map, interactive, compare, scatter, annual-compare)
- `--method` flag to select effective precipitation method
- `--percentage` flag for fixed_percentage method
- `--probability` flag for dependable_rainfall method
- `--awc-asset` flag for USDA-SCS method AWC GEE asset
- `--awc-band` flag for AWC band name
- `--eto-asset` flag for USDA-SCS method ETo GEE asset
- `--eto-band` flag for ETo band name
- `--eto-is-daily` flag for daily ETo aggregation to monthly
- `--rooting-depth` flag for crop rooting depth (USDA-SCS method)
- `--list-methods` to display available methods
- `--version` flag to display version
- Legacy mode support for backwards compatibility

### 📚 Documentation
- Added comprehensive MkDocs documentation with GitHub Pages deployment
- Added anomaly, climatology, and trend map visualization examples
- Added Arizona USDA-SCS example comparing U.S. vs Global datasets
- Added disk space requirements to installation guide
- Fixed image paths for GitHub README rendering

### 📁 Examples
- **South America (Rio de la Plata)**: Complete workflow with ERA5-Land and TerraClimate comparison
- **Arizona (USDA-SCS)**: U.S.-focused workflow with GridMET/PRISM precipitation, SSURGO AWC, and OpenET ETo
- **New Mexico**: 8-method comparison workflow with PRISM precipitation
- Example outputs (~32 GB) are generated locally by running the scripts

### 📦 New Dependencies
- `scipy>=1.9.0` - Statistical functions
- `matplotlib>=3.5.0` - Visualization
- `rasterstats>=0.18.0` - Zonal statistics
- `pandas>=1.4.0` - Data manipulation

### 📦 Optional Dependencies
- `leafmap>=0.30.0` - Interactive maps (optional)
- `folium>=0.14.0` - Alternative interactive maps (optional)

### 📁 New Files
- `pycropwat/methods.py` - Effective precipitation calculation methods
- `pycropwat/analysis.py` - Temporal aggregation, statistics, visualization
- `Examples/arizona_example.py` - Arizona workflow (8 methods, excludes PCML)
- `Examples/south_america_example.py` - Rio de la Plata workflow (8 methods, excludes PCML)
- `Examples/new_mexico_example.py` - New Mexico workflow (8 methods, excludes PCML)

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
- Christopher Pearson (Desert Research Institute)
- Soheil Nozari (Colorado State University)
- Justin Huntington (Desert Research Institute)
- Ryan Smith (Colorado State University)

### 🙏 Acknowledgments

This work was supported by the U.S. Army Corps of Engineers (Grant W912HZ25C0016) for the project *"Improved Characterization of Groundwater Resources in Transboundary Watersheds using Satellite Data and Integrated Models."*

### 📄 References

- Smith, M. (1992). *CROPWAT: A computer program for irrigation planning and management* (FAO Irrigation and Drainage Paper No. 46). Food and Agriculture Organization of the United Nations. https://www.fao.org/3/t7202e/t7202e00.htm
- Muratoglu, A., Bilgen, G. K., Angin, I., & Kodal, S. (2023). Performance analyses of effective rainfall estimation methods for accurate quantification of agricultural water footprint. *Water Research*, *238*, 120011. https://doi.org/10.1016/j.watres.2023.120011

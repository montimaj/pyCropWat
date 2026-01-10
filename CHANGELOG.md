# Changelog

All notable changes to pyCropWat will be documented in this file.

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
- Added Zenodo DOI (`10.5281/zenodo.18201620`) to citation

### 🔧 Fixes
- Fixed Git clone URLs to use correct repository path (`montimaj/pyCropWat`)

---

## [1.0.0] - 2026-01-08

### ✨ New Features

#### Multiple Effective Precipitation Methods
- **CROPWAT (default)**: USDA SCS method as implemented in FAO CROPWAT
- **FAO/AGLW**: FAO/AGLW formula from FAO Irrigation Paper No. 33
- **Fixed Percentage**: Configurable percentage method (default 70%)
- **Dependable Rainfall**: FAO method at specified probability levels (50-90%)
- **FarmWest**: Washington State University's simple empirical formula: `Peff = (P - 5) × 0.75`
- **USDA-SCS with AWC**: Site-specific method using Available Water Capacity and Reference ET from GEE assets

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
- Added pre-computed example outputs for both regions

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
- `Examples/arizona_usda_scs_example.py` - Arizona USDA-SCS workflow
- `Examples/south_america_cropwat_example.py` - Rio de la Plata CROPWAT workflow

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

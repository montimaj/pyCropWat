# Changelog

All notable changes to pyCropWat will be documented in this file.

## [1.3.0] - 2026-09-04

### ✨ Added
- **Local Precipitation Input**: effective precipitation can now be computed from **your own**
  precipitation rasters instead of Google Earth Engine. New module `pycropwat/local.py`:
  - `LocalPrecipitationSource` - indexed, lazily-read source for a directory of monthly rasters, a
    glob string, or one/many NetCDF files. Properties: `kind`, `files`, `crs`, `bounds`,
    `native_bounds`, `year_range`, `resolution`, `shape`; methods: `available_months()`,
    `has_month()`, `get_month()`, `close()`; also a context manager with `__len__`/`__repr__`
  - `open_local_precipitation(path, **kwargs)` - convenience factory
  - `parse_year_month(name, date_regex=None)` - dates a file name using the `YYYY_MM`, `YYYY-MM`,
    `YYYYMM` and `YYYY.MM` layouts, or a custom regex with named `year`/`month` groups
  - Monthly grids are returned as 2-D `('y', 'x')` `float32` DataArrays in mm, north-up, nodata as
    NaN, with the source CRS attached
- **New `EffectivePrecipitation` arguments** (all appended after the existing ones, so positional
  and keyword calls stay backward compatible): `local_precip`, `local_precip_pattern` (default
  `'*.tif'`), `local_precip_variable`, `local_precip_nodata`, `local_precip_crs`,
  `local_precip_date_regex`, `clip_to_geometry` (default `True`)
- **New CLI flags** on `pycropwat process`: `--local-precip`, `--local-precip-pattern`,
  `--local-precip-variable`, `--local-precip-nodata`, `--local-precip-date-regex` and
  `--local-precip-crs`, plus three new examples in the `process` help epilog (local rasters with no
  Earth Engine, local rasters with GEE AWC/ETo, and a local NetCDF file)
- **`--eto-scale-factor`** on `pycropwat process` (default `1.0`), the ETo counterpart to the
  existing `--awc-scale-factor`. It is passed to the `usda_scs`, `suet` and `ensemble` methods.
  Previously `eto_scale_factor` could only be set through the Python API's `method_params`, so any
  documented CLI command using TerraClimate `pet` (stored in 0.1 mm units) silently ran USDA-SCS
  with an ETo 10x too large - measured at 1.70x the mean and 5.3x the maximum effective
  precipitation. The affected commands in the docs now pass `--eto-scale-factor 0.1`
- **Lazy Earth Engine initialization**: `ee.Initialize()` is only called when something actually
  needs it. Local precipitation with a precipitation-only method (`cropwat`, `fao_aglw`,
  `fixed_percentage`, `dependable_rainfall`, `farmwest`) needs **no Earth Engine at all** and runs
  offline; `self.geometry` and `self.collection` are then `None`
- **Year inference and month filtering**: `start_year`/`end_year` are inferred from the local files
  when omitted, a supplied range is clamped to the available years (with a warning), and requested
  months with no file are skipped and reported instead of failing
- **Geometry is optional with local precipitation**: with no geometry the extent of the local files
  is used; a local vector `geometry_path` clips each monthly raster (`clip_to_geometry=True`)
- **`local_precip_crs` / `--local-precip-crs` is an override**: it replaces whatever CRS the local
  files declare, and supplies one when they declare none. Replacing a declared CRS is logged at
  INFO. It relabels the grid in place - values and transform are untouched - and never reprojects
- **AWC/ETo regridding for local grids**: GEE-sourced AWC and ETo are reprojected onto the local
  precipitation grid with `reproject_match()` (falling back to the previous zoom behaviour on
  error), and the local grid resolution is used as the GEE download scale. Pixels of the local grid
  that the download never covered stay **NaN** and propagate to NaN in the `usda_scs`, `ensemble`
  and `suet` outputs; the run logs a warning naming the uncovered fraction, once per field per
  process (so two lines for `usda_scs`/`ensemble`: one for AWC, one for ETo). pyCropWat does not
  substitute a filler value for a region it never downloaded
- **Works under `dask.distributed`**: `EffectivePrecipitation` is now picklable, so
  `process()` runs under a `Client`/`LocalCluster` with worker **processes**. Earth Engine is
  re-initialized inside each worker (its state is process-local), and the local source handle and
  its lock are rebuilt on unpickle rather than being serialized
- **Grid validation for local files**: CRS, shape and transform are checked across the whole
  directory when the index is built (metadata only - indexing 264 rasters takes ~0.12 s). A
  mismatched file raises an error naming it and how it differs, instead of silently producing a
  time series on inconsistent grids. With `local_precip_crs` set, a directory whose files *declare*
  different CRSs still warns
- **2-D curvilinear coordinates**: WRF-style NetCDF files are read whether their 2-D `XLAT`/`XLONG`
  arrays are promoted by a CF `coordinates` attribute or left as plain data variables, and are
  collapsed to 1-D when the grid is regular enough; a genuinely curvilinear grid raises rather than
  returning a wrong one
- `Examples/LocalPrecip/` - `local_precip_quickstart.py`, `wrf_south_america_local_example.py`
  (the local-file counterpart of `Examples/wrf/wrf_south_america_example.py`) and
  `local_netcdf_example.py`, with their own README
- `docs/user-guide/local-data.md` guide and `docs/api/local.md` API page, both added to the
  `mkdocs.yml` nav
- `tests/test_local.py` - tests for the local precipitation source

### 🔧 Changed
- `asset_id` and `precip_band` are now `Optional[str] = None`; they are required only when
  `local_precip` is not given and the method is not `pcml`. Omitting a precipitation source raises
  `ValueError("No precipitation source provided...")`
- CLI: `--start-year`/`--end-year` are no longer `required=True` (inferred from the files with
  `--local-precip`), and the `--asset`/`--band` check now reads
  `--asset and --band are required (unless --local-precip is given, or --method pcml is used)`
- New exports from `pycropwat`: `LocalPrecipitationSource`, `open_local_precipitation`,
  `parse_year_month`
- Output rasters are written with the CRS of the precipitation grid instead of a hard-coded
  `EPSG:4326`, so local grids in any projection round-trip correctly (no change for GEE runs, which
  are always EPSG:4326)
- `--local-precip` is rejected with `--method pcml` (PCML is a pre-computed Earth Engine product)
- Version bumped to 1.3.0 (`pycropwat/__init__.py`, `pyproject.toml`, `pycropwat --version`)

### 🐛 Fixes
- **NaN precipitation no longer reports a zero effective precipitation fraction**: nodata pixels now
  propagate as NaN to both `effective_precip_*.tif` and `effective_precip_fraction_*.tif`, the
  arithmetic runs under `np.errstate(divide='ignore', invalid='ignore')`, and NaN is declared as the
  raster nodata value. GEE downloads never produce NaN, so GEE outputs are unchanged
- Removed duplicate `'usda_scs'` and `'suet'` entries from the method list in the `pycropwat`
  package docstring
- Added the missing `suet_effective_precip`, `pcml_effective_precip` and `ensemble_effective_precip`
  exports to `pycropwat.__all__` (previously importable only from `pycropwat.methods`)
- CLI: the `--local-precip` summary no longer announces the source settings before it has managed to
  open the source. A source that cannot be read is now silent (debug only) instead of printing five
  INFO lines and a warning that merely duplicated the real error raised moments later
- **`TemporalAggregator` no longer turns nodata into 0 mm.** `sum`/`mean`/`min`/`max`/`std` use
  xarray's `skipna=True`, so a pixel that was NaN in *every* contributing month aggregated to
  exactly `0.0` and the result was written with no nodata value - making ocean and outside-domain
  pixels indistinguishable from genuinely dry land in annual, seasonal, growing-season, custom and
  climatology products, and via `pycropwat aggregate`. Such pixels now stay NaN and every aggregate
  declares `nodata=NaN`. A pixel missing only *some* months still aggregates over what exists, as
  before. This was unreachable before 1.3.0 because GEE downloads use `defaultValue=0` and never
  produced NaN; pixel values for NaN-free inputs are unchanged, only the nodata tag is new
- `StatisticalAnalyzer.calculate_anomaly` and `.calculate_trend` now declare `nodata=NaN` on their
  outputs too, matching the aggregator instead of writing untagged NaN rasters
- An Earth Engine authentication or initialization failure is no longer swallowed by the AWC/ETo
  fallback and silently reported as `0.15` / `100 mm`. It now propagates; the constant fallback is
  kept only for its original case - an asset that genuinely has no data over the region - and still
  logs at WARNING
- Suppressed two `RuntimeWarning`s raised by all-NaN pixels that had no way to occur before local
  nodata existed: `Mean of empty slice` from the `ensemble` method, and
  `Degrees of freedom <= 0 for slice` from the `std` aggregation and standardized anomalies
- Constructing `EffectivePrecipitation` without a usable precipitation source now reports which
  piece is missing - a missing `precip_band` names the asset and suggests a band for it, rather
  than claiming no source was provided at all
- Docs: `zonal_statistics(zones_path=...)` corrected to `geometry_path=...`,
  `multi_year_climatology()` no longer shown with a non-existent `method=` argument, USDA-SCS
  keywords moved inside `method_params` where several examples had them at top level, and the
  unresolvable `LICENSE` link and dead New Mexico anchor in `docs/index.md` fixed. All pre-existing
  since before 1.3.0

### 📚 Documentation
- New "Local Precipitation Data" user guide and API reference page
- Quick Start, CLI Reference, Python API, Examples and `Examples/README.md` updated for local
  precipitation input
- `--asset`/`--band` help now says they are not required with `--local-precip` either, not just for
  `--method pcml`
- Documented that `precip_scale_factor` is a single multiplier and so cannot convert a rate such as
  kg m⁻² s⁻¹ to a monthly total (the factor depends on the length of each month); `docs/examples.md`
  no longer suggests `86400` for that case

### 🏗️ Packaging
- **Minimum Python is now 3.10** (was 3.9). `requires-python`, the PyPI classifiers, the Black
  `target-version`, `environment.yml` and the documented requirements were updated together.
  Python 3.9 reached end of life in October 2025. No source change was needed - the codebase
  already avoided 3.10-only syntax - so 1.3.0 still runs on 3.9 if installed from source; the
  bound simply stops pip resolving it there
- `.gitignore`: `Examples/*` excluded the `Examples/LocalPrecip/` directory itself, and git cannot
  re-include a file whose parent directory is excluded, so the new example scripts and their README
  could never be committed. The directory is now re-included before its `*.py`/`*.md` contents,
  while the generated output directories (Peff/AWC/ETo GeoTIFFs, NetCDF stacks) stay ignored
- `.gitignore`: `Examples/wrf/` was ignored wholesale, which also hid the example scripts the docs
  link to. Its `*.py` files are now tracked (except the untitled `test.py` scratch script) while the
  ~7 GB of WRF/ERA5/AgERA5 inputs and `peff_ensemble/` outputs stay ignored

---

## [1.2.1] - 2026-02-19

### 🔧 Improvements
- **AWC Scale Factor**: Added `awc_scale_factor` parameter to handle FAO HWSD AWC unit conversion (mm/m → volumetric fraction with `awc_scale_factor=0.001`). SSURGO AWC requires no conversion (default 1.0).
- CLI: Added `--awc-scale-factor` argument for `usda_scs` and `ensemble` methods.
---

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

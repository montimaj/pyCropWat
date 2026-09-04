# CLI Reference

The `pycropwat` command-line interface provides access to all functionality through subcommands.

## Synopsis

```bash
pycropwat <command> [OPTIONS]
pycropwat --help
pycropwat --version
pycropwat --list-methods
```

## Commands

| Command | Description |
|---------|-------------|
| `process` | Calculate effective precipitation from GEE or local climate data |
| `aggregate` | Temporal aggregation (annual, seasonal, growing season) |
| `analyze` | Statistical analysis (anomaly, trend, zonal statistics) |
| `export` | Export to NetCDF or Cloud-Optimized GeoTIFF |
| `plot` | Create visualizations (time series, climatology, maps) |

---

## Process Command

Calculate effective precipitation from Google Earth Engine climate data, or from precipitation
rasters and NetCDF files you already have on disk (`--local-precip`).

```bash
pycropwat process [OPTIONS]
```

### Required Options

| Option | Description |
|--------|-------------|
| `--asset`, `-a` | GEE ImageCollection asset ID (not required with `--local-precip`, or for `--method pcml`) |
| `--band`, `-b` | Precipitation band name in the asset (not required with `--local-precip`, or for `--method pcml`) |
| `--start-year`, `-s` | Start year for processing (inclusive). Not required with `--local-precip` - inferred from the files |
| `--end-year`, `-e` | End year for processing (inclusive). Not required with `--local-precip` - inferred from the files |
| `--output`, `-o` | Output directory for GeoTIFF files |

!!! note "`--local-precip` relaxes the required options"
    When `--local-precip` is given, precipitation comes from disk instead of Earth Engine, so
    `--asset` and `--band` are not needed (if you pass them anyway they are simply unused).
    `--start-year`/`--end-year` are optional too: they are inferred from the available files, and a
    range you do supply is clamped to the years that exist on disk (with a warning). Only
    `--output` stays mandatory.

    Omitting `--asset`/`--band` **without** `--local-precip` (and without `--method pcml`) exits
    with `--asset and --band are required (unless --local-precip is given, or --method pcml is used)`.

### Geometry Options (one required)

| Option | Description |
|--------|-------------|
| `--geometry`, `-g` | Path to local shapefile or GeoJSON |
| `--gee-geometry`, `-G` | GEE FeatureCollection asset ID |

A geometry is **optional** with `--local-precip` (and for `--method pcml`). With local files and no
geometry, the extent of the precipitation files is used. With `--geometry` pointing at a local
vector file, each monthly raster is clipped to it.

### Local Precipitation Options

Use these to read monthly precipitation from disk instead of downloading it from Earth Engine.

| Option | Default | Description |
|--------|---------|-------------|
| `--local-precip` | None | Use local precipitation instead of GEE: a directory of monthly rasters, a NetCDF file, or a glob string. AWC/ETo still come from GEE; precipitation-only methods need no GEE |
| `--local-precip-pattern` | `*.tif` | Glob used when `--local-precip` is a directory (default: `'*.tif'`) |
| `--local-precip-variable` | None | NetCDF variable holding precipitation (auto-detected by default) |
| `--local-precip-nodata` | None | Extra nodata sentinel for local files, e.g. -9999. Applied on top of the value stored in the file metadata |
| `--local-precip-crs` | None | CRS **override** for local files, e.g. `'EPSG:4326'`. Replaces whatever CRS the files declare (logged at INFO) and supplies one when they declare none. Relabels the grid only - nothing is reprojected |
| `--local-precip-date-regex` | None | Regex with named groups `'year'` and `'month'` for dating local files, e.g. `'(?P<month>[0-9]{2})_(?P<year>[0-9]{4})'` |

`--scale-factor` is applied to local files as well, so use it to convert your units to millimetres
per month (e.g. `--scale-factor 1000` for metres). Being a single multiplier, it cannot turn a rate
such as kg m⁻² s⁻¹ into a monthly total - accumulate those to monthly totals beforehand.

!!! warning "Earth Engine is still used for AWC and ETo"
    `--method usda_scs`, `--method ensemble` and `--method suet` download AWC and/or ETo from Earth
    Engine even with `--local-precip`, so they still require GEE credentials. The
    precipitation-only methods (`cropwat`, `fao_aglw`, `fixed_percentage`, `dependable_rainfall`,
    `farmwest`) need no Earth Engine at all - it is never initialized.

    `--local-precip` cannot be combined with `--method pcml`; the run exits with
    `--local-precip cannot be used with --method pcml (PCML is a pre-computed Earth Engine product)`.

    For those three methods, AWC and ETo are downloaded only for the geometry you supply. Any part
    of the local precipitation grid the download did not cover stays **NaN** in the outputs -
    pyCropWat does not backfill it - and the run logs a WARNING naming the uncovered pixel count,
    once per field per process (so `--method usda_scs` and `--method ensemble` log it twice: once
    for AWC, once for ETo):

    ```
    AWC covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
    ```

    The `clip_to_geometry=True` the message suggests is the Python-API name for what the CLI
    already does: passing `--geometry` with a local vector file clips the rasters to the downloaded
    region, which gives full coverage (there is no `--clip-to-geometry` flag, and no way to turn
    the clipping off from the CLI). A `--gee-geometry` FeatureCollection never clips them, so make
    sure it spans the whole grid.

### Optional Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `--scale`, `-r` | Native | Output resolution in meters |
| `--scale-factor`, `-f` | 1.0 | Factor to convert precipitation to mm |
| `--months`, `-m` | All | Specific months to process (1-12) |
| `--workers`, `-w` | 4 | Number of parallel workers |
| `--project`, `-p` | None | GEE project ID for authentication |
| `--method` | ensemble | Effective precipitation method |
| `--percentage` | 0.7 | Percentage for fixed_percentage method |
| `--probability` | 0.75 | Probability for dependable_rainfall method |
| `--awc-asset` | None | GEE AWC asset for usda_scs method |
| `--awc-band` | AWC | AWC band name |
| `--eto-asset` | None | GEE ETo asset for usda_scs method |
| `--eto-band` | eto | ETo band name |
| `--eto-is-daily` | False | Set if ETo asset is daily (aggregates to monthly) |
| `--eto-scale-factor` | 1.0 | Scale factor for ETo, applied after aggregation. TerraClimate `pet` is stored in 0.1 mm units: use `0.1`. GridMET/AgERA5 ETo is already in mm: `1.0` |
| `--rooting-depth` | 1.0 | Crop rooting depth in meters for usda_scs |
| `--mad-factor` | 0.5 | Management Allowed Depletion factor (0-1) for usda_scs |
| `--sequential` | False | Process sequentially instead of parallel |

### Available Methods

| Method | Description |
|--------|-------------|
| `cropwat` | CROPWAT method from FAO |
| `fao_aglw` | FAO/AGLW Dependable Rainfall (80% exceedance) |
| `fixed_percentage` | Fixed percentage of rainfall |
| `dependable_rainfall` | Dependable rainfall approach |
| `farmwest` | FarmWest method: Peff = (P - 5) × 0.75 |
| `usda_scs` | USDA-SCS with AWC and ETo (requires additional assets) |
| `suet` | TAGEM-SuET method: P - ETo with 75mm threshold (requires ETo asset) |
| `pcml` | Physics-Constrained ML (Western U.S. only, uses default GEE asset) |
| `ensemble` | Mean of 6 methods - default (excludes TAGEM-SuET and PCML; requires AWC and ETo) |

### Examples

#### Basic Usage with Local Geometry

```bash
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry ./data/study_area.shp \
    --start-year 2020 \
    --end-year 2023 \
    --scale-factor 1000 \
    --output ./outputs
```

#### Using GEE Vector Asset

```bash
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --gee-geometry projects/my-project/assets/watersheds \
    --start-year 2015 \
    --end-year 2020 \
    --scale-factor 1000 \
    --output ./era5_outputs
```

#### Using Different Peff Methods

```bash
# FAO AGLW method
pycropwat process \
    --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson \
    --start-year 2015 --end-year 2020 \
    --method fao_aglw --output ./output

# Fixed percentage (80%)
pycropwat process \
    --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson \
    --start-year 2015 --end-year 2020 \
    --method fixed_percentage --percentage 0.8 --output ./output

# Dependable rainfall (80% probability)
pycropwat process \
    --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson \
    --start-year 2015 --end-year 2020 \
    --method dependable_rainfall --probability 0.8 --output ./output

# FarmWest method
pycropwat process \
    --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \
    --geometry roi.geojson \
    --start-year 2015 --end-year 2020 \
    --method farmwest --output ./output
```

#### Local Precipitation Files (No Earth Engine)

Point `--local-precip` at a directory of monthly rasters. With a precipitation-only method
(`cropwat`, `fao_aglw`, `fixed_percentage`, `dependable_rainfall`, `farmwest`) Earth Engine is never
initialized, so this runs offline and without credentials:

```bash
# 264 monthly GeoTIFFs named Precip_YYYY_MM.tif, nodata -9999, already in mm
pycropwat process \
    --local-precip ../pyCropWat_Data/Precip \
    --local-precip-pattern 'Precip_*.tif' \
    --local-precip-nodata -9999 \
    --method cropwat \
    --output ./out_local_cropwat
```

`--start-year`/`--end-year` were omitted, so the CLI logs
`Date range: inferred from the local precipitation files` and then the resolved
`Processing years: 2000 - 2021`. Add them to restrict the run, clip to a region with `--geometry`,
and select months as usual:

```bash
pycropwat process \
    --local-precip ../pyCropWat_Data/Precip \
    --local-precip-pattern 'Precip_*.tif' \
    --local-precip-nodata -9999 \
    --geometry roi.geojson \
    --start-year 2005 --end-year 2010 \
    --months 10 11 12 1 2 3 \
    --method farmwest \
    --workers 8 \
    --output ./out_local_farmwest
```

A NetCDF stack works the same way; name the variable when it cannot be auto-detected, and use
`--local-precip-crs` to supply a CRS the file lacks or to override one it declares wrongly (it
relabels the grid, it never reprojects):

```bash
pycropwat process \
    --local-precip ./wrf_precip.nc \
    --local-precip-variable RAINNC \
    --local-precip-crs EPSG:4326 \
    --method fao_aglw \
    --output ./out_local_netcdf
```

If your file names use a layout the built-in parser does not recognise (`YYYY_MM`, `YYYY-MM`,
`YYYYMM`, `YYYY.MM`), supply your own regex with `year` and `month` named groups:

```bash
pycropwat process \
    --local-precip ./rain \
    --local-precip-pattern 'rain_*.tif' \
    --local-precip-date-regex '(?P<month>[0-9]{2})_(?P<year>[0-9]{4})' \
    --method cropwat \
    --output ./out_local_regex
```

#### Local Precipitation with GEE AWC and ETo (Ensemble Method)

`--method ensemble` (and `usda_scs`, `suet`) still pulls AWC and ETo from Earth Engine, so this is a
hybrid run: precipitation from disk, soil and evaporative demand from GEE. The AWC/ETo grids are
reprojected onto the local precipitation grid, and any pixel the download did not reach is left NaN
(see the warning under [Local Precipitation Options](#local-precipitation-options)). The
`--geometry` below clips the rasters to exactly the downloaded region, so every output pixel has
precipitation, AWC and ETo:

```bash
# WRF South America precipitation on disk + FAO HWSD AWC + TerraClimate ETo from GEE
pycropwat process \
    --local-precip ../pyCropWat_Data/Precip \
    --local-precip-pattern 'Precip_*.tif' \
    --local-precip-nodata -9999 \
    --geometry riodelaplata.geojson \
    --start-year 2005 --end-year 2010 \
    --method ensemble \
    --awc-asset projects/sat-io/open-datasets/FAO/HWSD_V2_SMU \
    --awc-band AWC --awc-scale-factor 0.001 \
    --eto-asset IDAHO_EPSCOR/TERRACLIMATE --eto-band pet --eto-scale-factor 0.1 \
    --rooting-depth 2.0 --mad-factor 1.0 \
    --project your-gee-project \
    --output ./out_local_ensemble
```

!!! warning "ETo unit scaling"
    ETo assets are not all stored in millimetres. TerraClimate `pet` is in 0.1 mm units, so it needs
    `--eto-scale-factor 0.1`; omitting it feeds USDA-SCS an ETo **10x too large** and produces
    plausible-looking but wrong effective precipitation. GridMET and AgERA5 ETo are already in mm,
    so they keep the default `1.0`. The example above sets it explicitly for that reason.

#### USDA-SCS Method (U.S. Data)

The USDA-SCS method requires Available Water Capacity (AWC) and Reference Evapotranspiration (ETo) datasets:

```bash
# U.S. with SSURGO AWC and GridMET ETo
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry arizona.geojson \
    --start-year 2015 --end-year 2020 \
    --scale-factor 1000 \
    --method usda_scs \
    --awc-asset projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite \
    --eto-asset projects/openet/assets/reference_et/conus/gridmet/monthly/v1 \
    --rooting-depth 1.0 --mad-factor 0.5 \
    --output ./output_usda_scs
```

#### USDA-SCS Method (Global Data)

For global applications, use FAO HWSD AWC and AgERA5 ETo:

```bash
# Global with FAO HWSD AWC and AgERA5 ETo (daily - requires aggregation)
pycropwat process \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry study_area.geojson \
    --start-year 2015 --end-year 2020 \
    --scale-factor 1000 \
    --method usda_scs \
    --awc-asset projects/sat-io/open-datasets/FAO/HWSD_V2_SMU \
    --awc-band AWC \
    --eto-asset projects/climate-engine-pro/assets/ce-ag-era5-v2/daily \
    --eto-band ReferenceET_PenmanMonteith_FAO56 \
    --eto-is-daily \
    --rooting-depth 1.0 --mad-factor 0.5 \
    --output ./output_usda_scs_global
```

#### PCML Method (Western U.S. Only)

The Physics-Constrained Machine Learning (PCML) method uses a pre-computed GEE asset covering the Western United States. **No geometry is required** - the PCML asset's built-in geometry automatically covers the 17 western states.

```bash
# PCML for Western U.S. - no geometry needed
pycropwat process \
    --start-year 2015 --end-year 2020 \
    --method pcml \
    --output ./output_pcml
```

**PCML Method Details:**

| Attribute | Value |
|-----------|-------|
| **Coverage** | Western U.S. (17 states) |
| **Temporal** | January 2000 - September 2024 |
| **Resolution** | ~2 km (retrieved dynamically from asset) |
| **Outputs** | PCML Peff + annual (water year, Oct-Sep) fraction rasters (from GEE asset, band format: `bYYYY`) |

**PCML Geometry Options:**

- **No geometry provided**: Downloads the entire PCML asset (full Western U.S. - 17 states)
- **User provides geometry**: PCML data is clipped/subsetted to that geometry (e.g., just Pacific Northwest states), reducing download size and processing time

```bash
# PCML for a specific region (e.g., Pacific Northwest)
pycropwat process \
    --geometry pacific_northwest.geojson \
    --start-year 2015 --end-year 2020 \
    --method pcml \
    --output ./output_pcml_pnw
```

---

## Aggregate Command

Temporal aggregation of effective precipitation rasters.

```bash
pycropwat aggregate [OPTIONS]
```

### Options

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory with monthly rasters (required) |
| `--type`, `-t` | Aggregation type: annual, seasonal, growing-season, custom, climatology (required) |
| `--output`, `-o` | Output file or directory (required) |
| `--year`, `-y` | Year to aggregate |
| `--start-year` | Start year (for climatology) |
| `--end-year` | End year (for climatology) |
| `--season` | Season code: DJF, MAM, JJA, SON |
| `--start-month` | Growing season start month (default: 4) |
| `--end-month` | Growing season end month (default: 10) |
| `--months`, `-m` | Specific months for custom aggregation |
| `--method` | Aggregation method: sum, mean, min, max, std (default: sum) |
| `--pattern` | File glob pattern (default: effective_precip_[0-9]*.tif, excludes fraction files) |

### Examples

```bash
# Annual total
pycropwat aggregate --input ./output --type annual --year 2020 --output ./annual_2020.tif

# Seasonal (summer)
pycropwat aggregate --input ./output --type seasonal --year 2020 --season JJA --output ./summer_2020.tif

# Growing season (April-October)
pycropwat aggregate --input ./output --type growing-season --year 2020 \
    --start-month 4 --end-month 10 --output ./growing_2020.tif

# Multi-year climatology
pycropwat aggregate --input ./output --type climatology \
    --start-year 2000 --end-year 2020 --output ./climatology/
```

---

## Analyze Command

Statistical analysis of effective precipitation data.

```bash
pycropwat analyze <analysis_type> [OPTIONS]
```

### Anomaly Analysis

Calculate anomaly relative to climatology.

```bash
pycropwat analyze anomaly --input ./output --year 2020 --month 6 \
    --clim-start 1990 --clim-end 2020 --output ./anomaly_2020_06.tif
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory (required) |
| `--year`, `-y` | Target year (required) |
| `--month`, `-m` | Target month 1-12 (required) |
| `--clim-start` | Climatology start year (required) |
| `--clim-end` | Climatology end year (required) |
| `--anomaly-type` | Type: absolute, percent, standardized (default: absolute) |
| `--output`, `-o` | Output file (required) |

### Trend Analysis

Calculate temporal trend.

```bash
pycropwat analyze trend --input ./output --start-year 2000 --end-year 2020 \
    --trend-method sen --output ./trend/
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory (required) |
| `--start-year` | Start year (required) |
| `--end-year` | End year (required) |
| `--month`, `-m` | Specific month (or annual if omitted) |
| `--trend-method` | Method: linear, sen (default: linear) |
| `--output`, `-o` | Output directory (required) |

### Zonal Statistics

Calculate statistics by zone.

```bash
pycropwat analyze zonal --input ./output --zones ./regions.shp \
    --start-year 2000 --end-year 2020 --output ./zonal_stats.csv
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory (required) |
| `--zones`, `-z` | Path to zones shapefile/GeoJSON (required) |
| `--start-year` | Start year (required) |
| `--end-year` | End year (required) |
| `--months`, `-m` | Specific months |
| `--stats` | Statistics: mean,sum,min,max,std (default: all) |
| `--output`, `-o` | Output CSV file (required) |

---

## Export Command

Export data to different formats.

```bash
pycropwat export <format> [OPTIONS]
```

### NetCDF Export

```bash
pycropwat export netcdf --input ./output --output ./data.nc
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory (required) |
| `--output`, `-o` | Output NetCDF file (required) |
| `--pattern` | File glob pattern (default: effective_precip_[0-9]*.tif, excludes fraction files) |
| `--variable` | Variable name in NetCDF |
| `--no-compression` | Disable compression |

### Cloud-Optimized GeoTIFF

```bash
pycropwat export cog --input ./effective_precip_2020_06.tif --output ./cog_2020_06.tif
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input GeoTIFF file (required) |
| `--output`, `-o` | Output COG file (required) |

---

## Plot Command

Create visualizations.

```bash
pycropwat plot <plot_type> [OPTIONS]
```

### Time Series Plot

```bash
pycropwat plot timeseries --input ./output --start-year 2000 --end-year 2020 \
    --output ./timeseries.png
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory (required) |
| `--start-year` | Start year (required) |
| `--end-year` | End year (required) |
| `--months`, `-m` | Specific months |
| `--stat` | Spatial statistic: mean, sum, min, max (default: mean) |
| `--title` | Plot title |
| `--width`, `--height` | Figure size in inches (default: 12x6) |
| `--output`, `-o` | Output image file (required) |

### Climatology Bar Chart

```bash
pycropwat plot climatology --input ./output --start-year 2000 --end-year 2020 \
    --output ./climatology.png
```

### Raster Map

```bash
pycropwat plot map --input ./output --year 2020 --month 6 --output ./map_2020_06.png
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory (required) |
| `--year`, `-y` | Year (required) |
| `--month`, `-m` | Month 1-12 (required) |
| `--cmap` | Colormap name (default: YlGnBu) |
| `--vmin`, `--vmax` | Color scale limits |
| `--title` | Plot title |
| `--output`, `-o` | Output image file (required) |

### Interactive Map

Create an interactive HTML map with raster overlay (requires `leafmap` or `folium`).

```bash
pycropwat plot interactive --input ./output --year 2020 --month 6 --output ./map.html
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | Input directory (required) |
| `--year`, `-y` | Year (required) |
| `--month`, `-m` | Month 1-12 (required) |
| `--cmap` | Colormap name (default: YlGnBu) |
| `--opacity` | Layer opacity 0-1 (default: 0.7) |
| `--basemap` | Basemap type (default: OpenStreetMap) |
| `--output`, `-o` | Output HTML file (required) |

### Compare Datasets

Create side-by-side comparison with difference map.

```bash
pycropwat plot compare --input ./era5_output --other-input ./terraclimate_output \
    --year 2020 --month 6 --label1 ERA5 --label2 TerraClimate --output ./comparison.png
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | First input directory (required) |
| `--other-input` | Second input directory (required) |
| `--year`, `-y` | Year (required) |
| `--month`, `-m` | Month 1-12 (required) |
| `--label1`, `--label2` | Labels for the datasets |
| `--cmap` | Colormap for values (default: YlGnBu) |
| `--diff-cmap` | Colormap for difference (default: RdBu) |
| `--output`, `-o` | Output image file (required) |

### Scatter Comparison

Create scatter plot with R², RMSE, and bias statistics.

```bash
pycropwat plot scatter --input ./era5_output --other-input ./terraclimate_output \
    --start-year 2000 --end-year 2020 --output ./scatter.png
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | First input directory (required) |
| `--other-input` | Second input directory (required) |
| `--start-year` | Start year (required) |
| `--end-year` | End year (required) |
| `--months`, `-m` | Specific months (default: all) |
| `--sample-size` | Max points to plot (default: 10000) |
| `--label1`, `--label2` | Axis labels for datasets |
| `--output`, `-o` | Output image file (required) |

### Annual Comparison

Create bar chart comparing annual totals between datasets.

```bash
pycropwat plot annual-compare --input ./era5_output --other-input ./terraclimate_output \
    --start-year 2000 --end-year 2020 --output ./annual_compare.png
```

| Option | Description |
|--------|-------------|
| `--input`, `-i` | First input directory (required) |
| `--other-input` | Second input directory (required) |
| `--start-year` | Start year (required) |
| `--end-year` | End year (required) |
| `--label1`, `--label2` | Legend labels for datasets |
| `--output`, `-o` | Output image file (required) |

---

## Common GEE Assets

<table>
<thead>
<tr>
<th>Dataset</th>
<th>Asset ID</th>
<th>Precipitation Band</th>
<th>Scale Factor</th>
<th>Resolution</th>
</tr>
</thead>
<tbody>
<tr>
<td>ERA5-Land Monthly</td>
<td><code>ECMWF/ERA5_LAND/MONTHLY_AGGR</code></td>
<td><code>total_precipitation_sum</code></td>
<td>1000</td>
<td>~11 km</td>
</tr>
<tr>
<td>TerraClimate</td>
<td><code>IDAHO_EPSCOR/TERRACLIMATE</code></td>
<td><code>pr</code></td>
<td>1</td>
<td>~4 km</td>
</tr>
<tr>
<td>CHIRPS Daily</td>
<td><code>UCSB-CHG/CHIRPS/DAILY</code></td>
<td><code>precipitation</code></td>
<td>1</td>
<td>~5.5 km</td>
</tr>
<tr>
<td>CHIRPS Pentad</td>
<td><code>UCSB-CHG/CHIRPS/PENTAD</code></td>
<td><code>precipitation</code></td>
<td>1</td>
<td>~5.5 km</td>
</tr>
<tr>
<td>GPM IMERG</td>
<td><code>NASA/GPM_L3/IMERG_MONTHLY_V06</code></td>
<td><code>precipitation</code></td>
<td>1</td>
<td>~10 km</td>
</tr>
<tr>
<td>AgERA5</td>
<td><code>projects/climate-engine-pro/assets/ce-ag-era5-v2/daily</code></td>
<td><code>Precipitation_Flux</code></td>
<td>1</td>
<td>~9 km</td>
</tr>
<tr>
<td>GridMET</td>
<td><code>IDAHO_EPSCOR/GRIDMET</code></td>
<td><code>pr</code></td>
<td>1</td>
<td>~4 km (U.S. only)</td>
</tr>
<tr>
<td>PRISM</td>
<td><code>OREGONSTATE/PRISM/AN81m</code></td>
<td><code>ppt</code></td>
<td>1</td>
<td>~4 km (U.S. only)</td>
</tr>
</tbody>
</table>

## AWC and ETo Assets for USDA-SCS Method

### Available Water Capacity (AWC)

| Region | Asset ID | Band | Units |
|--------|----------|------|-------|
| U.S. (SSURGO) | `projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite` | `AWC` | cm/cm |
| Global (FAO HWSD) | `projects/sat-io/open-datasets/FAO/HWSD_V2_SMU` | `AWC` | cm/cm |

### Reference Evapotranspiration (ETo)

| Region | Asset ID | Band | Temporal |
|--------|----------|------|----------|
| U.S. (GridMET) | `projects/openet/assets/reference_et/conus/gridmet/monthly/v1` | `eto` | Monthly |
| Global (AgERA5) | `projects/climate-engine-pro/assets/ce-ag-era5-v2/daily` | `ReferenceET_PenmanMonteith_FAO56` | Daily (use `--eto-is-daily`) |

## Exit Codes

| Code | Description |
|------|-------------|
| 0 | Success |
| 1 | Error (invalid arguments, processing failure, etc.) |
| 130 | Interrupted (Ctrl+C) |

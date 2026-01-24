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
| `process` | Calculate effective precipitation from GEE climate data |
| `aggregate` | Temporal aggregation (annual, seasonal, growing season) |
| `analyze` | Statistical analysis (anomaly, trend, zonal statistics) |
| `export` | Export to NetCDF or Cloud-Optimized GeoTIFF |
| `plot` | Create visualizations (time series, climatology, maps) |

---

## Process Command

Calculate effective precipitation from Google Earth Engine climate data.

```bash
pycropwat process [OPTIONS]
```

### Required Options

| Option | Description |
|--------|-------------|
| `--asset`, `-a` | GEE ImageCollection asset ID |
| `--band`, `-b` | Precipitation band name in the asset |
| `--start-year`, `-s` | Start year for processing (inclusive) |
| `--end-year`, `-e` | End year for processing (inclusive) |
| `--output`, `-o` | Output directory for GeoTIFF files |

### Geometry Options (one required)

| Option | Description |
|--------|-------------|
| `--geometry`, `-g` | Path to local shapefile or GeoJSON |
| `--gee-geometry`, `-G` | GEE FeatureCollection asset ID |

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

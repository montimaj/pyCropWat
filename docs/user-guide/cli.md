# CLI Reference

The `pycropwat` command-line interface provides access to all functionality.

## Synopsis

```bash
pycropwat [OPTIONS]
```

## Options

### Required Options

| Option | Description |
|--------|-------------|
| `--asset` | GEE ImageCollection asset ID |
| `--band` | Precipitation band name in the asset |
| `--start-year` | Start year for processing (inclusive) |
| `--end-year` | End year for processing (inclusive) |
| `--output`, `-o` | Output directory for GeoTIFF files |

### Geometry Options (one required)

| Option | Description |
|--------|-------------|
| `--geometry`, `-g` | Path to local shapefile or GeoJSON |
| `--gee-geometry`, `-G` | GEE FeatureCollection asset ID |

### Optional Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `--scale` | Native | Output resolution in meters |
| `--scale-factor` | 1.0 | Factor to convert precipitation to mm |
| `--months` | All | Specific months to process (1-12) |
| `--workers` | 4 | Number of parallel workers |
| `--gee-project` | None | GEE project ID for authentication |

## Examples

### Basic Usage with Local Geometry

```bash
pycropwat \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry ./data/study_area.shp \
    --start-year 2020 \
    --end-year 2023 \
    --scale-factor 1000 \
    --output ./outputs
```

### Using GEE Vector Asset

```bash
pycropwat \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --gee-geometry projects/my-project/assets/watersheds \
    --start-year 2015 \
    --end-year 2020 \
    --scale-factor 1000 \
    --output ./era5_outputs
```

### TerraClimate Dataset

```bash
pycropwat \
    --asset IDAHO_EPSCOR/TERRACLIMATE \
    --band pr \
    --gee-geometry projects/my-project/assets/study_area \
    --start-year 2000 \
    --end-year 2023 \
    --output ./terraclimate_outputs
```

!!! note
    TerraClimate precipitation is already in mm, so no `--scale-factor` is needed.

### CHIRPS Daily Dataset

```bash
pycropwat \
    --asset UCSB-CHG/CHIRPS/DAILY \
    --band precipitation \
    --geometry ./africa_basin.geojson \
    --start-year 2010 \
    --end-year 2020 \
    --output ./chirps_outputs
```

!!! note
    Daily data is automatically summed to monthly totals.

### Custom Resolution

```bash
pycropwat \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry ./study_area.geojson \
    --start-year 2020 \
    --end-year 2023 \
    --scale 5000 \
    --scale-factor 1000 \
    --output ./outputs_5km
```

### Process Specific Months

```bash
pycropwat \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry ./study_area.geojson \
    --start-year 2020 \
    --end-year 2023 \
    --months 6 7 8 \
    --scale-factor 1000 \
    --output ./summer_outputs
```

### Increase Parallel Workers

```bash
pycropwat \
    --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
    --band total_precipitation_sum \
    --geometry ./study_area.geojson \
    --start-year 2000 \
    --end-year 2023 \
    --scale-factor 1000 \
    --workers 8 \
    --output ./outputs
```

## Common GEE Assets

| Dataset | Asset ID | Precipitation Band | Scale Factor | Resolution |
|---------|----------|-------------------|--------------|------------|
| ERA5-Land Monthly | "ECMWF/ERA5_LAND/MONTHLY_AGGR" | total_precipitation_sum | 1000 | ~11 km |
| TerraClimate | "IDAHO_EPSCOR/TERRACLIMATE" | pr | 1 | ~4 km |
| CHIRPS Daily | "UCSB-CHG/CHIRPS/DAILY" | precipitation | 1 | ~5.5 km |
| CHIRPS Pentad | "UCSB-CHG/CHIRPS/PENTAD" | precipitation | 1 | ~5.5 km |
| GPM IMERG | "NASA/GPM_L3/IMERG_MONTHLY_V06" | precipitation | 1 | ~10 km |
| AgERA5 | "projects/climate-engine-pro/assets/ce-ag-era5-v2/daily" | Precipitation_Flux | 1 | ~9 km |

## Exit Codes

| Code | Description |
|------|-------------|
| 0 | Success |
| 1 | Error (invalid arguments, processing failure, etc.) |
| 130 | Interrupted (Ctrl+C) |

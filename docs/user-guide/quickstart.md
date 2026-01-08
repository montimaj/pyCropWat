# Quick Start

This guide will help you get started with pyCropWat in minutes.

## Basic Usage

### Using the Command Line

The simplest way to use pyCropWat is through the CLI:

```bash
pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
          --band total_precipitation_sum \
          --geometry study_area.geojson \
          --start-year 2020 \
          --end-year 2023 \
          --scale-factor 1000 \
          --output ./outputs
```

### Using the Python API

```python
from pycropwat import EffectivePrecipitation

# Initialize with a local geometry file
ep = EffectivePrecipitation(
    asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    precip_band='total_precipitation_sum',
    geometry_path='study_area.geojson',
    start_year=2020,
    end_year=2023,
    precip_scale_factor=1000  # Convert ERA5 meters to mm
)

# Process all months
results = ep.process(output_dir='./outputs', n_workers=4)
```

## Using GEE Vector Assets

You can use a GEE FeatureCollection asset instead of a local file:

=== "CLI"
    ```bash
    pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
              --band total_precipitation_sum \
              --gee-geometry projects/my-project/assets/study_area \
              --start-year 2020 \
              --end-year 2023 \
              --scale-factor 1000 \
              --output ./outputs
    ```

=== "Python"
    ```python
    ep = EffectivePrecipitation(
        asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
        precip_band='total_precipitation_sum',
        gee_geometry_asset='projects/my-project/assets/study_area',
        start_year=2020,
        end_year=2023,
        precip_scale_factor=1000
    )
    ```

## Processing Specific Months

To process only certain months (e.g., growing season):

=== "CLI"
    ```bash
    pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \
              --band total_precipitation_sum \
              --geometry study_area.geojson \
              --start-year 2020 \
              --end-year 2023 \
              --months 4 5 6 7 8 9 \
              --scale-factor 1000 \
              --output ./outputs
    ```

=== "Python"
    ```python
    results = ep.process(
        output_dir='./outputs',
        n_workers=4,
        months=[4, 5, 6, 7, 8, 9]  # April through September
    )
    ```

## Custom Resolution

By default, pyCropWat uses the native resolution of the dataset. You can specify a custom resolution:

```bash
pycropwat --asset IDAHO_EPSCOR/TERRACLIMATE \
          --band pr \
          --geometry study_area.geojson \
          --start-year 2020 \
          --end-year 2023 \
          --scale 1000 \
          --output ./outputs
```

The `--scale` parameter is in meters.

## Output Files

pyCropWat generates two GeoTIFF files per month:

| File | Description | Units |
|------|-------------|-------|
| `effective_precip_YYYY_MM.tif` | Effective precipitation | mm |
| `effective_precip_fraction_YYYY_MM.tif` | Effective/Total precipitation ratio | fraction (0-1) |

Example output directory:

```
outputs/
├── effective_precip_2020_01.tif
├── effective_precip_2020_02.tif
├── ...
├── effective_precip_fraction_2020_01.tif
├── effective_precip_fraction_2020_02.tif
└── ...
```

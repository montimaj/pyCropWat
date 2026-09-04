# Local Precipitation Examples

These examples drive pyCropWat with **your own precipitation rasters** instead of downloading precipitation from Google Earth Engine. They are built on real WRF regional climate model output for South America, but nothing here is specific to that dataset - swap in your own files with `--precip-dir` and, if needed, `--pattern` / `--nodata`.

The headline: **precipitation-only methods need no Earth Engine at all.** No `ee.Initialize()`, no authentication, no network.

---

## The Example Dataset

WRF (Weather Research and Forecasting) regional climate model monthly precipitation totals over the **Rio de la Plata basin domain**, South America.

| Attribute | Value |
|-----------|-------|
| **Files** | 264 GeoTIFFs, `Precip_YYYY_MM.tif` |
| **Period** | January 2000 - December 2021 (monthly) |
| **CRS** | EPSG:4326 |
| **Grid** | 689 rows x 799 cols |
| **Resolution** | 0.03974775585088471 deg (~4.4 km) |
| **Bands / dtype** | 1 band, float32 |
| **Nodata** | -9999.0 |
| **Units** | millimetres per month (no scaling needed) |
| **Bounds** | left `-71.30183990477435`, bottom `-40.350990224291216`, right `-39.54338297991747`, top `-12.964786443031647` |

### Where it lives

The data sits **next to** the repository, not inside it:

```
<parent>/
├── pyCropWat/                     # this repository
│   └── Examples/LocalPrecip/      # you are here: the three scripts and this README
└── pyCropWat_Data/
    └── Precip/                    # the 264 monthly GeoTIFFs
        ├── Precip_2000_01.tif
        ├── ...
        └── Precip_2021_12.tif
```

`Examples/LocalPrecip/` is two levels below the repository root, so the same directory has a different relative spelling depending on where you stand - and getting that wrong is the classic way for these commands to fail:

| Written relative to | Path to the data |
|---------------------|------------------|
| `Examples/LocalPrecip/` - **what every command on this page uses** | `../../../pyCropWat_Data/Precip` |
| `Examples/` | `../../pyCropWat_Data/Precip` |
| the repository root | `../pyCropWat_Data/Precip` |

(`../..` gets you to the repository root; the data is one level above *that*, hence three.)

The scripts sidestep the question entirely by resolving an absolute path from `__file__`, so their *defaults* are correct from any working directory:

```python
REPO_ROOT = Path(__file__).resolve().parents[2]          # .../pyCropWat
DEFAULT_PRECIP_DIR = REPO_ROOT.parent / 'pyCropWat_Data' / 'Precip'
```

If the directory is missing, each script exits with a short message telling you to pass `--precip-dir` - no traceback. Point `--precip-dir` at any directory of monthly rasters to run these examples on your own data.

### Where to run the commands on this page

**Everything below - scripts, CLI commands and Python snippets - is written to be run from `Examples/LocalPrecip/`,** the directory holding this README, and every relative path in them is relative to *that* directory:

```bash
cd pyCropWat/Examples/LocalPrecip
```

That is why the data is spelled `../../../pyCropWat_Data/Precip` here, and why the default output directories (`./LocalPrecip_Quickstart`, `./WRF_SA_Local`, `./LocalPrecip_NetCDF`) land beside the scripts. From anywhere else, pass an absolute `--precip-dir` / `--output` instead of adjusting the `../../..` prefix by hand. Those generated directories are `.gitignore`d.

The directory also contains a `.geeup-state.json` companion file. `LocalPrecipitationSource` ignores anything that is not a recognised raster/NetCDF suffix, so stray sidecars never break indexing.

---

## The Scripts

| Script | What it shows | Earth Engine | Runtime at defaults |
|--------|---------------|--------------|---------------------|
| `local_precip_quickstart.py` | Smallest working example: local GeoTIFFs + a precipitation-only method | **Not used** | ~3 s (2 months) |
| `wrf_south_america_local_example.py` | Flagship: local precipitation + GEE AWC/ETo, USDA-SCS then all 8 methods | Required (AWC + ETo) | ~30 s (24 months) |
| `local_netcdf_example.py` | GeoTIFFs -> NetCDF -> Peff, and proof the two paths agree | **Not used** | ~3 s (12 months, both paths) |

Those are wall-clock times for a cold run with no cached output (Apple M2 Max, data on a local SSD, `-w 4`). The two offline scripts are CPU- and disk-bound and should land in the same ballpark anywhere; **the WRF figure is dominated by the Earth Engine AWC/ETo download** (Step 1 is ~90% of it; Step 2's eight methods take about 2 s) and will move with your network and GEE quota. Re-running the WRF script without `-f` skips whatever is already on disk and finishes in a couple of seconds. Scaling is linear in months: the full 2000-2021 record is 11x the default period.

All scripts assume pyCropWat is importable - install it first (`pip install -e .` from the repository root, or `pip install pycropwat`).

---

## 1. Quickstart (`local_precip_quickstart.py`)

The smallest thing that works: read monthly GeoTIFFs, calculate effective precipitation, write rasters. Zero Earth Engine involvement.

```bash
# Defaults: 2005, months 7 and 8, method cropwat
python local_precip_quickstart.py

# A different precipitation-only method
python local_precip_quickstart.py --method farmwest

# A longer period, whole summer, 8 workers
python local_precip_quickstart.py --start-year 2005 --end-year 2010 \
                                  --months 6 7 8 -w 8

# Your own data
python local_precip_quickstart.py --precip-dir /data/my_precip \
                                  --pattern '*.tif' --nodata -9999 \
                                  --output ./my_output
```

| Option | Default | Description |
|--------|---------|-------------|
| `--precip-dir` | `../../../pyCropWat_Data/Precip` | Directory of monthly rasters |
| `--pattern` | `Precip_*.tif` | Glob applied inside `--precip-dir` |
| `--nodata` | `-9999` | Extra nodata sentinel, on top of the file metadata |
| `--start-year` / `--end-year` | `2005` / `2005` | Period to process |
| `--months` | `7 8` | Months to process (1-12) |
| `--method` | `cropwat` | One of `cropwat`, `fao_aglw`, `fixed_percentage`, `dependable_rainfall`, `farmwest` |
| `--output` | `./LocalPrecip_Quickstart` | Output directory |
| `-w`, `--workers` | `4` | Parallel workers; `1` runs sequentially |

It prints a summary of the source it opened, then the statistics of one output raster:

```
Local precipitation source
  Path          : /.../pyCropWat_Data/Precip
  Kind          : raster
  Files         : 264
  Months        : 264 (2000-01 to 2021-12)
  Year range    : 2000 - 2021
  CRS           : EPSG:4326
  Grid shape    : (689, 799) (rows, cols)
  Resolution    : 0.039748, 0.039748 (CRS units)
  Bounds (WGS84): -71.3018, -40.3510, -39.5434, -12.9648
...
Sample output: effective_precip_2005_07.tif
  Shape       : (689, 799)
  CRS         : EPSG:4326
  Valid pixels: 485202 of 550511 (65309 nodata/NaN)
  Peff min    : 0.000 mm
  Peff mean   : 39.240 mm
  Peff max    : 197.317 mm
```

### Output

```
LocalPrecip_Quickstart/
├── effective_precip_2005_07.tif           # Peff (mm)
├── effective_precip_2005_08.tif
├── effective_precip_fraction_2005_07.tif  # Peff / P (0-1)
└── effective_precip_fraction_2005_08.tif
```

---

## 2. WRF South America, local precipitation (`wrf_south_america_local_example.py`)

The local-file twin of [`../wrf/wrf_south_america_example.py`](../wrf/wrf_south_america_example.py). Same science, same AWC/ETo assets, same rooting depth and MAD factor - only the precipitation source changes, from the GEE asset `projects/azhydro/assets/WRF-ET-SA` to the monthly GeoTIFFs on disk.

| Input | Source | Earth Engine |
|-------|--------|--------------|
| Precipitation | Local `Precip_YYYY_MM.tif` | No |
| AWC | `projects/sat-io/open-datasets/FAO/HWSD_V2_SMU`, band `AWC`, scale 0.001 | **Yes** |
| ETo | `IDAHO_EPSCOR/TERRACLIMATE`, band `pet`, scale 0.1, monthly | **Yes** |
| Rooting depth / MAD | 2.0 m / 1.0 | - |

**Workflow**

1. **Step 1** - Run `EffectivePrecipitation` with `method='usda_scs'` and `save_inputs=True`. AWC and ETo are downloaded from GEE and reprojected onto the local precipitation grid (`reproject_match`), so `awc.tif` and `eto_YYYY_MM.tif` land on exactly the same transform as `precip_YYYY_MM.tif`.
2. **Step 2** - Calculate the other 7 methods locally from those saved rasters, one output directory per method. No further GEE calls.

Both steps skip work that is already on disk unless you pass `-f`.

```bash
# Default period 2000-2001 (24 months) so the first run is quick
python wrf_south_america_local_example.py

# The full WRF record: 2000-2021, 264 months
python wrf_south_america_local_example.py --start-year 2000 --end-year 2021

# More workers for the AWC/ETo downloads
python wrf_south_america_local_example.py -w 8

# Named GEE project
python wrf_south_america_local_example.py --gee-project your-project-id

# Clip to an area of interest: your own local vector file, or a GEE FeatureCollection asset
python wrf_south_america_local_example.py --geometry /path/to/basin.geojson
python wrf_south_america_local_example.py --geometry projects/my-project/assets/roi

# Only Step 1, only Step 2, or force a full redo
python wrf_south_america_local_example.py --download-only
python wrf_south_america_local_example.py --calc-only
python wrf_south_america_local_example.py -f
```

| Option | Default | Description |
|--------|---------|-------------|
| `--precip-dir` | `../../../pyCropWat_Data/Precip` | Directory of monthly rasters |
| `--pattern` | `Precip_*.tif` | Glob applied inside `--precip-dir` |
| `--nodata` | `-9999` | Extra nodata sentinel |
| `--geometry` | `None` | Local vector file **or** GEE FeatureCollection asset. Default: full raster extent |
| `--gee-project` | `None` | GEE project ID for the AWC/ETo downloads |
| `--start-year` / `--end-year` | `2000` / `2001` | Period to process (widen `--end-year` to 2021 for the full record) |
| `-w`, `--workers` | `4` | Parallel workers |
| `--download-only` / `--calc-only` | off | Run only Step 1 / only Step 2 |
| `-f`, `--force` | off | Recompute even if output exists |
| `--output` | `./WRF_SA_Local` | Output directory |

**A note on `--geometry`:** a *local* vector file (`.shp`, `.geojson`, `.json`) clips the local precipitation rasters (`clip_to_geometry=True` is the default) and is also used for the AWC/ETo requests. A *GEE asset ID* is used for the AWC/ETo requests only - local rasters are not clipped to it. With no `--geometry`, the full extent of the local rasters is used.

Wherever the precipitation grid is left wider than the AWC/ETo download region, the pixels the download never reached have no AWC or ETo to pair with. They stay **NaN** in the `usda_scs`, `ensemble` and `suet` outputs - nothing is invented to fill them - and the run logs a WARNING naming the covered percentage and the uncovered pixel count. It is logged once per field per process, so `usda_scs` (and `ensemble`) log it twice - once for AWC, once for ETo:

```
AWC covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
ETo covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
```

Clip to the download region, or widen the geometry to span the whole grid, for full coverage.

### Output

```
WRF_SA_Local/
├── usda_scs/                              # Step 1
│   ├── effective_precip_YYYY_MM.tif
│   └── effective_precip_fraction_YYYY_MM.tif
├── analysis_inputs/                       # Step 1, save_inputs=True
│   ├── precip_YYYY_MM.tif                 #   the local precipitation actually used
│   ├── awc.tif                            #   FAO HWSD, regridded from GEE
│   └── eto_YYYY_MM.tif                    #   TerraClimate, regridded from GEE
└── peff_by_method/                        # Step 2 (8 methods)
    ├── usda_scs/effective_precip_YYYY_MM.tif
    ├── ensemble/effective_precip_YYYY_MM.tif
    ├── cropwat/effective_precip_YYYY_MM.tif
    ├── fao_aglw/effective_precip_YYYY_MM.tif
    ├── fixed_percentage/effective_precip_YYYY_MM.tif
    ├── dependable_rainfall/effective_precip_YYYY_MM.tif
    ├── farmwest/effective_precip_YYYY_MM.tif
    └── suet/effective_precip_YYYY_MM.tif
```

---

## 3. NetCDF input (`local_netcdf_example.py`)

Local precipitation does not have to be GeoTIFF. This script builds a NetCDF from the monthly rasters, runs pyCropWat on it, and then proves the NetCDF and GeoTIFF paths give the same answer.

1. **Step 1** - Stack twelve monthly GeoTIFFs on a `time` dimension with xarray, name the variable `precip`, tag it `units = 'mm'` with a CF-ish attribute set, and write it compressed.
2. **Step 2** - Open that NetCDF through `LocalPrecipitationSource` / `EffectivePrecipitation` and compute Peff.
3. **Step 3** - Compute the same months straight from the GeoTIFFs and assert the two agree to within floating point tolerance (1e-4 mm), including matching NaN masks.

```bash
# Defaults: build precip_2005.nc, compute all 12 months, method cropwat
python local_netcdf_example.py

# Two months only - fastest way to see the whole round trip
python local_netcdf_example.py --months 7 8

# Another year / method
python local_netcdf_example.py --year 2010 --method fao_aglw

# Reuse the NetCDF from a previous run
python local_netcdf_example.py --skip-build
```

| Option | Default | Description |
|--------|---------|-------------|
| `--precip-dir` | `../../../pyCropWat_Data/Precip` | Directory of monthly GeoTIFFs to stack |
| `--year` | `2005` | Year to stack into the NetCDF (12 months) |
| `--output` | `./LocalPrecip_NetCDF` | Output directory |
| `--method` | `cropwat` | Precipitation-only Peff method |
| `--nodata` | `-9999` | Nodata sentinel for both the GeoTIFFs and the NetCDF `_FillValue` |
| `--skip-build` | off | Reuse an existing NetCDF instead of rebuilding it |
| `--months` | all 12 | Months to compute and compare (the NetCDF always holds the full year) |

Expected tail of a run:

```
  Month | Peff max |diff| | Fraction max |diff| | NaN match
2005-07 |       0.000e+00 |           0.000e+00 | yes
2005-08 |       0.000e+00 |           0.000e+00 | yes
Assertion passed: NetCDF and GeoTIFF results agree for 2 month(s) within 0.0001 mm.
```

**Compression:** the script writes NetCDF4 with zlib when `netCDF4` or `h5netcdf` is installed. With neither, it falls back to xarray's built-in `scipy` backend, which writes classic uncompressed NetCDF3 and logs a warning - the example still runs, the file is just larger. `pip install netCDF4` to get compression.

**Writing your own NetCDF?** One gotcha the script guards against: `to_netcdf(encoding=...)` **replaces** a variable's encoding, so a hand-built encoding dict silently drops the `grid_mapping` link that `rio.write_crs()` set, and pyCropWat then falls back to "assuming EPSG:4326". Start from `dict(ds[var].encoding)` and update it, or pass `--local-precip-crs` when reading.

### Output

```
LocalPrecip_NetCDF/
├── precip_2005.nc          # 12 months, variable 'precip', dims (time, y, x)
├── from_netcdf/
│   ├── effective_precip_2005_MM.tif
│   └── effective_precip_fraction_2005_MM.tif
└── from_geotiff/
    ├── effective_precip_2005_MM.tif
    └── effective_precip_fraction_2005_MM.tif
```

---

## What needs Earth Engine, and what does not

Precipitation always comes from disk in these examples. Whether Earth Engine is touched depends entirely on the **method**, because some methods need AWC and/or ETo, and those are still GEE datasets.

| Method | Extra inputs | Earth Engine required |
|--------|--------------|-----------------------|
| `cropwat` | - | No |
| `fao_aglw` | - | No |
| `fixed_percentage` | - | No |
| `dependable_rainfall` | - | No |
| `farmwest` | - | No |
| `suet` | ETo | **Yes** |
| `usda_scs` | AWC + ETo | **Yes** |
| `ensemble` | AWC + ETo | **Yes** |
| `pcml` | - | Incompatible with local precipitation (pre-computed Western U.S. GEE product) |

A GEE `FeatureCollection` geometry also requires Earth Engine, whichever method you choose.

When nothing needs it, pyCropWat says so and never calls `ee.Initialize()`:

```
INFO - Earth Engine not required: local precipitation with method 'cropwat' needs no GEE assets. Skipping initialization.
```

---

## Command line equivalents

Two of the three scripts have a faithful `pycropwat` CLI equivalent. **The third - the USDA-SCS
workflow - does not, and cannot;** see the warning below rather than assembling one yourself. All of
these are run from `Examples/LocalPrecip/`.

### 1. Quickstart, on the command line

The CLI form of `local_precip_quickstart.py` at its defaults (2005, months 7 and 8, `cropwat`):

```bash
pycropwat process --local-precip ../../../pyCropWat_Data/Precip \
                  --local-precip-pattern 'Precip_*.tif' \
                  --local-precip-nodata -9999 \
                  --method cropwat \
                  --start-year 2005 --end-year 2005 --months 7 8 \
                  --output ./LocalPrecip_Quickstart_cli
```

`--start-year`/`--end-year` are what make it *equivalent*. They are optional with `--local-precip` -
omit them and the range is inferred from the files, which for this dataset is 2000-2021, so
`--months 7 8` would process **44 months instead of 2**. Handy for a full-record run, surprising if
you expected the quickstart.

### 2. NetCDF input, on the command line

The CLI form of Step 2 of `local_netcdf_example.py`. Run the script once first so the NetCDF exists:

```bash
python local_netcdf_example.py --months 7 8      # writes LocalPrecip_NetCDF/precip_2005.nc

pycropwat process --local-precip ./LocalPrecip_NetCDF/precip_2005.nc \
                  --local-precip-variable precip \
                  --local-precip-nodata -9999 \
                  --method cropwat \
                  --start-year 2005 --end-year 2005 --months 7 8 \
                  --output ./LocalPrecip_NetCDF/from_cli
```

Step 1 (building the NetCDF) and Step 3 (the GeoTIFF-vs-NetCDF assertion) are the script's own code;
the CLI only does the middle part.

### 3. WRF South America - Step 1 on the CLI

Step 1 of `wrf_south_america_local_example.py` (the USDA-SCS download) maps onto `pycropwat process`.
Run this from `Examples/LocalPrecip/`:

```bash
pycropwat process --local-precip ../../../pyCropWat_Data/Precip \
                  --local-precip-pattern 'Precip_*.tif' \
                  --local-precip-nodata -9999 \
                  --start-year 2000 --end-year 2001 \
                  --method usda_scs \
                  --awc-asset projects/sat-io/open-datasets/FAO/HWSD_V2_SMU \
                  --awc-band AWC --awc-scale-factor 0.001 \
                  --eto-asset IDAHO_EPSCOR/TERRACLIMATE \
                  --eto-band pet --eto-scale-factor 0.1 \
                  --rooting-depth 2.0 --mad-factor 1.0 \
                  --workers 4 \
                  --output ./WRF_SA_Local/usda_scs
```

> [!WARNING]
> **`--eto-scale-factor 0.1` is not optional here.** TerraClimate stores `pet` in units of 0.1 mm, so
> omitting it feeds USDA-SCS an ETo 10x too large. Measured on 2000-01 over this domain:
>
> - with `--eto-scale-factor 0.1`: Peff **mean 74.791 mm, max 219.300 mm**
> - without it: Peff **mean 127.105 mm, max 1169.000 mm**
>
> That is 1.70x the mean and 5.3x the maximum - rasters that look entirely plausible and are wrong.
> AWC has the same trap: HWSD is in mm/m, hence `--awc-scale-factor 0.001`.

Step 2 (computing the other 7 methods from the saved rasters) has no CLI equivalent - it is the
script's own code. The Python form of Step 1, which is what the script actually runs:

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    local_precip='../../../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    start_year=2000,
    end_year=2001,
    method='usda_scs',
    method_params={
        'awc_asset': 'projects/sat-io/open-datasets/FAO/HWSD_V2_SMU',
        'awc_band': 'AWC',
        'awc_scale_factor': 0.001,             # HWSD AWC is mm/m -> volumetric fraction
        'eto_asset': 'IDAHO_EPSCOR/TERRACLIMATE',
        'eto_band': 'pet',
        'eto_scale_factor': 0.1,               # TerraClimate 'pet' is 0.1 mm (CLI: --eto-scale-factor)
        'eto_is_daily': False,
        'rooting_depth': 2.0,
        'mad_factor': 1.0,
    },
)
ep.process(
    output_dir='./WRF_SA_Local/usda_scs',
    n_workers=4,
    save_inputs=True,
    input_dir='./WRF_SA_Local/analysis_inputs',
)
```

The CLI is fine for `usda_scs`/`suet`/`ensemble` when the ETo asset is **already in mm per month** -
GridMET `eto`, for instance - because then there is no ETo scaling to lose.

---

## Bring your own data

`--local-precip` (CLI) / `local_precip=` (API) accepts four layouts:

| Layout | Pass | Notes |
|--------|------|-------|
| **Directory of GeoTIFFs** | the directory, plus `--local-precip-pattern` | One file per month, dated by file name. Also accepts `.vrt`, `.img`, `.asc`, `.jp2` and friends |
| **Single NetCDF** | the `.nc` file | Needs a `time` coordinate covering the months, or a datable file name |
| **Multi-file NetCDF** | a directory with `--local-precip-pattern '*.nc'`, or a glob | Combined with `open_mfdataset(combine='by_coords')` when the files carry a time axis; otherwise each file is one month, dated by name |
| **Glob string** | e.g. `'/data/precip/**/Precip_*.tif'` | Expanded recursively |

Files whose suffix is not a recognised raster/NetCDF type are ignored, so sidecars (`.aux.xml`, `.json`, `.csv`) are harmless.

### File naming that just works

`parse_year_month()` reads the year and month out of the file name. Layouts are tried most-specific-first and the **last** match in the string wins, so prefixes and directory names never confuse it.

| Layout | Example | Parsed |
|--------|---------|--------|
| `YYYY_MM` | `Precip_2005_07.tif` | (2005, 7) |
| `YYYY-MM` | `precip-2005-07.nc` | (2005, 7) |
| `YYYYMM` | `pr200507.tif` | (2005, 7) |
| `YYYY.MM` | `pr.2005.07.tif` | (2005, 7) |

### When it does not

| Situation | Fix |
|-----------|-----|
| Names in another layout, e.g. `rain_07_2005.tif` | `--local-precip-date-regex '(?P<month>[0-9]{2})_(?P<year>[0-9]{4})'` - any regex with named groups `year` and `month` |
| NetCDF with several data variables, or an unusual variable name | `--local-precip-variable RAINNC` (auto-detection covers `precip`, `precipitation`, `pr`, `prcp`, `ppt`, `tp`, `rain`, `RAINNC`, ...) |
| Files carry no CRS, or the stored one is wrong | `--local-precip-crs EPSG:4326` (any CRS `rasterio` understands). It is an **override**: it replaces a CRS the file already declares, not just a missing one, and says so at INFO - `Precip_2000_01.tif declares CRS EPSG:32719; overriding it with the requested EPSG:4326 (relabelled only - no reprojection)`. Nothing is reprojected: the grid is relabelled where it already sits |
| Extra nodata sentinel not in the file metadata | `--local-precip-nodata -9999` (applied *in addition to* the file's own nodata) |
| Units are not mm per month | `--scale-factor` / `precip_scale_factor=` (e.g. `1000` for metres, `0.1` for mm x 10) |
| One 12-band raster per year | Not supported - split into one file per month. pyCropWat detects this layout and says so |

Missing months are simply skipped, with a log line naming them; a requested year range is clamped to what actually exists on disk, and omitting the range altogether infers it from the files.

### Minimal API version

Runnable as-is from `Examples/LocalPrecip/`; swap the path and pattern for your own data.

```python
from pycropwat import EffectivePrecipitation, open_local_precipitation

# Inspect first
with open_local_precipitation('../../../pyCropWat_Data/Precip',
                              pattern='Precip_*.tif', nodata=-9999) as src:
    print(src.kind, len(src), src.year_range, src.crs, src.shape)

# Then run - no Earth Engine for a precipitation-only method
ep = EffectivePrecipitation(
    local_precip='../../../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    start_year=2005,
    end_year=2005,
    method='cropwat',
)
ep.process(output_dir='./output', n_workers=4, months=[7, 8])
```

---

## See Also

- [Examples README](../README.md) - the GEE-based example workflows
- [`../wrf/wrf_south_america_example.py`](../wrf/wrf_south_america_example.py) - the GEE twin of the flagship script
- [pyCropWat Documentation](https://montimaj.github.io/pyCropWat/)
- [API Reference](https://montimaj.github.io/pyCropWat/api/core/)

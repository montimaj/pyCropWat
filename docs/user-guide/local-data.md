# Local Precipitation Data

pyCropWat can read monthly precipitation from **your own files** instead of downloading it from
Google Earth Engine. Point `local_precip` (Python) or `--local-precip` (CLI) at a directory of
monthly GeoTIFFs, a NetCDF file, or a glob, and every effective precipitation method that only
needs precipitation runs entirely on your data - with no Earth Engine account, no authentication
and no network access.

Everything else keeps working exactly as before: available water capacity (AWC), reference
evapotranspiration (ETo) and GEE `FeatureCollection` geometries are still read from Earth Engine,
so you can mix a locally-produced precipitation field with global GEE ancillary data.

## Why bring your own precipitation

| Situation | What local input gives you |
|-----------|----------------------------|
| Regional climate model output (WRF, RegCM, MPAS, a CORDEX RCM) | Run CROPWAT on the exact model grid, at the model's native resolution and projection |
| Gauge-interpolated or national products (e.g. a met-service gridded analysis) | Use the authoritative dataset for your country instead of a global reanalysis |
| Bias-corrected or downscaled precipitation | Feed the corrected field straight into the effective precipitation calculation |
| Climate scenarios / ensemble members | Process each member from disk; no upload to Earth Engine needed |
| Restricted or air-gapped environments | Precipitation-only methods never call Earth Engine, so the run works fully offline |
| Reproducibility | The precipitation input is a file you control and can archive with the results |

## Supported layouts

| Layout | `local_precip` value | Notes |
|--------|----------------------|-------|
| Directory of monthly rasters | `'./Precip'` | One file per month, dated by file name. Globbed with `local_precip_pattern` (default `'*.tif'`) |
| Single NetCDF with a time axis | `'./precip_2000_2021.nc'` | Months come from the `time` coordinate |
| Several NetCDF files | `'./nc_dir'` with `local_precip_pattern='*.nc'` | Combined with `xarray.open_mfdataset(..., combine='by_coords')` |
| Explicit glob | `'./Precip/Precip_*.tif'` or `'./nc/*.nc'` | Anything `glob.glob` accepts, including `**` |
| Single raster | `'./Precip/Precip_2005_07.tif'` | A one-month run |

Recognised raster suffixes are `.tif`, `.tiff`, `.vrt`, `.img`, `.bil`, `.bsq`, `.bip`, `.asc`,
`.jp2`, `.grd`, `.dat`, `.hgt`; NetCDF/HDF suffixes are `.nc`, `.nc4`, `.cdf`, `.netcdf`, `.h5`,
`.hdf5`. Anything else in the directory - `.aux.xml` sidecars, `metadata.csv`, a stray
`.geeup-state.json` - is ignored, so a broad pattern such as `'*'` is safe.

**Any CRS and any resolution work.** The precipitation grid is never resampled: it is read as-is,
and effective precipitation is written on the same grid, in the same CRS, with the same transform.
When AWC or ETo are needed they are downloaded from Earth Engine and reprojected *onto* your grid
with `reproject_match`, so all inputs land on identical pixels. Wherever that download does not
reach - a precipitation grid wider than the geometry AWC/ETo were requested for - the pixels stay
NaN rather than being filled with a stand-in value; see
[Geometry and clipping](#geometry-and-clipping).

!!! warning "Two layouts that are not supported"
    * **One multi-band raster per year** (12 bands = 12 months) is rejected with a clear error.
      Split it into one file per month first.
    * **Mixing rasters and NetCDF** in the same selection raises an error. Narrow the selection
      with `local_precip_pattern` (`'*.tif'` or `'*.nc'`).

## Quick start

**Python:**

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    method='cropwat',
    start_year=2005,
    end_year=2005
)

results = ep.process(output_dir='./outputs', n_workers=4, months=[7, 8])
```

**CLI:**

```bash
pycropwat process \
    --local-precip ../pyCropWat_Data/Precip \
    --local-precip-pattern 'Precip_*.tif' \
    --local-precip-nodata -9999 \
    --method cropwat \
    --start-year 2005 --end-year 2005 \
    --months 7 8 \
    --output ./outputs
```

Both produce the usual pair of rasters per month:

```
outputs/
├── effective_precip_2005_07.tif
├── effective_precip_2005_08.tif
├── effective_precip_fraction_2005_07.tif
└── effective_precip_fraction_2005_08.tif
```

`--asset`/`--band` (`asset_id`/`precip_band`) are not required when `--local-precip` is given, and
neither are `--start-year`/`--end-year` - they are inferred from the files.

## Do I need Earth Engine?

This is the whole matrix. `cropwat`, `fao_aglw`, `fixed_percentage`, `dependable_rainfall` and
`farmwest` are pure functions of precipitation; `usda_scs` needs AWC **and** ETo, `suet` needs ETo,
and `ensemble` is the mean of six methods including USDA-SCS, so it needs both.

| Precipitation | Method | Geometry | Earth Engine? |
|---------------|--------|----------|---------------|
| Local | `cropwat`, `fao_aglw`, `fixed_percentage`, `dependable_rainfall`, `farmwest` | none, or a local `.shp`/`.geojson` | **No** - runs fully offline |
| Local | `cropwat`, `fao_aglw`, `fixed_percentage`, `dependable_rainfall`, `farmwest` | GEE `FeatureCollection` asset | **Yes** - for the geometry only |
| Local | `suet` | any | **Yes** - ETo is downloaded from GEE |
| Local | `usda_scs`, `ensemble` | any | **Yes** - AWC and ETo are downloaded from GEE |
| Local | `pcml` | any | **Not supported** - raises `ValueError` |
| GEE (`asset_id` + `precip_band`) | any | required (except `pcml`, which brings its own) | **Yes** |

When Earth Engine is not required, `ee.Initialize()` is never called, `ep.geometry` and
`ep.collection` are `None`, and the log line reads:

```
Earth Engine not required: local precipitation with method 'cropwat' needs no GEE assets. Skipping initialization.
```

!!! note "PCML cannot use local precipitation"
    `pcml` is a pre-computed Earth Engine product for the Western U.S., not a formula applied to
    precipitation. Combining it with local files raises
    `ValueError: local_precip cannot be combined with method='pcml'...`, and the CLI exits with
    `--local-precip cannot be used with --method pcml`.

## File naming and dating

In raster mode (and for NetCDF files without a time coordinate) the month comes from the file name.
Layouts are tried most-specific-first and the **last** match in the string wins, so prefixes such
as `Precip_` and directory names such as `/data/2019/` never confuse the parser.

| Layout | Example file name | Parsed as |
|--------|-------------------|-----------|
| `YYYY_MM` | `Precip_2005_07.tif` | `(2005, 7)` |
| `YYYY-MM` | `precip-2005-07.nc` | `(2005, 7)` |
| `YYYYMM` | `pr200507.tif` | `(2005, 7)` |
| `YYYY.MM` | `pr.2005.07.tif` | `(2005, 7)` |

Months may be one or two digits (`2005_7` and `2005_07` both work), years must fall between 1700
and 2200, and the month must be 1-12.

Check a name before you commit to a run:

```python
from pycropwat import parse_year_month

parse_year_month('Precip_2005_07.tif')       # (2005, 7)
parse_year_month('/data/2019/x-2005-07.nc')  # (2005, 7)
parse_year_month('rain_july.tif')            # None
```

### Anything else: a custom regex

Supply a regular expression with named groups `year` and `month`:

**Python:**

```python
ep = EffectivePrecipitation(
    local_precip='./rain',
    local_precip_pattern='rain_*.tif',       # rain_07_2005.tif  (MM_YYYY)
    local_precip_date_regex=r'(?P<month>\d{2})_(?P<year>\d{4})',
    method='cropwat'
)
```

**CLI:**

```bash
pycropwat process \
    --local-precip ./rain \
    --local-precip-pattern 'rain_*.tif' \
    --local-precip-date-regex '(?P<month>[0-9]{2})_(?P<year>[0-9]{4})' \
    --method cropwat --output ./outputs
```

!!! tip "Quote the regex in the shell"
    Wrap the pattern in single quotes so the shell does not eat `(`, `)` or `?`, and prefer
    `[0-9]` over `\d` on the command line to avoid backslash mangling. In Python use a raw
    string (`r'...'`).

When a custom regex is given, the built-in layouts are not tried. A regex missing the named groups
raises `ValueError: date_regex must define named groups 'year' and 'month'...`.

Duplicate months are not an error: the first file (in sorted order) wins and the other is logged as
`Duplicate month 2005-07: keeping Precip_2005_07.tif, ignoring rain_2005_07.tif`. Files that cannot
be dated are skipped with a warning; if *nothing* can be dated, the run fails immediately with a
message telling you to pass a `date_regex`.

## NetCDF specifics

### Variable selection

The precipitation variable is auto-detected by looking (case-insensitively) for
`precip`, `precipitation`, `pr`, `prcp`, `ppt`, `tp`, `total_precipitation`, `rain`, `rainfall`,
`RAINNC`, `PRECIP`. If none of those exist but the file has exactly one 2-D-or-higher data
variable, that one is used. Otherwise you get an error listing the candidates.

```python
ep = EffectivePrecipitation(
    local_precip='./wrf_precip.nc',
    local_precip_variable='RAINNC',   # skip auto-detection
    local_precip_crs='EPSG:4326',
    method='fao_aglw'
)
```

```bash
pycropwat process --local-precip ./wrf_precip.nc \
    --local-precip-variable RAINNC \
    --local-precip-crs EPSG:4326 \
    --method fao_aglw --output ./outputs
```

### Dimensions

X is detected from `x`, `lon`, `longitude`, `west_east`, `XLONG`; Y from `y`, `lat`, `latitude`,
`south_north`, `XLAT`; time from `time`, `Time`, `valid_time`, `t`. Matching is
case-insensitive, 2-D curvilinear `XLAT`/`XLONG` coordinates are collapsed to 1-D when possible,
and singleton dimensions (`band`, `level`, a length-1 `time`) are dropped.

For dimension names outside those lists, use
[`LocalPrecipitationSource`](../api/local.md) directly, which exposes `x_dim`, `y_dim` and
`time_dim`:

```python
from pycropwat import LocalPrecipitationSource

with LocalPrecipitationSource('./odd_grid.nc',
                              variable='RAINNC',
                              x_dim='easting',
                              y_dim='northing',
                              time_dim='valid_time',
                              crs='EPSG:32721') as src:
    print(src.available_months()[:3])
    da = src.get_month(2005, 7)
```

!!! note "`x_dim`, `y_dim` and `time_dim` are class-level options"
    They are not exposed on `EffectivePrecipitation` or the CLI. If auto-detection fails for a
    processing run, rename the dimensions once and write a normalised NetCDF:

    ```python
    import xarray as xr

    ds = xr.open_dataset('odd_grid.nc')
    ds = ds.rename({'easting': 'x', 'northing': 'y', 'valid_time': 'time'})
    ds.to_netcdf('normalised.nc')
    ```

### Time and `_FillValue`

Times are decoded by xarray (`numpy.datetime64`, `datetime` and `cftime` calendars are all
handled) and reduced to `(year, month)` pairs - a file with daily or sub-monthly steps inside a
month is not aggregated for you, so supply monthly totals. If the time dimension has **no
coordinate values**, pyCropWat falls back to dating each NetCDF file by its name and logs
`No time coordinate found in ...; dating NetCDF files by name`.

`_FillValue` and `missing_value`, whether xarray already applied them or left them in `attrs`, are
masked to `NaN` automatically. `local_precip_nodata` adds a sentinel on top of them.

Latitude ascending? No action needed - arrays are flipped to north-up (descending `y`) before they
are written.

## Units, scaling and nodata

`precip_scale_factor` (`--scale-factor`) is a plain multiplier applied to every value read from
disk. Values must end up as **millimetres per month**.

| Source units | `precip_scale_factor` |
|--------------|-----------------------|
| mm per month | `1.0` (default) |
| metres per month | `1000` |
| cm per month | `10` |
| inches per month | `25.4` |

```python
ep = EffectivePrecipitation(
    local_precip='./era5_local.nc',
    local_precip_variable='tp',
    precip_scale_factor=1000,   # metres -> mm
    method='cropwat'
)
```

!!! warning "Rate units need pre-processing"
    A single multiplier cannot convert a rate such as `kg m-2 s-1` to a monthly total, because the
    factor depends on the number of days in each month. Accumulate to monthly totals before
    handing the data to pyCropWat.

### Nodata

Three sources of nodata are honoured, in this order:

1. the nodata value stored in the raster metadata (or `_FillValue`/`missing_value` in NetCDF),
2. any extra sentinel you pass in `local_precip_nodata` / `--local-precip-nodata`,
3. `NaN` already present in the data.

Everything masked becomes `NaN`, and **`NaN` propagates to both outputs** - the effective
precipitation raster *and* the fraction raster - which are written with `nodata = NaN`. Pixels
never silently become zero.

```python
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,   # sentinel not declared in the file metadata
    method='cropwat'
)
```

## CRS

Rasters and NetCDF files that carry a CRS are used as-is, in their own CRS. `local_precip_crs`
supplies one when the file declares none, and **overrides** the declared one when a file carries a
CRS that is wrong or unusable:

```bash
pycropwat process --local-precip ./precip_nc --local-precip-pattern '*.nc' \
    --local-precip-crs EPSG:4326 --method cropwat --output ./outputs
```

!!! warning "`local_precip_crs` relabels, it never reprojects"
    `local_precip_crs` is an **override**, not a fallback. Whatever the file declares, the grid is
    relabelled with the CRS you pass - a file that already carries one does **not** keep it.
    Replacing a declared CRS is logged at INFO, naming the file, so it is never silent:

    ```
    Precip_2000_01.tif declares CRS EPSG:32719; overriding it with the requested EPSG:4326 (relabelled only - no reprojection)
    ```

    A file that declares nothing logs the shorter
    `Precip_2000_01.tif declares no CRS; assigning EPSG:4326` instead.

    Relabelling changes the label only. Pixel values and the geotransform are left exactly as they
    are, so nothing is resampled and no coordinate moves. Pass a CRS that the grid's existing
    coordinates really are in: naming one whose units differ from the file's (`EPSG:4326` for a
    grid whose transform is in metres, say) silently puts the data in the wrong place. To genuinely
    move a grid between projections, reproject the files before pyCropWat reads them - with
    `rioxarray`'s `reproject()` or `gdalwarp` - and leave `local_precip_crs` unset.

    If nothing declares a CRS and none is supplied, EPSG:4326 is assumed and a warning is logged:

    ```
    Precip_2000_01.tif declares no CRS and none was supplied; assuming EPSG:4326. Pass local_precip_crs=... if this is wrong.
    ```

The outputs inherit the precipitation grid's CRS, so a projected local grid (UTM, Lambert
conformal, a WRF Lambert grid) round-trips without ever passing through EPSG:4326.

## Geometry and clipping

Geometry is **optional** with local precipitation. There are three cases.

=== "No geometry"

    The extent of the local files is used. Nothing is clipped and, for a precipitation-only
    method, Earth Engine is never touched.

    ```python
    ep = EffectivePrecipitation(
        local_precip='../pyCropWat_Data/Precip',
        local_precip_pattern='Precip_*.tif',
        local_precip_nodata=-9999,
        method='cropwat'
    )
    print(ep.bounds)   # ring built from the files' WGS84 bounds
    ```

    ```
    No geometry provided; using local precipitation extent (-71.3018, -40.3510, -39.5434, -12.9648)
    ```

=== "Local vector file"

    A `.shp`/`.geojson` is read with GeoPandas, reprojected to the raster CRS and used to clip
    each month (`all_touched=True`, empty rows/columns dropped). Still no Earth Engine for a
    precipitation-only method.

    ```python
    ep = EffectivePrecipitation(
        local_precip='../pyCropWat_Data/Precip',
        local_precip_pattern='Precip_*.tif',
        local_precip_nodata=-9999,
        geometry_path='basin.geojson',
        method='farmwest'
    )
    ```

    ```bash
    pycropwat process --local-precip ../pyCropWat_Data/Precip \
        --local-precip-pattern 'Precip_*.tif' --local-precip-nodata -9999 \
        --geometry basin.geojson --method farmwest --output ./outputs
    ```

=== "GEE FeatureCollection"

    A `--gee-geometry` asset (or a `geometry_path` that is an asset ID) is loaded from Earth
    Engine, so **this forces Earth Engine initialization** even for a precipitation-only method.
    Local rasters are *not* clipped in this case - the asset defines the region used for the
    AWC/ETo downloads.

    ```python
    ep = EffectivePrecipitation(
        local_precip='../pyCropWat_Data/Precip',
        local_precip_pattern='Precip_*.tif',
        gee_geometry_asset='projects/my-project/assets/basin',
        gee_project='my-gee-project',
        method='cropwat'
    )
    ```

    !!! warning "Uncovered pixels become NaN for AWC/ETo methods"
        Because the local rasters keep their full extent while AWC and ETo are downloaded only
        for the asset's region, any part of the precipitation grid outside that region is a
        pixel pyCropWat has no AWC or ETo for. Those pixels stay **NaN** and the NaN propagates
        into `usda_scs`, `ensemble` and `suet` outputs. The run logs a WARNING naming the
        covered percentage and the uncovered pixel count, once per field per process (twice for
        `usda_scs` and `ensemble`: once for AWC, once for ETo):

        ```
        AWC covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
        ```

        A GEE `FeatureCollection` geometry never clips the local rasters, so
        `clip_to_geometry=True` cannot help here: use a FeatureCollection that covers the whole
        grid if you need full coverage. The precipitation-only methods (`cropwat`, `fao_aglw`,
        `fixed_percentage`, `dependable_rainfall`, `farmwest`) download nothing and are
        unaffected.

To keep the full grid while still using a vector for the AWC/ETo region, set
`clip_to_geometry=False`:

```python
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    geometry_path='basin.geojson',
    clip_to_geometry=False,      # keep the full 689 x 799 grid
    method='cropwat'
)
```

!!! warning "`clip_to_geometry=False` leaves AWC/ETo methods NaN outside the geometry"
    This keeps the precipitation grid whole, but AWC and ETo are still downloaded only for
    `basin.geojson`. For `usda_scs`, `ensemble` and `suet` the pixels outside the basin have no
    AWC/ETo to pair with the precipitation, so they stay **NaN** in the outputs rather than being
    filled in - pyCropWat will not invent values for a region it never downloaded. A WARNING names
    how much of the grid is uncovered, once per field per process (twice for `usda_scs` and
    `ensemble`: once for AWC, once for ETo):

    ```
    ETo covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
    ```

    The example above uses `method='cropwat'`, a precipitation-only method that downloads no AWC
    or ETo at all, so it is unaffected. If you need an AWC/ETo method to cover the whole grid, get
    the download region to cover it: leave `clip_to_geometry=True` (the default) so the outputs
    are clipped to the region that was actually downloaded, or pass a `geometry_path` wide enough
    to span the full precipitation extent.

`clip_to_geometry` has no CLI flag - the CLI always uses the default `True`, so passing
`--geometry` with `--local-precip` clips. It also has no effect on GEE precipitation, which is
always clipped server-side, nor when the geometry is a GEE asset. If clipping fails (for example
a geometry that does not overlap the rasters) a warning is logged and the unclipped array is used.

## Years and missing months

| Situation | Behaviour |
|-----------|-----------|
| `start_year`/`end_year` omitted | Inferred from the files: `start_year not given; inferred 2000 from local files` |
| Range wider than the data | Clamped, with `Requested years 1995-2030 clamped to 2000-2021 (local data covers 2000-2021)` |
| Range with no overlap | `ValueError: Requested years 1990-1995 do not overlap the local precipitation data (2000-2021).` |
| Individual months missing | Skipped, not failed: `No local precipitation file for 2 requested month(s): 2005-03, 2005-04` |

Skipped months are dropped from the work list entirely, so `len(results)` reflects the months that
actually existed rather than a run of `(None, None)` entries.

```python
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    method='cropwat'
)
print(ep.start_year, ep.end_year)   # 2000 2021
```

## Inspecting a dataset before processing

`LocalPrecipitationSource` is the reader underneath, and it is worth a look before a long run.
It indexes file names up front (cheap) and only reads pixels in `get_month()`.

```python
from pycropwat import open_local_precipitation

src = open_local_precipitation('../pyCropWat_Data/Precip',
                               pattern='Precip_*.tif',
                               nodata=-9999)

print(src)
# LocalPrecipitationSource(kind='raster', files=264, months=264, years=2000-2021,
#                          crs='EPSG:4326', shape=(689, 799))

print(src.kind)                    # 'raster'
print(len(src))                    # 264
print(src.year_range)              # (2000, 2021)
print(src.crs)                     # EPSG:4326
print(src.shape)                   # (689, 799)
print(src.resolution)              # (0.03974775585088471, 0.03974775585088471)
print(src.bounds)                  # WGS84 (minx, miny, maxx, maxy)
print(src.native_bounds)           # same tuple in the source CRS
print(src.available_months()[:2])  # [(2000, 1), (2000, 2)]
print(src.has_month(2005, 7))      # True

da = src.get_month(2005, 7)        # xr.DataArray ('y', 'x'), float32, mm
print(da.dims, da.shape)           # ('y', 'x') (689, 799)
print(da.attrs['units'])           # 'mm'

src.close()
```

It is also a context manager:

```python
with open_local_precipitation('./precip_2000_2021.nc', variable='precip') as src:
    monthly_means = {(y, m): float(src.get_month(y, m).mean())
                     for y, m in src.available_months()[:12]}
```

Every array returned by `get_month()` has dims exactly `('y', 'x')`, dtype `float32`, values in mm,
nodata as `NaN`, a descending (north-up) `y` coordinate, the source CRS attached via `rioxarray`,
and attrs `units`, `long_name`, `year`, `month`, `source`.

## Worked example: WRF over South America

The bundled example uses monthly precipitation from a WRF regional climate simulation of the Rio de
la Plata basin:

| Property | Value |
|----------|-------|
| Files | 264 GeoTIFFs, `Precip_YYYY_MM.tif` |
| Period | 2000-01 to 2021-12 |
| CRS | EPSG:4326 |
| Grid | 689 rows x 799 columns |
| Resolution | 0.0397478° (~4.4 km) |
| Units | mm per month (`precip_scale_factor=1.0`) |
| Nodata | `-9999` |
| Extent | -71.302, -40.351 to -39.543, -12.965 |

### Offline: CROPWAT, no Earth Engine

**Python:**

```python
from pycropwat import EffectivePrecipitation

ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    method='cropwat',
    start_year=2005,
    end_year=2005
)

results = ep.process(output_dir='./wrf_outputs', n_workers=4, months=[7, 8])
```

**CLI:**

```bash
pycropwat process \
    --local-precip ../pyCropWat_Data/Precip \
    --local-precip-pattern 'Precip_*.tif' \
    --local-precip-nodata -9999 \
    --method cropwat \
    --start-year 2005 --end-year 2005 --months 7 8 \
    --workers 4 \
    --output ./wrf_outputs
```

The outputs keep the WRF grid exactly - 689 x 799, EPSG:4326, `nodata = NaN` - and the 65,309
nodata pixels of the input mask are `NaN` in both the effective precipitation and the fraction
raster.

### Hybrid: ensemble method with AWC and ETo from Earth Engine

The `ensemble` method needs AWC and ETo, which still come from GEE. They are downloaded for the
region and reprojected onto the WRF grid, so `awc.tif`, `eto_YYYY_MM.tif` and `precip_YYYY_MM.tif`
share one transform.

The run below passes `geometry_path='basin.geojson'` and leaves `clip_to_geometry` at its default
`True`, so the precipitation is clipped to exactly the region AWC and ETo were downloaded for and
every output pixel has all three inputs. That is the arrangement to aim for: **an AWC/ETo method
only produces values where the download covered the precipitation grid, and anywhere it did not is
left as NaN** (see [Geometry and clipping](#geometry-and-clipping)).

**Python:**

```python
ep = EffectivePrecipitation(
    local_precip='../pyCropWat_Data/Precip',
    local_precip_pattern='Precip_*.tif',
    local_precip_nodata=-9999,
    geometry_path='basin.geojson',
    start_year=2005,
    end_year=2005,
    method='ensemble',
    gee_project='my-gee-project',
    method_params={
        'awc_asset': 'projects/sat-io/open-datasets/FAO/HWSD_V2_SMU',
        'awc_band': 'AWC',
        'awc_scale_factor': 0.001,   # HWSD mm/m -> volumetric fraction
        'eto_asset': 'IDAHO_EPSCOR/TERRACLIMATE',
        'eto_band': 'pet',
        'eto_scale_factor': 0.1,     # TerraClimate PET is stored in 0.1 mm units
        'eto_is_daily': False,
        'rooting_depth': 2.0,
        'mad_factor': 1.0
    }
)

results = ep.process_sequential(
    output_dir='./wrf_ensemble',
    months=[7],
    input_dir='./wrf_inputs',
    save_inputs=True      # writes precip_2005_07.tif, awc.tif, eto_2005_07.tif
)
```

**CLI:**

```bash
pycropwat process \
    --local-precip ../pyCropWat_Data/Precip \
    --local-precip-pattern 'Precip_*.tif' \
    --local-precip-nodata -9999 \
    --geometry basin.geojson \
    --method ensemble \
    --awc-asset projects/sat-io/open-datasets/FAO/HWSD_V2_SMU \
    --awc-band AWC --awc-scale-factor 0.001 \
    --eto-asset IDAHO_EPSCOR/TERRACLIMATE --eto-band pet \
    --rooting-depth 2.0 --mad-factor 1.0 \
    --project my-gee-project \
    --start-year 2005 --end-year 2005 --months 7 \
    --output ./wrf_ensemble
```

!!! note "`save_inputs` is a Python-only option"
    `input_dir` and `save_inputs` are arguments of `process()` and `process_sequential()`; the
    `process` subcommand has no equivalent flags. Saving the inputs is the quickest way to confirm
    that AWC and ETo landed on your precipitation grid.

Runnable scripts live in **`Examples/LocalPrecip/`**; `Examples/wrf/` shows the same WRF dataset
driven through Earth Engine for comparison.

## Composing with the rest of pyCropWat

Outputs from a local run are ordinary GeoTIFFs with the standard names, so every analysis, export
and plotting tool works on them unchanged:

```python
from pycropwat.analysis import (
    TemporalAggregator, StatisticalAnalyzer, Visualizer, export_to_netcdf
)

agg = TemporalAggregator('./wrf_outputs')
annual = agg.annual_aggregate(2005, method='sum', output_path='./annual_2005.tif')
summer = agg.seasonal_aggregate(2005, 'JJA', method='sum')

# Southern Hemisphere growing season (Oct 2005 - Mar 2006)
growing = agg.growing_season_aggregate(2005, start_month=10, end_month=3)

stats = StatisticalAnalyzer('./wrf_outputs')
slope, pvalue = stats.calculate_trend(start_year=2000, end_year=2021, method='sen')

viz = Visualizer('./wrf_outputs')
viz.plot_raster(2005, 7, output_path='./map_2005_07.png')

export_to_netcdf('./wrf_outputs', './wrf_peff.nc',
                 pattern='effective_precip_[0-9]*.tif')
```

The same commands work on the CLI:

```bash
pycropwat aggregate --input ./wrf_outputs --type annual --year 2005 --output ./annual_2005.tif
pycropwat analyze trend --input ./wrf_outputs --start-year 2000 --end-year 2021 \
    --trend-method sen --output ./trend/
pycropwat plot map --input ./wrf_outputs --year 2005 --month 7 --output ./map.png
pycropwat export netcdf --input ./wrf_outputs --output ./wrf_peff.nc
```

Because outputs carry the local grid's CRS, comparing a local run against a GEE run of the same
region means reprojecting one onto the other first (`rioxarray`'s `reproject_match`) unless both
already share a grid.

## Troubleshooting

### `FileNotFoundError: Local precipitation path not found: ...`

The path does not exist and is not a glob. Check for a typo, and remember that relative paths are
resolved from the working directory, not from the output directory.

### `ValueError: No files matching pattern '*.nc' in directory ...`

The directory exists but the glob matched nothing. Adjust `local_precip_pattern` - `'*.tif'`,
`'Precip_*.tif'`, `'*.nc'`. List the directory first to confirm the extension and its case.

### `ValueError: Glob '...' matched no files.`

Quote the glob so the shell does not expand it before pyCropWat sees it:
`--local-precip './Precip/Precip_*.tif'`.

### `ValueError: No dated months could be resolved from ...`

No file name matched a known layout. Test one name first:

```python
from pycropwat import parse_year_month

parse_year_month('your_file.tif')   # None means the parser could not date it
```

If it returns `None`, pass `local_precip_date_regex` with named `year`/`month` groups, or rename
the files to `something_YYYY_MM.tif`.

### `ValueError: ... looks like a single multi-band raster for year 2005 (12 bands).`

One 12-band raster per year is not supported. Split it into monthly files first:

```python
import rioxarray

da = rioxarray.open_rasterio('Precip_2005.tif')
for i in range(da.sizes['band']):
    da.isel(band=i, drop=True).rio.to_raster(f'Precip_2005_{i + 1:02d}.tif')
```

### The wrong NetCDF variable was chosen, or `Could not determine which variable holds precipitation`

Pass `local_precip_variable='<name>'` (`--local-precip-variable`). The error lists the candidates,
and `xarray.open_dataset(path)` shows every variable in the file.

### `ValueError: Could not identify the spatial dimensions ...`

The dimension names are outside the recognised lists. Use `LocalPrecipitationSource` with
`x_dim=`/`y_dim=`, or rename the dimensions and write a normalised NetCDF - see
[Dimensions](#dimensions).

### `ValueError: Mixed raster and NetCDF inputs are not supported.`

The selection picked up both kinds. Narrow `local_precip_pattern` to one of them (`'*.tif'` or
`'*.nc'`).

### The output is entirely `NaN`

The whole grid was masked. Check that `local_precip_nodata` is not masking real data (passing `0`
when zero precipitation is valid, for instance) and that the file's declared nodata is correct.
Reading one month straight from the source tells you immediately whether the input or the
calculation is at fault:

```python
from pycropwat import open_local_precipitation
import numpy as np

with open_local_precipitation('./Precip', pattern='Precip_*.tif', nodata=-9999) as src:
    da = src.get_month(2005, 7)
    print(float(np.nanmean(da.values)), int(np.isnan(da.values).sum()), da.size)
```

### Values look about 1000x too small or too large

`precip_scale_factor` is wrong. Compare `float(src.get_month(y, m).max())` against a plausible
monthly total in millimetres and adjust - `1000` for metres, `10` for centimetres, `25.4` for
inches.

### Effective precipitation is zero over part of the domain

Either the clip landed somewhere unexpected or the precipitation really is zero. Compare
`src.bounds` with the vector's bounds: if they do not overlap, clipping either fails (a warning is
logged and the unclipped array is used) or yields an empty region.

### Part of a `usda_scs`, `ensemble` or `suet` output is NaN

The AWC/ETo download did not cover the whole precipitation grid. pyCropWat leaves those pixels NaN
instead of substituting a filler value, and says so at WARNING level - once per field per process,
so `usda_scs` and `ensemble` log it twice (AWC and ETo), `suet` once (ETo):

```
AWC covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
ETo covers only 3.0% of the precipitation grid; 534255 pixel(s) will be NaN in the output. Clip the precipitation to the geometry (clip_to_geometry=True) or widen the geometry.
```

Grep for `covers only` if the run's log has scrolled past.

This happens when the precipitation grid is larger than the geometry the AWC/ETo were requested
for - a GEE `FeatureCollection` geometry (which never clips the local rasters), or a local vector
with `clip_to_geometry=False`. Either clip to the download region (`clip_to_geometry=True`, the
default) or widen the geometry until it spans the whole precipitation extent. Precipitation-only
methods never hit this, because they download nothing.

### AWC or ETo look misaligned with precipitation

Run with `save_inputs=True` and compare `awc.tif`, `eto_YYYY_MM.tif` and `precip_YYYY_MM.tif` -
the CRS and transform should be identical, because AWC and ETo are reprojected onto the
precipitation grid. If the log shows a warning about the reprojection, the older resampling
fallback was used; supplying a geometry makes the downloaded region smaller and better
conditioned.

### `--scale` seems to do nothing

With local precipitation it does not change the precipitation grid, which is always used exactly as
stored. It only sets the resolution at which AWC and ETo are requested from Earth Engine.

### Memory pressure with large NetCDF stacks

Only one month is held in memory at a time, but a many-file `open_mfdataset` still builds a large
index. Lower `n_workers`, process fewer years per run with `start_year`/`end_year`, or split the
stack into per-year files. For very long records a directory of monthly rasters is the lighter
option: each month opens and closes its own file.

### Threading errors when reading NetCDF in parallel

This is handled for you: NetCDF datasets stay open and shared, so reads are serialised with an
internal lock, while raster reads open their own handle per month and run concurrently. If an
exotic NetCDF backend still raises file-handle errors, fall back to `process_sequential()` or
`process(..., n_workers=1)`.

## Option reference

| Python argument | CLI flag | Default | Purpose |
|-----------------|----------|---------|---------|
| `local_precip` | `--local-precip` | `None` | Directory, NetCDF file, or glob of local precipitation |
| `local_precip_pattern` | `--local-precip-pattern` | `'*.tif'` | Glob used when the path is a directory |
| `local_precip_variable` | `--local-precip-variable` | `None` | NetCDF variable name (auto-detected when omitted) |
| `local_precip_nodata` | `--local-precip-nodata` | `None` | Extra nodata sentinel, e.g. `-9999` |
| `local_precip_crs` | `--local-precip-crs` | `None` | CRS override: relabels the grid, replacing any CRS the file declares. Never reprojects |
| `local_precip_date_regex` | `--local-precip-date-regex` | `None` | Regex with named `year`/`month` groups |
| `clip_to_geometry` | *(none - always `True`)* | `True` | Clip local rasters to a local vector geometry |
| `precip_scale_factor` | `--scale-factor` | `1.0` | Multiplier to convert source units to mm |
| `start_year` / `end_year` | `--start-year` / `--end-year` | `None` | Optional; inferred from the files |
| `asset_id` / `precip_band` | `--asset` / `--band` | `None` | Not required when local precipitation is used |

## Next Steps

- [`pycropwat.local` API reference](../api/local.md) - full `LocalPrecipitationSource` documentation
- [Quick Start](quickstart.md) - the same workflow with Earth Engine precipitation
- [CLI Reference](cli.md) - every `process` flag
- [Python API](api.md) - aggregation, statistics and plotting in depth

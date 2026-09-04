"""
pyCropWat Local Precipitation from NetCDF
=========================================

pyCropWat can read local precipitation from NetCDF as easily as from a directory
of GeoTIFFs. This example shows the whole round trip and proves the two paths
give identical answers:

    Step 1 - Build a NetCDF stack from the monthly WRF GeoTIFFs (one year,
             12 months, variable 'precip', dim 'time', mm). Written compressed
             when netCDF4 or h5netcdf is installed, uncompressed NetCDF3 otherwise.
    Step 2 - Compute effective precipitation straight from that NetCDF.
    Step 3 - Compute the same months from the original GeoTIFFs and assert the
             two sets of rasters agree to within floating point tolerance.

Everything here uses a precipitation-only method, so no Earth Engine is involved.

Data:
    WRF regional climate model monthly precipitation for South America
    (Rio de la Plata basin domain), ``Precip_YYYY_MM.tif``, EPSG:4326,
    689 x 799 @ ~0.0397 deg, float32, nodata -9999, mm/month.
    Expected at ``<repo>/../pyCropWat_Data/Precip``; override with --precip-dir.

Outputs (under --output, default ./LocalPrecip_NetCDF):
    precip_YYYY.nc            The NetCDF built in Step 1
    from_netcdf/*.tif         Peff + fraction computed from the NetCDF
    from_geotiff/*.tif        Peff + fraction computed from the GeoTIFFs

Usage:
    python local_netcdf_example.py
    python local_netcdf_example.py --year 2010 --method farmwest
    python local_netcdf_example.py --months 7 8          # quick run
    python local_netcdf_example.py --skip-build          # reuse an existing .nc
    python local_netcdf_example.py --precip-dir /data/precip --output ./nc_test
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from pycropwat import EffectivePrecipitation, open_local_precipitation


# =============================================================================
# Configuration
# =============================================================================

# Examples/LocalPrecip/local_netcdf_example.py -> repository root
REPO_ROOT = Path(__file__).resolve().parents[2]

# The example dataset lives next to the repository, not inside it
DEFAULT_PRECIP_DIR = REPO_ROOT.parent / 'pyCropWat_Data' / 'Precip'

DEFAULT_PATTERN = 'Precip_*.tif'
DEFAULT_NODATA = -9999.0
DEFAULT_OUTPUT = Path('./LocalPrecip_NetCDF')
DEFAULT_YEAR = 2005

# NetCDF layout written by Step 1
NC_VARIABLE = 'precip'
NC_TIME_DIM = 'time'

# Methods that need precipitation only - these require no Earth Engine at all
PRECIP_ONLY_METHODS = [
    'cropwat',
    'fao_aglw',
    'fixed_percentage',
    'dependable_rainfall',
    'farmwest',
]

# Agreement tolerance between the NetCDF and GeoTIFF results (mm)
TOLERANCE = 1e-4


# =============================================================================
# Helpers
# =============================================================================

def check_precip_dir(precip_dir: Path) -> None:
    """
    Exit with a friendly message when the precipitation data is missing.

    Parameters
    ----------

    precip_dir : Path
        Directory that should contain the monthly precipitation rasters.
    """
    if precip_dir.exists():
        return

    logger.error("Precipitation directory not found: %s", precip_dir)
    logger.error(
        "This example expects the WRF South America monthly GeoTIFFs at "
        "'<repo>/../pyCropWat_Data/Precip'."
    )
    logger.error(
        "Pass --precip-dir /path/to/your/monthly/rasters to point at your own data "
        "(and --nodata if your nodata value differs)."
    )
    sys.exit(1)


def validate_months(months) -> None:
    """
    Exit with a friendly message when ``--months`` is outside 1-12.

    Mirrors :func:`check_precip_dir`: a request that cannot be satisfied is
    reported in plain language and exits 1, never as a traceback - and before
    the NetCDF is built, so no work is wasted on a request that can produce
    nothing.

    Parameters
    ----------

    months : list of int or None
        Months requested on the command line; None means all twelve.
    """
    if not months:
        return

    bad = sorted({month for month in months if not 1 <= month <= 12})
    if not bad:
        return

    logger.error(
        "Invalid --months value(s): %s.",
        ' '.join(str(month) for month in bad)
    )
    logger.error(
        "Months must be between 1 and 12. The NetCDF holds one calendar year - "
        "pass --year to choose a different year."
    )
    sys.exit(1)


def pick_netcdf_engine():
    """
    Choose the NetCDF writer, preferring one that supports compression.

    ``netCDF4`` or ``h5netcdf`` write compressed NetCDF4/HDF5. The ``scipy``
    backend (always available with xarray) only writes classic NetCDF3, which has
    no compression - usable, but the file is several times larger.

    Returns
    -------
    tuple
        ``(engine_name, supports_compression)``.
    """
    engines = set(xr.backends.list_engines())
    for name in ('netcdf4', 'h5netcdf'):
        if name in engines:
            return name, True
    if 'scipy' in engines:
        logger.warning(
            "Neither netCDF4 nor h5netcdf is installed; falling back to the scipy "
            "backend (classic NetCDF3, no compression). Install netCDF4 "
            "(pip install netCDF4) for compressed NetCDF4 output."
        )
        return 'scipy', False
    logger.error("No NetCDF backend available. Install netCDF4: pip install netCDF4")
    sys.exit(1)


# =============================================================================
# Step 1: Build a NetCDF from the monthly GeoTIFFs
# =============================================================================

def build_netcdf(precip_dir: Path, year: int, nodata: float, nc_path: Path) -> Path:
    """
    Stack twelve monthly GeoTIFFs into one compressed, CF-ish NetCDF file.

    Parameters
    ----------

    precip_dir : Path
        Directory holding ``Precip_YYYY_MM.tif``.

    year : int
        Year to stack.

    nodata : float
        Nodata sentinel written to the NetCDF as ``_FillValue``.

    nc_path : Path
        Destination NetCDF file.

    Returns
    -------
    Path
        The written NetCDF file.
    """
    import rioxarray  # noqa: F401  (registers the .rio accessor)

    logger.info(f"Building {nc_path.name} from monthly GeoTIFFs for {year}...")

    slices = []
    times = []
    for month in range(1, 13):
        tif = precip_dir / f'Precip_{year}_{month:02d}.tif'
        if not tif.exists():
            logger.warning(f"  Missing {tif.name}, skipping")
            continue
        da = rioxarray.open_rasterio(tif, masked=True).squeeze('band', drop=True)
        slices.append(da.astype('float32'))
        times.append(pd.Timestamp(year=year, month=month, day=1))

    if not slices:
        logger.error(f"No monthly GeoTIFFs found for {year} in {precip_dir}")
        logger.error("Pass --year with a year that exists, or --precip-dir with "
                     "the directory holding your rasters.")
        sys.exit(1)

    crs = slices[0].rio.crs

    stacked = xr.concat(slices, dim=pd.Index(times, name=NC_TIME_DIM))
    stacked.name = NC_VARIABLE
    stacked.attrs = {
        'units': 'mm',
        'long_name': 'monthly total precipitation',
        'standard_name': 'precipitation_amount',
        'cell_methods': 'time: sum',
    }

    dataset = stacked.to_dataset()
    dataset = dataset.rio.write_crs(crs)

    if crs is not None and crs.is_geographic:
        dataset['x'].attrs = {'standard_name': 'longitude', 'units': 'degrees_east',
                              'axis': 'X'}
        dataset['y'].attrs = {'standard_name': 'latitude', 'units': 'degrees_north',
                              'axis': 'Y'}
    else:
        dataset['x'].attrs = {'standard_name': 'projection_x_coordinate', 'axis': 'X'}
        dataset['y'].attrs = {'standard_name': 'projection_y_coordinate', 'axis': 'Y'}

    dataset.attrs = {
        'Conventions': 'CF-1.8',
        'title': f'WRF South America monthly precipitation, {year}',
        'summary': 'Monthly precipitation totals stacked from pyCropWat example GeoTIFFs',
        'source': str(precip_dir),
        'history': 'Created by pyCropWat Examples/LocalPrecip/local_netcdf_example.py',
    }

    engine, compress = pick_netcdf_engine()

    # Start from the encoding rioxarray already set (it carries
    # grid_mapping='spatial_ref'). xarray REPLACES a variable's encoding with the
    # dict passed to to_netcdf(), so building it from scratch would silently drop
    # the CRS link and pyCropWat would fall back to "assuming EPSG:4326".
    var_encoding = dict(dataset[NC_VARIABLE].encoding)
    var_encoding.update({
        'dtype': 'float32',
        '_FillValue': np.float32(nodata),
    })
    encoding = {
        NC_VARIABLE: var_encoding,
        NC_TIME_DIM: {
            'units': 'days since 1900-01-01',
            'calendar': 'standard',
        },
    }
    if compress:
        var_encoding.update({'zlib': True, 'complevel': 4})
    else:
        # Classic NetCDF3 has no int64; the CRS is carried as a small int32 variable
        if 'spatial_ref' in dataset.coords:
            dataset['spatial_ref'] = dataset['spatial_ref'].astype('int32')

    nc_path.parent.mkdir(parents=True, exist_ok=True)
    if nc_path.exists():
        nc_path.unlink()
    dataset.to_netcdf(nc_path, engine=engine, encoding=encoding)
    dataset.close()
    for da in slices:
        da.close()

    size_mb = nc_path.stat().st_size / 1e6
    logger.info(f"  Wrote {nc_path} ({len(times)} months, {size_mb:.1f} MB, "
                f"engine={engine}, compressed={compress})")
    return nc_path


def summarize_netcdf(nc_path: Path, nodata: float) -> None:
    """
    Open the NetCDF through pyCropWat and print what it found.

    Parameters
    ----------

    nc_path : Path
        NetCDF file to inspect.

    nodata : float
        Extra nodata sentinel.
    """
    with open_local_precipitation(nc_path, variable=NC_VARIABLE, nodata=nodata) as src:
        min_lon, min_lat, max_lon, max_lat = src.bounds
        logger.info("-" * 70)
        logger.info(f"  NetCDF source : {src.kind}, {len(src)} months, "
                    f"{src.year_range[0]}-{src.year_range[1]}")
        logger.info(f"  Variable      : {NC_VARIABLE}")
        logger.info(f"  Grid          : {src.shape} @ {src.resolution[0]:.6f}, "
                    f"CRS {src.crs}")
        logger.info(f"  Bounds (WGS84): {min_lon:.4f}, {min_lat:.4f}, "
                    f"{max_lon:.4f}, {max_lat:.4f}")
        sample = src.get_month(*src.available_months()[0])
        logger.info(f"  First month   : {sample.attrs['year']}-"
                    f"{sample.attrs['month']:02d}, "
                    f"mean {float(np.nanmean(sample.values)):.2f} "
                    f"{sample.attrs['units']}")
        logger.info("-" * 70)


# =============================================================================
# Step 2 / 3: Compute Peff from each source and compare
# =============================================================================

def compute_peff(local_precip, output_dir: Path, year: int, method: str,
                 nodata: float, months, variable=None, pattern=DEFAULT_PATTERN) -> list:
    """
    Run pyCropWat over a local precipitation source.

    Parameters
    ----------

    local_precip : str or Path
        Directory of rasters, or a NetCDF file.

    output_dir : Path
        Where to write the Peff rasters.

    year : int
        Year to process.

    method : str
        Precipitation-only Peff method.

    nodata : float
        Extra nodata sentinel.

    months : list of int or None
        Months to process; None means all twelve.

    variable : str, optional
        NetCDF variable name (ignored for rasters).

    pattern : str, optional
        Glob for raster directories (ignored for NetCDF).

    Returns
    -------
    list
        List of ``(peff_path, fraction_path)`` tuples.
    """
    try:
        ep = EffectivePrecipitation(
            local_precip=str(local_precip),
            local_precip_pattern=pattern,
            local_precip_variable=variable,
            local_precip_nodata=nodata,
            start_year=year,
            end_year=year,
            method=method,
        )
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Could not read precipitation for {year} from "
                     f"{local_precip}: {e}")
        logger.error("Pass --year with a year that exists in --precip-dir "
                     "(and in the NetCDF, when using --skip-build).")
        sys.exit(1)
    return ep.process_sequential(output_dir=str(output_dir), months=months)


def compare_results(nc_dir: Path, tif_dir: Path, year: int, months) -> bool:
    """
    Assert the NetCDF-derived and GeoTIFF-derived rasters agree.

    Parameters
    ----------

    nc_dir : Path
        Directory of rasters computed from the NetCDF.

    tif_dir : Path
        Directory of rasters computed from the GeoTIFFs.

    year : int
        Year that was processed.

    months : list of int or None
        Months that were processed; None means all twelve.

    Returns
    -------
    bool
        True when every compared pair agrees within :data:`TOLERANCE`.
    """
    import rioxarray

    month_list = months if months else list(range(1, 13))
    logger.info("-" * 70)
    logger.info(f"{'Month':>7} | {'Peff max |diff|':>15} | {'Fraction max |diff|':>19} | "
                f"NaN match")
    logger.info("-" * 70)

    all_ok = True
    compared = 0

    for month in month_list:
        row_ok = True
        diffs = []
        for prefix in ('effective_precip', 'effective_precip_fraction'):
            name = f'{prefix}_{year}_{month:02d}.tif'
            nc_file = nc_dir / name
            tif_file = tif_dir / name
            if not (nc_file.exists() and tif_file.exists()):
                diffs.append(None)
                continue
            nc_values = rioxarray.open_rasterio(
                nc_file, masked=True).squeeze('band', drop=True).values
            tif_values = rioxarray.open_rasterio(
                tif_file, masked=True).squeeze('band', drop=True).values

            nan_match = bool(np.array_equal(np.isnan(nc_values), np.isnan(tif_values)))
            valid = ~np.isnan(nc_values) & ~np.isnan(tif_values)
            max_diff = float(np.max(np.abs(nc_values[valid] - tif_values[valid]))) \
                if valid.any() else 0.0
            diffs.append((max_diff, nan_match))
            if max_diff > TOLERANCE or not nan_match:
                row_ok = False

        if any(item is None for item in diffs):
            logger.warning(f"{year}-{month:02d} : output missing, not compared")
            continue

        compared += 1
        all_ok = all_ok and row_ok
        logger.info(f"{year}-{month:02d} | {diffs[0][0]:15.3e} | {diffs[1][0]:19.3e} | "
                    f"{'yes' if diffs[0][1] and diffs[1][1] else 'NO'}")

    logger.info("-" * 70)
    if compared == 0:
        logger.error("Nothing was compared - no matching output rasters were found.")
        return False
    if all_ok:
        logger.info(f"Assertion passed: NetCDF and GeoTIFF results agree for "
                    f"{compared} month(s) within {TOLERANCE:g} mm.")
    else:
        logger.error(f"Assertion FAILED: differences above {TOLERANCE:g} mm "
                     f"(or mismatched NaN masks).")
    return all_ok


# =============================================================================
# Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            'Build a NetCDF from monthly precipitation GeoTIFFs, compute effective '
            'precipitation from it, and verify it matches the GeoTIFF result.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python local_netcdf_example.py
  python local_netcdf_example.py --year 2010 --method fao_aglw
  python local_netcdf_example.py --months 7 8
  python local_netcdf_example.py --skip-build --output ./LocalPrecip_NetCDF
        """
    )
    parser.add_argument('--precip-dir', type=Path, default=DEFAULT_PRECIP_DIR,
                        help=f'Directory of monthly precipitation GeoTIFFs '
                             f'(default: {DEFAULT_PRECIP_DIR})')
    parser.add_argument('--year', type=int, default=DEFAULT_YEAR,
                        help=f'Year to stack into the NetCDF (default: {DEFAULT_YEAR})')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT,
                        help=f'Output directory (default: {DEFAULT_OUTPUT})')
    parser.add_argument('--method', type=str, default='cropwat',
                        choices=PRECIP_ONLY_METHODS,
                        help='Precipitation-only Peff method (default: cropwat)')
    parser.add_argument('--nodata', type=float, default=DEFAULT_NODATA,
                        help=f'Nodata sentinel for both the GeoTIFFs and the NetCDF '
                             f'(default: {DEFAULT_NODATA})')
    parser.add_argument('--skip-build', action='store_true',
                        help='Reuse an existing NetCDF instead of rebuilding it')
    parser.add_argument('--months', type=int, nargs='+', default=None,
                        help='Months to compute and compare, 1-12 (default: all 12). '
                             'The NetCDF always holds the full year')
    return parser.parse_args()


def main():
    args = parse_args()

    precip_dir = Path(args.precip_dir).expanduser()
    output_dir = Path(args.output)
    nc_path = output_dir / f'precip_{args.year}.nc'
    nc_out = output_dir / 'from_netcdf'
    tif_out = output_dir / 'from_geotiff'

    logger.info("=" * 70)
    logger.info("pyCropWat - Local Precipitation from NetCDF")
    logger.info("=" * 70)
    logger.info(f"  Precipitation : {precip_dir}")
    logger.info(f"  Year          : {args.year}")
    logger.info(f"  Method        : {args.method}")
    logger.info(f"  Months        : {args.months or 'all 12'}")
    logger.info(f"  Output        : {output_dir.resolve()}")
    logger.info("=" * 70)

    check_precip_dir(precip_dir)
    validate_months(args.months)

    # ---------------- Step 1: GeoTIFFs -> NetCDF ----------------
    if args.skip_build:
        if not nc_path.exists():
            logger.error(f"--skip-build was given but {nc_path} does not exist.")
            logger.error("Run once without --skip-build to create it.")
            sys.exit(1)
        logger.info(f"Step 1: Reusing existing NetCDF {nc_path}")
    else:
        logger.info("Step 1: Building a NetCDF from the monthly GeoTIFFs...")
        build_netcdf(precip_dir, args.year, args.nodata, nc_path)

    summarize_netcdf(nc_path, args.nodata)

    # ---------------- Step 2: Peff from the NetCDF ----------------
    logger.info("Step 2: Computing effective precipitation from the NetCDF...")
    nc_results = compute_peff(nc_path, nc_out, args.year, args.method,
                              args.nodata, args.months, variable=NC_VARIABLE)
    logger.info(f"  Wrote {sum(1 for peff, _ in nc_results if peff)} raster pair(s) "
                f"to {nc_out}")

    # ---------------- Step 3: Peff from the GeoTIFFs, then compare ----------------
    logger.info("Step 3: Computing the same months from the GeoTIFFs and comparing...")
    tif_results = compute_peff(precip_dir, tif_out, args.year, args.method,
                               args.nodata, args.months)
    logger.info(f"  Wrote {sum(1 for peff, _ in tif_results if peff)} raster pair(s) "
                f"to {tif_out}")

    ok = compare_results(nc_out, tif_out, args.year, args.months)

    logger.info("=" * 70)
    logger.info("Done.")
    logger.info(f"  NetCDF          : {nc_path}")
    logger.info(f"  Peff (NetCDF)   : {nc_out}")
    logger.info(f"  Peff (GeoTIFF)  : {tif_out}")
    logger.info("=" * 70)

    if not ok:
        sys.exit(1)


if __name__ == '__main__':
    main()

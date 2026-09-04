"""
pyCropWat Local Precipitation Quickstart (no Earth Engine)
==========================================================

The smallest useful pyCropWat workflow: point the package at a directory of
monthly precipitation GeoTIFFs on disk and calculate effective precipitation
with a precipitation-only method. Nothing is downloaded, ``ee.Initialize()`` is
never called, and the script runs offline.

Precipitation-only methods (no AWC, no ETo, no Earth Engine):
    cropwat, fao_aglw, fixed_percentage, dependable_rainfall, farmwest

The methods ``usda_scs``, ``ensemble`` and ``suet`` also work with local
precipitation, but they need AWC and/or ETo, which pyCropWat still reads from
Google Earth Engine - see ``wrf_south_america_local_example.py`` for that.

Data:
    WRF regional climate model monthly precipitation for South America
    (Rio de la Plata basin domain), 264 GeoTIFFs named ``Precip_YYYY_MM.tif``
    covering 2000-01 to 2021-12, EPSG:4326, 689 x 799 @ ~0.0397 deg (~4.4 km),
    float32, nodata -9999, units mm/month.

    Expected at ``<repo>/../pyCropWat_Data/Precip``; override with --precip-dir.

Outputs (per processed month, in --output):
    effective_precip_YYYY_MM.tif           Effective precipitation (mm)
    effective_precip_fraction_YYYY_MM.tif  Peff / P (0-1)

Usage:
    python local_precip_quickstart.py
    python local_precip_quickstart.py --method farmwest
    python local_precip_quickstart.py --start-year 2005 --end-year 2010 --months 6 7 8
    python local_precip_quickstart.py --precip-dir /data/my_precip --pattern '*.tif'
    python local_precip_quickstart.py -w 8 --output ./my_output
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

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

# Examples/LocalPrecip/local_precip_quickstart.py -> repository root
REPO_ROOT = Path(__file__).resolve().parents[2]

# The example dataset lives next to the repository, not inside it
DEFAULT_PRECIP_DIR = REPO_ROOT.parent / 'pyCropWat_Data' / 'Precip'

DEFAULT_PATTERN = 'Precip_*.tif'
DEFAULT_NODATA = -9999.0
DEFAULT_OUTPUT = Path('./LocalPrecip_Quickstart')

# Methods that need precipitation only - these require no Earth Engine at all
PRECIP_ONLY_METHODS = [
    'cropwat',
    'fao_aglw',
    'fixed_percentage',
    'dependable_rainfall',
    'farmwest',
]


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
        "(and --pattern / --nodata if the naming or nodata value differs)."
    )
    sys.exit(1)


def summarize_source(precip_dir: Path, pattern: str, nodata: float) -> list:
    """
    Open the local precipitation source and print a short summary.

    Parameters
    ----------

    precip_dir : Path
        Directory of monthly rasters.

    pattern : str
        Glob applied inside ``precip_dir``.

    nodata : float
        Extra nodata sentinel, applied on top of the file metadata.

    Returns
    -------
    list of tuple
        Sorted ``(year, month)`` pairs found in the source.
    """
    try:
        with open_local_precipitation(precip_dir, pattern=pattern, nodata=nodata) as src:
            min_lon, min_lat, max_lon, max_lat = src.bounds
            first_year, last_year = src.year_range
            months = src.available_months()

            logger.info("-" * 70)
            logger.info("Local precipitation source")
            logger.info("-" * 70)
            logger.info(f"  Path          : {src.path}")
            logger.info(f"  Kind          : {src.kind}")
            logger.info(f"  Files         : {len(src.files)}")
            logger.info(f"  Months        : {len(src)} "
                        f"({months[0][0]}-{months[0][1]:02d} to "
                        f"{months[-1][0]}-{months[-1][1]:02d})")
            logger.info(f"  Year range    : {first_year} - {last_year}")
            logger.info(f"  CRS           : {src.crs}")
            logger.info(f"  Grid shape    : {src.shape} (rows, cols)")
            logger.info(f"  Resolution    : {src.resolution[0]:.6f}, "
                        f"{src.resolution[1]:.6f} (CRS units)")
            logger.info(f"  Bounds (WGS84): {min_lon:.4f}, {min_lat:.4f}, "
                        f"{max_lon:.4f}, {max_lat:.4f}")
            logger.info("-" * 70)
            return months
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Could not open the local precipitation data: {e}")
        logger.error("Check --precip-dir, --pattern and the file naming "
                     "(pyCropWat expects e.g. 'Precip_2005_07.tif').")
        sys.exit(1)


def select_months(available: list, start_year: int, end_year: int,
                  months: list) -> list:
    """
    Narrow the months on disk to the requested period, or exit with a message.

    Mirrors :func:`check_precip_dir`: a request that cannot be satisfied is
    reported in plain language and exits 1, never as a traceback.

    Parameters
    ----------

    available : list of tuple
        Sorted ``(year, month)`` pairs present on disk.

    start_year : int
        First year requested.

    end_year : int
        Last year requested.

    months : list of int
        Months requested (1-12).

    Returns
    -------
    list of tuple
        The ``(year, month)`` pairs that will actually be processed.
    """
    selected = [
        (year, month) for year, month in available
        if start_year <= year <= end_year and month in months
    ]
    if selected:
        return selected

    logger.error(
        "Nothing to process: no file matches --start-year %d --end-year %d "
        "--months %s.",
        start_year, end_year, ' '.join(str(month) for month in months)
    )
    if available:
        logger.error(
            "The data on disk covers %d-%02d to %d-%02d (years %d-%d).",
            available[0][0], available[0][1], available[-1][0], available[-1][1],
            available[0][0], available[-1][0]
        )
        logger.error(
            "Choose --start-year/--end-year inside that period and --months in 1-12, "
            "or pass --precip-dir to point at different data."
        )
    else:
        logger.error(
            "No datable rasters were found at all - check --precip-dir and --pattern."
        )
    sys.exit(1)


def report_output(peff_path: Path) -> None:
    """
    Print min/mean/max of one effective precipitation raster.

    Parameters
    ----------

    peff_path : Path
        Path to an ``effective_precip_YYYY_MM.tif`` file.
    """
    import rioxarray

    da = rioxarray.open_rasterio(peff_path, masked=True).squeeze('band', drop=True)
    values = da.values
    valid = int(np.isfinite(values).sum())

    logger.info("-" * 70)
    logger.info(f"Sample output: {peff_path.name}")
    logger.info(f"  Shape       : {values.shape}")
    logger.info(f"  CRS         : {da.rio.crs}")
    logger.info(f"  Valid pixels: {valid} of {values.size} "
                f"({values.size - valid} nodata/NaN)")
    logger.info(f"  Peff min    : {np.nanmin(values):.3f} mm")
    logger.info(f"  Peff mean   : {np.nanmean(values):.3f} mm")
    logger.info(f"  Peff max    : {np.nanmax(values):.3f} mm")
    logger.info("-" * 70)


# =============================================================================
# Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            'Quickstart: effective precipitation from local monthly GeoTIFFs. '
            'Uses a precipitation-only method, so no Earth Engine is required.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python local_precip_quickstart.py
  python local_precip_quickstart.py --method farmwest --months 1 2 3
  python local_precip_quickstart.py --start-year 2005 --end-year 2010 -w 8
  python local_precip_quickstart.py --precip-dir /data/precip --pattern '*.tif' \\
                                    --nodata -9999
        """
    )
    parser.add_argument('--precip-dir', type=Path, default=DEFAULT_PRECIP_DIR,
                        help=f'Directory of monthly precipitation rasters '
                             f'(default: {DEFAULT_PRECIP_DIR})')
    parser.add_argument('--pattern', type=str, default=DEFAULT_PATTERN,
                        help=f"Glob for the monthly rasters (default: '{DEFAULT_PATTERN}')")
    parser.add_argument('--nodata', type=float, default=DEFAULT_NODATA,
                        help=f'Extra nodata sentinel, applied on top of the file '
                             f'metadata (default: {DEFAULT_NODATA})')
    parser.add_argument('--start-year', type=int, default=2005,
                        help='First year to process (default: 2005)')
    parser.add_argument('--end-year', type=int, default=2005,
                        help='Last year to process (default: 2005)')
    parser.add_argument('--months', type=int, nargs='+', default=[7, 8],
                        help='Months to process, 1-12 (default: 7 8)')
    parser.add_argument('--method', type=str, default='cropwat',
                        choices=PRECIP_ONLY_METHODS,
                        help='Precipitation-only Peff method (default: cropwat)')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT,
                        help=f'Output directory (default: {DEFAULT_OUTPUT})')
    parser.add_argument('-w', '--workers', type=int, default=4,
                        help='Parallel workers; 1 runs sequentially (default: 4)')
    return parser.parse_args()


def main():
    args = parse_args()

    precip_dir = Path(args.precip_dir).expanduser()
    output_dir = Path(args.output)

    logger.info("=" * 70)
    logger.info("pyCropWat - Local Precipitation Quickstart (no Earth Engine)")
    logger.info("=" * 70)
    logger.info(f"  Precipitation : {precip_dir}")
    logger.info(f"  Pattern       : {args.pattern}")
    logger.info(f"  Method        : {args.method}")
    logger.info(f"  Period        : {args.start_year}-{args.end_year}, "
                f"months {args.months}")
    logger.info(f"  Output        : {output_dir.resolve()}")
    logger.info("=" * 70)

    check_precip_dir(precip_dir)
    available = summarize_source(precip_dir, args.pattern, args.nodata)
    selected = select_months(available, args.start_year, args.end_year, args.months)
    logger.info(f"Months to process: {len(selected)}")

    # ------------------------------------------------------------------
    # The whole workflow: local files in, Peff rasters out
    # ------------------------------------------------------------------
    ep = EffectivePrecipitation(
        local_precip=str(precip_dir),
        local_precip_pattern=args.pattern,
        local_precip_nodata=args.nodata,
        start_year=args.start_year,
        end_year=args.end_year,
        method=args.method,
    )

    if args.workers > 1:
        results = ep.process(
            output_dir=str(output_dir),
            n_workers=args.workers,
            months=args.months,
        )
    else:
        results = ep.process_sequential(
            output_dir=str(output_dir),
            months=args.months,
        )

    written = [peff for peff, _ in results if peff is not None]
    logger.info(f"Wrote {len(written)} effective precipitation raster(s) to {output_dir}")

    if not written:
        logger.error("No output rasters were produced. Check the requested "
                     "years/months against the months available on disk.")
        sys.exit(1)

    report_output(Path(written[0]))

    logger.info("Done.")
    logger.info(f"  Peff rasters     : {output_dir}/effective_precip_YYYY_MM.tif")
    logger.info(f"  Fraction rasters : {output_dir}/effective_precip_fraction_YYYY_MM.tif")


if __name__ == '__main__':
    main()

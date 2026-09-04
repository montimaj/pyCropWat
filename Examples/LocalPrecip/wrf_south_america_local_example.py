"""
pyCropWat WRF South America Effective Precipitation - LOCAL Precipitation
=========================================================================

This is the local-file twin of ``Examples/wrf/wrf_south_america_example.py``.
The science is identical - same WRF South America precipitation, same global AWC
and ETo assets, same rooting depth and MAD factor - but precipitation is read
from monthly GeoTIFFs **on disk** instead of the Earth Engine asset
``projects/azhydro/assets/WRF-ET-SA``.

What comes from where:
    Precipitation : local monthly rasters (``Precip_YYYY_MM.tif``)   - no GEE
    AWC           : FAO HWSD v2 (GEE)                                - needs GEE
    ETo           : TerraClimate (GEE)                               - needs GEE

Because this example uses ``method='usda_scs'``, Earth Engine credentials ARE
required (AWC and ETo are downloaded and regridded onto the local precipitation
grid). If you only want precipitation-only methods, see
``local_precip_quickstart.py``, which runs completely offline.

Workflow:
    1. Run USDA-SCS with ``save_inputs=True``: writes Peff + fraction rasters and
       saves the precipitation, AWC and ETo grids that were used.
    2. Calculate the other 7 Peff methods locally from those saved rasters
       (no additional GEE calls needed).

Methods computed (8 total, excludes PCML - PCML is a pre-computed Western U.S.
GEE product and cannot be combined with local precipitation):
    1. USDA-SCS     - Local precip + GEE AWC + GEE ETo
    2. Ensemble     - Mean of 6 methods (calculated locally)
    3. CROPWAT      - FAO CROPWAT method
    4. FAO/AGLW     - FAO Dependable Rainfall (80% exceedance)
    5. Fixed %      - Simple 70% fixed percentage
    6. Dependable   - FAO Dependable Rainfall (75% probability)
    7. FarmWest     - WSU irrigation scheduling formula
    8. TAGEM-SuET   - Turkish irrigation method (P - ETo, if P > 75 mm)

Local Precipitation Data:
    WRF regional climate model, Rio de la Plata basin domain.
    264 GeoTIFFs ``Precip_YYYY_MM.tif`` (2000-01 to 2021-12), EPSG:4326,
    689 x 799 @ ~0.0397 deg (~4.4 km), float32, nodata -9999, mm/month.
    Expected at ``<repo>/../pyCropWat_Data/Precip``; override with --precip-dir.

AWC and ETo Assets (Global, from GEE):
    AWC: FAO HWSD v2 (projects/sat-io/open-datasets/FAO/HWSD_V2_SMU, band: AWC)
         Units: mm/m -> awc_scale_factor=0.001 to convert to volumetric fraction
    ETo: TerraClimate (IDAHO_EPSCOR/TERRACLIMATE, band: pet, scale factor 0.1)

Study Area:
    By default the FULL extent of the local rasters is used - no geometry needed.
    Pass --geometry to clip: either a local vector file (.shp/.geojson), which
    clips the local rasters, or a GEE FeatureCollection asset ID, which is used
    for the AWC/ETo requests (local rasters are not clipped in that case).

Output Structure (under --output, default ./WRF_SA_Local):
    ./WRF_SA_Local/usda_scs/                  - USDA-SCS Peff + fraction rasters
    ./WRF_SA_Local/analysis_inputs/           - Saved precip, AWC and ETo inputs
    ./WRF_SA_Local/peff_by_method/{method}/   - All 8 methods' Peff rasters

Usage:
    # Default: 2000-2001 (24 months) so the first run is quick
    python wrf_south_america_local_example.py

    # The full WRF record (2000-2021, 264 months)
    python wrf_south_america_local_example.py --start-year 2000 --end-year 2021

    python wrf_south_america_local_example.py -w 8
    python wrf_south_america_local_example.py --gee-project your-project-id
    python wrf_south_america_local_example.py --geometry basin.geojson
    python wrf_south_america_local_example.py --download-only
    python wrf_south_america_local_example.py --calc-only
    python wrf_south_america_local_example.py -f   # force, ignore existing output
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path

import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from pycropwat import EffectivePrecipitation, open_local_precipitation
from pycropwat.methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip,
    farmwest_effective_precip,
    usda_scs_effective_precip,
    suet_effective_precip,
    ensemble_effective_precip,
)


# =============================================================================
# Configuration
# =============================================================================

# Examples/LocalPrecip/wrf_south_america_local_example.py -> repository root
REPO_ROOT = Path(__file__).resolve().parents[2]

# The example dataset lives next to the repository, not inside it
DEFAULT_PRECIP_DIR = REPO_ROOT.parent / 'pyCropWat_Data' / 'Precip'

# Local WRF precipitation
DEFAULT_PATTERN = 'Precip_*.tif'
DEFAULT_NODATA = -9999.0
WRF_SCALE_FACTOR = 1.0  # Already in mm

# Global AWC and ETo datasets (still downloaded from GEE)
AWC_ASSET = "projects/sat-io/open-datasets/FAO/HWSD_V2_SMU"
AWC_BAND = "AWC"
AWC_SCALE_FACTOR = 0.001  # FAO HWSD AWC is in mm/m -> convert to volumetric fraction
ETO_ASSET = "IDAHO_EPSCOR/TERRACLIMATE"
ETO_BAND = "pet"
ETO_SCALE_FACTOR = 0.1
ETO_IS_DAILY = False

# Physical parameters
ROOTING_DEPTH = 2  # meters
MAD_FACTOR = 1

# Default output base (kept small by the default 2000-2001 period)
DEFAULT_OUTPUT = Path('./WRF_SA_Local')

# All methods written by Step 2
ALL_METHODS = [
    'usda_scs',
    'cropwat',
    'fao_aglw',
    'fixed_percentage',
    'dependable_rainfall',
    'farmwest',
    'suet',
    'ensemble',
]


# =============================================================================
# Helpers
# =============================================================================

def output_dirs(output_base: Path) -> dict:
    """
    Build the output directory layout.

    Parameters
    ----------

    output_base : Path
        Root output directory.

    Returns
    -------
    dict
        Mapping with keys ``'usda_scs'``, ``'inputs'`` and ``'by_method'``.
    """
    return {
        'usda_scs': output_base / 'usda_scs',
        'inputs': output_base / 'analysis_inputs',
        'by_method': output_base / 'peff_by_method',
    }


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


def list_months(precip_dir: Path, pattern: str, nodata: float,
                start_year: int, end_year: int) -> list:
    """
    List the months available on disk within the requested year range.

    Also logs a short summary of the local precipitation source.

    Parameters
    ----------

    precip_dir : Path
        Directory of monthly rasters.

    pattern : str
        Glob applied inside ``precip_dir``.

    nodata : float
        Extra nodata sentinel.

    start_year : int
        First year of interest.

    end_year : int
        Last year of interest.

    Returns
    -------
    list of tuple
        Sorted ``(year, month)`` pairs present on disk in the range.
    """
    try:
        with open_local_precipitation(precip_dir, pattern=pattern, nodata=nodata) as src:
            min_lon, min_lat, max_lon, max_lat = src.bounds
            logger.info(f"  Local files   : {len(src.files)} ({src.kind}), "
                        f"{src.year_range[0]}-{src.year_range[1]}")
            logger.info(f"  Grid          : {src.shape} @ {src.resolution[0]:.6f} deg, "
                        f"CRS {src.crs}")
            logger.info(f"  Bounds (WGS84): {min_lon:.4f}, {min_lat:.4f}, "
                        f"{max_lon:.4f}, {max_lat:.4f}")
            return [
                (year, month) for year, month in src.available_months()
                if start_year <= year <= end_year
            ]
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Could not open the local precipitation data: {e}")
        logger.error("Check --precip-dir, --pattern and the file naming "
                     "(pyCropWat expects e.g. 'Precip_2005_07.tif').")
        sys.exit(1)


def resample_to(values: np.ndarray, target: np.ndarray) -> np.ndarray:
    """
    Resample an array onto the shape of another array if they differ.

    AWC and ETo are regridded onto the precipitation grid by pyCropWat, so this
    is only a defensive fallback for hand-assembled input directories.

    Parameters
    ----------

    values : np.ndarray
        Array to resample.

    target : np.ndarray
        Array whose shape is the target.

    Returns
    -------
    np.ndarray
        ``values`` resampled to ``target.shape``.
    """
    if values.shape == target.shape:
        return values
    from scipy.ndimage import zoom
    factors = (target.shape[0] / values.shape[0], target.shape[1] / values.shape[1])
    return zoom(values, factors, order=1)


# =============================================================================
# Step 1: Compute USDA-SCS from local precipitation + GEE AWC/ETo
# =============================================================================

def run_usda_scs(args, dirs: dict, months: list) -> None:
    """
    Calculate USDA-SCS effective precipitation and save all inputs.

    Precipitation is read from the local rasters; AWC and ETo are downloaded
    from Google Earth Engine and reprojected onto the local precipitation grid.

    This writes:
        - Effective precipitation (USDA-SCS)
        - Effective precipitation fraction (Peff / P)
        - The precipitation actually used (``precip_YYYY_MM.tif``)
        - AWC (``awc.tif``) and ETo (``eto_YYYY_MM.tif``)

    Parameters
    ----------

    args : argparse.Namespace
        Parsed command line arguments.

    dirs : dict
        Output directory layout from :func:`output_dirs`.

    months : list of tuple
        ``(year, month)`` pairs available on disk in the requested range.
    """
    expected_files = len(months)

    if not args.force and dirs['usda_scs'].exists():
        existing = list(dirs['usda_scs'].glob('effective_precip_[0-9]*.tif'))
        if len(existing) >= expected_files:
            logger.info(
                f"USDA-SCS output already complete "
                f"({len(existing)}/{expected_files} files); skipping. Use -f to redo."
            )
            return

    logger.info("Calculating USDA-SCS from LOCAL precipitation + GEE AWC/ETo...")
    logger.info(f"  Precipitation : {args.precip_dir} (pattern: {args.pattern}, "
                f"nodata: {args.nodata})")
    logger.info(f"  AWC           : {AWC_ASSET} (band: {AWC_BAND}, "
                f"scale={AWC_SCALE_FACTOR})")
    logger.info(f"  ETo           : {ETO_ASSET} (band: {ETO_BAND}, "
                f"scale={ETO_SCALE_FACTOR}, daily={ETO_IS_DAILY})")
    logger.info(f"  Rooting Depth : {ROOTING_DEPTH} m, MAD Factor: {MAD_FACTOR}")

    method_params = {
        'awc_asset': AWC_ASSET,
        'awc_band': AWC_BAND,
        'awc_scale_factor': AWC_SCALE_FACTOR,
        'eto_asset': ETO_ASSET,
        'eto_band': ETO_BAND,
        'eto_scale_factor': ETO_SCALE_FACTOR,
        'eto_is_daily': ETO_IS_DAILY,
        'rooting_depth': ROOTING_DEPTH,
        'mad_factor': MAD_FACTOR,
    }

    try:
        ep = EffectivePrecipitation(
            # Precipitation comes from disk - no asset_id / precip_band needed
            local_precip=str(args.precip_dir),
            local_precip_pattern=args.pattern,
            local_precip_nodata=args.nodata,
            # Optional clip: local vector file, or a GEE FeatureCollection asset ID
            geometry_path=args.geometry,
            start_year=args.start_year,
            end_year=args.end_year,
            precip_scale_factor=WRF_SCALE_FACTOR,
            gee_project=args.gee_project,
            method='usda_scs',
            method_params=method_params,
        )
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Could not set up the USDA-SCS run: {e}")
        logger.error("Adjust --start-year/--end-year to a period covered by "
                     "--precip-dir, or point --precip-dir at different data.")
        sys.exit(1)

    ep.process(
        output_dir=str(dirs['usda_scs']),
        n_workers=args.workers,
        save_inputs=True,
        input_dir=str(dirs['inputs']),
    )

    logger.info("USDA-SCS step complete")


# =============================================================================
# Step 2: Calculate All Methods from Saved Rasters
# =============================================================================

def calculate_all_methods(args, dirs: dict, months: list) -> None:
    """
    Calculate Peff for all 8 methods from the rasters saved by Step 1.

    Loads precipitation, AWC and ETo directly from the saved input files
    (written by ``save_inputs=True``), then calculates each method locally.
    No Earth Engine calls are made here.

    Parameters
    ----------

    args : argparse.Namespace
        Parsed command line arguments.

    dirs : dict
        Output directory layout from :func:`output_dirs`.

    months : list of tuple
        ``(year, month)`` pairs available on disk in the requested range.
    """
    import rioxarray

    input_dir = dirs['inputs']
    usda_scs_dir = dirs['usda_scs']

    # Verify input data exists
    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        logger.error("Run Step 1 first (python wrf_south_america_local_example.py)")
        return
    if not usda_scs_dir.exists():
        logger.error(f"USDA-SCS directory not found: {usda_scs_dir}")
        logger.error("Run Step 1 first (python wrf_south_america_local_example.py)")
        return

    logger.info("Calculating all 8 Peff methods from saved rasters...")

    # Load AWC data (static, loaded once)
    awc_file = input_dir / 'awc.tif'
    if awc_file.exists():
        da_awc = rioxarray.open_rasterio(awc_file, masked=True).squeeze('band', drop=True)
        awc_data = da_awc.values
        logger.info(f"  Loaded AWC: shape={awc_data.shape}, "
                    f"mean={np.nanmean(awc_data):.4f}")
    else:
        logger.warning(f"  AWC file not found: {awc_file}, using default 0.15")
        awc_data = None

    processed_count = 0
    missing_precip = 0

    for year, month in months:
        precip_file = input_dir / f'precip_{year}_{month:02d}.tif'
        peff_file = usda_scs_dir / f'effective_precip_{year}_{month:02d}.tif'

        if not precip_file.exists():
            missing_precip += 1
            continue

        # Skip if all method outputs already exist for this month
        if not args.force:
            all_exist = all(
                (dirs['by_method'] / name /
                 f'effective_precip_{year}_{month:02d}.tif').exists()
                for name in ALL_METHODS
            )
            if all_exist:
                processed_count += 1
                continue

        da_precip = rioxarray.open_rasterio(
            precip_file, masked=True).squeeze('band', drop=True)
        precip = da_precip.values
        da_template = da_precip

        # USDA-SCS Peff was already computed in Step 1; reuse it when present
        if peff_file.exists():
            da_peff = rioxarray.open_rasterio(
                peff_file, masked=True).squeeze('band', drop=True)
            peff_usda = da_peff.values
        else:
            peff_usda = None

        # Load ETo for this month
        eto_file = input_dir / f'eto_{year}_{month:02d}.tif'
        if eto_file.exists():
            da_eto = rioxarray.open_rasterio(
                eto_file, masked=True).squeeze('band', drop=True)
            eto = resample_to(da_eto.values, precip)
        else:
            logger.warning(f"  ETo missing for {year}-{month:02d}, using 120 mm default")
            eto = np.full_like(precip, 120.0)

        # Resample AWC if shapes differ
        if awc_data is not None:
            awc = resample_to(awc_data, precip)
        else:
            awc = np.full_like(precip, 0.15)

        # Calculate all methods (NaN precipitation propagates as NaN).
        # The ensemble averages its 6 component methods with np.nanmean, which
        # raises "Mean of empty slice" on every all-NaN pixel - here, the WRF
        # -9999 model-domain border. That is expected, so silence exactly that
        # message; every other RuntimeWarning still surfaces.
        with np.errstate(divide='ignore', invalid='ignore'), warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Mean of empty slice',
                                    category=RuntimeWarning)
            method_results = {
                'usda_scs': peff_usda if peff_usda is not None else usda_scs_effective_precip(
                    precip, eto, awc, ROOTING_DEPTH, MAD_FACTOR),
                'cropwat': cropwat_effective_precip(precip),
                'fao_aglw': fao_aglw_effective_precip(precip),
                'fixed_percentage': fixed_percentage_effective_precip(precip, 0.7),
                'dependable_rainfall': dependable_rainfall_effective_precip(precip, 0.75),
                'farmwest': farmwest_effective_precip(precip),
                'suet': suet_effective_precip(precip, eto),
                'ensemble': ensemble_effective_precip(
                    precip, eto, awc, ROOTING_DEPTH, 0.7, 0.75
                ),
            }

        # Save each method
        for method_name, peff_values in method_results.items():
            method_dir = dirs['by_method'] / method_name
            method_dir.mkdir(parents=True, exist_ok=True)

            output_file = method_dir / f'effective_precip_{year}_{month:02d}.tif'

            da_result = da_template.copy(data=peff_values)
            da_result = da_result.where(~np.isnan(precip))
            da_result.attrs = {
                'units': 'mm',
                'long_name': f'effective_precipitation_{method_name}',
                'method': method_name,
                'year': year,
                'month': month,
            }
            # Keep the CRS of the local precipitation grid (not necessarily EPSG:4326)
            da_result = da_result.rio.write_crs(da_template.rio.crs or "EPSG:4326")
            da_result = da_result.rio.write_nodata(np.nan, encoded=False)
            da_result.rio.to_raster(output_file)

        processed_count += 1
        if month == 1:
            logger.info(f"  Processed {year}...")

    if missing_precip:
        logger.warning(
            f"  {missing_precip} month(s) had no saved precip_YYYY_MM.tif and were "
            f"skipped. Re-run Step 1 (save_inputs=True) to produce them."
        )

    logger.info(f"Completed all 8 methods ({processed_count} months)")


# =============================================================================
# Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            'WRF effective precipitation for South America from LOCAL monthly '
            'rasters. Step 1 runs USDA-SCS (AWC/ETo from GEE), Step 2 calculates '
            'all other methods locally.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python wrf_south_america_local_example.py
  python wrf_south_america_local_example.py --start-year 2000 --end-year 2021
  python wrf_south_america_local_example.py --geometry basin.geojson -w 8
  python wrf_south_america_local_example.py --gee-project your-project-id
  python wrf_south_america_local_example.py --calc-only
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
    parser.add_argument('--geometry', type=str, default=None,
                        help='Optional area of interest: a local vector file '
                             '(.shp/.geojson), which clips the local rasters, or a GEE '
                             'FeatureCollection asset ID. Default: the full raster extent')
    parser.add_argument('--gee-project', type=str, default=None,
                        help='GEE project ID for authentication (AWC/ETo downloads)')
    parser.add_argument('--start-year', type=int, default=2000,
                        help='First year to process (default: 2000)')
    parser.add_argument('--end-year', type=int, default=2001,
                        help='Last year to process (default: 2001; use 2021 for the '
                             'full WRF record)')
    parser.add_argument('-w', '--workers', type=int, default=4,
                        help='Parallel workers for the AWC/ETo downloads (default: 4)')
    parser.add_argument('--download-only', action='store_true',
                        help='Only run Step 1 (USDA-SCS), skip the local calculation')
    parser.add_argument('--calc-only', action='store_true',
                        help='Only run Step 2 (assumes Step 1 is done)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Recompute even if output already exists')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT,
                        help=f'Output directory (default: {DEFAULT_OUTPUT})')
    return parser.parse_args()


def main():
    args = parse_args()

    args.precip_dir = Path(args.precip_dir).expanduser()
    dirs = output_dirs(Path(args.output))

    logger.info("=" * 70)
    logger.info("pyCropWat - WRF South America Effective Precipitation (local files)")
    logger.info("=" * 70)
    logger.info(f"  Precipitation : {args.precip_dir} (local, no GEE)")
    logger.info(f"  Study Area    : {args.geometry or 'full local raster extent'}")
    logger.info(f"  Period        : {args.start_year}-{args.end_year}")
    logger.info(f"  Output        : {Path(args.output).resolve()}")

    check_precip_dir(args.precip_dir)
    months = list_months(args.precip_dir, args.pattern, args.nodata,
                         args.start_year, args.end_year)
    logger.info(f"  Months to run : {len(months)}")
    logger.info("=" * 70)

    if not months:
        logger.error(
            f"No local precipitation files fall inside {args.start_year}-{args.end_year}. "
            f"Adjust --start-year/--end-year or --precip-dir."
        )
        sys.exit(1)

    if not args.calc_only:
        logger.info("Step 1: USDA-SCS from local precipitation + GEE AWC/ETo...")
        run_usda_scs(args, dirs, months)

    if not args.download_only:
        logger.info("Step 2: Calculating all 8 Peff methods from saved rasters...")
        calculate_all_methods(args, dirs, months)

    logger.info("=" * 70)
    logger.info("Done.")
    logger.info(f"  USDA-SCS rasters       : {dirs['usda_scs']}")
    logger.info(f"  Precip / AWC / ETo     : {dirs['inputs']}")
    logger.info(f"  All methods            : {dirs['by_method']}")
    logger.info("=" * 70)


if __name__ == '__main__':
    import atexit

    def _silence_excepthook():
        """Suppress 'Error in sys.excepthook' noise during interpreter shutdown."""
        sys.excepthook = lambda *_: None
        try:
            sys.stderr.close()
        except Exception:
            pass

    atexit.register(_silence_excepthook)
    main()

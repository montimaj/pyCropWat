"""
pyCropWat ERA5 South America Effective Precipitation Example
============================================================

This script demonstrates an efficient workflow for computing effective
precipitation using ERA5 data for a South American study region.

Workflow:
    1. Download rasters using USDA-SCS method (saves Peff, fraction, AWC, ETo)
    2. Recover total precipitation from Peff and fraction files
    3. Calculate all other 7 Peff methods locally from saved rasters
       (no additional GEE calls needed)

Methods computed (8 total, excludes PCML):
    1. USDA-SCS     - Downloaded from GEE (site-specific, requires AWC + ETo)
    2. Ensemble     - Mean of 6 methods (calculated locally)
    3. CROPWAT      - FAO CROPWAT method
    4. FAO/AGLW     - FAO Dependable Rainfall (80% exceedance)
    5. Fixed %      - Simple 70% fixed percentage
    6. Dependable   - FAO Dependable Rainfall (75% probability)
    7. FarmWest     - WSU irrigation scheduling formula
    8. TAGEM-SuET   - Turkish irrigation method (P - ETo, if P > 75 mm)

ERA5 Precipitation Asset:
    ECMWF/ERA5_LAND/MONTHLY_AGGR  (monthly, band: total_precipitation_sum, 01/2000-12/2021)

AWC and ETo Default Assets (Global):
    AWC: FAO HWSD v2 (projects/sat-io/open-datasets/FAO/HWSD_V2_SMU, band: AWC)
         Units: mm/m -> awc_scale_factor=0.001 to convert to volumetric fraction
    ETo: ECMWF/ERA5_LAND/MONTHLY_AGGR (projects/climate-engine-pro/assets/ce-ag-era5-v2/daily,
         band: potential_evaporation_sum, scale=1000 for mm, daily=False)

Study Area:
    Rio de la Plata Basin (projects/ssebop-471916/assets/Riodelaplata)

Output Structure:
    ./ERA5_SA/usda_scs/                  - USDA-SCS Peff + fraction rasters
    ./ERA5_SA/analysis_inputs/           - Saved AWC and ETo inputs
    ./ERA5_SA/peff_by_method/{method}/   - All 8 methods' Peff rasters

Usage:
    python era5_south_america_example.py
    python era5_south_america_example.py -w 8
    python era5_south_america_example.py --gee-project your-project-id
    python era5_south_america_example.py --download-only
    python era5_south_america_example.py --calc-only
    python era5_south_america_example.py -f   # force re-download
"""

import argparse
import logging
from pathlib import Path

import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from pycropwat import EffectivePrecipitation
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

# GEE Configuration
GEE_PROJECT = None  # Set to your GEE project ID, or None to use default
STUDY_AREA_ASSET = "projects/ssebop-471916/assets/riodelaplata_guarani"

# ERA5 precipitation asset
ERA5_ASSET = "ECMWF/ERA5_LAND/MONTHLY_AGGR"
ERA5_BAND = "total_precipitation_sum"
ERA5_SCALE_FACTOR = 1000  # ERA5-Land precipitation is in m -> convert to mm

# Global AWC and ETo datasets
AWC_ASSET = "projects/sat-io/open-datasets/FAO/HWSD_V2_SMU"
AWC_BAND = "AWC"
AWC_SCALE_FACTOR = 0.001  # FAO HWSD AWC is in mm/m -> convert to volumetric fraction
ETO_ASSET = "ECMWF/ERA5_LAND/MONTHLY_AGGR"
ETO_BAND = "potential_evaporation_sum"
# ERA5-Land ETo is in m -> convert to mm, ERA5 evaporation is negative, so use negative scale factor to convert to positive mm
ETO_SCALE_FACTOR = -1000
ETO_IS_DAILY = False

# Physical parameters
ROOTING_DEPTH = 2  # meters
MAD_FACTOR = 1

# Time period
START_YEAR = 2000
END_YEAR = 2021

# Output directories
OUTPUT_BASE = Path('./ERA5_SA')
USDA_SCS_DIR = OUTPUT_BASE / 'usda_scs'
INPUT_DIR = OUTPUT_BASE / 'analysis_inputs'
PEFF_BY_METHOD_DIR = OUTPUT_BASE / 'peff_by_method'


# =============================================================================
# Step 1: Download USDA-SCS Rasters from GEE
# =============================================================================

def download_usda_scs(n_workers: int = 10, skip_if_exists: bool = True):
    """
    Download rasters using USDA-SCS method.

    This downloads and saves:
        - Effective precipitation (USDA-SCS)
        - Effective precipitation fraction (Peff / P)
        - AWC data (FAO HWSD)
        - ETo data (ERA5-Land)

    Other methods are then calculated locally from these saved files.
    """
    expected_files = (END_YEAR - START_YEAR + 1) * 12

    if skip_if_exists and USDA_SCS_DIR.exists():
        existing = list(USDA_SCS_DIR.glob('effective_precip_[0-9]*.tif'))
        if len(existing) >= expected_files:
            logger.info(
                f"USDA-SCS download already complete ({len(existing)}/{expected_files} files)"
            )
            return

    logger.info("Downloading ERA5 rasters using USDA-SCS method...")
    logger.info(f"  Precipitation : {ERA5_ASSET} (band: {ERA5_BAND})")
    logger.info(f"  AWC           : {AWC_ASSET} (band: {AWC_BAND}, scale={AWC_SCALE_FACTOR})")
    logger.info(f"  ETo           : {ETO_ASSET} (band: {ETO_BAND}, scale={ETO_SCALE_FACTOR}, daily={ETO_IS_DAILY})")
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

    ep = EffectivePrecipitation(
        asset_id=ERA5_ASSET,
        precip_band=ERA5_BAND,
        gee_geometry_asset=STUDY_AREA_ASSET,
        start_year=START_YEAR,
        end_year=END_YEAR,
        precip_scale_factor=ERA5_SCALE_FACTOR,
        gee_project=GEE_PROJECT,
        method='usda_scs',
        method_params=method_params,
    )

    ep.process(
        output_dir=str(USDA_SCS_DIR),
        n_workers=n_workers,
        save_inputs=True,
        input_dir=str(INPUT_DIR),
    )

    logger.info("USDA-SCS download complete")


# =============================================================================
# Step 2: Calculate All Methods from Saved Rasters
# =============================================================================

def calculate_all_methods():
    """
    Calculate Peff for all 8 methods from saved USDA-SCS rasters.

    Loads precipitation, AWC, and ETo directly from saved input files
    (written by save_inputs=True), then calculates each method locally.
    """
    import rioxarray

    # Verify input data exists
    if not INPUT_DIR.exists():
        logger.error(f"Input directory not found: {INPUT_DIR}")
        logger.error("Run download step first with save_inputs=True")
        return
    if not USDA_SCS_DIR.exists():
        logger.error(f"USDA-SCS directory not found: {USDA_SCS_DIR}")
        logger.error("Run download step first (python era5_south_america_example.py)")
        return

    logger.info("Calculating all 8 Peff methods from saved rasters...")

    # Load AWC data (static, loaded once)
    awc_file = INPUT_DIR / 'awc.tif'
    if awc_file.exists():
        da_awc = rioxarray.open_rasterio(awc_file).squeeze('band', drop=True)
        awc_data = da_awc.values
        logger.info(f"  Loaded AWC: shape={awc_data.shape}, mean={np.nanmean(awc_data):.4f}")
    else:
        logger.warning(f"  AWC file not found: {awc_file}, using default 0.15")
        awc_data = None

    processed_count = 0
    for year in range(START_YEAR, END_YEAR + 1):
        for month in range(1, 13):
            # Load saved precipitation directly from input files
            precip_file = INPUT_DIR / f'precip_{year}_{month:02d}.tif'
            peff_file = USDA_SCS_DIR / f'effective_precip_{year}_{month:02d}.tif'

            if not precip_file.exists():
                if not peff_file.exists():
                    continue
                # Fallback: recover precip from Peff / fraction
                logger.warning(f"  Precip file missing for {year}-{month:02d}, "
                               "recovering from Peff/fraction")
                da_peff = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
                peff_usda = da_peff.values
                fraction_file = USDA_SCS_DIR / f'effective_precip_fraction_{year}_{month:02d}.tif'
                if fraction_file.exists():
                    da_frac = rioxarray.open_rasterio(fraction_file).squeeze('band', drop=True)
                    frac = da_frac.values.copy()
                    frac[frac == 0] = 1
                    precip = peff_usda / frac
                else:
                    precip = peff_usda * 1.5
            else:
                da_precip = rioxarray.open_rasterio(precip_file).squeeze('band', drop=True)
                precip = da_precip.values
                # Load USDA-SCS Peff (already computed)
                if peff_file.exists():
                    da_peff = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
                    peff_usda = da_peff.values
                else:
                    peff_usda = None

            # Use precip DataArray as template for output
            if precip_file.exists():
                da_template = da_precip
            else:
                da_template = da_peff

            # Skip if all method outputs already exist for this month
            all_exist = all(
                (PEFF_BY_METHOD_DIR / m / f'effective_precip_{year}_{month:02d}.tif').exists()
                for m in ['usda_scs', 'cropwat', 'fao_aglw', 'fixed_percentage',
                          'dependable_rainfall', 'farmwest', 'suet', 'ensemble']
            )
            if all_exist:
                processed_count += 1
                continue

            # Load ETo for this month
            eto_file = INPUT_DIR / f'eto_{year}_{month:02d}.tif'
            if eto_file.exists():
                da_eto = rioxarray.open_rasterio(eto_file).squeeze('band', drop=True)
                eto = da_eto.values
                if eto.shape != precip.shape:
                    from scipy.ndimage import zoom
                    factors = (precip.shape[0] / eto.shape[0],
                               precip.shape[1] / eto.shape[1])
                    eto = zoom(eto, factors, order=1)
            else:
                logger.warning(f"  ETo missing for {year}-{month:02d}, using 120 mm default")
                eto = np.full_like(precip, 120.0)

            # Resample AWC if shapes differ
            if awc_data is not None:
                awc = awc_data
                if awc_data.shape != precip.shape:
                    from scipy.ndimage import zoom
                    factors = (precip.shape[0] / awc_data.shape[0],
                               precip.shape[1] / awc_data.shape[1])
                    awc = zoom(awc_data, factors, order=1)
            else:
                awc = np.full_like(precip, 0.15)

            # Calculate all methods
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
                method_dir = PEFF_BY_METHOD_DIR / method_name
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
                da_result = da_result.rio.write_crs("EPSG:4326")
                da_result.rio.to_raster(output_file)

            processed_count += 1
            if month == 1:
                logger.info(f"  Processed {year}...")

    logger.info(f"Completed all 8 methods ({processed_count} months)")


# =============================================================================
# Main
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            'ERA5 effective precipitation for South America. '
            'Downloads USDA-SCS from GEE, then calculates all other methods locally.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-w', '--workers', type=int, default=10,
                        help='Number of parallel workers for GEE download (default: 10)')
    parser.add_argument('--gee-project', type=str, default=None,
                        help='GEE project ID for authentication')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Force re-download even if output exists')
    parser.add_argument('--download-only', action='store_true',
                        help='Only download USDA-SCS rasters, skip local calculation')
    parser.add_argument('--calc-only', action='store_true',
                        help='Only calculate methods locally (assumes download is done)')
    return parser.parse_args()


def main():
    args = parse_args()

    global GEE_PROJECT
    if args.gee_project:
        GEE_PROJECT = args.gee_project

    logger.info("=" * 70)
    logger.info("pyCropWat - ERA5 South America Effective Precipitation")
    logger.info("=" * 70)
    logger.info(f"  ERA5 Asset    : {ERA5_ASSET} (band: {ERA5_BAND})")
    logger.info(f"  Study Area    : {STUDY_AREA_ASSET}")
    logger.info(f"  Period        : {START_YEAR}-{END_YEAR}")
    logger.info(f"  Output        : {OUTPUT_BASE.resolve()}")
    logger.info("=" * 70)

    if not args.calc_only:
        logger.info("Step 1: Downloading USDA-SCS rasters from GEE...")
        download_usda_scs(
            n_workers=args.workers,
            skip_if_exists=not args.force,
        )

    if not args.download_only:
        logger.info("Step 2: Calculating all 8 Peff methods from saved rasters...")
        calculate_all_methods()

    logger.info("=" * 70)
    logger.info("Done.")
    logger.info(f"  USDA-SCS rasters : {USDA_SCS_DIR}")
    logger.info(f"  AWC / ETo inputs : {INPUT_DIR}")
    logger.info(f"  All methods      : {PEFF_BY_METHOD_DIR}")
    logger.info("=" * 70)


if __name__ == '__main__':
    import sys
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

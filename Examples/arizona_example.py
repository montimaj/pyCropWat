"""
pyCropWat Arizona USDA-SCS Workflow Example
============================================

This script demonstrates effective precipitation calculation for Arizona, USA
comparing U.S.-based and global climate datasets with all 9 available methods.

The workflow uses an efficient approach:
1. Download rasters once using USDA-SCS method (saves P, AWC, ETo)
2. Calculate all 9 Peff methods locally from saved rasters
3. Perform temporal aggregation, statistical analysis, and visualization
4. Compare U.S. vs Global datasets and all Peff methods

This approach is more efficient than downloading for each method separately.

Precipitation Data Sources:
  U.S. Datasets:
    - GridMET (IDAHO_EPSCOR/GRIDMET) - ~4km daily
    - PRISM (projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY) - ~800m monthly
  Global Datasets:
    - ERA5-Land (ECMWF/ERA5_LAND/MONTHLY_AGGR) - ~11km monthly
    - TerraClimate (IDAHO_EPSCOR/TERRACLIMATE) - ~4km monthly

USDA-SCS Required Data (downloaded and saved for local calculation):
  U.S. Datasets:
    - AWC: SSURGO (projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite)
    - ETo: GridMET Monthly (projects/openet/assets/reference_et/conus/gridmet/monthly/v1)
  Global Datasets:
    - AWC: FAO HWSD v2 (projects/sat-io/open-datasets/FAO/HWSD_V2_SMU, band: AWC)
    - ETo: AgERA5 Daily (projects/climate-engine-pro/assets/ce-ag-era5-v2/daily, 
           band: ReferenceET_PenmanMonteith_FAO56)

Effective Precipitation Methods Compared (8 total):
- CROPWAT - Method from FAO CROPWAT
- FAO/AGLW - FAO Dependable Rainfall (80% exceedance)
- Fixed Percentage (70%) - Simple empirical method
- Dependable Rainfall (75% probability) - Statistical approach
- FarmWest - WSU irrigation scheduling formula
- USDA-SCS - Site-specific method with AWC and ETo
- TAGEM-SuET - Turkish Irrigation Management System (P - ETo, if P > 75mm)
- Ensemble - Mean of 6 methods (default, excludes TAGEM-SuET and PCML)

Study Area:
- Arizona (users/montimajumdar/AZ)

Output Structure:
- ./Arizona/AZ_{Dataset}_USDA_SCS/ - USDA-SCS Peff downloads
- ./Arizona/analysis_inputs/AZ_{Dataset}_USDA_SCS/ - AWC and ETo inputs
- ./Arizona/analysis_outputs/peff_by_method/{Dataset}/{Method}/ - All 8 methods

Requirements:
- Google Earth Engine account with access to the study area asset
- pycropwat package installed

Usage:
    python arizona_usda_scs_example.py # default 4 workers
    python arizona_usda_scs_example.py -f -w 32 # 32 workers, force reprocess
    python arizona_usda_scs_example.py --analysis-only
    python arizona_usda_scs_example.py --gee-project your-project-id
"""

import os
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import pyCropWat modules
from pycropwat import EffectivePrecipitation
from pycropwat.analysis import (
    TemporalAggregator,
    StatisticalAnalyzer,
    Visualizer,
    export_to_netcdf
)


# =============================================================================
# Configuration
# =============================================================================

# GEE Configuration
GEE_PROJECT = None  # Set to your GEE project ID, or None to use default
STUDY_AREA_GEOJSON = "./AZ.geojson"  # Local GeoJSON file for Arizona boundary

# USDA-SCS Required Assets
# U.S. datasets (high-resolution SSURGO and GridMET)
US_AWC_ASSET = "projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite"
US_AWC_BAND = None  # Single band asset, no band name needed
US_ETO_ASSET = "projects/openet/assets/reference_et/conus/gridmet/monthly/v1"
US_ETO_BAND = "eto"

# Global datasets (FAO HWSD AWC and AgERA5 ETo for global coverage)
GLOBAL_AWC_ASSET = "projects/sat-io/open-datasets/FAO/HWSD_V2_SMU"  # FAO Harmonized World Soil Database
GLOBAL_AWC_BAND = "AWC"  # Available Water Capacity
GLOBAL_ETO_ASSET = "projects/climate-engine-pro/assets/ce-ag-era5-v2/daily"  # AgERA5 daily
GLOBAL_ETO_BAND = "ReferenceET_PenmanMonteith_FAO56"  # FAO-56 Penman-Monteith ETo
GLOBAL_ETO_IS_DAILY = True  # AgERA5 is daily data, needs monthly aggregation

ROOTING_DEPTH = 1.0  # meters

# Region-specific output directory
REGION_DIR = Path('./Arizona')

# Precipitation datasets configuration
# U.S. Datasets (high-resolution, use U.S. AWC and ETo)
US_DATASETS = {
    'GridMET': {
        'asset_id': 'IDAHO_EPSCOR/GRIDMET',
        'precip_band': 'pr',
        'scale_factor': 1.0,  # Already in mm
        'output_dir': './Arizona/AZ_GridMET_USDA_SCS',
        'description': 'University of Idaho GridMET ~4km daily aggregated to monthly',
        'region': 'US',
        'awc_asset': US_AWC_ASSET,
        'awc_band': US_AWC_BAND,
        'eto_asset': US_ETO_ASSET,
        'eto_band': US_ETO_BAND,
        'eto_is_daily': False,  # GridMET ETo is monthly
        'eto_scale_factor': 1.0  # Already in mm
    },
    'PRISM': {
        'asset_id': 'projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY',
        'precip_band': 'ppt',
        'scale_factor': 1.0,  # Already in mm
        'output_dir': './Arizona/AZ_PRISM_USDA_SCS',
        'description': 'Oregon State PRISM 800m monthly precipitation',
        'region': 'US',
        'awc_asset': US_AWC_ASSET,
        'awc_band': US_AWC_BAND,
        'eto_asset': US_ETO_ASSET,
        'eto_band': US_ETO_BAND,
        'eto_is_daily': False,  # GridMET ETo is monthly
        'eto_scale_factor': 1.0  # Already in mm
    }
}

# Global Datasets (use FAO HWSD AWC and AgERA5 ETo)
GLOBAL_DATASETS = {
    'ERA5Land': {
        'asset_id': 'ECMWF/ERA5_LAND/MONTHLY_AGGR',
        'precip_band': 'total_precipitation_sum',
        'scale_factor': 1000.0,  # Convert m to mm
        'output_dir': './Arizona/AZ_ERA5Land_USDA_SCS',
        'description': 'ECMWF ERA5-Land ~11km monthly reanalysis',
        'region': 'Global',
        'awc_asset': GLOBAL_AWC_ASSET,
        'awc_band': GLOBAL_AWC_BAND,
        'eto_asset': GLOBAL_ETO_ASSET,
        'eto_band': GLOBAL_ETO_BAND,
        'eto_is_daily': GLOBAL_ETO_IS_DAILY,  # AgERA5 is daily
        'eto_scale_factor': 1.0  # AgERA5 ETo is in mm/day
    },
    'TerraClimate': {
        'asset_id': 'IDAHO_EPSCOR/TERRACLIMATE',
        'precip_band': 'pr',
        'scale_factor': 1.0,  # Already in mm
        'output_dir': './Arizona/AZ_TerraClimate_USDA_SCS',
        'description': 'University of Idaho TerraClimate ~4km monthly',
        'region': 'Global',
        'awc_asset': GLOBAL_AWC_ASSET,
        'awc_band': GLOBAL_AWC_BAND,
        'eto_asset': GLOBAL_ETO_ASSET,
        'eto_band': GLOBAL_ETO_BAND,
        'eto_is_daily': GLOBAL_ETO_IS_DAILY,  # AgERA5 is daily
        'eto_scale_factor': 1.0  # AgERA5 ETo is in mm/day
    }
}

# Combined datasets dictionary
DATASETS = {**US_DATASETS, **GLOBAL_DATASETS}

# Effective precipitation methods to compare (8 methods - excludes PCML)
PEFF_METHODS = ['cropwat', 'fao_aglw', 'fixed_percentage', 'dependable_rainfall', 'farmwest', 'usda_scs', 'suet', 'ensemble']

# Time period
START_YEAR = 1985
END_YEAR = 2025

# Climatology period (for anomaly calculations)
CLIM_START = 1985
CLIM_END = 2020

# Output scale (resolution in meters)
OUTPUT_SCALE = None  # Use native resolution

# Analysis output directories (within region folder)
ANALYSIS_DIR = REGION_DIR / 'analysis_outputs'
ANNUAL_DIR = ANALYSIS_DIR / 'annual'
CLIMATOLOGY_DIR = ANALYSIS_DIR / 'climatology'
ANOMALY_DIR = ANALYSIS_DIR / 'anomalies'
TREND_DIR = ANALYSIS_DIR / 'trend'
FIGURES_DIR = ANALYSIS_DIR / 'figures'
COMPARISON_DIR = ANALYSIS_DIR / 'comparisons'
ZONAL_DIR = ANALYSIS_DIR / 'zonal_stats'
METHOD_COMPARISON_DIR = ANALYSIS_DIR / 'method_comparison'
US_VS_GLOBAL_DIR = ANALYSIS_DIR / 'us_vs_global'

# Zone file for zonal statistics (Arizona counties or irrigation districts)
ZONE_FILE = ANALYSIS_DIR / 'az_zones.geojson'


def create_output_directories():
    """Create all necessary output directories."""
    directories = [
        REGION_DIR,
        ANALYSIS_DIR,
        ANNUAL_DIR,
        CLIMATOLOGY_DIR,
        ANOMALY_DIR,
        TREND_DIR,
        FIGURES_DIR,
        COMPARISON_DIR,
        ZONAL_DIR,
        METHOD_COMPARISON_DIR,
        US_VS_GLOBAL_DIR
    ]
    
    # Create method-specific directories for each dataset
    for dataset_name in DATASETS.keys():
        for method_name in PEFF_METHODS:
            directories.append(ANALYSIS_DIR / 'peff_by_method' / dataset_name / method_name)
            directories.append(ANNUAL_DIR / dataset_name / method_name)
            directories.append(CLIMATOLOGY_DIR / dataset_name / method_name)
            directories.append(ANOMALY_DIR / dataset_name / method_name)
            directories.append(TREND_DIR / dataset_name / method_name)
        # Also create dataset-specific input directories
        directories.append(REGION_DIR / 'analysis_inputs' / Path(DATASETS[dataset_name]['output_dir']).name)
    
    for d in directories:
        d.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directories created in {REGION_DIR}")


def create_sample_zones():
    """
    Create sample zone polygons for zonal statistics demonstration.
    
    Creates three zones representing different Arizona regions:
    - Central Arizona (Phoenix metropolitan area)
    - Southern Arizona (Tucson/Santa Cruz Valley)
    - Northern Arizona (Flagstaff/Colorado Plateau)
    """
    import json
    
    if ZONE_FILE.exists():
        logger.info(f"Zone file already exists: {ZONE_FILE}")
        return
    
    logger.info("Creating sample zone polygons...")
    
    # Define three zones within Arizona
    zones_geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 1,
                    "name": "Central AZ",
                    "description": "Central Arizona - Phoenix Metropolitan / Salt River Valley"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-113.0, 33.0],
                        [-113.0, 34.0],
                        [-111.0, 34.0],
                        [-111.0, 33.0],
                        [-113.0, 33.0]
                    ]]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 2,
                    "name": "Southern AZ",
                    "description": "Southern Arizona - Tucson / Santa Cruz Valley"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-112.0, 31.5],
                        [-112.0, 32.5],
                        [-110.0, 32.5],
                        [-110.0, 31.5],
                        [-112.0, 31.5]
                    ]]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 3,
                    "name": "Northern AZ",
                    "description": "Northern Arizona - Flagstaff / Colorado Plateau"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-112.5, 34.5],
                        [-112.5, 36.0],
                        [-110.5, 36.0],
                        [-110.5, 34.5],
                        [-112.5, 34.5]
                    ]]
                }
            }
        ]
    }
    
    with open(ZONE_FILE, 'w') as f:
        json.dump(zones_geojson, f, indent=2)
    
    logger.info(f"Created zone file: {ZONE_FILE}")


# =============================================================================
# Step 1: Download Rasters Using USDA-SCS Method
# =============================================================================

def download_usda_scs_rasters(
    dataset_name: str,
    skip_if_exists: bool = True,
    n_workers: int = 4
):
    """
    Download rasters using USDA-SCS method for a given dataset.
    
    This downloads and saves:
    - Effective precipitation (USDA-SCS method)
    - Effective precipitation fraction (for back-calculating total P)
    - AWC data (SSURGO for US, FAO HWSD for Global)
    - ETo data (gridMET for US, AgERA5 for Global)
    
    Other methods will be calculated locally from these saved rasters.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate')
    skip_if_exists : bool
        Skip processing if output files already exist
    n_workers : int
        Number of parallel workers for processing
    """
    config = DATASETS[dataset_name]
    output_dir = Path(config['output_dir'])
    
    # Check if data already exists
    if skip_if_exists and output_dir.exists():
        existing_files = list(output_dir.glob('effective_precip_[0-9]*.tif'))
        expected_files = (END_YEAR - START_YEAR + 1) * 12
        if len(existing_files) >= expected_files * 0.9:  # Allow 10% tolerance
            logger.info(f"Skipping {dataset_name} download - data already exists ({len(existing_files)} files)")
            return
    
    # Get dataset-specific AWC and ETo assets
    awc_asset = config.get('awc_asset', US_AWC_ASSET)
    awc_band = config.get('awc_band', US_AWC_BAND)
    eto_asset = config.get('eto_asset', US_ETO_ASSET)
    eto_band = config.get('eto_band', US_ETO_BAND)
    eto_is_daily = config.get('eto_is_daily', False)
    eto_scale_factor = config.get('eto_scale_factor', 1.0)
    region = config.get('region', 'US')
    
    logger.info(f"Downloading {dataset_name} ({region}) rasters using USDA-SCS method...")
    logger.info(f"  This will save P, AWC, and ETo data for local method calculation")
    logger.info(f"  Precipitation: {config['asset_id']} (band: {config['precip_band']})")
    logger.info(f"  AWC: {awc_asset}" + (f" (band: {awc_band})" if awc_band else ""))
    logger.info(f"  ETo: {eto_asset} (band: {eto_band}, daily={eto_is_daily}, scale={eto_scale_factor})")
    logger.info(f"  Rooting Depth: {ROOTING_DEPTH} m")
    
    # USDA-SCS specific parameters
    method_params = {
        'awc_asset': awc_asset,
        'awc_band': awc_band,
        'eto_asset': eto_asset,
        'eto_band': eto_band,
        'eto_is_daily': eto_is_daily,
        'eto_scale_factor': eto_scale_factor,
        'rooting_depth': ROOTING_DEPTH
    }
    
    ep = EffectivePrecipitation(
        asset_id=config['asset_id'],
        precip_band=config['precip_band'],
        geometry_path=STUDY_AREA_GEOJSON,
        start_year=START_YEAR,
        end_year=END_YEAR,
        precip_scale_factor=config['scale_factor'],
        scale=OUTPUT_SCALE,
        gee_project=GEE_PROJECT,
        method='usda_scs',
        method_params=method_params
    )
    
    # Input data will always be saved to region-specific analysis_inputs folder
    # e.g., ./Arizona/analysis_inputs/AZ_GridMET_USDA_SCS/
    input_dir = REGION_DIR / 'analysis_inputs' / output_dir.name
    
    ep.process(
        output_dir=str(output_dir),
        n_workers=n_workers,
        save_inputs=True,  # Always save inputs for local method calculations
        input_dir=str(input_dir)
    )
    
    logger.info(f"Completed downloading {dataset_name} rasters")


def calculate_all_methods_from_rasters(dataset_name: str):
    """
    Calculate effective precipitation using all 8 methods from saved rasters.
    
    Loads the precipitation, AWC, and ETo data saved during USDA-SCS processing
    and calculates Peff for each method locally (without additional GEE calls).
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate')
    """
    import rioxarray
    import numpy as np
    import xarray as xr
    from pycropwat.methods import (
        cropwat_effective_precip,
        fao_aglw_effective_precip,
        fixed_percentage_effective_precip,
        dependable_rainfall_effective_precip,
        farmwest_effective_precip,
        usda_scs_effective_precip,
        suet_effective_precip,
        ensemble_effective_precip
    )
    
    config = DATASETS[dataset_name]
    peff_dir = Path(config['output_dir'])
    input_dir = REGION_DIR / 'analysis_inputs' / peff_dir.name
    
    # Check if input data exists
    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        logger.error(f"Run download_usda_scs_rasters('{dataset_name}') first")
        return
    
    logger.info(f"Calculating all 9 Peff methods for {dataset_name} from saved rasters...")
    
    # Load AWC data (static, only once)
    awc_file = input_dir / 'awc.tif'
    if awc_file.exists():
        da_awc = rioxarray.open_rasterio(awc_file).squeeze('band', drop=True)
        awc_data = da_awc.values
        logger.info(f"  Loaded AWC data: shape={awc_data.shape}, mean={np.nanmean(awc_data):.2f} mm/m")
    else:
        logger.warning(f"  AWC file not found: {awc_file}, using representative value")
        awc_data = None
    
    # Process each month
    processed_count = 0
    for year in range(START_YEAR, END_YEAR + 1):
        for month in range(1, 13):
            # Load USDA-SCS Peff and fraction files
            peff_file = peff_dir / f'effective_precip_{year}_{month:02d}.tif'
            fraction_file = peff_dir / f'effective_precip_fraction_{year}_{month:02d}.tif'
            
            if not peff_file.exists():
                continue
            
            # Load USDA-SCS Peff (already computed)
            da_peff = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
            
            # Calculate total precipitation from fraction
            if fraction_file.exists():
                da_fraction = rioxarray.open_rasterio(fraction_file).squeeze('band', drop=True)
                # Replace zero fractions to avoid division by zero (zero fraction means zero precip)
                fraction_vals = da_fraction.values
                fraction_vals[fraction_vals == 0] = 1
                precip = da_peff.values / fraction_vals
            else:
                # Approximate (shouldn't happen if USDA-SCS processed correctly)
                logger.warning(f"Fraction file not found for {year}-{month:02d}, estimating precipitation")
                precip = da_peff.values * 1.5
            
            # Load ETo for this month
            eto_file = input_dir / f'eto_{year}_{month:02d}.tif'
            if eto_file.exists():
                da_eto = rioxarray.open_rasterio(eto_file).squeeze('band', drop=True)
                eto_data = da_eto.values
                # Resample if shapes don't match
                if eto_data.shape != precip.shape:
                    from scipy.ndimage import zoom
                    zoom_factors = (precip.shape[0] / eto_data.shape[0],
                                   precip.shape[1] / eto_data.shape[1])
                    eto_data = zoom(eto_data, zoom_factors, order=1)
            else:
                logger.warning(f"  ETo file not found for {year}-{month:02d}, using regional mean")
                eto_data = np.full_like(precip, 180.0)  # Representative for Arizona
            
            # Resample AWC if needed
            if awc_data is not None:
                awc_resampled = awc_data
                if awc_data.shape != precip.shape:
                    from scipy.ndimage import zoom
                    zoom_factors = (precip.shape[0] / awc_data.shape[0],
                                   precip.shape[1] / awc_data.shape[1])
                    awc_resampled = zoom(awc_data, zoom_factors, order=1)
            else:
                awc_resampled = np.full_like(precip, 150.0)  # Representative value
            
            # Calculate Peff for each method
            method_results = {
                'cropwat': cropwat_effective_precip(precip),
                'fao_aglw': fao_aglw_effective_precip(precip),
                'fixed_percentage': fixed_percentage_effective_precip(precip, 0.7),
                'dependable_rainfall': dependable_rainfall_effective_precip(precip, 0.75),
                'farmwest': farmwest_effective_precip(precip),
                'usda_scs': da_peff.values,  # Already computed
                'suet': suet_effective_precip(precip, eto_data),
                'ensemble': ensemble_effective_precip(precip, eto_data, awc_resampled, ROOTING_DEPTH, 0.7, 0.75)
            }
            
            # Save each method's output
            for method_name, peff_values in method_results.items():
                method_dir = ANALYSIS_DIR / 'peff_by_method' / dataset_name / method_name
                method_dir.mkdir(parents=True, exist_ok=True)
                
                output_file = method_dir / f'effective_precip_{year}_{month:02d}.tif'
                
                # Create DataArray with same coordinates as template
                da_result = da_peff.copy(data=peff_values)
                da_result = da_result.where(~np.isnan(da_peff.values))
                da_result.attrs = {
                    'units': 'mm',
                    'long_name': f'effective_precipitation_{method_name}',
                    'method': method_name,
                    'dataset': dataset_name,
                    'year': year,
                    'month': month
                }
                da_result = da_result.rio.write_crs("EPSG:4326")
                da_result.rio.to_raster(output_file)
            
            processed_count += 1
            if month == 1:
                logger.info(f"  Processed {year}...")
    
    logger.info(f"Completed calculation of all 8 methods for {dataset_name} ({processed_count} months)")


# =============================================================================
# Step 2: Temporal Aggregation
# =============================================================================

def run_temporal_aggregation(dataset_name: str):
    """
    Run temporal aggregation for a dataset.
    
    Creates:
    - Annual totals
    - Monsoon season aggregations (July-September)
    - Winter season aggregations (January-February)
    - Monthly climatology
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate').
        
    Returns
    -------
    TemporalAggregator
        The aggregator instance used for the dataset.
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    logger.info(f"Running temporal aggregation for {dataset_name}...")
    
    # Use pattern that excludes fraction files
    agg = TemporalAggregator(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    # Create dataset-specific output subdirectories
    dataset_annual_dir = ANNUAL_DIR / dataset_name
    dataset_clim_dir = CLIMATOLOGY_DIR / dataset_name
    dataset_annual_dir.mkdir(exist_ok=True)
    dataset_clim_dir.mkdir(exist_ok=True)
    
    # Annual totals
    logger.info("Computing annual totals...")
    for year in range(START_YEAR, END_YEAR + 1):
        output_path = dataset_annual_dir / f'annual_{year}.tif'
        if not output_path.exists():
            try:
                agg.annual_aggregate(
                    year=year,
                    method='sum',
                    output_path=str(output_path)
                )
            except Exception as e:
                logger.warning(f"Could not aggregate year {year}: {e}")
    
    # Monthly climatology
    logger.info("Computing monthly climatology...")
    agg.multi_year_climatology(
        start_year=CLIM_START,
        end_year=CLIM_END,
        output_dir=str(dataset_clim_dir)
    )
    
    # Monsoon season aggregations (July-September - Arizona monsoon)
    logger.info("Computing monsoon season totals...")
    monsoon_dir = ANALYSIS_DIR / 'monsoon_season' / dataset_name
    monsoon_dir.mkdir(parents=True, exist_ok=True)
    
    for year in range(START_YEAR, END_YEAR + 1):
        output_path = monsoon_dir / f'monsoon_{year}.tif'
        if not output_path.exists():
            try:
                agg.custom_aggregate(
                    year=year,
                    months=[7, 8, 9],  # Jul-Sep (Arizona Monsoon)
                    method='sum',
                    output_path=str(output_path)
                )
            except Exception as e:
                logger.warning(f"Could not aggregate monsoon season {year}: {e}")
    
    # Winter season (Dec-Feb - Pacific storms)
    logger.info("Computing winter season totals...")
    winter_dir = ANALYSIS_DIR / 'winter_season' / dataset_name
    winter_dir.mkdir(parents=True, exist_ok=True)
    
    for year in range(START_YEAR + 1, END_YEAR + 1):  # Start from year+1 for Dec of prev year
        output_path = winter_dir / f'winter_{year}.tif'
        if not output_path.exists():
            try:
                # Note: For winter (Dec-Feb), we'd need cross-year handling
                # Here we use Jan-Feb only for simplicity
                agg.custom_aggregate(
                    year=year,
                    months=[1, 2],  # Jan-Feb
                    method='sum',
                    output_path=str(output_path)
                )
            except Exception as e:
                logger.warning(f"Could not aggregate winter season {year}: {e}")
    
    logger.info(f"Completed temporal aggregation for {dataset_name}")
    return agg


# =============================================================================
# Step 3: Statistical Analysis
# =============================================================================

def run_statistical_analysis(dataset_name: str):
    """
    Run statistical analysis for a dataset.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate').
        
    Returns
    -------
    StatisticalAnalyzer
        The analyzer instance used for the dataset.
        
    Computes:
    - Anomalies for recent years
    - Long-term trends
    - Zonal statistics by region
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    logger.info(f"Running statistical analysis for {dataset_name}...")
    
    # Use pattern that excludes fraction files
    stats = StatisticalAnalyzer(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    # Dataset-specific output directories
    dataset_anomaly_dir = ANOMALY_DIR / dataset_name
    dataset_trend_dir = TREND_DIR / dataset_name
    dataset_anomaly_dir.mkdir(exist_ok=True)
    dataset_trend_dir.mkdir(exist_ok=True)
    
    # Calculate anomalies for recent years
    logger.info("Computing anomalies...")
    for year in range(2021, min(END_YEAR + 1, 2025)):
        for month in range(1, 13):
            output_path = dataset_anomaly_dir / f'anomaly_{year}_{month:02d}.tif'
            if not output_path.exists():
                try:
                    stats.calculate_anomaly(
                        year=year,
                        month=month,
                        clim_start=CLIM_START,
                        clim_end=CLIM_END,
                        anomaly_type='percent',
                        output_path=str(output_path)
                    )
                except Exception as e:
                    logger.debug(f"Could not calculate anomaly for {year}-{month:02d}: {e}")
    
    # Calculate trend
    logger.info("Computing trends...")
    try:
        stats.calculate_trend(
            start_year=START_YEAR,
            end_year=END_YEAR,
            method='sen',
            output_dir=str(dataset_trend_dir)
        )
    except Exception as e:
        logger.warning(f"Could not calculate trend: {e}")
    
    # Calculate zonal statistics
    logger.info("Computing zonal statistics...")
    if ZONE_FILE.exists():
        dataset_zonal_dir = ZONAL_DIR / dataset_name
        dataset_zonal_dir.mkdir(exist_ok=True)
        
        try:
            zonal_df = stats.zonal_statistics(
                geometry_path=str(ZONE_FILE),
                start_year=START_YEAR,
                end_year=END_YEAR,
                stats=['mean', 'sum', 'min', 'max', 'std', 'count'],
                output_path=str(dataset_zonal_dir / 'zonal_stats.csv')
            )
            logger.info(f"Zonal statistics saved: {dataset_zonal_dir / 'zonal_stats.csv'}")
            
            # Create zonal summary plot
            create_zonal_summary_plot(zonal_df, dataset_name, dataset_zonal_dir)
        except Exception as e:
            logger.warning(f"Could not calculate zonal statistics: {e}")
    else:
        logger.warning("Zone file not found. Run create_sample_zones() first.")
    
    logger.info(f"Completed statistical analysis for {dataset_name}")
    return stats


def create_zonal_summary_plot(zonal_df, dataset_name: str, output_dir: Path):
    """
    Create summary plots for zonal statistics.
    
    Generates a figure with two panels:
    - Time series of mean effective precipitation by zone
    - Monthly climatology bar chart by zone
    
    Parameters
    ----------
    zonal_df : pandas.DataFrame
        DataFrame containing zonal statistics with columns:
        zone_id, year, month, mean, etc.
    dataset_name : str
        Name of the dataset for plot titles.
    output_dir : Path
        Directory to save the output plot.
    """
    import matplotlib.pyplot as plt
    
    try:
        # Aggregate by zone and year
        annual_by_zone = zonal_df.groupby(['zone_id', 'year'])['mean'].mean().reset_index()
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        zone_names = {1: 'Central AZ', 2: 'Southern AZ', 3: 'Northern AZ'}
        colors = {1: '#e41a1c', 2: '#377eb8', 3: '#4daf4a'}
        
        # Plot 1: Time series by zone
        for zone_id in annual_by_zone['zone_id'].unique():
            zone_data = annual_by_zone[annual_by_zone['zone_id'] == zone_id]
            zone_name = zone_names.get(zone_id, f'Zone {zone_id}')
            axes[0].plot(zone_data['year'], zone_data['mean'], 
                        marker='o', label=zone_name, 
                        color=colors.get(zone_id, 'gray'), linewidth=2)
        
        axes[0].set_xlabel('Year')
        axes[0].set_ylabel('Mean Effective Precipitation [mm]')
        axes[0].set_title(f'{dataset_name} USDA-SCS - Zonal Mean Time Series')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Plot 2: Monthly climatology by zone
        monthly_by_zone = zonal_df.groupby(['zone_id', 'month'])['mean'].mean().reset_index()
        
        month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        
        width = 0.25
        x = range(1, 13)
        
        zone_ids = sorted(monthly_by_zone['zone_id'].unique())
        for i, zone_id in enumerate(zone_ids):
            zone_data = monthly_by_zone[monthly_by_zone['zone_id'] == zone_id]
            zone_name = zone_names.get(zone_id, f'Zone {zone_id}')
            offset = (i - 1) * width
            axes[1].bar([m + offset for m in zone_data['month']], 
                       zone_data['mean'], width=width, 
                       label=zone_name, color=colors.get(zone_id, 'gray'))
        
        axes[1].set_xlabel('Month')
        axes[1].set_ylabel('Mean Effective Precipitation [mm]')
        axes[1].set_title(f'{dataset_name} USDA-SCS - Monthly Zonal Climatology')
        axes[1].set_xticks(range(1, 13))
        axes[1].set_xticklabels(month_names)
        axes[1].legend()
        axes[1].grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        output_path = output_dir / 'zonal_summary.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved zonal summary plot: {output_path}")
        
    except Exception as e:
        logger.warning(f"Could not create zonal summary plot: {e}")


# =============================================================================
# Step 4: Visualization
# =============================================================================

def run_visualization(dataset_name: str):
    """
    Create visualizations for a dataset.
    
    Creates:
    - Time series plots
    - Monthly climatology bar charts  
    - Static raster maps for different seasons
    - Interactive HTML maps
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate').
        
    Returns
    -------
    Visualizer
        The visualizer instance used for the dataset.
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    logger.info(f"Creating visualizations for {dataset_name}...")
    
    # Use pattern that excludes fraction files
    viz = Visualizer(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    # Dataset-specific figure directory
    dataset_fig_dir = FIGURES_DIR / dataset_name
    dataset_fig_dir.mkdir(exist_ok=True)
    
    # Time series plot
    logger.info("Creating time series plot...")
    try:
        viz.plot_time_series(
            start_year=START_YEAR,
            end_year=END_YEAR,
            output_path=str(dataset_fig_dir / 'time_series.png'),
            title=f'{dataset_name} USDA-SCS Effective Precipitation - Arizona'
        )
    except Exception as e:
        logger.warning(f"Could not create time series plot: {e}")
    
    # Monthly climatology plot
    logger.info("Creating climatology plot...")
    try:
        viz.plot_monthly_climatology(
            start_year=CLIM_START,
            end_year=CLIM_END,
            output_path=str(dataset_fig_dir / 'monthly_climatology.png'),
            title=f'{dataset_name} Monthly Climatology ({CLIM_START}-{CLIM_END})'
        )
    except Exception as e:
        logger.warning(f"Could not create climatology plot: {e}")
    
    # Static maps for different seasons
    logger.info("Creating static maps...")
    months_to_plot = [
        (2023, 1, 'January (Winter)'),
        (2023, 4, 'April (Spring/Dry)'),
        (2023, 8, 'August (Monsoon Peak)'),
        (2023, 10, 'October (Fall)'),
    ]
    for year, month, desc in months_to_plot:
        try:
            viz.plot_raster(
                year=year,
                month=month,
                output_path=str(dataset_fig_dir / f'map_{year}_{month:02d}.png'),
                title=f'{dataset_name} USDA-SCS Peff - Arizona ({desc} {year})'
            )
        except Exception as e:
            logger.warning(f"Could not create static map for {year}-{month:02d}: {e}")
    
    # Interactive maps
    logger.info("Creating interactive maps...")
    for year, month, desc in months_to_plot:
        try:
            viz.plot_interactive_map(
                year=year,
                month=month,
                output_path=str(dataset_fig_dir / f'interactive_map_{year}_{month:02d}.html'),
                title=f'{dataset_name} USDA-SCS Peff - Arizona ({desc} {year})',
                zoom_start=6,
                opacity=0.8
            )
        except Exception as e:
            logger.warning(f"Could not create interactive map for {year}-{month:02d}: {e}")
    
    # Drought year comparison
    logger.info("Creating maps for drought/wet years...")
    notable_years = [
        (2002, 8, 'Drought Aug 2002'),
        (2021, 8, 'Wet Monsoon Aug 2021'),
    ]
    for year, month, desc in notable_years:
        try:
            viz.plot_raster(
                year=year,
                month=month,
                output_path=str(dataset_fig_dir / f'map_notable_{year}_{month:02d}.png'),
                title=f'{dataset_name} USDA-SCS - {desc}'
            )
        except Exception as e:
            logger.warning(f"Could not create map for {desc}: {e}")
    
    logger.info(f"Completed visualizations for {dataset_name}")
    return viz


def plot_anomaly_maps(dataset_name: str):
    """
    Create anomaly maps for recent years.
    
    Plots percent anomalies relative to climatology for years 2021-2023.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate').
    """
    """
    Create anomaly map visualizations for a dataset using Visualizer class.
    
    Plots percent anomalies for selected months using a diverging colormap.
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    dataset_anomaly_dir = ANOMALY_DIR / dataset_name
    dataset_fig_dir = FIGURES_DIR / dataset_name
    dataset_fig_dir.mkdir(exist_ok=True)
    
    if not dataset_anomaly_dir.exists():
        logger.warning(f"Anomaly directory not found: {dataset_anomaly_dir}")
        return
    
    logger.info(f"Creating anomaly maps for {dataset_name}...")
    
    viz = Visualizer(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    # Months to plot anomalies
    months_to_plot = [
        (2023, 8, 'August 2023 (Monsoon)'),
        (2022, 8, 'August 2022'),
        (2021, 8, 'August 2021'),
    ]
    
    for year, month, desc in months_to_plot:
        anomaly_file = dataset_anomaly_dir / f'anomaly_{year}_{month:02d}.tif'
        if not anomaly_file.exists():
            logger.debug(f"Anomaly file not found: {anomaly_file}")
            continue
        
        try:
            viz.plot_anomaly_map(
                anomaly_path=anomaly_file,
                title=f'{dataset_name} USDA-SCS Peff Anomaly - {desc}',
                output_path=dataset_fig_dir / f'anomaly_map_{year}_{month:02d}.png'
            )
        except Exception as e:
            logger.warning(f"Could not create anomaly map for {year}-{month:02d}: {e}")


def plot_climatology_maps(dataset_name: str):
    """
    Create monthly climatology maps.
    
    Plots long-term mean effective precipitation for each month.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate').
    """
    """
    Create monthly climatology map visualizations for a dataset using Visualizer class.
    
    Plots climatology maps for selected months (e.g., typical wet/dry season months).
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    dataset_clim_dir = CLIMATOLOGY_DIR / dataset_name
    dataset_fig_dir = FIGURES_DIR / dataset_name
    dataset_fig_dir.mkdir(exist_ok=True)
    
    if not dataset_clim_dir.exists():
        logger.warning(f"Climatology directory not found: {dataset_clim_dir}")
        return
    
    logger.info(f"Creating climatology maps for {dataset_name}...")
    
    viz = Visualizer(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    # Key months representing different seasons
    months_to_plot = [
        (1, 'January (Winter)'),
        (4, 'April (Spring/Dry)'),
        (8, 'August (Monsoon Peak)'),
        (10, 'October (Fall)'),
    ]
    
    for month, desc in months_to_plot:
        # Try different naming patterns for climatology files
        clim_file = dataset_clim_dir / f'climatology_{CLIM_START}_{CLIM_END}_month_{month:02d}.tif'
        if not clim_file.exists():
            # Fallback to simpler naming pattern
            clim_file = dataset_clim_dir / f'climatology_{month:02d}.tif'
        if not clim_file.exists():
            logger.debug(f"Climatology file not found for month {month:02d}")
            continue
        
        try:
            viz.plot_climatology_map(
                climatology_path=clim_file,
                title=f'{dataset_name} USDA-SCS Peff Climatology - {desc}',
                output_path=dataset_fig_dir / f'climatology_map_{month:02d}.png'
            )
        except Exception as e:
            logger.warning(f"Could not create climatology map for month {month:02d}: {e}")


def plot_trend_maps(dataset_name: str):
    """
    Create trend maps showing Sen's slope and significance.
    
    Plots the Sen's slope (mm/year) and Mann-Kendall p-values.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate').
    """
    """
    Create trend map visualizations for a dataset using Visualizer class.
    
    Plots:
    - Combined slope and p-value panel
    - Slope map with significance overlay
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    dataset_trend_dir = TREND_DIR / dataset_name
    dataset_fig_dir = FIGURES_DIR / dataset_name
    dataset_fig_dir.mkdir(exist_ok=True)
    
    if not dataset_trend_dir.exists():
        logger.warning(f"Trend directory not found: {dataset_trend_dir}")
        return
    
    logger.info(f"Creating trend maps for {dataset_name}...")
    
    viz = Visualizer(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    # Try different naming patterns for trend files
    slope_file = dataset_trend_dir / f'trend_slope_{START_YEAR}_{END_YEAR}_annual.tif'
    pvalue_file = dataset_trend_dir / f'trend_pvalue_{START_YEAR}_{END_YEAR}_annual.tif'
    
    # Fallback to simpler naming pattern
    if not slope_file.exists():
        slope_file = dataset_trend_dir / 'slope.tif'
        pvalue_file = dataset_trend_dir / 'pvalue.tif'
    
    if not slope_file.exists():
        logger.warning(f"Slope file not found in {dataset_trend_dir}")
        return
    
    # Create combined trend panel (slope + p-value)
    try:
        if pvalue_file.exists():
            viz.plot_trend_panel(
                slope_path=slope_file,
                pvalue_path=pvalue_file,
                title=f'{dataset_name} USDA-SCS Peff Trend ({START_YEAR}-{END_YEAR})',
                output_path=dataset_fig_dir / 'trend_maps.png'
            )
    except Exception as e:
        logger.warning(f"Could not create trend panel: {e}")
    
    # Create single trend map with significance stippling
    try:
        viz.plot_trend_map(
            slope_path=slope_file,
            pvalue_path=pvalue_file if pvalue_file.exists() else None,
            title=f'{dataset_name} USDA-SCS Peff Trend with Significance ({START_YEAR}-{END_YEAR})',
            output_path=dataset_fig_dir / 'trend_map_with_significance.png',
            show_significance=True
        )
    except Exception as e:
        logger.warning(f"Could not create trend map: {e}")


# =============================================================================
# Step 5: Dataset Comparison (GridMET vs PRISM)
# =============================================================================

def compare_datasets():
    """
    Compare GridMET and PRISM effective precipitation.
    
    Creates:
    - Side-by-side comparison plots
    - Scatter comparison plots
    - Annual comparison bar charts
    """
    logger.info("Comparing GridMET and PRISM datasets...")
    
    # Initialize visualizers
    gridmet_dir = Path(DATASETS['GridMET']['output_dir'])
    prism_dir = Path(DATASETS['PRISM']['output_dir'])
    
    if not gridmet_dir.exists() or not prism_dir.exists():
        logger.warning("Both datasets must be processed for comparison")
        return
    
    # Pattern to exclude fraction files
    pattern = 'effective_precip_[0-9]*.tif'
    
    viz_gridmet = Visualizer(str(gridmet_dir), pattern=pattern)
    
    # Side-by-side comparison for monsoon month
    logger.info("Creating side-by-side comparison...")
    try:
        viz_gridmet.plot_comparison(
            year=2023,
            month=8,  # August - peak monsoon
            other_dir=str(prism_dir),
            other_pattern=pattern,
            labels=('GridMET', 'PRISM'),
            title='USDA-SCS Effective Precipitation Comparison - Arizona (August 2023)',
            output_path=str(COMPARISON_DIR / 'comparison_2023_08.png')
        )
    except Exception as e:
        logger.warning(f"Could not create comparison plot: {e}")
    
    # Winter month comparison
    try:
        viz_gridmet.plot_comparison(
            year=2023,
            month=1,  # January
            other_dir=str(prism_dir),
            other_pattern=pattern,
            labels=('GridMET', 'PRISM'),
            title='USDA-SCS Effective Precipitation Comparison - Arizona (January 2023)',
            output_path=str(COMPARISON_DIR / 'comparison_2023_01.png')
        )
    except Exception as e:
        logger.warning(f"Could not create winter comparison plot: {e}")
    
    # Scatter comparison
    logger.info("Creating scatter comparison...")
    try:
        viz_gridmet.plot_scatter_comparison(
            start_year=START_YEAR,
            end_year=END_YEAR,
            other_dir=str(prism_dir),
            other_pattern=pattern,
            labels=('GridMET', 'PRISM'),
            output_path=str(COMPARISON_DIR / 'scatter_comparison.png')
        )
    except Exception as e:
        logger.warning(f"Could not create scatter plot: {e}")
    
    # Annual comparison
    logger.info("Creating annual comparison...")
    try:
        viz_gridmet.plot_annual_comparison(
            start_year=START_YEAR,
            end_year=END_YEAR,
            other_dir=str(prism_dir),
            other_pattern=pattern,
            labels=('GridMET', 'PRISM'),
            output_path=str(COMPARISON_DIR / 'annual_comparison.png')
        )
    except Exception as e:
        logger.warning(f"Could not create annual comparison: {e}")
    
    # Zonal comparison between datasets
    logger.info("Creating zonal statistics comparison...")
    compare_zonal_statistics()
    
    logger.info("Completed dataset comparison")


def compare_zonal_statistics():
    """
    Compare zonal statistics between GridMET and PRISM.
    
    Creates comparison plots showing time series and monthly climatology
    for each Arizona zone (Central, Southern, Northern).
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    
    gridmet_zonal = ZONAL_DIR / 'GridMET' / 'zonal_stats.csv'
    prism_zonal = ZONAL_DIR / 'PRISM' / 'zonal_stats.csv'
    
    if not gridmet_zonal.exists() or not prism_zonal.exists():
        logger.warning("Zonal statistics files not found for comparison")
        return
    
    try:
        df_gridmet = pd.read_csv(gridmet_zonal)
        df_prism = pd.read_csv(prism_zonal)
        
        # Add dataset identifier
        df_gridmet['dataset'] = 'GridMET'
        df_prism['dataset'] = 'PRISM'
        
        # Combine datasets
        df_combined = pd.concat([df_gridmet, df_prism], ignore_index=True)
        
        # Create comparison figure
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        zone_names = {1: 'Central AZ', 2: 'Southern AZ', 3: 'Northern AZ'}
        colors = {'GridMET': '#1f77b4', 'PRISM': '#ff7f0e'}
        
        # Plot 1: Annual time series comparison - Central AZ
        ax = axes[0, 0]
        for dataset in ['GridMET', 'PRISM']:
            data = df_combined[(df_combined['dataset'] == dataset) & (df_combined['zone_id'] == 1)]
            annual = data.groupby('year')['mean'].mean()
            ax.plot(annual.index, annual.values, marker='o', 
                   label=dataset, color=colors[dataset], linewidth=2)
        ax.set_xlabel('Year')
        ax.set_ylabel('Mean Effective Precipitation [mm]')
        ax.set_title('Central AZ - Annual Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Monthly climatology - All zones comparison
        ax = axes[0, 1]
        for dataset in ['GridMET', 'PRISM']:
            data = df_combined[df_combined['dataset'] == dataset]
            monthly = data.groupby('month')['mean'].mean()
            ax.plot(monthly.index, monthly.values, marker='o', 
                   label=dataset, color=colors[dataset], linewidth=2)
        ax.set_xlabel('Month')
        ax.set_ylabel('Mean Effective Precipitation [mm]')
        ax.set_title('All Arizona - Monthly Climatology Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
        
        # Plot 3: Monsoon (Jul-Sep) comparison by zone
        ax = axes[1, 0]
        monsoon_data = df_combined[df_combined['month'].isin([7, 8, 9])]
        monsoon_by_zone = monsoon_data.groupby(['dataset', 'zone_id'])['mean'].mean().unstack()
        
        x = range(len(zone_names))
        width = 0.35
        
        for i, dataset in enumerate(['GridMET', 'PRISM']):
            if dataset in monsoon_by_zone.index:
                offset = (i - 0.5) * width
                values = [monsoon_by_zone.loc[dataset, z] for z in sorted(zone_names.keys())]
                ax.bar([xi + offset for xi in x], values, width=width, 
                      label=dataset, color=colors[dataset])
        
        ax.set_xlabel('Zone')
        ax.set_ylabel('Mean Monsoon Peff [mm]')
        ax.set_title('Monsoon Season (Jul-Sep) by Zone')
        ax.set_xticks(x)
        ax.set_xticklabels([zone_names[z] for z in sorted(zone_names.keys())])
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Plot 4: Interannual variability (coefficient of variation)
        ax = axes[1, 1]
        cv_data = df_combined.groupby(['dataset', 'zone_id']).apply(
            lambda g: g['mean'].std() / g['mean'].mean() * 100
        ).unstack()
        
        for i, dataset in enumerate(['GridMET', 'PRISM']):
            if dataset in cv_data.index:
                offset = (i - 0.5) * width
                values = [cv_data.loc[dataset, z] for z in sorted(zone_names.keys())]
                ax.bar([xi + offset for xi in x], values, width=width, 
                      label=dataset, color=colors[dataset])
        
        ax.set_xlabel('Zone')
        ax.set_ylabel('Coefficient of Variation [%]')
        ax.set_title('Interannual Variability by Zone')
        ax.set_xticks(x)
        ax.set_xticklabels([zone_names[z] for z in sorted(zone_names.keys())])
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.suptitle('GridMET vs PRISM USDA-SCS Comparison - Arizona', fontsize=14, y=1.02)
        plt.tight_layout()
        
        output_path = COMPARISON_DIR / 'zonal_comparison.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved zonal comparison: {output_path}")
        
    except Exception as e:
        logger.warning(f"Could not create zonal comparison: {e}")


# =============================================================================
# Step 5b: U.S. vs Global Dataset Comparison
# =============================================================================

def compare_us_vs_global():
    """
    Compare U.S. and Global datasets for Arizona.
    
    Creates comparison plots between U.S. datasets (GridMET, PRISM) and
    global datasets (ERA5-Land, TerraClimate) including:
    - Side-by-side spatial maps
    - Scatter plot comparisons
    - Time series comparisons
    """
    """
    Compare U.S.-specific datasets (GridMET, PRISM) with global datasets (ERA5-Land, TerraClimate).
    
    Creates:
    - Side-by-side comparison of U.S. vs Global datasets
    - Scatter plots between datasets
    - Summary statistics comparison
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    logger.info("=" * 40)
    logger.info("Comparing U.S. vs Global Datasets")
    logger.info("=" * 40)
    
    # Check which datasets are available
    available_us = [name for name, config in US_DATASETS.items() 
                    if Path(config['output_dir']).exists()]
    available_global = [name for name, config in GLOBAL_DATASETS.items() 
                        if Path(config['output_dir']).exists()]
    
    if not available_us or not available_global:
        logger.warning("Need at least one U.S. and one Global dataset for comparison")
        return
    
    logger.info(f"Available U.S. datasets: {available_us}")
    logger.info(f"Available Global datasets: {available_global}")
    
    pattern = 'effective_precip_[0-9]*.tif'
    
    # Compare GridMET (U.S.) vs ERA5-Land (Global) if both available
    if 'GridMET' in available_us and 'ERA5Land' in available_global:
        logger.info("Comparing GridMET (U.S.) vs ERA5-Land (Global)...")
        gridmet_dir = Path(US_DATASETS['GridMET']['output_dir'])
        era5_dir = Path(GLOBAL_DATASETS['ERA5Land']['output_dir'])
        
        viz = Visualizer(str(gridmet_dir), pattern=pattern)
        
        try:
            viz.plot_comparison(
                year=2023, month=8,
                other_dir=str(era5_dir),
                other_pattern=pattern,
                labels=('GridMET (U.S.)', 'ERA5-Land (Global)'),
                title='U.S. vs Global: GridMET vs ERA5-Land - Arizona (August 2023)',
                output_path=str(US_VS_GLOBAL_DIR / 'gridmet_vs_era5land_2023_08.png')
            )
        except Exception as e:
            logger.warning(f"Could not create GridMET vs ERA5-Land comparison: {e}")
        
        try:
            viz.plot_scatter_comparison(
                start_year=START_YEAR, end_year=END_YEAR,
                other_dir=str(era5_dir),
                other_pattern=pattern,
                labels=('GridMET (U.S.)', 'ERA5-Land (Global)'),
                output_path=str(US_VS_GLOBAL_DIR / 'scatter_gridmet_vs_era5land.png')
            )
        except Exception as e:
            logger.warning(f"Could not create scatter plot: {e}")
    
    # Compare PRISM (U.S.) vs TerraClimate (Global) if both available
    if 'PRISM' in available_us and 'TerraClimate' in available_global:
        logger.info("Comparing PRISM (U.S.) vs TerraClimate (Global)...")
        prism_dir = Path(US_DATASETS['PRISM']['output_dir'])
        terra_dir = Path(GLOBAL_DATASETS['TerraClimate']['output_dir'])
        
        viz = Visualizer(str(prism_dir), pattern=pattern)
        
        try:
            viz.plot_comparison(
                year=2023, month=8,
                other_dir=str(terra_dir),
                other_pattern=pattern,
                labels=('PRISM (U.S.)', 'TerraClimate (Global)'),
                title='U.S. vs Global: PRISM vs TerraClimate - Arizona (August 2023)',
                output_path=str(US_VS_GLOBAL_DIR / 'prism_vs_terraclimate_2023_08.png')
            )
        except Exception as e:
            logger.warning(f"Could not create PRISM vs TerraClimate comparison: {e}")
    
    # Create multi-dataset summary comparison
    create_multi_dataset_comparison(available_us + available_global)
    
    logger.info("Completed U.S. vs Global comparison")


def create_multi_dataset_comparison(dataset_names: list):
    """
    Create comprehensive multi-dataset comparison.
    
    Generates a multi-panel figure comparing all available datasets
    for a monsoon month (August), including spatial maps and statistics.
    
    Parameters
    ----------
    dataset_names : list of str
        List of dataset names to compare.
    """
    """
    Create summary comparison across all available datasets.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import rioxarray
    
    logger.info("Creating multi-dataset summary comparison...")
    
    results = []
    sample_years = [2010, 2015, 2020]
    sample_months = [1, 4, 7, 10]  # Quarterly
    
    for dataset_name in dataset_names:
        config = DATASETS[dataset_name]
        input_dir = Path(config['output_dir'])
        region = config.get('region', 'Unknown')
        
        for year in sample_years:
            for month in sample_months:
                peff_file = input_dir / f'effective_precip_{year}_{month:02d}.tif'
                if not peff_file.exists():
                    continue
                
                try:
                    da = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
                    da = da.where(da != 0)
                    vals = da.values[~np.isnan(da.values)]
                    
                    results.append({
                        'dataset': dataset_name,
                        'region': region,
                        'year': year,
                        'month': month,
                        'mean_peff': np.nanmean(vals),
                        'std_peff': np.nanstd(vals),
                        'median_peff': np.nanmedian(vals)
                    })
                except Exception as e:
                    logger.debug(f"Error reading {peff_file}: {e}")
    
    if not results:
        logger.warning("No data available for multi-dataset comparison")
        return
    
    df = pd.DataFrame(results)
    
    # Create comparison figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Color scheme: U.S. datasets in blues, Global in oranges/reds
    colors = {
        'GridMET': '#1f77b4',
        'PRISM': '#17becf',
        'ERA5Land': '#ff7f0e',
        'TerraClimate': '#d62728'
    }
    
    # Plot 1: Box plot by dataset
    ax = axes[0, 0]
    datasets_available = df['dataset'].unique()
    box_data = [df[df['dataset'] == d]['mean_peff'].values for d in datasets_available]
    bp = ax.boxplot(box_data, labels=datasets_available, patch_artist=True)
    for patch, dataset in zip(bp['boxes'], datasets_available):
        patch.set_facecolor(colors.get(dataset, 'gray'))
        patch.set_alpha(0.7)
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Distribution by Dataset')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Monthly climatology by dataset
    ax = axes[0, 1]
    for dataset in datasets_available:
        data = df[df['dataset'] == dataset]
        monthly = data.groupby('month')['mean_peff'].mean()
        ax.plot(monthly.index, monthly.values, marker='o',
               label=dataset, color=colors.get(dataset, 'gray'), linewidth=2)
    ax.set_xlabel('Month')
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Monthly Climatology by Dataset')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xticks([1, 4, 7, 10])
    ax.set_xticklabels(['Jan', 'Apr', 'Jul', 'Oct'])
    
    # Plot 3: U.S. vs Global bar comparison
    ax = axes[1, 0]
    us_mean = df[df['region'] == 'US']['mean_peff'].mean()
    global_mean = df[df['region'] == 'Global']['mean_peff'].mean()
    us_std = df[df['region'] == 'US']['mean_peff'].std()
    global_std = df[df['region'] == 'Global']['mean_peff'].std()
    
    x = [0, 1]
    ax.bar(x, [us_mean, global_mean], yerr=[us_std, global_std],
           color=['#1f77b4', '#ff7f0e'], alpha=0.8, capsize=5)
    ax.set_xticks(x)
    ax.set_xticklabels(['U.S. Datasets\n(GridMET, PRISM)', 'Global Datasets\n(ERA5-Land, TerraClimate)'])
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('U.S. vs Global Dataset Comparison')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Annual time series by region
    ax = axes[1, 1]
    for region, color in [('US', '#1f77b4'), ('Global', '#ff7f0e')]:
        data = df[df['region'] == region]
        annual = data.groupby('year')['mean_peff'].mean()
        ax.plot(annual.index, annual.values, marker='o',
               label=f'{region} Datasets', color=color, linewidth=2)
    ax.set_xlabel('Year')
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Annual Mean by Dataset Region')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('U.S. vs Global Precipitation Dataset Comparison - Arizona', fontsize=14, y=1.02)
    plt.tight_layout()
    
    output_path = US_VS_GLOBAL_DIR / 'multi_dataset_summary.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved multi-dataset summary: {output_path}")
    
    # Save statistics to CSV
    summary_stats = df.groupby(['dataset', 'region']).agg({
        'mean_peff': ['mean', 'std', 'min', 'max'],
        'median_peff': 'mean'
    }).round(2)
    summary_stats.to_csv(US_VS_GLOBAL_DIR / 'dataset_statistics.csv')
    logger.info(f"Saved dataset statistics: {US_VS_GLOBAL_DIR / 'dataset_statistics.csv'}")


# =============================================================================
# Step 5c: Compare Effective Precipitation Methods
# =============================================================================

def compare_peff_methods():
    """
    Compare different effective precipitation methods for Arizona.
    
    Compares CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall,
    FarmWest, USDA-SCS, and TAGEM-SuET methods using the same input precipitation data.
    
    Note: USDA-SCS method requires AWC and ETo. For demonstration,
    representative values are used: AWC=150 mm/m, ETo=150 mm/month (arid climate).
    TAGEM-SuET requires ETo data.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    logger.info("=" * 40)
    logger.info("Comparing Effective Precipitation Methods")
    logger.info("=" * 40)
    
    # Method display names and colors
    method_info = {
        'cropwat': {'name': 'CROPWAT', 'color': '#1f77b4'},
        'fao_aglw': {'name': 'FAO/AGLW', 'color': '#ff7f0e'},
        'fixed_percentage': {'name': 'Fixed 70%', 'color': '#2ca02c'},
        'dependable_rainfall': {'name': 'Dependable Rain', 'color': '#d62728'},
        'farmwest': {'name': 'FarmWest', 'color': '#9467bd'},
        'usda_scs': {'name': 'USDA-SCS', 'color': '#8c564b'},
        'suet': {'name': 'TAGEM-SuET', 'color': '#e377c2'},
        'ensemble': {'name': 'Ensemble', 'color': '#17becf'}
    }
    
    # Check available datasets
    available_datasets = [name for name, config in DATASETS.items() 
                         if Path(config['output_dir']).exists()]
    
    if not available_datasets:
        logger.warning("No processed datasets available for method comparison")
        return
    
    # Create method comparison for each available dataset
    for dataset_name in available_datasets:
        logger.info(f"\nComparing methods for {dataset_name}...")
        
        try:
            compare_methods_for_dataset(dataset_name, method_info)
        except Exception as e:
            logger.warning(f"Could not compare methods for {dataset_name}: {e}")
    
    # Create theoretical curve comparison
    try:
        plot_method_curves(method_info)
    except Exception as e:
        logger.warning(f"Could not create method curves: {e}")
    
    # Create cross-dataset method comparison
    try:
        compare_methods_across_datasets(available_datasets, method_info)
    except Exception as e:
        logger.warning(f"Could not create cross-dataset comparison: {e}")
    
    logger.info("Completed method comparison")


def compare_methods_for_dataset(dataset_name: str, method_info: dict):
    """
    Compare effective precipitation methods for a specific dataset.
    
    Creates a 2x3 plot showing the spatial distribution of effective
    precipitation calculated by each method for a sample monsoon month.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset to use as base precipitation.
    method_info : dict
        Dictionary with method names as keys and sub-dicts containing
        'name' (display name) and 'color' for plotting.
    """
    import matplotlib.pyplot as plt
    import rioxarray
    import numpy as np
    from pycropwat.methods import (
        cropwat_effective_precip,
        fao_aglw_effective_precip,
        fixed_percentage_effective_precip,
        dependable_rainfall_effective_precip,
        farmwest_effective_precip,
        usda_scs_effective_precip,
        suet_effective_precip,
        ensemble_effective_precip
    )
    
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    # Input directory where AWC and ETo were saved
    analysis_inputs_dir = REGION_DIR / 'analysis_inputs' / input_dir.name
    
    # Use a sample month
    sample_year = 2020
    sample_month = 8  # August - monsoon season
    
    # Find the effective precipitation and fraction files
    peff_file = input_dir / f'effective_precip_{sample_year}_{sample_month:02d}.tif'
    fraction_file = input_dir / f'effective_precip_fraction_{sample_year}_{sample_month:02d}.tif'
    
    if not peff_file.exists():
        logger.warning(f"File not found: {peff_file}")
        return
    
    # Load the default Peff
    da_default = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
    
    # Calculate total precipitation from fraction
    if fraction_file.exists():
        da_fraction = rioxarray.open_rasterio(fraction_file).squeeze('band', drop=True)
        fraction_values = da_fraction.values
        fraction_values[fraction_values == 0] = 1  # Avoid division by zero
        precip_total = da_default.values / fraction_values
    else:
        # Approximate inverse for USDA-SCS (rough estimate)
        precip_total = da_default.values * 1.5  # Rough approximation
    
    # Load actual AWC data (SSURGO for US, FAO for global)
    awc_file = analysis_inputs_dir / 'awc.tif'
    if awc_file.exists():
        da_awc = rioxarray.open_rasterio(awc_file).squeeze('band', drop=True)
        awc_data = da_awc.values
        # Resample if shapes don't match
        if awc_data.shape != precip_total.shape:
            from scipy.ndimage import zoom
            zoom_factors = (precip_total.shape[0] / awc_data.shape[0],
                           precip_total.shape[1] / awc_data.shape[1])
            awc_data = zoom(awc_data, zoom_factors, order=1)
        logger.info(f"Loaded actual AWC data from {awc_file}")
    else:
        logger.warning(f"AWC file not found: {awc_file}, using representative value")
        awc_data = np.full_like(precip_total, 150.0)  # SSURGO representative
    
    # Load actual ETo data (gridMET for US, AgERA5 for global)
    eto_file = analysis_inputs_dir / f'eto_{sample_year}_{sample_month:02d}.tif'
    if eto_file.exists():
        da_eto = rioxarray.open_rasterio(eto_file).squeeze('band', drop=True)
        eto_data = da_eto.values
        # Resample if shapes don't match
        if eto_data.shape != precip_total.shape:
            from scipy.ndimage import zoom
            zoom_factors = (precip_total.shape[0] / eto_data.shape[0],
                           precip_total.shape[1] / eto_data.shape[1])
            eto_data = zoom(eto_data, zoom_factors, order=1)
        logger.info(f"Loaded actual ETo data from {eto_file}")
    else:
        logger.warning(f"ETo file not found: {eto_file}, using representative value")
        eto_data = np.full_like(precip_total, 180.0)  # gridMET representative for Arizona
    
    # Calculate Peff using each method with actual data
    peff_methods = {
        'cropwat': cropwat_effective_precip(precip_total),
        'fao_aglw': fao_aglw_effective_precip(precip_total),
        'fixed_percentage': fixed_percentage_effective_precip(precip_total, 0.7),
        'dependable_rainfall': dependable_rainfall_effective_precip(precip_total, 0.75),
        'farmwest': farmwest_effective_precip(precip_total),
        'usda_scs': usda_scs_effective_precip(precip_total, eto_data, awc_data),
        'suet': suet_effective_precip(precip_total, eto_data),
        'ensemble': ensemble_effective_precip(precip_total, eto_data, awc_data, ROOTING_DEPTH, 0.7, 0.75)
    }
    
    # Create 3x3 comparison plot (9 methods)
    fig, axes = plt.subplots(3, 3, figsize=(24, 18))
    axes = axes.flatten()
    
    # Common color scale
    vmin = min(np.nanmin(v) for v in peff_methods.values())
    vmax = max(np.nanmax(v) for v in peff_methods.values())
    
    for idx, (method, peff_vals) in enumerate(peff_methods.items()):
        ax = axes[idx]
        
        # Create xarray for plotting
        da_method = da_default.copy(data=peff_vals)
        da_method = da_method.where(~np.isnan(da_default.values))
        
        im = da_method.plot(
            ax=ax,
            cmap='YlGnBu',
            vmin=vmin,
            vmax=vmax,
            add_colorbar=False
        )
        
        ax.set_title(f"{method_info[method]['name']}", fontsize=12, fontweight='bold')
        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
    
    # Hide unused axes (if any)
    for idx in range(len(peff_methods), len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle(
        f'{dataset_name} - Method Comparison (August {sample_year})',
        fontsize=14, fontweight='bold'
    )
    
    # Add shared colorbar
    plt.tight_layout()
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label='Effective Precipitation [mm]')
    
    output_path = METHOD_COMPARISON_DIR / f'{dataset_name}_method_maps_{sample_year}_{sample_month:02d}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method comparison maps: {output_path}")


def plot_method_curves(method_info: dict):
    """
    Plot theoretical effective precipitation curves for all methods.
    
    Creates a figure showing how each method transforms total precipitation
    to effective precipitation across a range of 0-500 mm.
    
    Parameters
    ----------
    method_info : dict
        Dictionary with method names as keys and sub-dicts containing
        'name' (display name) and 'color' for plotting.
    """
    """
    Plot theoretical Peff curves for each method.
    
    Shows how each method transforms precipitation to effective precipitation,
    with Arizona-specific parameters for USDA-SCS.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from pycropwat.methods import (
        cropwat_effective_precip,
        fao_aglw_effective_precip,
        fixed_percentage_effective_precip,
        dependable_rainfall_effective_precip,
        farmwest_effective_precip,
        usda_scs_effective_precip,
        suet_effective_precip,
        ensemble_effective_precip
    )
    
    # Create precipitation range (Arizona typically 0-150mm/month)
    precip = np.linspace(0, 200, 1000)
    
    # For theoretical curves, use regional mean AWC and ETo values
    # These are illustrative - actual spatial data is used in compare_methods_for_dataset()
    # AWC: SSURGO regional mean (~150 mm/m for Arizona desert soils)
    # ETo: gridMET regional mean (~180 mm/month high evaporative demand)
    awc_ssurgo = np.full_like(precip, 150.0)  # mm/m (SSURGO regional mean)
    eto_gridmet = np.full_like(precip, 180.0)  # mm/month (gridMET regional mean)
    
    # Calculate Peff for each method
    methods = {
        'cropwat': cropwat_effective_precip(precip),
        'fao_aglw': fao_aglw_effective_precip(precip),
        'fixed_percentage': fixed_percentage_effective_precip(precip, 0.7),
        'dependable_rainfall': dependable_rainfall_effective_precip(precip, 0.75),
        'farmwest': farmwest_effective_precip(precip),
        'usda_scs': usda_scs_effective_precip(precip, eto_gridmet, awc_ssurgo),
        'suet': suet_effective_precip(precip, eto_gridmet),
        'ensemble': ensemble_effective_precip(precip, eto_gridmet, awc_ssurgo, ROOTING_DEPTH, 0.7, 0.75)
    }
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: All methods on same plot
    ax = axes[0]
    for method, peff in methods.items():
        ax.plot(precip, peff, 
               label=method_info[method]['name'],
               color=method_info[method]['color'],
               linewidth=2)
    
    ax.plot(precip, precip, 'k--', alpha=0.3, label='1:1 line')
    ax.set_xlabel('Total Precipitation [mm]')
    ax.set_ylabel('Effective Precipitation [mm]')
    ax.set_title('Peff Methods Comparison (Arizona Parameters)')
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 150)
    
    # Shade typical Arizona monsoon range
    ax.axvspan(20, 80, alpha=0.1, color='blue', label='Typical Monsoon Range')
    
    # Plot 2: Efficiency ratio (Peff/P) vs precipitation
    ax = axes[1]
    for method, peff in methods.items():
        efficiency = np.divide(peff, precip, where=precip > 0, out=np.zeros_like(peff))
        ax.plot(precip[precip > 0], efficiency[precip > 0],
               label=method_info[method]['name'],
               color=method_info[method]['color'],
               linewidth=2)
    
    ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.3, label='100%')
    ax.set_xlabel('Total Precipitation [mm]')
    ax.set_ylabel('Efficiency (Peff/P)')
    ax.set_title('Rainfall Effectiveness by Method')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 1.1)
    
    plt.suptitle('Theoretical Method Curves - Arizona Conditions', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = METHOD_COMPARISON_DIR / 'method_curves_arizona.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method curves: {output_path}")


def compare_methods_across_datasets(dataset_names: list, method_info: dict):
    """
    Compare effective precipitation methods across multiple datasets.
    
    Creates summary statistics and comparison plots showing how different
    methods perform across U.S. and Global datasets.
    
    Parameters
    ----------
    dataset_names : list of str
        List of dataset names to include in comparison.
    method_info : dict
        Dictionary with method names as keys and sub-dicts containing
        'name' (display name) and 'color' for plotting.
    """
    """
    Create summary comparison of methods across all datasets.
    
    Shows how method choice affects results for each precipitation dataset.
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import rioxarray
    from pycropwat.methods import (
        cropwat_effective_precip,
        fao_aglw_effective_precip,
        fixed_percentage_effective_precip,
        dependable_rainfall_effective_precip,
        farmwest_effective_precip,
        usda_scs_effective_precip,
        suet_effective_precip,
        ensemble_effective_precip
    )
    
    results = []
    
    for dataset_name in dataset_names:
        config = DATASETS[dataset_name]
        input_dir = Path(config['output_dir'])
        region = config.get('region', 'Unknown')
        
        # Input directory where AWC and ETo were saved
        analysis_inputs_dir = REGION_DIR / 'analysis_inputs' / input_dir.name
        
        # Load AWC data once per dataset (static)
        awc_file = analysis_inputs_dir / 'awc.tif'
        awc_loaded = None
        if awc_file.exists():
            da_awc = rioxarray.open_rasterio(awc_file).squeeze('band', drop=True)
            awc_loaded = da_awc.values
            logger.debug(f"Loaded AWC data from {awc_file}")
        
        # Sample months
        sample_years = [2010, 2015, 2020]
        sample_months = [1, 4, 7, 10]
        
        for year in sample_years:
            for month in sample_months:
                peff_file = input_dir / f'effective_precip_{year}_{month:02d}.tif'
                fraction_file = input_dir / f'effective_precip_fraction_{year}_{month:02d}.tif'
                
                if not peff_file.exists():
                    continue
                
                try:
                    da = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
                    da = da.where(da != 0)
                    peff_vals = da.values[~np.isnan(da.values)]
                    
                    # Get total precipitation
                    if fraction_file.exists():
                        da_frac = rioxarray.open_rasterio(fraction_file).squeeze('band', drop=True)
                        da_frac = da_frac.where(da_frac > 0)
                        precip_approx = peff_vals / da_frac.values[~np.isnan(da.values)]
                    else:
                        precip_approx = peff_vals * 1.5
                    
                    precip_approx = np.clip(precip_approx, 0, 500)
                    
                    # Load actual AWC data or use representative value
                    if awc_loaded is not None:
                        awc_flat = awc_loaded[~np.isnan(da.values)]
                        if len(awc_flat) != len(precip_approx):
                            awc_data = np.full_like(precip_approx, np.nanmean(awc_loaded))
                        else:
                            awc_data = awc_flat
                    else:
                        awc_data = np.full_like(precip_approx, 150.0)  # SSURGO representative
                    
                    # Load actual ETo data or use representative value
                    eto_file = analysis_inputs_dir / f'eto_{year}_{month:02d}.tif'
                    if eto_file.exists():
                        da_eto = rioxarray.open_rasterio(eto_file).squeeze('band', drop=True)
                        eto_flat = da_eto.values[~np.isnan(da.values)]
                        if len(eto_flat) != len(precip_approx):
                            eto_data = np.full_like(precip_approx, np.nanmean(da_eto.values))
                        else:
                            eto_data = eto_flat
                    else:
                        eto_data = np.full_like(precip_approx, 180.0)  # gridMET representative
                    
                    # Calculate mean Peff for each method using actual data
                    for method_name, func in [
                        ('cropwat', lambda p: cropwat_effective_precip(p)),
                        ('fao_aglw', lambda p: fao_aglw_effective_precip(p)),
                        ('fixed_percentage', lambda p: fixed_percentage_effective_precip(p, 0.7)),
                        ('dependable_rainfall', lambda p: dependable_rainfall_effective_precip(p, 0.75)),
                        ('farmwest', lambda p: farmwest_effective_precip(p)),
                        ('usda_scs', lambda p: usda_scs_effective_precip(p, eto_data, awc_data)),
                        ('suet', lambda p: suet_effective_precip(p, eto_data)),
                        ('ensemble', lambda p: ensemble_effective_precip(p, eto_data, awc_data, ROOTING_DEPTH, 0.7, 0.75))
                    ]:
                        peff_method = func(precip_approx)
                        results.append({
                            'dataset': dataset_name,
                            'region': region,
                            'year': year,
                            'month': month,
                            'method': method_name,
                            'mean_peff': np.nanmean(peff_method),
                            'mean_precip': np.nanmean(precip_approx)
                        })
                except Exception as e:
                    logger.debug(f"Error processing {peff_file}: {e}")
    
    if not results:
        logger.warning("No data available for cross-dataset method comparison")
        return
    
    df = pd.DataFrame(results)
    
    # Create comparison figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Box plot of Peff by method (all datasets)
    ax = axes[0, 0]
    methods_list = list(method_info.keys())
    box_data = [df[df['method'] == m]['mean_peff'].values for m in methods_list]
    bp = ax.boxplot(box_data, labels=[method_info[m]['name'] for m in methods_list],
                   patch_artist=True)
    for patch, method in zip(bp['boxes'], methods_list):
        patch.set_facecolor(method_info[method]['color'])
        patch.set_alpha(0.7)
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Method Distribution (All Datasets)')
    ax.grid(True, alpha=0.3, axis='y')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=15)
    
    # Plot 2: Method comparison by dataset
    ax = axes[0, 1]
    method_means = df.groupby(['dataset', 'method'])['mean_peff'].mean().unstack()
    x = np.arange(len(methods_list))
    width = 0.15
    
    for i, dataset in enumerate(dataset_names):
        if dataset in method_means.index:
            offset = (i - len(dataset_names)/2 + 0.5) * width
            values = [method_means.loc[dataset, m] if m in method_means.columns else 0 
                     for m in methods_list]
            ax.bar(x + offset, values, width, label=dataset, alpha=0.8)
    
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Method × Dataset Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels([method_info[m]['name'] for m in methods_list], rotation=15)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: U.S. vs Global by method
    ax = axes[1, 0]
    region_method = df.groupby(['region', 'method'])['mean_peff'].mean().unstack()
    
    width = 0.35
    for i, region in enumerate(['US', 'Global']):
        if region in region_method.index:
            offset = (i - 0.5) * width
            values = [region_method.loc[region, m] if m in region_method.columns else 0 
                     for m in methods_list]
            color = '#1f77b4' if region == 'US' else '#ff7f0e'
            ax.bar(x + offset, values, width, label=f'{region} Datasets', color=color, alpha=0.8)
    
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('U.S. vs Global Datasets by Method')
    ax.set_xticks(x)
    ax.set_xticklabels([method_info[m]['name'] for m in methods_list], rotation=15)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Efficiency by method and dataset
    ax = axes[1, 1]
    df['efficiency'] = df['mean_peff'] / df['mean_precip']
    eff_by_method = df.groupby(['dataset', 'method'])['efficiency'].mean().unstack()
    
    for i, dataset in enumerate(dataset_names):
        if dataset in eff_by_method.index:
            offset = (i - len(dataset_names)/2 + 0.5) * width
            values = [eff_by_method.loc[dataset, m] if m in eff_by_method.columns else 0 
                     for m in methods_list]
            ax.bar(x + offset, values, width, label=dataset, alpha=0.8)
    
    ax.set_ylabel('Mean Efficiency (Peff/P)')
    ax.set_title('Rainfall Effectiveness by Method & Dataset')
    ax.set_xticks(x)
    ax.set_xticklabels([method_info[m]['name'] for m in methods_list], rotation=15)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1)
    
    plt.suptitle('Effective Precipitation Method Comparison - Arizona', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = METHOD_COMPARISON_DIR / 'method_comparison_summary.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method comparison summary: {output_path}")
    
    # Save statistics to CSV
    summary_stats = df.groupby(['dataset', 'region', 'method']).agg({
        'mean_peff': ['mean', 'std', 'min', 'max'],
        'mean_precip': 'mean'
    }).round(2)
    summary_stats.to_csv(METHOD_COMPARISON_DIR / 'method_statistics.csv')
    logger.info(f"Saved method statistics: {METHOD_COMPARISON_DIR / 'method_statistics.csv'}")


# =============================================================================
# Step 6: Export to NetCDF
# =============================================================================

def export_netcdf(dataset_name: str):
    """
    Export dataset to NetCDF format.
    
    Combines all monthly effective precipitation rasters into a single
    NetCDF file with CF-compliant metadata.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('GridMET', 'PRISM', 'ERA5Land', or 'TerraClimate').
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    output_file = ANALYSIS_DIR / f'{dataset_name}_usda_scs_peff_{START_YEAR}_{END_YEAR}.nc'
    
    if output_file.exists():
        logger.info(f"NetCDF export already exists: {output_file}")
        return
    
    logger.info(f"Exporting {dataset_name} to NetCDF...")
    
    try:
        export_to_netcdf(
            input_dir=str(input_dir),
            output_path=str(output_file),
            pattern='effective_precip_[0-9]*.tif',
            variable_name='effective_precipitation'
        )
        logger.info(f"NetCDF exported: {output_file}")
    except Exception as e:
        logger.warning(f"Could not export to NetCDF: {e}")


# =============================================================================
# Main Workflow
# =============================================================================

def run_full_workflow(skip_processing: bool = True, n_workers: int = 4):
    """
    Run the complete Arizona workflow with U.S. vs Global comparison and method comparison.
    
    The workflow uses an efficient approach:
    1. Download rasters once using USDA-SCS method (saves P, AWC, ETo)
    2. Calculate all 9 Peff methods locally from saved rasters
    3. Perform temporal aggregation, statistics, and visualization
    
    Parameters
    ----------
    skip_processing : bool
        If True, skip GEE processing if data already exists.
    n_workers : int
        Number of parallel workers for GEE processing.
    """
    logger.info("=" * 60)
    logger.info("pyCropWat Arizona Workflow - U.S. vs Global & Method Comparison")
    logger.info("Study Area: Arizona")
    logger.info(f"Period: {START_YEAR}-{END_YEAR}")
    logger.info(f"Method: Download with USDA-SCS, then calculate all 8 methods locally")
    logger.info("=" * 60)
    logger.info(f"U.S. Datasets: {list(US_DATASETS.keys())}")
    logger.info(f"Global Datasets: {list(GLOBAL_DATASETS.keys())}")
    logger.info(f"Peff Methods: {PEFF_METHODS}")
    logger.info("=" * 60)
    
    # Create output directories
    create_output_directories()
    
    # Create sample zones for zonal statistics
    create_sample_zones()
    
    # Process ALL datasets (U.S. and Global)
    all_datasets = list(DATASETS.keys())
    
    for dataset_name in all_datasets:
        logger.info(f"\n{'=' * 40}")
        config = DATASETS[dataset_name]
        logger.info(f"Processing {dataset_name} ({config.get('region', 'Unknown')})")
        logger.info("=" * 40)
        
        # Step 1a: Download rasters using USDA-SCS method (saves P, AWC, ETo)
        download_usda_scs_rasters(
            dataset_name, 
            skip_if_exists=skip_processing, 
            n_workers=n_workers
        )
        
        # Step 1b: Calculate all 9 Peff methods locally from saved rasters
        calculate_all_methods_from_rasters(dataset_name)
        
        # Step 2: Temporal aggregation
        agg = run_temporal_aggregation(dataset_name)
        
        # Step 3: Statistical analysis (includes zonal stats)
        stats = run_statistical_analysis(dataset_name)
        
        # Step 4: Visualization
        viz = run_visualization(dataset_name)
        
        # Step 4b: Anomaly, climatology, and trend maps
        plot_anomaly_maps(dataset_name)
        plot_climatology_maps(dataset_name)
        plot_trend_maps(dataset_name)
        
        # Export to NetCDF
        export_netcdf(dataset_name)
    
    
    # Step 5a: U.S. Dataset comparison (GridMET vs PRISM)
    logger.info(f"\n{'=' * 40}")
    logger.info("U.S. Dataset Comparison (GridMET vs PRISM)")
    logger.info("=" * 40)
    compare_datasets()
    
    # Step 5b: U.S. vs Global dataset comparison
    logger.info(f"\n{'=' * 40}")
    logger.info("U.S. vs Global Dataset Comparison")
    logger.info("=" * 40)
    compare_us_vs_global()
    
    # Step 5c: Effective precipitation method comparison
    logger.info(f"\n{'=' * 40}")
    logger.info("Effective Precipitation Method Comparison")
    logger.info("=" * 40)
    compare_peff_methods()
    
    logger.info("\n" + "=" * 60)
    logger.info("Arizona Workflow Complete!")
    logger.info("=" * 60)
    logger.info(f"Analysis outputs saved to: {ANALYSIS_DIR.absolute()}")
    logger.info(f"  - U.S. vs Global comparison: {US_VS_GLOBAL_DIR}")
    logger.info(f"  - Method comparison: {METHOD_COMPARISON_DIR}")


def run_analysis_only():
    """
    Run analysis workflow only (without GEE processing).
    
    Use this if you already have processed effective precipitation data
    from a previous run. Calculates all methods locally and performs 
    temporal aggregation, statistical analysis, visualization, and comparisons.
    """
    logger.info("Running analysis-only workflow...")
    
    create_output_directories()
    create_sample_zones()
    
    # Process all available datasets
    for dataset_name in DATASETS.keys():
        config = DATASETS[dataset_name]
        input_dir = Path(config['output_dir'])
        
        if not input_dir.exists():
            logger.warning(f"Input directory not found: {input_dir}")
            logger.warning(f"Run full workflow or process {dataset_name} data first")
            continue
        
        logger.info(f"\nAnalyzing {dataset_name} ({config.get('region', 'Unknown')})...")
        
        # Calculate all methods from rasters (if input data exists)
        calculate_all_methods_from_rasters(dataset_name)
        
        agg = run_temporal_aggregation(dataset_name)
        stats = run_statistical_analysis(dataset_name)
        viz = run_visualization(dataset_name)
        
        # Anomaly, climatology, and trend maps
        plot_anomaly_maps(dataset_name)
        plot_climatology_maps(dataset_name)
        plot_trend_maps(dataset_name)
        
        export_netcdf(dataset_name)
    
    # Run all comparisons
    compare_datasets()
    compare_us_vs_global()
    compare_peff_methods()
    
    logger.info("\nAnalysis workflow complete!")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='pyCropWat Arizona USDA-SCS Workflow Example',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--analysis-only', '-a',
        action='store_true',
        help='Run analysis only (skip GEE processing)'
    )
    parser.add_argument(
        '--force-reprocess', '-f',
        action='store_true',
        help='Force reprocessing even if data exists'
    )
    parser.add_argument(
        '--gee-project', '-p',
        type=str,
        help='GEE project ID for authentication'
    )
    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=4,
        help='Number of parallel workers for GEE processing (default: 4)'
    )
    
    args = parser.parse_args()
    
    # Update GEE project if provided
    if args.gee_project:
        GEE_PROJECT = args.gee_project
    
    # Run appropriate workflow
    if args.analysis_only:
        run_analysis_only()
    else:
        run_full_workflow(skip_processing=not args.force_reprocess, n_workers=args.workers)

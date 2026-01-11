"""
pyCropWat New Mexico Method Comparison Example
===============================================

This script demonstrates comparison of all effective precipitation methods
for New Mexico using PRISM precipitation data.

The workflow includes:
1. Download rasters using USDA-SCS method (saves P, AWC, ETo)
2. Calculate all other methods locally from saved rasters
3. Temporal aggregation (annual, seasonal)
4. Statistical analysis (anomalies, trends)
5. Visualization of method differences

This approach is more efficient than downloading for each method separately.

Precipitation Data Source:
- PRISM (projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY) - ~800m monthly

USDA-SCS Required Data (US):
- AWC: SSURGO (projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite)
- ETo: gridMET Monthly (projects/openet/assets/reference_et/conus/gridmet/monthly/v1)

Effective Precipitation Methods Compared:
- CROPWAT - Default method from FAO CROPWAT
- FAO/AGLW - FAO Dependable Rainfall (80% exceedance)
- Fixed Percentage (70%) - Simple empirical method
- Dependable Rainfall (80% probability) - Statistical approach
- FarmWest - WSU irrigation scheduling formula
- USDA-SCS - Site-specific method with AWC and ETo (SSURGO + gridMET)
- TAGEM-SuET - Turkish Irrigation Management System based on P - ETo
- Ensemble - Mean of 6 methods (excludes TAGEM-SuET)

Study Area:
- New Mexico (NM.geojson)

Requirements:
- Google Earth Engine account
- pycropwat package installed

Usage:
    python new_mexico_example.py
    python new_mexico_example.py -f -w 8  # Force reprocess with 8 workers
    python new_mexico_example.py --analysis-only
    python new_mexico_example.py --gee-project your-project-id
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
STUDY_AREA_GEOJSON = "./NM.geojson"  # Local GeoJSON file for New Mexico boundary

# PRISM Dataset Configuration
PRISM_CONFIG = {
    'asset_id': 'projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY',
    'precip_band': 'ppt',
    'scale_factor': 1.0,  # Already in mm
    'description': 'Oregon State PRISM 800m monthly precipitation'
}

# USDA-SCS Required Assets (for usda_scs method)
AWC_ASSET = "projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite"
AWC_BAND = None  # Single band asset
ETO_ASSET = "projects/openet/assets/reference_et/conus/gridmet/monthly/v1"
ETO_BAND = "eto"
ROOTING_DEPTH = 1.0  # meters

# Effective precipitation methods to compare
PEFF_METHODS = {
    'cropwat': {
        'name': 'CROPWAT',
        'color': '#1f77b4',
        'description': 'CROPWAT method from FAO',
        'params': {}
    },
    'fao_aglw': {
        'name': 'FAO/AGLW',
        'color': '#ff7f0e',
        'description': 'FAO Dependable Rainfall (80% exceedance)',
        'params': {}
    },
    'fixed_percentage': {
        'name': 'Fixed 70%',
        'color': '#2ca02c',
        'description': 'Simple 70% of total precipitation',
        'params': {'percentage': 0.7}
    },
    'dependable_rainfall': {
        'name': 'Dependable Rain (75%)',
        'color': '#d62728',
        'description': 'FAO dependable rainfall at 75% probability',
        'params': {'probability': 0.75}
    },
    'farmwest': {
        'name': 'FarmWest',
        'color': '#9467bd',
        'description': 'WSU FarmWest formula: (P-5) × 0.75',
        'params': {}
    },
    'usda_scs': {
        'name': 'USDA-SCS',
        'color': '#8c564b',
        'description': 'Site-specific method with AWC and ETo',
        'params': {
            'awc_asset': AWC_ASSET,
            'awc_band': AWC_BAND,
            'eto_asset': ETO_ASSET,
            'eto_band': ETO_BAND,
            'eto_is_daily': False,
            'eto_scale_factor': 1.0,
            'rooting_depth': ROOTING_DEPTH
        }
    },
    'suet': {
        'name': 'TAGEM-SuET',
        'color': '#e377c2',
        'description': 'Turkish Irrigation Management System (P - ETo)',
        'params': {
            'eto_asset': ETO_ASSET,
            'eto_band': ETO_BAND,
            'eto_is_daily': False,
            'eto_scale_factor': 1.0
        }
    },
    'ensemble': {
        'name': 'Ensemble',
        'color': '#17becf',
        'description': 'Mean of 6 methods (excludes TAGEM-SuET)',
        'params': {
            'awc_asset': AWC_ASSET,
            'awc_band': AWC_BAND,
            'eto_asset': ETO_ASSET,
            'eto_band': ETO_BAND,
            'eto_is_daily': False,
            'eto_scale_factor': 1.0,
            'rooting_depth': ROOTING_DEPTH
        }
    }
}

# Region-specific output directory
REGION_DIR = Path('./NewMexico')

# Primary output directory (USDA-SCS downloads go here)
PRIMARY_OUTPUT_DIR = REGION_DIR / 'NM_PRISM_USDA_SCS'
PRIMARY_INPUT_DIR = REGION_DIR / 'analysis_inputs'

# Time period
START_YEAR = 1986
END_YEAR = 2025

# Climatology period (for anomaly calculations)
CLIM_START = 1986
CLIM_END = 2020

# Output scale (resolution in meters)
OUTPUT_SCALE = None  # Use native resolution

# Analysis output directories
ANALYSIS_DIR = REGION_DIR / 'analysis_outputs'
ANNUAL_DIR = ANALYSIS_DIR / 'annual'
CLIMATOLOGY_DIR = ANALYSIS_DIR / 'climatology'
ANOMALY_DIR = ANALYSIS_DIR / 'anomalies'
TREND_DIR = ANALYSIS_DIR / 'trend'
FIGURES_DIR = ANALYSIS_DIR / 'figures'
METHOD_COMPARISON_DIR = ANALYSIS_DIR / 'method_comparison'
ZONAL_DIR = ANALYSIS_DIR / 'zonal_stats'

# Zone file for zonal statistics
ZONE_FILE = ANALYSIS_DIR / 'nm_zones.geojson'


def create_output_directories():
    """Create all necessary output directories."""
    directories = [
        REGION_DIR,
        PRIMARY_OUTPUT_DIR,
        PRIMARY_INPUT_DIR,
        ANALYSIS_DIR,
        ANNUAL_DIR,
        CLIMATOLOGY_DIR,
        ANOMALY_DIR,
        TREND_DIR,
        FIGURES_DIR,
        METHOD_COMPARISON_DIR,
        ZONAL_DIR
    ]
    
    # Create method-specific directories for analysis outputs
    for method in PEFF_METHODS.keys():
        directories.append(ANNUAL_DIR / method)
        directories.append(CLIMATOLOGY_DIR / method)
        directories.append(ANOMALY_DIR / method)
        directories.append(TREND_DIR / method)
    
    for d in directories:
        d.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Output directories created in {REGION_DIR}")


def create_sample_zones():
    """
    Create sample zone polygons for zonal statistics demonstration.
    
    Creates zones representing different New Mexico regions:
    - Northern NM (Santa Fe/Taos - Rio Grande Valley)
    - Central NM (Albuquerque - Middle Rio Grande)
    - Southern NM (Las Cruces - Mesilla Valley)
    - Eastern NM (Clovis/Portales - Llano Estacado)
    """
    import json
    
    if ZONE_FILE.exists():
        logger.info(f"Zone file already exists: {ZONE_FILE}")
        return
    
    logger.info("Creating sample zone polygons...")
    
    zones_geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 1,
                    "name": "Northern NM",
                    "description": "Northern New Mexico - Santa Fe/Taos/Rio Grande Valley"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-106.5, 35.5],
                        [-106.5, 37.0],
                        [-105.0, 37.0],
                        [-105.0, 35.5],
                        [-106.5, 35.5]
                    ]]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 2,
                    "name": "Central NM",
                    "description": "Central New Mexico - Albuquerque/Middle Rio Grande"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-107.0, 34.5],
                        [-107.0, 35.5],
                        [-106.0, 35.5],
                        [-106.0, 34.5],
                        [-107.0, 34.5]
                    ]]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 3,
                    "name": "Southern NM",
                    "description": "Southern New Mexico - Las Cruces/Mesilla Valley"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-107.5, 32.0],
                        [-107.5, 33.0],
                        [-106.0, 33.0],
                        [-106.0, 32.0],
                        [-107.5, 32.0]
                    ]]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 4,
                    "name": "Eastern NM",
                    "description": "Eastern New Mexico - Clovis/Portales/Llano Estacado"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-104.5, 33.5],
                        [-104.5, 35.0],
                        [-103.0, 35.0],
                        [-103.0, 33.5],
                        [-104.5, 33.5]
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
    skip_if_exists: bool = True,
    n_workers: int = 4
):
    """
    Download rasters using USDA-SCS method.
    
    This downloads:
    - Effective precipitation (USDA-SCS method)
    - Total precipitation (via fraction file)
    - AWC data (SSURGO)
    - ETo data (gridMET)
    
    Other methods will be calculated locally from these saved rasters.
    
    Parameters
    ----------
    skip_if_exists : bool
        Skip processing if output files already exist
    n_workers : int
        Number of parallel workers for processing
    """
    # Check if data already exists
    if skip_if_exists and PRIMARY_OUTPUT_DIR.exists():
        existing_files = list(PRIMARY_OUTPUT_DIR.glob('effective_precip_[0-9]*.tif'))
        expected_files = (END_YEAR - START_YEAR + 1) * 12
        if len(existing_files) >= expected_files * 0.9:
            logger.info(f"Skipping USDA-SCS download - data already exists ({len(existing_files)} files)")
            return
    
    logger.info("Downloading rasters using USDA-SCS method...")
    logger.info("  This will save P, AWC (SSURGO), and ETo (gridMET) data")
    
    # USDA-SCS method parameters
    method_params = PEFF_METHODS['usda_scs']['params']
    
    ep = EffectivePrecipitation(
        asset_id=PRISM_CONFIG['asset_id'],
        precip_band=PRISM_CONFIG['precip_band'],
        geometry_path=STUDY_AREA_GEOJSON,
        start_year=START_YEAR,
        end_year=END_YEAR,
        precip_scale_factor=PRISM_CONFIG['scale_factor'],
        scale=OUTPUT_SCALE,
        gee_project=GEE_PROJECT,
        method='usda_scs',
        method_params=method_params
    )
    
    ep.process(
        output_dir=str(PRIMARY_OUTPUT_DIR),
        n_workers=n_workers,
        save_inputs=True,
        input_dir=str(PRIMARY_INPUT_DIR)
    )
    
    logger.info("Completed USDA-SCS raster download")


def calculate_all_methods_from_rasters():
    """
    Calculate effective precipitation using all methods from saved rasters.
    
    Loads the precipitation, AWC, and ETo data saved during USDA-SCS processing
    and calculates Peff for each method locally (without additional GEE calls).
    
    Outputs are saved as method-specific GeoTIFF files for analysis.
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
        suet_effective_precip,
        ensemble_effective_precip
    )
    
    logger.info("Calculating all methods from saved rasters...")
    
    # Load AWC data (static, only once)
    # SSURGO AWC is in inches/inch (volumetric fraction)
    awc_file = PRIMARY_INPUT_DIR / 'awc.tif'
    if awc_file.exists():
        da_awc = rioxarray.open_rasterio(awc_file).squeeze('band', drop=True)
        awc_data = da_awc.values
        logger.info(f"Loaded AWC data: shape={awc_data.shape}, mean={np.nanmean(awc_data):.2f} in/in")
    else:
        logger.error(f"AWC file not found: {awc_file}")
        logger.error("Run download_usda_scs_rasters() first")
        return
    
    # Process each month
    for year in range(START_YEAR, END_YEAR + 1):
        
        for month in range(1, 13):
            # Load precipitation from fraction file and Peff
            peff_file = PRIMARY_OUTPUT_DIR / f'effective_precip_{year}_{month:02d}.tif'
            fraction_file = PRIMARY_OUTPUT_DIR / f'effective_precip_fraction_{year}_{month:02d}.tif'
            
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
            eto_file = PRIMARY_INPUT_DIR / f'eto_{year}_{month:02d}.tif'
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
                logger.warning(f"ETo file not found for {year}-{month:02d}, using regional mean")
                eto_data = np.full_like(precip, 160.0)
            
            # Resample AWC if needed
            awc_resampled = awc_data
            if awc_data.shape != precip.shape:
                from scipy.ndimage import zoom
                zoom_factors = (precip.shape[0] / awc_data.shape[0],
                               precip.shape[1] / awc_data.shape[1])
                awc_resampled = zoom(awc_data, zoom_factors, order=1)
            
            # Calculate Peff for each method (suppress warnings for masked data)
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
                method_dir = ANALYSIS_DIR / 'peff_by_method' / method_name
                method_dir.mkdir(parents=True, exist_ok=True)
                
                output_file = method_dir / f'effective_precip_{year}_{month:02d}.tif'
                
                # Create DataArray with same coordinates as template
                da_result = da_peff.copy(data=peff_values)
                da_result = da_result.where(~np.isnan(da_peff.values))
                
                da_result.attrs = {
                    'units': 'mm',
                    'long_name': f'effective_precipitation_{method_name}',
                    'method': method_name,
                    'year': year,
                    'month': month
                }
                da_result = da_result.rio.write_crs("EPSG:4326")
                da_result.rio.to_raster(output_file)
            
            if month == 1:
                logger.info(f"  Processed {year}...")
    
    logger.info("Completed calculation of all methods from rasters")


# =============================================================================
# Step 2: Temporal Aggregation
# =============================================================================

def run_temporal_aggregation(method_name: str):
    """
    Run temporal aggregation for a method's outputs.
    
    Parameters
    ----------
    method_name : str
        Name of the method.
        
    Returns
    -------
    TemporalAggregator
        The aggregator instance.
    """
    # Use the calculated method outputs
    input_dir = ANALYSIS_DIR / 'peff_by_method' / method_name
    
    if not input_dir.exists():
        logger.warning(f"Input directory not found: {input_dir}")
        return None
    
    logger.info(f"Running temporal aggregation for {method_name}...")
    
    agg = TemporalAggregator(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    # Dataset-specific output subdirectories
    method_annual_dir = ANNUAL_DIR / method_name
    method_clim_dir = CLIMATOLOGY_DIR / method_name
    
    # Annual totals
    logger.info("Computing annual totals...")
    for year in range(START_YEAR, END_YEAR + 1):
        output_path = method_annual_dir / f'annual_{year}.tif'
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
        output_dir=str(method_clim_dir)
    )
    
    logger.info(f"Completed temporal aggregation for {method_name}")
    return agg


# =============================================================================
# Step 3: Statistical Analysis
# =============================================================================

def run_statistical_analysis(method_name: str):
    """
    Run statistical analysis for a method's outputs.
    
    Parameters
    ----------
    method_name : str
        Name of the method.
        
    Returns
    -------
    StatisticalAnalyzer
        The analyzer instance.
    """
    # Use the calculated method outputs
    input_dir = ANALYSIS_DIR / 'peff_by_method' / method_name
    
    if not input_dir.exists():
        logger.warning(f"Input directory not found: {input_dir}")
        return None
    
    logger.info(f"Running statistical analysis for {method_name}...")
    
    stats = StatisticalAnalyzer(str(input_dir), pattern='effective_precip_[0-9]*.tif')
    
    method_anomaly_dir = ANOMALY_DIR / method_name
    method_trend_dir = TREND_DIR / method_name
    
    # Calculate anomalies for recent years
    logger.info("Computing anomalies...")
    for year in range(2021, min(END_YEAR + 1, 2025)):
        for month in range(1, 13):
            output_path = method_anomaly_dir / f'anomaly_{year}_{month:02d}.tif'
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
            output_dir=str(method_trend_dir)
        )
    except Exception as e:
        logger.warning(f"Could not calculate trend: {e}")
    
    logger.info(f"Completed statistical analysis for {method_name}")
    return stats


# =============================================================================
# Step 4: Method Comparison Visualization
# =============================================================================

def plot_method_maps(sample_year: int = 2020, sample_month: int = 8):
    """
    Create side-by-side maps comparing all methods for a sample month.
    
    Parameters
    ----------
    sample_year : int
        Year to plot
    sample_month : int
        Month to plot (default: August - monsoon season)
    """
    import matplotlib.pyplot as plt
    import rioxarray
    import numpy as np
    
    logger.info(f"Creating method comparison maps for {sample_year}-{sample_month:02d}...")
    
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    axes = axes.flatten()
    
    all_data = {}
    vmin, vmax = float('inf'), float('-inf')
    
    # Load all method data
    for method_name in PEFF_METHODS.keys():
        input_dir = ANALYSIS_DIR / 'peff_by_method' / method_name
        peff_file = input_dir / f'effective_precip_{sample_year}_{sample_month:02d}.tif'
        
        if peff_file.exists():
            da = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
            da = da.where(da != 0)
            all_data[method_name] = da
            
            valid_data = da.values[~np.isnan(da.values)]
            if len(valid_data) > 0:
                vmin = min(vmin, np.nanpercentile(valid_data, 2))
                vmax = max(vmax, np.nanpercentile(valid_data, 98))
    
    if not all_data:
        logger.warning("No data found for method comparison maps")
        return
    
    # Create plots
    for idx, (method_name, da) in enumerate(all_data.items()):
        ax = axes[idx]
        method_config = PEFF_METHODS[method_name]
        
        im = da.plot(
            ax=ax,
            cmap='YlGnBu',
            vmin=vmin,
            vmax=vmax,
            add_colorbar=False
        )
        
        ax.set_title(f"{method_config['name']}", fontsize=16, fontweight='bold')
        ax.set_xlabel('Longitude [°]', fontsize=16)
        ax.set_ylabel('Latitude [°]', fontsize=16)
    
    # Hide empty axes
    for idx in range(len(all_data), len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle(
        f'PRISM Effective Precipitation - Method Comparison\nNew Mexico, August {sample_year}',
        fontsize=16, fontweight='bold'
    )
    
    plt.tight_layout()
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label='Effective Precipitation [mm]')
    
    output_path = METHOD_COMPARISON_DIR / f'method_maps_{sample_year}_{sample_month:02d}.png'
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method comparison maps: {output_path}")


def compare_mean_annual_maps():
    """
    Compare mean annual (1986-2025) effective precipitation maps for each method.
    
    Computes the long-term mean of annual totals for each method and creates
    side-by-side comparison maps. Shows 6 methods (excluding SuET and Ensemble),
    the pre-computed Ensemble method, and the coefficient of variation.
    """
    import matplotlib.pyplot as plt
    import rioxarray
    import numpy as np
    import xarray as xr
    
    logger.info(f"Computing mean annual effective precipitation maps ({START_YEAR}-{END_YEAR})...")
    
    mean_annual_data = {}
    
    # Exclude 'suet' and 'ensemble' from first 6 panels - ensemble shown in panel 7
    methods_to_plot = [m for m in PEFF_METHODS.keys() if m not in ('suet', 'ensemble')]
    
    for method_name in PEFF_METHODS.keys():
        method_annual_dir = ANNUAL_DIR / method_name
        
        if not method_annual_dir.exists():
            logger.warning(f"Annual directory not found for {method_name}: {method_annual_dir}")
            continue
        
        # Collect all annual files
        annual_arrays = []
        for year in range(START_YEAR, END_YEAR + 1):
            annual_file = method_annual_dir / f'annual_{year}.tif'
            if annual_file.exists():
                try:
                    da = rioxarray.open_rasterio(annual_file).squeeze('band', drop=True)
                    da = da.where(da != 0)
                    annual_arrays.append(da)
                except Exception as e:
                    logger.debug(f"Could not read {annual_file}: {e}")
        
        if annual_arrays:
            # Stack and compute mean
            stacked = xr.concat(annual_arrays, dim='year')
            mean_annual = stacked.mean(dim='year')
            mean_annual_data[method_name] = mean_annual
            logger.info(f"  {method_name}: computed mean from {len(annual_arrays)} years")
    
    if not mean_annual_data:
        logger.warning("No data available for mean annual comparison")
        return
    
    # Create comparison figure
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    axes = axes.flatten()
    
    # Compute common color scale for the method maps (excluding suet)
    vmin, vmax = float('inf'), float('-inf')
    for method_name, da in mean_annual_data.items():
        if method_name in methods_to_plot:
            valid_data = da.values[~np.isnan(da.values)]
            if len(valid_data) > 0:
                vmin = min(vmin, np.nanpercentile(valid_data, 2))
                vmax = max(vmax, np.nanpercentile(valid_data, 98))
    
    # Plot each method (first 6 panels, excluding suet)
    im = None
    panel_idx = 0
    for method_name in methods_to_plot:
        if method_name not in mean_annual_data:
            continue
        da = mean_annual_data[method_name]
        ax = axes[panel_idx]
        method_config = PEFF_METHODS[method_name]
        
        im = da.plot(
            ax=ax,
            cmap='YlGnBu',
            vmin=vmin,
            vmax=vmax,
            add_colorbar=False
        )
        
        # Compute statistics for title
        mean_val = np.nanmean(da.values)
        std_val = np.nanstd(da.values)
        ax.set_title(f"{method_config['name']}\nMean: {mean_val:.0f} ± {std_val:.0f} mm/yr", 
                    fontsize=16, fontweight='bold')
        ax.set_xlabel('Longitude [°]', fontsize=16)
        ax.set_ylabel('Latitude [°]', fontsize=16)
        panel_idx += 1
    
    # 7th panel: Ensemble Method (use pre-computed ensemble rasters)
    ax_mean = axes[6]
    
    # Use ensemble data directly instead of computing mean
    if 'ensemble' in mean_annual_data:
        ensemble_data = mean_annual_data['ensemble']
    else:
        # Fallback: compute mean if ensemble not available
        logger.warning("Ensemble rasters not found, computing mean of methods")
        methods_for_mean = [mean_annual_data[m] for m in methods_to_plot if m in mean_annual_data]
        method_stack = xr.concat(methods_for_mean, dim='method')
        ensemble_data = method_stack.mean(dim='method')
    
    # Plot ensemble with same colormap
    im_mean = ensemble_data.plot(
        ax=ax_mean,
        cmap='YlGnBu',
        vmin=vmin,
        vmax=vmax,
        add_colorbar=False
    )
    
    mean_val = np.nanmean(ensemble_data.values)
    std_val = np.nanstd(ensemble_data.values)
    ax_mean.set_title(f"Ensemble\nMean: {mean_val:.0f} ± {std_val:.0f} mm/yr", 
                     fontsize=16, fontweight='bold')
    ax_mean.set_xlabel('Longitude [°]', fontsize=16)
    ax_mean.set_ylabel('Latitude [°]', fontsize=16)
    
    # 8th panel: Coefficient of Variation (CV) across methods (excluding suet and ensemble)
    ax_cv = axes[7]
    
    # Compute CV from method data (excluding suet and ensemble)
    methods_for_cv = [m for m in methods_to_plot if m in mean_annual_data and m != 'ensemble']
    cv_stack = xr.concat([mean_annual_data[m] for m in methods_for_cv], dim='method')
    method_mean_cv = cv_stack.mean(dim='method')
    method_std = cv_stack.std(dim='method')
    cv = (method_std / method_mean_cv) * 100  # CV in percent
    
    # Plot CV with its own colormap and colorbar
    im_cv = cv.plot(
        ax=ax_cv,
        cmap='YlOrRd',
        vmin=0,
        vmax=np.nanpercentile(cv.values, 98),
        add_colorbar=False
    )
    
    cv_mean = np.nanmean(cv.values)
    cv_std = np.nanstd(cv.values)
    ax_cv.set_title(f"Coefficient of Variation\nCV: {cv_mean:.1f} ± {cv_std:.1f} %", 
                   fontsize=16, fontweight='bold')
    ax_cv.set_xlabel('Longitude [°]', fontsize=16)
    ax_cv.set_ylabel('Latitude [°]', fontsize=16)
    
    # Add CV colorbar below the CV panel
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax_cv)
    cax_cv = divider.append_axes("bottom", size="5%", pad=0.6)
    fig.colorbar(im_cv, cax=cax_cv, orientation='horizontal', label='CV [%]')
    
    plt.suptitle(
        f'Mean Annual Effective Precipitation by Method\nNew Mexico, {START_YEAR}-{END_YEAR}',
        fontsize=16, fontweight='bold'
    )
    
    plt.tight_layout()
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label='Mean Annual Effective Precipitation [mm/yr]')
    
    output_path = METHOD_COMPARISON_DIR / f'mean_annual_comparison_{START_YEAR}_{END_YEAR}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved mean annual comparison: {output_path}")
    
    # Create difference maps (relative to CROPWAT as reference)
    if 'cropwat' in mean_annual_data:
        fig, axes = plt.subplots(2, 4, figsize=(24, 12))
        axes = axes.flatten()
        
        reference = mean_annual_data['cropwat']
        
        # Plot methods (excluding suet)
        panel_idx = 0
        for method_name in methods_to_plot:
            if method_name not in mean_annual_data:
                continue
            da = mean_annual_data[method_name]
            ax = axes[panel_idx]
            method_config = PEFF_METHODS[method_name]
            
            if method_name == 'cropwat':
                # Show reference map
                im_ref = da.plot(
                    ax=ax,
                    cmap='YlGnBu',
                    add_colorbar=False
                )
                mean_val = np.nanmean(da.values)
                std_val = np.nanstd(da.values)
                ax.set_title(f"{method_config['name']} (Reference)\nMean: {mean_val:.0f} ± {std_val:.0f} mm/yr", 
                            fontsize=16, fontweight='bold')
            else:
                # Show difference
                diff = da - reference
                diff_mean = np.nanmean(diff.values)
                diff_std = np.nanstd(diff.values)
                
                # Use diverging colormap
                vmax_diff = max(abs(np.nanpercentile(diff.values, 5)), 
                               abs(np.nanpercentile(diff.values, 95)))
                
                im_diff = diff.plot(
                    ax=ax,
                    cmap='RdBu',
                    vmin=-vmax_diff,
                    vmax=vmax_diff,
                    add_colorbar=True,
                    cbar_kwargs={'label': 'Difference [mm/yr]', 'shrink': 0.8}
                )
                ax.set_title(f"{method_config['name']} - CROPWAT\nMean diff: {diff_mean:+.0f} ± {diff_std:.0f} mm/yr", 
                            fontsize=16, fontweight='bold')
            
            ax.set_xlabel('Longitude [°]', fontsize=16)
            ax.set_ylabel('Latitude [°]', fontsize=16)
            panel_idx += 1
        
        # 7th panel: Ensemble - CROPWAT difference
        ax_mean_diff = axes[6]
        
        # Use ensemble data directly
        if 'ensemble' in mean_annual_data:
            ensemble_for_diff = mean_annual_data['ensemble']
        else:
            # Fallback: compute mean if ensemble not available
            methods_for_mean = [mean_annual_data[m] for m in methods_to_plot if m in mean_annual_data]
            ensemble_for_diff = xr.concat(methods_for_mean, dim='method').mean(dim='method')
        
        diff_mean_all = ensemble_for_diff - reference
        diff_mean_val = np.nanmean(diff_mean_all.values)
        diff_std_val = np.nanstd(diff_mean_all.values)
        
        vmax_diff = max(abs(np.nanpercentile(diff_mean_all.values, 5)), 
                       abs(np.nanpercentile(diff_mean_all.values, 95)))
        
        diff_mean_all.plot(
            ax=ax_mean_diff,
            cmap='RdBu',
            vmin=-vmax_diff,
            vmax=vmax_diff,
            add_colorbar=True,
            cbar_kwargs={'label': 'Difference [mm/yr]', 'shrink': 0.8}
        )
        ax_mean_diff.set_title(f"Ensemble - CROPWAT\nMean diff: {diff_mean_val:+.0f} ± {diff_std_val:.0f} mm/yr", 
                              fontsize=16, fontweight='bold')
        ax_mean_diff.set_xlabel('Longitude [°]', fontsize=16)
        ax_mean_diff.set_ylabel('Latitude [°]', fontsize=16)
        
        # Hide 8th axis
        axes[7].set_visible(False)
        
        plt.suptitle(
            f'Method Differences from CROPWAT Reference\nNew Mexico, {START_YEAR}-{END_YEAR}',
            fontsize=16, fontweight='bold'
        )
        
        plt.tight_layout()
        
        output_path = METHOD_COMPARISON_DIR / f'mean_annual_differences_{START_YEAR}_{END_YEAR}.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved method differences: {output_path}")
    
    # Save summary statistics to CSV
    import pandas as pd
    stats_list = []
    for method_name, da in mean_annual_data.items():
        valid_data = da.values[~np.isnan(da.values)]
        stats_list.append({
            'method': method_name,
            'method_name': PEFF_METHODS[method_name]['name'],
            'mean_mm_yr': np.mean(valid_data),
            'std_mm_yr': np.std(valid_data),
            'min_mm_yr': np.min(valid_data),
            'max_mm_yr': np.max(valid_data),
            'median_mm_yr': np.median(valid_data)
        })
    
    stats_df = pd.DataFrame(stats_list)
    stats_csv = METHOD_COMPARISON_DIR / f'mean_annual_statistics_{START_YEAR}_{END_YEAR}.csv'
    stats_df.to_csv(stats_csv, index=False)
    logger.info(f"Saved statistics: {stats_csv}")
    
    # Print summary
    logger.info("\nMean Annual Effective Precipitation Summary:")
    logger.info("-" * 50)
    for _, row in stats_df.iterrows():
        logger.info(f"  {row['method_name']:20s}: {row['mean_mm_yr']:.0f} ± {row['std_mm_yr']:.0f} mm/yr")


def plot_method_curves():
    """
    Plot theoretical effective precipitation curves for all methods.
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
        suet_effective_precip
    )
    
    logger.info("Creating theoretical method curves...")
    
    # Create precipitation range (NM typically 0-100mm/month)
    precip = np.linspace(0, 150, 1000)
    
    # For theoretical curves, use regional mean AWC and ETo values
    # These are illustrative - actual spatial data is used in GEE processing
    # AWC: SSURGO regional mean (~120 mm/m for New Mexico sandy loam)
    # ETo: gridMET regional mean (~160 mm/month high desert climate)
    awc_ssurgo = np.full_like(precip, 120.0)  # mm/m (SSURGO regional mean)
    eto_gridmet = np.full_like(precip, 160.0)  # mm/month (gridMET regional mean)
    
    # Calculate Peff for each method
    methods = {
        'cropwat': cropwat_effective_precip(precip),
        'fao_aglw': fao_aglw_effective_precip(precip),
        'fixed_percentage': fixed_percentage_effective_precip(precip, 0.7),
        'dependable_rainfall': dependable_rainfall_effective_precip(precip, 0.75),
        'farmwest': farmwest_effective_precip(precip),
        'usda_scs': usda_scs_effective_precip(precip, eto_gridmet, awc_ssurgo),
        'suet': suet_effective_precip(precip, eto_gridmet)
    }
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: All methods comparison
    ax = axes[0]
    for method_name, peff in methods.items():
        config = PEFF_METHODS[method_name]
        ax.plot(precip, peff, 
               label=config['name'],
               color=config['color'],
               linewidth=2)
    
    ax.plot(precip, precip, 'k--', alpha=0.3, label='1:1 line')
    ax.set_xlabel('Total Precipitation [mm]')
    ax.set_ylabel('Effective Precipitation [mm]')
    ax.set_title('Peff Methods Comparison (New Mexico Parameters)')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 150)
    ax.set_ylim(0, 120)
    
    # Shade typical NM monsoon range
    ax.axvspan(20, 60, alpha=0.1, color='blue', label='Typical Monsoon Range')
    
    # Plot 2: Efficiency ratio
    ax = axes[1]
    for method_name, peff in methods.items():
        config = PEFF_METHODS[method_name]
        efficiency = np.divide(peff, precip, where=precip > 0, out=np.zeros_like(peff))
        ax.plot(precip[precip > 0], efficiency[precip > 0],
               label=config['name'],
               color=config['color'],
               linewidth=2)
    
    ax.axhline(y=1.0, color='k', linestyle='--', alpha=0.3, label='100%')
    ax.set_xlabel('Total Precipitation [mm]')
    ax.set_ylabel('Efficiency (Peff/P)')
    ax.set_title('Rainfall Effectiveness by Method')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 150)
    ax.set_ylim(0, 1.1)
    
    plt.suptitle('Theoretical Method Curves - New Mexico Conditions', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = METHOD_COMPARISON_DIR / 'method_curves_new_mexico.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method curves: {output_path}")


def compare_methods_statistics():
    """
    Compare methods using statistical summaries across all processed data.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import rioxarray
    
    logger.info("Computing method comparison statistics...")
    
    results = []
    sample_years = range(START_YEAR, END_YEAR + 1)
    sample_months = range(1, 13)
    
    for method_name in PEFF_METHODS.keys():
        input_dir = ANALYSIS_DIR / 'peff_by_method' / method_name
        
        if not input_dir.exists():
            continue
        
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
                    
                    if len(peff_vals) == 0:
                        continue
                    
                    # Get total precipitation from fraction
                    if fraction_file.exists():
                        da_frac = rioxarray.open_rasterio(fraction_file).squeeze('band', drop=True)
                        da_frac = da_frac.where(da_frac > 0)
                        frac_vals = da_frac.values[~np.isnan(da.values)]
                        precip_vals = peff_vals / frac_vals
                        precip_mean = np.nanmean(precip_vals)
                    else:
                        precip_mean = np.nan
                    
                    results.append({
                        'method': method_name,
                        'year': year,
                        'month': month,
                        'mean_peff': np.nanmean(peff_vals),
                        'std_peff': np.nanstd(peff_vals),
                        'min_peff': np.nanmin(peff_vals),
                        'max_peff': np.nanmax(peff_vals),
                        'mean_precip': precip_mean
                    })
                except Exception as e:
                    logger.debug(f"Error processing {peff_file}: {e}")
    
    if not results:
        logger.warning("No data available for method comparison")
        return
    
    df = pd.DataFrame(results)
    
    # Create comparison figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Box plot of mean Peff by method
    ax = axes[0, 0]
    methods_list = list(PEFF_METHODS.keys())
    box_data = [df[df['method'] == m]['mean_peff'].values for m in methods_list 
                if m in df['method'].values]
    labels = [PEFF_METHODS[m]['name'] for m in methods_list if m in df['method'].values]
    colors = [PEFF_METHODS[m]['color'] for m in methods_list if m in df['method'].values]
    
    bp = ax.boxplot(box_data, labels=labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Monthly Peff Distribution by Method')
    ax.grid(True, alpha=0.3, axis='y')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=15)
    
    # Plot 2: Annual means by method
    ax = axes[0, 1]
    annual_means = df.groupby(['method', 'year'])['mean_peff'].sum().reset_index()
    
    for method_name in methods_list:
        if method_name in annual_means['method'].values:
            method_data = annual_means[annual_means['method'] == method_name]
            ax.plot(method_data['year'], method_data['mean_peff'],
                   label=PEFF_METHODS[method_name]['name'],
                   color=PEFF_METHODS[method_name]['color'],
                   linewidth=2, marker='o', markersize=3)
    
    ax.set_xlabel('Year')
    ax.set_ylabel('Annual Effective Precipitation [mm]')
    ax.set_title('Annual Peff Time Series by Method')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Monthly climatology by method
    ax = axes[1, 0]
    monthly_clim = df.groupby(['method', 'month'])['mean_peff'].mean().reset_index()
    
    for method_name in methods_list:
        if method_name in monthly_clim['method'].values:
            method_data = monthly_clim[monthly_clim['method'] == method_name]
            ax.plot(method_data['month'], method_data['mean_peff'],
                   label=PEFF_METHODS[method_name]['name'],
                   color=PEFF_METHODS[method_name]['color'],
                   linewidth=2, marker='o')
    
    ax.set_xlabel('Month')
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Monthly Climatology by Method')
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Method efficiency comparison
    ax = axes[1, 1]
    df_eff = df.dropna(subset=['mean_precip'])
    if len(df_eff) > 0:
        df_eff['efficiency'] = df_eff['mean_peff'] / df_eff['mean_precip']
        eff_by_method = df_eff.groupby('method')['efficiency'].mean()
        
        methods_with_data = [m for m in methods_list if m in eff_by_method.index]
        bars = ax.bar([PEFF_METHODS[m]['name'] for m in methods_with_data],
                     [eff_by_method[m] for m in methods_with_data],
                     color=[PEFF_METHODS[m]['color'] for m in methods_with_data],
                     alpha=0.7)
        ax.set_ylabel('Mean Efficiency (Peff/P)')
        ax.set_title('Average Rainfall Effectiveness')
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3, axis='y')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=15)
    else:
        ax.text(0.5, 0.5, 'Efficiency data not available', 
               ha='center', va='center', transform=ax.transAxes)
    
    plt.suptitle('Effective Precipitation Method Comparison - New Mexico (PRISM)', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = METHOD_COMPARISON_DIR / 'method_comparison_summary.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method comparison summary: {output_path}")
    
    # Save statistics to CSV
    summary_stats = df.groupby('method').agg({
        'mean_peff': ['mean', 'std', 'min', 'max'],
        'mean_precip': 'mean'
    }).round(2)
    summary_stats.to_csv(METHOD_COMPARISON_DIR / 'method_statistics.csv')
    logger.info(f"Saved method statistics: {METHOD_COMPARISON_DIR / 'method_statistics.csv'}")


# =============================================================================
# Step 5: Export to NetCDF
# =============================================================================

def export_netcdf(method_name: str):
    """
    Export method outputs to NetCDF format.
    
    Parameters
    ----------
    method_name : str
        Name of the method.
    """
    input_dir = ANALYSIS_DIR / 'peff_by_method' / method_name
    output_file = ANALYSIS_DIR / f'PRISM_{method_name}_peff_{START_YEAR}_{END_YEAR}.nc'
    
    if output_file.exists():
        logger.info(f"NetCDF export already exists: {output_file}")
        return
    
    if not input_dir.exists():
        logger.warning(f"Input directory not found: {input_dir}")
        return
    
    logger.info(f"Exporting {method_name} to NetCDF...")
    
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
    Run the complete New Mexico method comparison workflow.
    
    This workflow:
    1. Downloads precipitation, AWC, and ETo rasters using USDA-SCS method
    2. Calculates all 8 Peff methods locally from saved rasters
    3. Performs temporal aggregation and statistical analysis
    4. Creates comparison visualizations
    
    Parameters
    ----------
    skip_processing : bool
        If True, skip GEE processing if data already exists.
    n_workers : int
        Number of parallel workers for GEE processing.
    """
    logger.info("=" * 60)
    logger.info("pyCropWat New Mexico Method Comparison Workflow")
    logger.info("Data Source: PRISM 800m Monthly Precipitation")
    logger.info(f"Period: {START_YEAR}-{END_YEAR}")
    logger.info(f"Methods: {list(PEFF_METHODS.keys())}")
    logger.info("=" * 60)
    
    # Create output directories
    create_output_directories()
    
    # Create sample zones
    create_sample_zones()
    
    # Step 1: Download rasters using USDA-SCS method (gets P, AWC, ETo)
    logger.info(f"\n{'=' * 40}")
    logger.info("Step 1: Downloading Rasters via USDA-SCS Method")
    logger.info("=" * 40)
    download_usda_scs_rasters(skip_if_exists=skip_processing, n_workers=n_workers)
    
    # Step 2: Calculate all methods from saved rasters
    logger.info(f"\n{'=' * 40}")
    logger.info("Step 2: Calculating All Methods from Saved Rasters")
    logger.info("=" * 40)
    calculate_all_methods_from_rasters()
    
    # Step 3: Temporal aggregation and statistical analysis for each method
    for method_name in PEFF_METHODS.keys():
        logger.info(f"\n{'=' * 40}")
        logger.info(f"Processing {PEFF_METHODS[method_name]['name']}")
        logger.info("=" * 40)
        
        # Temporal aggregation
        run_temporal_aggregation(method_name)
        
        # Statistical analysis
        run_statistical_analysis(method_name)
        
        # Export to NetCDF
        export_netcdf(method_name)
    
    # Step 4: Method comparison visualizations
    logger.info(f"\n{'=' * 40}")
    logger.info("Method Comparison Analysis")
    logger.info("=" * 40)
    
    plot_method_curves()
    plot_method_maps()
    compare_mean_annual_maps()
    compare_methods_statistics()
    
    logger.info("\n" + "=" * 60)
    logger.info("New Mexico Method Comparison Complete!")
    logger.info("=" * 60)
    logger.info(f"Analysis outputs saved to: {ANALYSIS_DIR.absolute()}")
    logger.info(f"  - Method comparison: {METHOD_COMPARISON_DIR}")


def run_analysis_only(comparison_only: bool = False):
    """
    Run analysis workflow only (without GEE processing).
    
    Assumes rasters have already been downloaded and methods calculated.
    
    Parameters
    ----------
    comparison_only : bool
        If True, skip the per-method loop and only run method comparison.
    """
    logger.info("Running analysis-only workflow...")
    
    create_output_directories()
    create_sample_zones()
    
    if not comparison_only:
        for method_name in PEFF_METHODS.keys():
            input_dir = ANALYSIS_DIR / 'peff_by_method' / method_name
            
            if not input_dir.exists():
                logger.warning(f"Input directory not found: {input_dir}")
                continue
            
            logger.info(f"\nAnalyzing {method_name}...")
            
            run_temporal_aggregation(method_name)
            run_statistical_analysis(method_name)
            export_netcdf(method_name)
    else:
        logger.info("Skipping per-method analysis, running comparison only...")
    
    # Method comparison
    plot_method_curves()
    plot_method_maps()
    compare_mean_annual_maps()
    compare_methods_statistics()
    
    logger.info("\nAnalysis workflow complete!")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='pyCropWat New Mexico Method Comparison Example',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--analysis-only', '-a',
        action='store_true',
        help='Run analysis only (skip GEE processing)'
    )
    parser.add_argument(
        '--comparison-only', '-c',
        action='store_true',
        help='Run method comparison only (skip per-method analysis loop)'
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
    if args.comparison_only:
        run_analysis_only(comparison_only=True)
    elif args.analysis_only:
        run_analysis_only()
    else:
        run_full_workflow(skip_processing=not args.force_reprocess, n_workers=args.workers)

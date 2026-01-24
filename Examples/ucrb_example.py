"""
pyCropWat Upper Colorado River Basin (UCRB) Effective Precipitation Example
============================================================================

This script demonstrates effective precipitation calculation for field-level 
data in the Upper Colorado River Basin using existing precipitation volumes
stored in a GeoPackage file.

Key Features:
- Reads field-level precipitation volumes from a GeoPackage (PPT_VOLUME_ET_Demands layer)
- Uses existing precipitation columns (PPT_VOLUME_mm_yy_acft format) - no GEE download required
- Calculates 6 effective precipitation methods (excludes ensemble and SuET)
- For USDA-SCS method, downloads AWC data from GEE for the entire region
- Calculates mean AWC per field using zonal statistics
- Converts between acre-feet and mm for method application

Unit Conversions:
- Input: Precipitation volumes in acre-feet (acft)
- Divide by ACRES_FTR_GEOM to get depth in feet
- Convert feet to mm (1 ft = 304.8 mm) to apply pyCropWat methods
- Convert effective precipitation mm back to feet, then multiply by area for volume in acft

Effective Precipitation Methods Calculated (8 methods):
- CROPWAT - Method from FAO CROPWAT
- FAO/AGLW - FAO Dependable Rainfall (80% exceedance)
- Fixed Percentage (70%) - Simple empirical method
- Dependable Rainfall (75% probability) - Statistical approach
- FarmWest - WSU irrigation scheduling formula
- USDA-SCS (MAD=0.5) - Site-specific method with AWC, default MAD factor
- USDA-SCS (MAD=1.0) - Site-specific method with AWC, full soil storage
- USDA-SCS (MAD=1.0, RD=2m) - Site-specific method with AWC, full soil storage, 2m rooting depth

Data Requirements:
- GeoPackage: ucrb_field_effective_precip_intercomparison_geopackage.gpkg
- Layer: PPT_VOLUME_ET_Demands
- AWC: SSURGO (projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite) - downloaded from GEE

Output:
- GeoPackage with effective precipitation volumes (acft) for each method
- CSV summary statistics per field

Usage:
    python ucrb_example.py
    python ucrb_example.py --gee-project your-project-id
    python ucrb_example.py --skip-awc-download  # Use cached AWC if available
"""

import re
import logging
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

# Dask for parallelization
import dask
from dask import delayed, compute
from dask.diagnostics import ProgressBar

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import pyCropWat methods
from pycropwat.methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip,
    farmwest_effective_precip,
    usda_scs_effective_precip
)
from pycropwat.utils import initialize_gee

# =============================================================================
# Configuration
# =============================================================================

# Input GeoPackage
GPKG_PATH = Path('./ucrb_field_effective_precip_intercomparison_geopackage.gpkg')
PPT_LAYER_NAME = 'PPT_VOLUME_ET_Demands'
ETO_LAYER_NAME = 'ETDa_VOLUME_ET_Demands'  # ETo (ETc) data for USDA-SCS method

# GEE Configuration
GEE_PROJECT = None  # Set to your GEE project ID, or None to use default

# AWC Asset (SSURGO for US)
AWC_ASSET = "projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite"
AWC_BAND = None  # Single-band asset

# USDA-SCS Parameters
ROOTING_DEPTH = 1.0  # meters

# Output directory
OUTPUT_DIR = Path('./UCRB')

# Unit conversion constants
FT_TO_MM = 304.8  # 1 foot = 304.8 mm
MM_TO_FT = 1 / FT_TO_MM

# Effective precipitation methods to calculate (7 methods, excluding ensemble and SuET)
PEFF_METHODS = {
    'cropwat': {
        'name': 'CROPWAT',
        'function': cropwat_effective_precip,
        'description': 'CROPWAT method from FAO',
        'requires_awc': False,
        'requires_eto': False,
        'mad_factor': None
    },
    'fao_aglw': {
        'name': 'FAO_AGLW',
        'function': fao_aglw_effective_precip,
        'description': 'FAO Dependable Rainfall (80% exceedance)',
        'requires_awc': False,
        'requires_eto': False,
        'mad_factor': None
    },
    'fixed_percentage': {
        'name': 'Fixed_70pct',
        'function': lambda pr: fixed_percentage_effective_precip(pr, percentage=0.7),
        'description': 'Simple 70% of total precipitation',
        'requires_awc': False,
        'requires_eto': False,
        'mad_factor': None
    },
    'dependable_rainfall': {
        'name': 'Dependable_Rain',
        'function': lambda pr: dependable_rainfall_effective_precip(pr, probability=0.75),
        'description': 'FAO dependable rainfall at 75% probability',
        'requires_awc': False,
        'requires_eto': False,
        'mad_factor': None
    },
    'farmwest': {
        'name': 'FarmWest',
        'function': farmwest_effective_precip,
        'description': 'WSU FarmWest formula: (P-5) × 0.75',
        'requires_awc': False,
        'requires_eto': False,
        'mad_factor': None
    },
    'usda_scs': {
        'name': 'USDA_SCS',
        'function': usda_scs_effective_precip,
        'description': 'Site-specific method with AWC and ETo (MAD=0.5)',
        'requires_awc': True,
        'requires_eto': True,
        'mad_factor': 0.5,
        'rooting_depth': ROOTING_DEPTH
    },
    'usda_scs_mad1': {
        'name': 'USDA_SCS_MAD1',
        'function': usda_scs_effective_precip,
        'description': 'Site-specific method with AWC and ETo (MAD=1.0)',
        'requires_awc': True,
        'requires_eto': True,
        'mad_factor': 1.0,
        'rooting_depth': ROOTING_DEPTH
    },
    'usda_scs_mad1_rd2': {
        'name': 'USDA_SCS_MAD1_RD2',
        'function': usda_scs_effective_precip,
        'description': 'Site-specific method with AWC and ETo (MAD=1.0, RD=2m)',
        'requires_awc': True,
        'requires_eto': True,
        'mad_factor': 1.0,
        'rooting_depth': 2.0
    }
}


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Calculate effective precipitation for UCRB fields'
    )
    parser.add_argument(
        '--gee-project',
        type=str,
        default=None,
        help='Google Earth Engine project ID'
    )
    parser.add_argument(
        '--skip-awc-download',
        action='store_true',
        help='Skip AWC download if cached data exists'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=str(OUTPUT_DIR),
        help='Output directory'
    )
    parser.add_argument(
        '-w', '--n-workers',
        type=int,
        default=8,
        help='Number of parallel workers (default: 8)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    return parser.parse_args()


def create_output_directories(output_dir: Path):
    """Create all necessary output directories."""
    directories = [
        output_dir,
        output_dir / 'awc_data',
        output_dir / 'summaries'
    ]
    
    for d in directories:
        d.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Output directories created in {output_dir}")


def parse_ppt_column_name(col_name: str) -> Optional[Tuple[int, int]]:
    """
    Parse precipitation column name to extract month and year.
    
    Format: PPT_VOLUME_mm_yy_acft
    Example: PPT_VOLUME_11_90_acft -> (11, 1990)
             PPT_VOLUME_01_00_acft -> (1, 2000)
    
    Parameters
    ----------
    col_name : str
        Column name to parse
        
    Returns
    -------
    tuple or None
        (month, year) tuple or None if not a valid PPT column
    """
    pattern = r'PPT_VOLUME_(\d{2})_(\d{2})_acft'
    match = re.match(pattern, col_name)
    
    if match:
        month = int(match.group(1))
        yy = int(match.group(2))
        # Convert 2-digit year to 4-digit
        # Assume 90-99 = 1990-1999, 00-89 = 2000-2089
        year = 1900 + yy if yy >= 90 else 2000 + yy
        return (month, year)
    
    return None


def get_ppt_columns(gdf: gpd.GeoDataFrame) -> Dict[str, Tuple[int, int]]:
    """
    Get all precipitation columns and their month/year info.
    
    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Input GeoDataFrame
        
    Returns
    -------
    dict
        Dictionary mapping column names to (month, year) tuples
    """
    ppt_cols = {}
    for col in gdf.columns:
        result = parse_ppt_column_name(col)
        if result:
            ppt_cols[col] = result
    
    logger.info(f"Found {len(ppt_cols)} precipitation columns")
    
    # Log year range
    years = sorted(set(y for _, y in ppt_cols.values()))
    logger.info(f"Year range: {min(years)} - {max(years)}")
    
    return ppt_cols


def download_awc_zonal_stats(
    gdf: gpd.GeoDataFrame,
    awc_cache_path: Path,
    gee_project: Optional[str] = None,
    skip_if_exists: bool = True,
    n_workers: int = 8
) -> pd.Series:
    """
    Download AWC data from GEE and calculate mean AWC per field.
    
    Uses Google Earth Engine to download SSURGO AWC data and calculates
    zonal statistics (mean) for each field polygon. Uses parallel processing
    for efficiency with large numbers of fields.
    
    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame with field geometries (must have 'DRI_ID' column)
    awc_cache_path : Path
        Path to cache AWC values
    gee_project : str, optional
        GEE project ID
    skip_if_exists : bool
        Skip download if cache exists
    n_workers : int
        Number of parallel workers for batch processing
        
    Returns
    -------
    pd.Series
        Mean AWC values indexed by DRI_ID
    """
    import ee
    
    # Check cache
    if skip_if_exists and awc_cache_path.exists():
        logger.info(f"Loading cached AWC data from {awc_cache_path}")
        awc_df = pd.read_csv(awc_cache_path, index_col='DRI_ID')
        return awc_df['mean_awc']
    
    logger.info("Downloading AWC data from Google Earth Engine...")
    
    # Initialize GEE
    initialize_gee(project=gee_project)
    
    # Load AWC image
    awc_img = ee.Image(AWC_ASSET)
    if AWC_BAND:
        awc_img = awc_img.select(AWC_BAND)
    
    # Ensure geometries are in EPSG:4326 for GEE
    gdf_wgs84 = gdf.to_crs('EPSG:4326')
    
    # Calculate zonal stats in batches to avoid GEE limits
    batch_size = 500
    n_fields = len(gdf_wgs84)
    
    # Limit GEE workers to avoid connection pool exhaustion (default pool size is 10)
    gee_workers = min(n_workers, 10)
    logger.info(f"Calculating mean AWC for {n_fields} fields using {gee_workers} GEE workers...")
    
    def process_awc_batch(batch_info):
        """Process a single batch of fields for AWC zonal stats."""
        start_idx, end_idx, batch_gdf = batch_info
        batch_results = {}
        
        # Convert batch to FeatureCollection
        features = []
        for idx, row in batch_gdf.iterrows():
            geom = row.geometry
            # Convert to GeoJSON and create EE geometry
            if geom.geom_type == 'MultiPolygon':
                coords = [list(p.exterior.coords) for p in geom.geoms]
                ee_geom = ee.Geometry.MultiPolygon(coords)
            else:
                coords = list(geom.exterior.coords)
                ee_geom = ee.Geometry.Polygon(coords)
            
            feature = ee.Feature(ee_geom, {'DRI_ID': row['DRI_ID']})
            features.append(feature)
        
        fc = ee.FeatureCollection(features)
        
        # Calculate mean AWC for each feature
        try:
            result = awc_img.reduceRegions(
                collection=fc,
                reducer=ee.Reducer.mean(),
                scale=90  # SSURGO native resolution
            )
            
            # Get results
            result_list = result.getInfo()['features']
            
            for feat in result_list:
                dri_id = feat['properties']['DRI_ID']
                mean_awc = feat['properties'].get('mean', 0.15)  # Default 0.15 if no data
                if mean_awc == 0.15:
                    logger.debug(f"No AWC data for DRI_ID {dri_id}, using default 0.15")
                batch_results[dri_id] = mean_awc
                
        except Exception as e:
            logger.warning(f"Error in batch {start_idx}-{end_idx}: {e}")
            # Use default AWC for failed fields
            for idx, row in batch_gdf.iterrows():
                if row['DRI_ID'] not in batch_results:
                    batch_results[row['DRI_ID']] = 0.15
        
        return batch_results
    
    # Prepare batches
    batches = []
    for start_idx in range(0, n_fields, batch_size):
        end_idx = min(start_idx + batch_size, n_fields)
        batch_gdf = gdf_wgs84.iloc[start_idx:end_idx]
        batches.append((start_idx, end_idx, batch_gdf))
    
    # Process batches in parallel using ThreadPoolExecutor
    # Use gee_workers to limit concurrent connections to GEE
    awc_values = {}
    with ThreadPoolExecutor(max_workers=gee_workers) as executor:
        futures = {executor.submit(process_awc_batch, batch): batch[0] for batch in batches}
        
        with tqdm(total=len(batches), desc="AWC zonal stats") as pbar:
            for future in as_completed(futures):
                batch_results = future.result()
                awc_values.update(batch_results)
                pbar.update(1)
    
    # Create Series
    awc_series = pd.Series(awc_values, name='mean_awc')
    awc_series.index.name = 'DRI_ID'
    
    # Cache results
    awc_df = awc_series.reset_index()
    awc_df.to_csv(awc_cache_path, index=False)
    logger.info(f"Cached AWC data to {awc_cache_path}")
    
    # Log statistics
    logger.info(f"AWC statistics: min={awc_series.min():.4f}, "
                f"max={awc_series.max():.4f}, mean={awc_series.mean():.4f}")
    
    return awc_series


def parse_eto_column_name(col_name: str) -> Optional[Tuple[int, int]]:
    """
    Parse ETo column name to extract month and year.
    
    Format: ETDa_VOLUME_mm_yy_acft
    Example: ETDa_VOLUME_11_90_acft -> (11, 1990)
             ETDa_VOLUME_01_00_acft -> (1, 2000)
    
    Parameters
    ----------
    col_name : str
        Column name to parse
        
    Returns
    -------
    tuple or None
        (month, year) tuple or None if not a valid ETo column
    """
    pattern = r'ETDa_VOLUME_(\d{2})_(\d{2})_acft'
    match = re.match(pattern, col_name)
    
    if match:
        month = int(match.group(1))
        yy = int(match.group(2))
        # Convert 2-digit year to 4-digit
        # Assume 90-99 = 1990-1999, 00-89 = 2000-2089
        year = 1900 + yy if yy >= 90 else 2000 + yy
        return (month, year)
    
    return None


def get_eto_columns(gdf: gpd.GeoDataFrame) -> Dict[str, Tuple[int, int]]:
    """
    Get all ETo columns and their month/year info.
    
    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Input GeoDataFrame with ETo data
        
    Returns
    -------
    dict
        Dictionary mapping column names to (month, year) tuples
    """
    eto_cols = {}
    for col in gdf.columns:
        result = parse_eto_column_name(col)
        if result:
            eto_cols[col] = result
    
    logger.info(f"Found {len(eto_cols)} ETo columns")
    
    return eto_cols


def load_eto_data(gpkg_path: Path) -> Tuple[gpd.GeoDataFrame, Dict[str, Tuple[int, int]]]:
    """
    Load ETo data from the ETDa_VOLUME_ET_Demands layer.
    
    Parameters
    ----------
    gpkg_path : Path
        Path to GeoPackage file
        
    Returns
    -------
    tuple
        (GeoDataFrame with ETo data, dict mapping column names to (month, year))
    """
    logger.info(f"Loading ETo data from '{ETO_LAYER_NAME}' layer...")
    eto_gdf = gpd.read_file(gpkg_path, layer=ETO_LAYER_NAME)
    eto_columns = get_eto_columns(eto_gdf)
    
    return eto_gdf, eto_columns


def acft_to_mm(volume_acft: np.ndarray, area_acres: np.ndarray) -> np.ndarray:
    """
    Convert precipitation volume (acre-feet) to depth (mm).
    
    Parameters
    ----------
    volume_acft : np.ndarray
        Precipitation volume in acre-feet
    area_acres : np.ndarray
        Field area in acres
        
    Returns
    -------
    np.ndarray
        Precipitation depth in mm
    """
    # Volume (acft) / Area (acres) = Depth (ft)
    depth_ft = volume_acft / np.maximum(area_acres, 0.001)  # Avoid division by zero
    # Convert feet to mm
    depth_mm = depth_ft * FT_TO_MM
    return depth_mm


def mm_to_acft(depth_mm: np.ndarray, area_acres: np.ndarray) -> np.ndarray:
    """
    Convert precipitation depth (mm) to volume (acre-feet).
    
    Parameters
    ----------
    depth_mm : np.ndarray
        Precipitation depth in mm
    area_acres : np.ndarray
        Field area in acres
        
    Returns
    -------
    np.ndarray
        Precipitation volume in acre-feet
    """
    # Convert mm to feet
    depth_ft = depth_mm * MM_TO_FT
    # Depth (ft) × Area (acres) = Volume (acft)
    volume_acft = depth_ft * area_acres
    return volume_acft


def process_all_methods(
    gdf: gpd.GeoDataFrame,
    ppt_columns: Dict[str, Tuple[int, int]],
    awc_series: pd.Series,
    eto_gdf: gpd.GeoDataFrame,
    eto_columns: Dict[str, Tuple[int, int]],
    gpkg_path: Path,
    n_workers: int = 8
) -> Dict[str, gpd.GeoDataFrame]:
    """
    Calculate effective precipitation for all methods and all time periods.
    
    Uses Dask for parallel processing across methods and time periods.
    Results are saved as new layers in the input GeoPackage.
    
    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Input GeoDataFrame with precipitation volumes
    ppt_columns : dict
        Dictionary mapping column names to (month, year) tuples
    awc_series : pd.Series
        AWC values indexed by DRI_ID
    eto_gdf : gpd.GeoDataFrame
        GeoDataFrame with ETo volumes
    eto_columns : dict
        Dictionary mapping ETo column names to (month, year) tuples
    gpkg_path : Path
        Path to the input GeoPackage (layers will be added here)
    n_workers : int
        Number of parallel workers
        
    Returns
    -------
    dict
        Dictionary mapping method names to result GeoDataFrames
    """
    results = {}
    
    # Sort columns by date for consistent processing
    sorted_cols = sorted(ppt_columns.items(), key=lambda x: (x[1][1], x[1][0]))
    
    # Create a lookup for ETo columns by (month, year)
    eto_col_lookup = {v: k for k, v in eto_columns.items()}
    
    # Pre-compute AWC values array for faster access
    awc_dict = awc_series.to_dict()
    awc_values_arr = np.array([awc_dict.get(dri_id, 0.15) for dri_id in gdf['DRI_ID']])
    
    # Pre-compute area array
    area_acres = gdf['ACRES_FTR_GEOM'].values.astype(np.float64)
    
    # Pre-compute ETo mapping for all columns
    eto_mappings = {}
    for eto_col, (month, year) in eto_columns.items():
        eto_dict = dict(zip(eto_gdf['DRI_ID'], eto_gdf[eto_col]))
        eto_values = np.array([eto_dict.get(dri_id, 0.0) for dri_id in gdf['DRI_ID']], dtype=np.float64)
        eto_values = np.nan_to_num(eto_values, nan=0.0)
        eto_mappings[(month, year)] = eto_values
    
    # Configure Dask scheduler
    dask.config.set(scheduler='threads', num_workers=n_workers)
    
    logger.info(f"Processing with {n_workers} parallel workers using Dask")
    
    for method_key, method_config in PEFF_METHODS.items():
        logger.info(f"\nProcessing method: {method_config['name']}")
        
        # Create result GeoDataFrame with base columns (fid is auto-generated by GPKG)
        result_gdf = gdf[['DRI_ID', 'ACRES_FTR_GEOM', 'STATE', 'HUC8', 
                          'HUC8_all', 'geometry']].copy()
        
        # Define delayed computation function for a single column
        @delayed
        def compute_peff_column(ppt_col, month, year, method_key, method_func, 
                                requires_awc, requires_eto, mad_factor, rooting_depth):
            """Compute effective precipitation for a single column."""
            try:
                # Get precipitation values
                ppt_volume_acft = gdf[ppt_col].values.astype(np.float64)
                ppt_volume_acft = np.nan_to_num(ppt_volume_acft, nan=0.0)
                
                # Convert to mm depth
                ppt_mm = acft_to_mm(ppt_volume_acft, area_acres)
                
                if method_key in ('usda_scs', 'usda_scs_mad1', 'usda_scs_mad1_rd2'):
                    # Get ETo for this month/year
                    eto_volume_acft = eto_mappings.get((month, year))
                    if eto_volume_acft is None:
                        return None, None
                    eto_mm = acft_to_mm(eto_volume_acft, area_acres)
                    
                    # Calculate Peff with AWC, ETo, and MAD factor
                    peff_mm = method_func(ppt_mm, eto_mm, awc_values_arr, rooting_depth, mad_factor)
                else:
                    # Other methods only need precipitation
                    peff_mm = method_func(ppt_mm)
                
                # Convert back to acre-feet
                peff_acft = mm_to_acft(peff_mm, area_acres)
                
                # Create output column name
                yy = year % 100
                out_col = f"PEFF_{method_config['name']}_{month:02d}_{yy:02d}_acft"
                
                return out_col, peff_acft
                
            except Exception as e:
                logger.warning(f"Error processing {ppt_col} with {method_key}: {e}")
                return None, None
        
        # Create delayed tasks for all columns
        delayed_results = []
        for ppt_col, (month, year) in sorted_cols:
            # Check if ETo is available for USDA-SCS methods
            if method_config['requires_eto'] and (month, year) not in eto_mappings:
                continue
            
            task = compute_peff_column(
                ppt_col, month, year, method_key,
                method_config['function'],
                method_config['requires_awc'],
                method_config['requires_eto'],
                method_config['mad_factor'],
                method_config.get('rooting_depth', ROOTING_DEPTH)
            )
            delayed_results.append(task)
        
        # Execute all tasks in parallel with progress bar
        logger.info(f"  Computing {len(delayed_results)} time periods in parallel...")
        with ProgressBar():
            computed_results = compute(*delayed_results)
        
        # Collect all results into a dictionary first to avoid DataFrame fragmentation
        peff_columns = {}
        for out_col, peff_acft in computed_results:
            if out_col is not None and peff_acft is not None:
                peff_columns[out_col] = peff_acft
        
        # For USDA_SCS methods, include AWC column
        if method_key in ('usda_scs', 'usda_scs_mad1', 'usda_scs_mad1_rd2'):
            peff_columns['AWC'] = awc_values_arr
        
        # Create DataFrame from all Peff columns at once and concat with base columns
        peff_df = pd.DataFrame(peff_columns)
        result_gdf = pd.concat([result_gdf, peff_df], axis=1)
        result_gdf = gpd.GeoDataFrame(result_gdf, geometry='geometry', crs=gdf.crs)
        
        results[method_key] = result_gdf
        
        # Save results as a new layer in the input GeoPackage
        layer_name = f"Peff_{method_config['name']}"
        result_gdf.to_file(gpkg_path, layer=layer_name, driver='GPKG')
        logger.info(f"Added layer '{layer_name}' to {gpkg_path.name} ({len(result_gdf.columns)} columns)")
    
    return results


def create_summary_statistics(
    results: Dict[str, gpd.GeoDataFrame],
    ppt_columns: Dict[str, Tuple[int, int]],
    output_dir: Path
):
    """
    Create summary statistics comparing all methods.
    
    Parameters
    ----------
    results : dict
        Dictionary mapping method names to result GeoDataFrames
    ppt_columns : dict
        Dictionary mapping column names to (month, year) tuples
    output_dir : Path
        Output directory
    """
    summary_dir = output_dir / 'summaries'
    
    logger.info("\nCreating summary statistics...")
    
    # Get year range
    years = sorted(set(y for _, y in ppt_columns.values()))
    
    # Create annual totals comparison
    annual_summaries = []
    
    for year in years:
        year_cols = [col for col, (m, y) in ppt_columns.items() if y == year]
        
        if not year_cols:
            continue
            
        year_summary = {'year': year}
        
        for method_key, result_gdf in results.items():
            method_name = PEFF_METHODS[method_key]['name']
            
            # Get corresponding output columns for this year
            yy = year % 100
            peff_cols = [c for c in result_gdf.columns 
                        if c.startswith(f'PEFF_{method_name}') and f'_{yy:02d}_acft' in c]
            
            if peff_cols:
                # Calculate total annual Peff for all fields
                annual_total = result_gdf[peff_cols].sum().sum()
                year_summary[f'{method_name}_total_acft'] = annual_total
        
        annual_summaries.append(year_summary)
    
    # Save annual summary
    annual_df = pd.DataFrame(annual_summaries)
    annual_csv_path = summary_dir / 'annual_totals_by_method.csv'
    annual_df.to_csv(annual_csv_path, index=False)
    logger.info(f"Saved annual totals: {annual_csv_path}")
    
    # Create method comparison statistics (mean across all years and fields)
    method_stats = []
    for method_key, result_gdf in results.items():
        method_name = PEFF_METHODS[method_key]['name']
        
        # Get all Peff columns
        peff_cols = [c for c in result_gdf.columns if c.startswith(f'PEFF_{method_name}')]
        
        if peff_cols:
            all_values = result_gdf[peff_cols].values.flatten()
            all_values = all_values[~np.isnan(all_values)]
            
            method_stats.append({
                'method': method_name,
                'mean_peff_acft': np.mean(all_values),
                'std_peff_acft': np.std(all_values),
                'min_peff_acft': np.min(all_values),
                'max_peff_acft': np.max(all_values),
                'total_peff_acft': np.sum(all_values)
            })
    
    stats_df = pd.DataFrame(method_stats)
    stats_csv_path = summary_dir / 'method_comparison_statistics.csv'
    stats_df.to_csv(stats_csv_path, index=False)
    logger.info(f"Saved method comparison: {stats_csv_path}")
    
    # Print summary
    logger.info("\n" + "="*60)
    logger.info("Method Comparison Summary (all years)")
    logger.info("="*60)
    print(stats_df.to_string(index=False))


def main():
    """Main function to run the UCRB effective precipitation calculation."""
    args = parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    output_dir = Path(args.output_dir)
    gee_project = args.gee_project or GEE_PROJECT
    n_workers = args.n_workers
    
    logger.info("="*60)
    logger.info("UCRB Effective Precipitation Calculation")
    logger.info("="*60)
    logger.info(f"Using {n_workers} parallel workers")
    
    # Create output directories
    create_output_directories(output_dir)
    
    # Load GeoPackage - Precipitation data
    logger.info(f"\nLoading GeoPackage: {GPKG_PATH}")
    gdf = gpd.read_file(GPKG_PATH, layer=PPT_LAYER_NAME)
    logger.info(f"Loaded {len(gdf)} fields from '{PPT_LAYER_NAME}' layer")
    
    # Get precipitation columns
    ppt_columns = get_ppt_columns(gdf)
    
    # Load ETo data for USDA-SCS method
    eto_gdf, eto_columns = load_eto_data(GPKG_PATH)
    
    if not ppt_columns:
        logger.error("No precipitation columns found!")
        return
    
    # Download AWC data for USDA-SCS method
    # Check if existing AWC cache is available in Examples/UCRB/awc_data/
    existing_awc_cache = Path('./UCRB/awc_data/field_awc_values.csv')
    awc_cache_path = output_dir / 'awc_data' / 'field_awc_values.csv'
    
    # Use existing cache if available, otherwise use output_dir path
    if existing_awc_cache.exists():
        awc_cache_path = existing_awc_cache
        logger.info(f"Found existing AWC cache at {existing_awc_cache}")
    
    awc_series = download_awc_zonal_stats(
        gdf=gdf,
        awc_cache_path=awc_cache_path,
        gee_project=gee_project,
        skip_if_exists=True,  # Always skip if cache exists
        n_workers=n_workers
    )
    
    # Process all methods
    logger.info("\n" + "="*60)
    logger.info("Calculating Effective Precipitation for All Methods")
    logger.info("="*60)
    
    results = process_all_methods(
        gdf=gdf,
        ppt_columns=ppt_columns,
        awc_series=awc_series,
        eto_gdf=eto_gdf,
        eto_columns=eto_columns,
        gpkg_path=GPKG_PATH,
        n_workers=n_workers
    )
    
    # Create summary statistics
    create_summary_statistics(results, ppt_columns, output_dir)
    
    logger.info("\n" + "="*60)
    logger.info("Processing Complete!")
    logger.info("="*60)
    logger.info(f"Peff layers added to: {GPKG_PATH}")
    logger.info(f"Layer names: {['Peff_' + PEFF_METHODS[m]['name'] for m in PEFF_METHODS]}")
    logger.info(f"Summary statistics saved to: {output_dir / 'summaries'}")
    logger.info(f"Methods calculated: {list(PEFF_METHODS.keys())}")
    logger.info(f"Time periods: {len(ppt_columns)} months")
    logger.info(f"Fields processed: {len(gdf)}")


if __name__ == '__main__':
    main()

"""
pyCropWat South America Complete Workflow Example
===================================

This script demonstrates the complete pyCropWat workflow for effective
precipitation analysis using real data from the Rio de la Plata basin in South America.

The workflow includes:
1. Processing effective precipitation from GEE climate data
2. Temporal aggregation (annual, seasonal, climatology)
3. Statistical analysis (anomalies, trends, zonal statistics)
4. Visualization (time series, maps, comparisons)
5. Dataset comparison (ERA5-Land vs TerraClimate)
6. Effective precipitation method comparison (CROPWAT, FAO/AGLW, Fixed %, Dependable Rainfall, FarmWest, USDA-SCS)
7. Export to NetCDF format

Data Sources:
- ERA5-Land (ECMWF/ERA5_LAND/MONTHLY_AGGR)
- TerraClimate (IDAHO_EPSCOR/TERRACLIMATE)

Study Area:
- Rio de la Plata Basin (projects/ssebop-471916/assets/Riodelaplata)

Effective Precipitation Methods Compared:
- CROPWAT (USDA SCS) - Default method
- FAO/AGLW - FAO Land and Water Division formula
- Fixed Percentage (70%)
- Dependable Rainfall (75% probability)
- FarmWest - WSU irrigation scheduling formula
- USDA-SCS - Site-specific method with AWC and ETo

Requirements:
- Google Earth Engine account with access to the study area asset
- pycropwat package installed

Usage:
    python south_america_cropwat_example.py # default 4 workers
    python south_america_cropwat_example.py -f -w 32 # 32 workers, force reprocess
    python south_america_cropwat_example.py --analysis-only
    python south_america_cropwat_example.py --gee-project your-project-id
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
STUDY_AREA_ASSET = "projects/ssebop-471916/assets/Riodelaplata"

# Data sources configuration
DATASETS = {
    'ERA5Land': {
        'asset_id': 'ECMWF/ERA5_LAND/MONTHLY_AGGR',
        'precip_band': 'total_precipitation_sum',
        'scale_factor': 1000,  # Convert m to mm
        'output_dir': './RioDelaPlata/RDP_ERA5Land'
    },
    'TerraClimate': {
        'asset_id': 'IDAHO_EPSCOR/TERRACLIMATE',
        'precip_band': 'pr',
        'scale_factor': 1.0,  # Already in mm
        'output_dir': './RioDelaPlata/RDP_TerraClimate'
    }
}

# Time period
START_YEAR = 2000
END_YEAR = 2025

# Climatology period (for anomaly calculations)
CLIM_START = 2000
CLIM_END = 2020

# Region-specific output directory
REGION_DIR = Path('./RioDelaPlata')

# Analysis output directories
ANALYSIS_DIR = REGION_DIR / 'analysis_outputs'
ANNUAL_DIR = ANALYSIS_DIR / 'annual'
CLIMATOLOGY_DIR = ANALYSIS_DIR / 'climatology'
ANOMALY_DIR = ANALYSIS_DIR / 'anomalies'
TREND_DIR = ANALYSIS_DIR / 'trend'
FIGURES_DIR = ANALYSIS_DIR / 'figures'
COMPARISON_DIR = ANALYSIS_DIR / 'comparisons'
ZONAL_DIR = ANALYSIS_DIR / 'zonal_stats'
METHOD_COMPARISON_DIR = ANALYSIS_DIR / 'method_comparison'

# Zone file for zonal statistics
ZONE_FILE = ANALYSIS_DIR / 'rdp_zones.geojson'

# Effective precipitation methods to compare
PEFF_METHODS = ['cropwat', 'fao_aglw', 'fixed_percentage', 'dependable_rainfall', 'farmwest', 'usda_scs']


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
        METHOD_COMPARISON_DIR
    ]
    for d in directories:
        d.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directories created in {REGION_DIR}")


def create_sample_zones():
    """
    Create sample zone polygons for zonal statistics demonstration.
    
    Creates two zones:
    - Eastern RDP: Eastern portion of the Rio de la Plata basin (Uruguay/Brazil border)
    - Western RDP: Western portion (Argentina/Paraguay)
    """
    import json
    
    if ZONE_FILE.exists():
        logger.info(f"Zone file already exists: {ZONE_FILE}")
        return
    
    logger.info("Creating sample zone polygons...")
    
    # Define two zones within the RDP basin
    # These coordinates are approximate for the Rio de la Plata region
    zones_geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 1,
                    "name": "Eastern RDP",
                    "description": "Eastern Rio de la Plata - Uruguay/Southern Brazil"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-54.0, -35.0],
                        [-54.0, -28.0],
                        [-50.0, -28.0],
                        [-50.0, -35.0],
                        [-54.0, -35.0]
                    ]]
                }
            },
            {
                "type": "Feature",
                "properties": {
                    "zone_id": 2,
                    "name": "Western RDP",
                    "description": "Western Rio de la Plata - Argentina/Paraguay"
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [-62.0, -35.0],
                        [-62.0, -25.0],
                        [-54.0, -25.0],
                        [-54.0, -35.0],
                        [-62.0, -35.0]
                    ]]
                }
            }
        ]
    }
    
    with open(ZONE_FILE, 'w') as f:
        json.dump(zones_geojson, f, indent=2)
    
    logger.info(f"Created zone file: {ZONE_FILE}")


# =============================================================================
# Step 1: Process Effective Precipitation
# =============================================================================

def process_effective_precipitation(dataset_name: str, skip_if_exists: bool = True, n_workers: int = 4, save_inputs: bool = False):
    """
    Process effective precipitation for a given dataset.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('ERA5Land' or 'TerraClimate')
    skip_if_exists : bool
        Skip processing if output files already exist
    n_workers : int
        Number of parallel workers for processing
    save_inputs : bool
        Save downloaded input data (precipitation) to analysis_inputs folder
    """
    config = DATASETS[dataset_name]
    output_dir = Path(config['output_dir'])
    
    # Check if data already exists
    if skip_if_exists and output_dir.exists():
        existing_files = list(output_dir.glob('effective_precip_[0-9]*.tif'))
        expected_files = (END_YEAR - START_YEAR + 1) * 12
        if len(existing_files) >= expected_files * 0.9:  # Allow 10% tolerance
            logger.info(f"Skipping {dataset_name} - data already exists ({len(existing_files)} files)")
            return
    
    logger.info(f"Processing {dataset_name} effective precipitation...")
    
    ep = EffectivePrecipitation(
        asset_id=config['asset_id'],
        precip_band=config['precip_band'],
        gee_geometry_asset=STUDY_AREA_ASSET,
        start_year=START_YEAR,
        end_year=END_YEAR,
        precip_scale_factor=config['scale_factor'],
        gee_project=GEE_PROJECT
    )
    
    # Input data will be saved to region-specific analysis_inputs folder
    # e.g., ./RioDelaPlata/analysis_inputs/RDP_ERA5Land/
    input_dir = REGION_DIR / 'analysis_inputs' / output_dir.name if save_inputs else None
    
    ep.process(
        output_dir=str(output_dir),
        n_workers=n_workers,
        save_inputs=save_inputs,
        input_dir=str(input_dir) if input_dir else None
    )
    
    logger.info(f"Completed processing {dataset_name}")


# =============================================================================
# Step 2: Temporal Aggregation
# =============================================================================

def run_temporal_aggregation(dataset_name: str):
    """
    Run temporal aggregation for a dataset.
    
    Creates:
    - Annual totals
    - Growing season aggregations (April-September)
    - Monthly climatology
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('ERA5Land' or 'TerraClimate').
        
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
    
    # Seasonal aggregations (example: growing season Apr-Sep)
    logger.info("Computing growing season totals...")
    growing_season_dir = ANALYSIS_DIR / 'growing_season' / dataset_name
    growing_season_dir.mkdir(parents=True, exist_ok=True)
    
    for year in range(START_YEAR, END_YEAR + 1):
        output_path = growing_season_dir / f'growing_season_{year}.tif'
        if not output_path.exists():
            try:
                agg.custom_aggregate(
                    year=year,
                    months=[4, 5, 6, 7, 8, 9],  # Apr-Sep
                    method='sum',
                    output_path=str(output_path)
                )
            except Exception as e:
                logger.warning(f"Could not aggregate growing season {year}: {e}")
    
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
        Name of the dataset ('ERA5Land' or 'TerraClimate').
        
    Returns
    -------
    StatisticalAnalyzer
        The analyzer instance used for the dataset.
        
    Computes:
    - Anomalies for recent years
    - Long-term trends
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
    
    # Calculate anomalies for recent years (2021-2023)
    logger.info("Computing anomalies...")
    for year in range(2021, min(END_YEAR + 1, 2024)):
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
    import pandas as pd
    
    try:
        # Aggregate by zone and year
        annual_by_zone = zonal_df.groupby(['zone_id', 'year'])['mean'].mean().reset_index()
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Plot 1: Time series by zone
        for zone_id in annual_by_zone['zone_id'].unique():
            zone_data = annual_by_zone[annual_by_zone['zone_id'] == zone_id]
            zone_name = 'Eastern RDP' if zone_id == 0 else 'Western RDP'
            axes[0].plot(zone_data['year'], zone_data['mean'], 
                        marker='o', label=zone_name, linewidth=2)
        
        axes[0].set_xlabel('Year')
        axes[0].set_ylabel('Mean Effective Precipitation [mm]')
        axes[0].set_title(f'{dataset_name} - Zonal Mean Time Series')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Plot 2: Monthly climatology by zone
        monthly_by_zone = zonal_df.groupby(['zone_id', 'month'])['mean'].mean().reset_index()
        
        month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        
        width = 0.35
        x = range(1, 13)
        
        zone_ids = monthly_by_zone['zone_id'].unique()
        for i, zone_id in enumerate(zone_ids):
            zone_data = monthly_by_zone[monthly_by_zone['zone_id'] == zone_id]
            zone_name = 'Eastern RDP' if zone_id == 0 else 'Western RDP'
            offset = (i - 0.5) * width
            axes[1].bar([m + offset for m in zone_data['month']], 
                       zone_data['mean'], width=width, label=zone_name)
        
        axes[1].set_xlabel('Month')
        axes[1].set_ylabel('Mean Effective Precipitation [mm]')
        axes[1].set_title(f'{dataset_name} - Monthly Zonal Climatology')
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
    - Static raster maps (wet/dry seasons, notable events)
    - Interactive HTML maps
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('ERA5Land' or 'TerraClimate').
        
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
            title=f'{dataset_name} Effective Precipitation Time Series'
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
    
    # Static maps for multiple months (wet season and dry season)
    logger.info("Creating static maps for multiple months...")
    months_to_plot = [
        (2023, 1, 'January (Summer/Wet)'),
        (2023, 6, 'June (Winter/Dry)'),
        (2023, 10, 'October (Spring)'),
    ]
    for year, month, desc in months_to_plot:
        try:
            viz.plot_raster(
                year=year,
                month=month,
                output_path=str(dataset_fig_dir / f'map_{year}_{month:02d}.png'),
                title=f'{dataset_name} Effective Precipitation ({desc} {year})'
            )
        except Exception as e:
            logger.warning(f"Could not create static map for {year}-{month:02d}: {e}")
    
    # Interactive maps for multiple months
    logger.info("Creating interactive maps...")
    for year, month, desc in months_to_plot:
        try:
            viz.plot_interactive_map(
                year=year,
                month=month,
                output_path=str(dataset_fig_dir / f'interactive_map_{year}_{month:02d}.html'),
                title=f'{dataset_name} Effective Precipitation ({desc} {year})',
                zoom_start=5,
                opacity=0.8
            )
        except Exception as e:
            logger.warning(f"Could not create interactive map for {year}-{month:02d}: {e}")
    
    # Create maps for anomaly years (e.g., El Niño/La Niña)
    logger.info("Creating maps for notable years...")
    notable_years = [
        (2015, 12, 'El Niño Dec 2015'),
        (2022, 6, 'La Niña Jun 2022'),
    ]
    for year, month, desc in notable_years:
        try:
            viz.plot_raster(
                year=year,
                month=month,
                output_path=str(dataset_fig_dir / f'map_notable_{year}_{month:02d}.png'),
                title=f'{dataset_name} - {desc}'
            )
            viz.plot_interactive_map(
                year=year,
                month=month,
                output_path=str(dataset_fig_dir / f'interactive_notable_{year}_{month:02d}.html'),
                title=f'{dataset_name} - {desc}',
                zoom_start=5
            )
        except Exception as e:
            logger.warning(f"Could not create map for {desc}: {e}")
    
    logger.info(f"Completed visualizations for {dataset_name}")
    return viz


def plot_anomaly_maps(dataset_name: str):
    """
    Create anomaly map visualizations for a dataset using Visualizer class.
    
    Plots percent anomalies for selected months using a diverging colormap,
    including notable climate events (El Niño, La Niña).
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('ERA5Land' or 'TerraClimate').
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
    
    # Months to plot anomalies (Southern Hemisphere wet/dry seasons)
    months_to_plot = [
        (2023, 1, 'January 2023 (Summer/Wet)'),
        (2022, 6, 'June 2022 (Winter/Dry)'),
        (2015, 12, 'December 2015 (El Niño)'),
    ]
    
    for year, month, desc in months_to_plot:
        anomaly_file = dataset_anomaly_dir / f'anomaly_{year}_{month:02d}.tif'
        if not anomaly_file.exists():
            logger.debug(f"Anomaly file not found: {anomaly_file}")
            continue
        
        try:
            viz.plot_anomaly_map(
                anomaly_path=anomaly_file,
                title=f'{dataset_name} Peff Anomaly - {desc}',
                output_path=dataset_fig_dir / f'anomaly_map_{year}_{month:02d}.png'
            )
        except Exception as e:
            logger.warning(f"Could not create anomaly map for {year}-{month:02d}: {e}")


def plot_climatology_maps(dataset_name: str):
    """
    Create monthly climatology map visualizations for a dataset.
    
    Plots long-term mean effective precipitation for key months representing
    different seasons (Southern Hemisphere: Jan=summer, Jun=winter).
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('ERA5Land' or 'TerraClimate').
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
    
    # Key months representing different seasons (Southern Hemisphere)
    months_to_plot = [
        (1, 'January (Summer/Wet)'),
        (6, 'June (Winter/Dry)'),
        (10, 'October (Spring)'),
        (12, 'December (Late Spring)'),
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
                title=f'{dataset_name} Peff Climatology - {desc}',
                output_path=dataset_fig_dir / f'climatology_map_{month:02d}.png'
            )
        except Exception as e:
            logger.warning(f"Could not create climatology map for month {month:02d}: {e}")


def plot_trend_maps(dataset_name: str):
    """
    Create trend map visualizations for a dataset.
    
    Plots Sen's slope (mm/year) with Mann-Kendall significance testing:
    - Combined slope and p-value panel
    - Slope map with significance stippling overlay
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('ERA5Land' or 'TerraClimate').
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
                title=f'{dataset_name} Peff Trend ({START_YEAR}-{END_YEAR})',
                output_path=dataset_fig_dir / 'trend_maps.png'
            )
    except Exception as e:
        logger.warning(f"Could not create trend panel: {e}")
    
    # Create single trend map with significance stippling
    try:
        viz.plot_trend_map(
            slope_path=slope_file,
            pvalue_path=pvalue_file if pvalue_file.exists() else None,
            title=f'{dataset_name} Peff Trend with Significance ({START_YEAR}-{END_YEAR})',
            output_path=dataset_fig_dir / 'trend_map_with_significance.png',
            show_significance=True
        )
    except Exception as e:
        logger.warning(f"Could not create trend map: {e}")


# =============================================================================
# Step 5: Dataset Comparison
# =============================================================================

def compare_datasets():
    """
    Compare ERA5-Land and TerraClimate effective precipitation.
    
    Creates:
    - Side-by-side spatial comparison plots
    - Scatter plot showing pixel-by-pixel correlation
    - Annual mean time series comparison
    - Zonal statistics comparison
    """
    logger.info("Comparing ERA5-Land and TerraClimate datasets...")
    
    # Initialize visualizers
    era5_dir = Path(DATASETS['ERA5Land']['output_dir'])
    terra_dir = Path(DATASETS['TerraClimate']['output_dir'])
    
    if not era5_dir.exists() or not terra_dir.exists():
        logger.warning("Both datasets must be processed for comparison")
        return
    
    # Pattern to exclude fraction files
    pattern = 'effective_precip_[0-9]*.tif'
    
    # Use pattern that excludes fraction files
    viz_era5 = Visualizer(str(era5_dir), pattern=pattern)
    
    # Side-by-side comparison
    logger.info("Creating side-by-side comparison...")
    try:
        viz_era5.plot_comparison(
            year=2023,
            month=6,
            other_dir=str(terra_dir),
            other_pattern=pattern,
            labels=('ERA5-Land', 'TerraClimate'),
            title='Effective Precipitation Comparison (June 2023)',
            output_path=str(COMPARISON_DIR / 'comparison_2023_06.png')
        )
    except Exception as e:
        logger.warning(f"Could not create comparison plot: {e}")
    
    # Scatter comparison
    logger.info("Creating scatter comparison...")
    try:
        viz_era5.plot_scatter_comparison(
            start_year=START_YEAR,
            end_year=END_YEAR,
            other_dir=str(terra_dir),
            other_pattern=pattern,
            labels=('ERA5-Land', 'TerraClimate'),
            output_path=str(COMPARISON_DIR / 'scatter_comparison.png')
        )
    except Exception as e:
        logger.warning(f"Could not create scatter plot: {e}")
    
    # Annual comparison
    logger.info("Creating annual comparison...")
    try:
        viz_era5.plot_annual_comparison(
            start_year=START_YEAR,
            end_year=END_YEAR,
            other_dir=str(terra_dir),
            other_pattern=pattern,
            labels=('ERA5-Land', 'TerraClimate'),
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
    Compare zonal statistics between ERA5-Land and TerraClimate.
    
    Creates a 2x2 comparison figure showing:
    - Annual time series for Eastern RDP
    - Annual time series for Western RDP
    - Monthly climatology for Eastern RDP
    - Monthly climatology for Western RDP
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    
    era5_zonal = ZONAL_DIR / 'ERA5Land' / 'zonal_stats.csv'
    terra_zonal = ZONAL_DIR / 'TerraClimate' / 'zonal_stats.csv'
    
    if not era5_zonal.exists() or not terra_zonal.exists():
        logger.warning("Zonal statistics files not found for comparison")
        return
    
    try:
        df_era5 = pd.read_csv(era5_zonal)
        df_terra = pd.read_csv(terra_zonal)
        
        # Add dataset identifier
        df_era5['dataset'] = 'ERA5-Land'
        df_terra['dataset'] = 'TerraClimate'
        
        # Combine datasets
        df_combined = pd.concat([df_era5, df_terra], ignore_index=True)
        
        # Create comparison figure
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        zone_names = {0: 'Eastern RDP', 1: 'Western RDP'}
        colors = {'ERA5-Land': '#1f77b4', 'TerraClimate': '#ff7f0e'}
        
        # Plot 1: Annual time series comparison - Eastern RDP
        ax = axes[0, 0]
        for dataset in ['ERA5-Land', 'TerraClimate']:
            data = df_combined[(df_combined['dataset'] == dataset) & (df_combined['zone_id'] == 0)]
            annual = data.groupby('year')['mean'].mean()
            ax.plot(annual.index, annual.values, marker='o', 
                   label=dataset, color=colors[dataset], linewidth=2)
        ax.set_xlabel('Year')
        ax.set_ylabel('Mean Effective Precipitation [mm]')
        ax.set_title('Eastern RDP - Annual Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Annual time series comparison - Western RDP
        ax = axes[0, 1]
        for dataset in ['ERA5-Land', 'TerraClimate']:
            data = df_combined[(df_combined['dataset'] == dataset) & (df_combined['zone_id'] == 1)]
            annual = data.groupby('year')['mean'].mean()
            ax.plot(annual.index, annual.values, marker='o', 
                   label=dataset, color=colors[dataset], linewidth=2)
        ax.set_xlabel('Year')
        ax.set_ylabel('Mean Effective Precipitation [mm]')
        ax.set_title('Western RDP - Annual Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Monthly climatology comparison - Eastern RDP
        ax = axes[1, 0]
        month_names = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
        width = 0.35
        
        for i, dataset in enumerate(['ERA5-Land', 'TerraClimate']):
            data = df_combined[(df_combined['dataset'] == dataset) & (df_combined['zone_id'] == 0)]
            monthly = data.groupby('month')['mean'].mean()
            offset = (i - 0.5) * width
            ax.bar([m + offset for m in monthly.index], monthly.values, 
                  width=width, label=dataset, color=colors[dataset])
        ax.set_xlabel('Month')
        ax.set_ylabel('Mean Effective Precipitation [mm]')
        ax.set_title('Eastern RDP - Monthly Climatology')
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(month_names)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Plot 4: Monthly climatology comparison - Western RDP
        ax = axes[1, 1]
        for i, dataset in enumerate(['ERA5-Land', 'TerraClimate']):
            data = df_combined[(df_combined['dataset'] == dataset) & (df_combined['zone_id'] == 1)]
            monthly = data.groupby('month')['mean'].mean()
            offset = (i - 0.5) * width
            ax.bar([m + offset for m in monthly.index], monthly.values, 
                  width=width, label=dataset, color=colors[dataset])
        ax.set_xlabel('Month')
        ax.set_ylabel('Mean Effective Precipitation [mm]')
        ax.set_title('Western RDP - Monthly Climatology')
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(month_names)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.suptitle('ERA5-Land vs TerraClimate Zonal Comparison', fontsize=14, y=1.02)
        plt.tight_layout()
        
        output_path = COMPARISON_DIR / 'zonal_comparison.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved zonal comparison: {output_path}")
        
    except Exception as e:
        logger.warning(f"Could not create zonal comparison: {e}")


# =============================================================================
# Step 5b: Compare Effective Precipitation Methods
# =============================================================================

def compare_peff_methods():
    """
    Compare different effective precipitation methods for both datasets.
    
    Compares CROPWAT, FAO/AGLW, Fixed Percentage, Dependable Rainfall,
    FarmWest, and USDA-SCS methods using the same input precipitation data.
    
    Note: USDA-SCS method requires AWC (Available Water Capacity) and ETo
    (Reference Evapotranspiration). For demonstration purposes, representative
    values are used: AWC=150 mm/m (medium-textured soil), ETo=120 mm/month.
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
        'usda_scs': {'name': 'USDA-SCS', 'color': '#8c564b'}
    }
    
    # Create comparison for both datasets
    for dataset_name in ['ERA5Land', 'TerraClimate']:
        logger.info(f"\nComparing methods for {dataset_name}...")
        
        config = DATASETS[dataset_name]
        input_dir = Path(config['output_dir'])
        
        if not input_dir.exists():
            logger.warning(f"Dataset directory not found: {input_dir}")
            continue
        
        # Load a sample month for visual comparison
        sample_year = 2020
        sample_month = 1
        
        try:
            compare_methods_for_month(
                dataset_name, sample_year, sample_month, method_info
            )
        except Exception as e:
            logger.warning(f"Could not compare methods for {dataset_name}: {e}")
        
        # Create theoretical curve comparison
        try:
            plot_method_curves(dataset_name, method_info)
        except Exception as e:
            logger.warning(f"Could not create method curves: {e}")
    
    # Create cross-dataset method comparison
    try:
        compare_methods_across_datasets(method_info)
    except Exception as e:
        logger.warning(f"Could not create cross-dataset comparison: {e}")
    
    logger.info("Completed method comparison")


def compare_methods_for_month(
    dataset_name: str,
    year: int,
    month: int,
    method_info: dict
):
    """
    Compare effective precipitation methods for a specific month.
    
    Creates a 2x3 plot showing the spatial distribution of effective
    precipitation calculated by each of the six methods.
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset ('ERA5Land' or 'TerraClimate').
    year : int
        Year to analyze.
    month : int
        Month to analyze (1-12).
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
        usda_scs_effective_precip
    )
    
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    
    # Find the effective precipitation and fraction files
    peff_file = input_dir / f'effective_precip_{year}_{month:02d}.tif'
    fraction_file = input_dir / f'effective_precip_fraction_{year}_{month:02d}.tif'
    
    if not peff_file.exists():
        logger.warning(f"File not found: {peff_file}")
        return
    
    # Load the CROPWAT Peff (default)
    da_cropwat = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
    da_cropwat = da_cropwat.where(da_cropwat != 0)
    
    # Calculate total precipitation
    if fraction_file.exists():
        # Use fraction file: P = Peff / fraction
        da_fraction = rioxarray.open_rasterio(fraction_file).squeeze('band', drop=True)
        da_fraction = da_fraction.where(da_fraction > 0)  # Avoid division by zero
        precip_total = da_cropwat.values / da_fraction.values
        logger.info(f"Using fraction file to calculate total precipitation")
    else:
        # Fallback: approximate inverse of CROPWAT formula
        logger.warning(f"Fraction file not found, using approximate inverse formula")
        peff_vals = da_cropwat.values.copy()
        # For P <= 250: Peff = P * (125 - 0.2*P) / 125
        # Approximate inverse: P ≈ Peff / (1 - 0.0016 * Peff) for Peff < 125
        precip_total = np.where(
            peff_vals < 125,  # If Peff < 125, P was likely < 250
            peff_vals / (1 - 0.0016 * peff_vals),  # Approximate inverse
            (peff_vals - 125) / 0.1  # If Peff was from P > 250 formula
        )
    
    # Calculate Peff using each method
    # Note: USDA-SCS requires AWC (mm/m) and ETo (mm). Using representative values
    # for demonstration: AWC=150 mm/m (medium-textured soil), ETo=120 mm (typical monthly)
    awc_demo = np.full_like(precip_total, 150.0)  # mm/m
    eto_demo = np.full_like(precip_total, 120.0)  # mm/month
    
    peff_methods = {
        'cropwat': cropwat_effective_precip(precip_total),
        'fao_aglw': fao_aglw_effective_precip(precip_total),
        'fixed_percentage': fixed_percentage_effective_precip(precip_total, 0.7),
        'dependable_rainfall': dependable_rainfall_effective_precip(precip_total, 0.75),
        'farmwest': farmwest_effective_precip(precip_total),
        'usda_scs': usda_scs_effective_precip(precip_total, awc_demo, eto_demo)
    }
    
    # Create 2x3 comparison plot (6 methods)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    # Common color scale
    vmin = min(np.nanmin(v) for v in peff_methods.values())
    vmax = max(np.nanmax(v) for v in peff_methods.values())
    
    for idx, (method, peff_vals) in enumerate(peff_methods.items()):
        ax = axes[idx]
        
        # Create xarray for plotting
        da_method = da_cropwat.copy(data=peff_vals)
        da_method = da_method.where(~np.isnan(da_cropwat.values))
        
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
    
    plt.suptitle(
        f'{dataset_name} - Method Comparison ({year}/{month:02d})',
        fontsize=14, fontweight='bold'
    )
    
    # Add shared colorbar outside the figure on the right
    plt.tight_layout()
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    fig.colorbar(im, cax=cbar_ax, label='Effective Precipitation [mm]')
    
    output_path = METHOD_COMPARISON_DIR / f'{dataset_name}_method_maps_{year}_{month:02d}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method comparison maps: {output_path}")


def plot_method_curves(dataset_name: str, method_info: dict):
    """
    Plot theoretical effective precipitation curves for each method.
    
    Creates a two-panel figure showing:
    - Peff vs P curves for all methods
    - Efficiency ratio (Peff/P) vs P for all methods
    
    Parameters
    ----------
    dataset_name : str
        Name of the dataset (used for plot titles).
    method_info : dict
        Dictionary with method names as keys and sub-dicts containing
        'name' (display name) and 'color' for plotting.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from pycropwat.methods import (
        cropwat_effective_precip,
        fao_aglw_effective_precip,
        fixed_percentage_effective_precip,
        dependable_rainfall_effective_precip,
        farmwest_effective_precip,
        usda_scs_effective_precip
    )
    
    # Create precipitation range
    precip = np.linspace(0, 400, 1000)
    
    # For USDA-SCS: use representative values for demonstration
    # AWC=150 mm/m (medium-textured soil), ETo=120 mm (typical monthly)
    awc_demo = np.full_like(precip, 150.0)  # mm/m
    eto_demo = np.full_like(precip, 120.0)  # mm/month
    
    # Calculate Peff for each method
    methods = {
        'cropwat': cropwat_effective_precip(precip),
        'fao_aglw': fao_aglw_effective_precip(precip),
        'fixed_percentage': fixed_percentage_effective_precip(precip, 0.7),
        'dependable_rainfall': dependable_rainfall_effective_precip(precip, 0.75),
        'farmwest': farmwest_effective_precip(precip),
        'usda_scs': usda_scs_effective_precip(precip, awc_demo, eto_demo)
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
    ax.set_title('Peff Methods Comparison')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 400)
    ax.set_ylim(0, 300)
    
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
    ax.set_xlim(0, 400)
    ax.set_ylim(0, 1.1)
    
    plt.suptitle(f'{dataset_name} - Theoretical Method Curves', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = METHOD_COMPARISON_DIR / f'{dataset_name}_method_curves.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method curves: {output_path}")


def compare_methods_across_datasets(method_info: dict):
    """
    Create summary comparison of methods across both datasets.
    
    Samples multiple months from both ERA5-Land and TerraClimate to compute
    statistics showing how method choice affects results for each dataset.
    
    Parameters
    ----------
    method_info : dict
        Dictionary with method names as keys and sub-dicts containing
        'name' (display name) and 'color' for plotting.
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
        usda_scs_effective_precip
    )
    
    results = []
    
    for dataset_name in ['ERA5Land', 'TerraClimate']:
        config = DATASETS[dataset_name]
        input_dir = Path(config['output_dir'])
        
        if not input_dir.exists():
            continue
        
        # Sample a few months to get statistics
        sample_years = [2010, 2015, 2020]
        sample_months = [1, 4, 7, 10]  # Jan, Apr, Jul, Oct
        
        for year in sample_years:
            for month in sample_months:
                peff_file = input_dir / f'effective_precip_{year}_{month:02d}.tif'
                
                if not peff_file.exists():
                    continue
                
                try:
                    da = rioxarray.open_rasterio(peff_file).squeeze('band', drop=True)
                    da = da.where(da != 0)
                    peff_vals = da.values[~np.isnan(da.values)]
                    
                    # Approximate precipitation
                    precip_approx = np.where(
                        peff_vals < 125,
                        peff_vals / (1 - 0.0016 * peff_vals),
                        (peff_vals - 125) / 0.1
                    )
                    precip_approx = np.clip(precip_approx, 0, 500)
                    
                    # For USDA-SCS: use representative values
                    awc_demo = np.full_like(precip_approx, 150.0)  # mm/m
                    eto_demo = np.full_like(precip_approx, 120.0)  # mm/month
                    
                    # Calculate mean Peff for each method
                    for method_name, func in [
                        ('cropwat', lambda p: cropwat_effective_precip(p)),
                        ('fao_aglw', lambda p: fao_aglw_effective_precip(p)),
                        ('fixed_percentage', lambda p: fixed_percentage_effective_precip(p, 0.7)),
                        ('dependable_rainfall', lambda p: dependable_rainfall_effective_precip(p, 0.75)),
                        ('farmwest', lambda p: farmwest_effective_precip(p)),
                        ('usda_scs', lambda p: usda_scs_effective_precip(p, awc_demo, eto_demo))
                    ]:
                        peff_method = func(precip_approx)
                        results.append({
                            'dataset': dataset_name,
                            'year': year,
                            'month': month,
                            'method': method_name,
                            'mean_peff': np.nanmean(peff_method),
                            'mean_precip': np.nanmean(precip_approx)
                        })
                except Exception as e:
                    logger.debug(f"Error processing {peff_file}: {e}")
    
    if not results:
        logger.warning("No data available for cross-dataset comparison")
        return
    
    df = pd.DataFrame(results)
    
    # Create comparison figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Box plot of Peff by method - ERA5-Land
    ax = axes[0, 0]
    era5_data = df[df['dataset'] == 'ERA5Land']
    if not era5_data.empty:
        box_data = [era5_data[era5_data['method'] == m]['mean_peff'].values 
                   for m in PEFF_METHODS]
        bp = ax.boxplot(box_data, labels=[method_info[m]['name'] for m in PEFF_METHODS],
                       patch_artist=True)
        for patch, method in zip(bp['boxes'], PEFF_METHODS):
            patch.set_facecolor(method_info[method]['color'])
            patch.set_alpha(0.7)
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('ERA5-Land - Method Distribution')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Box plot of Peff by method - TerraClimate
    ax = axes[0, 1]
    terra_data = df[df['dataset'] == 'TerraClimate']
    if not terra_data.empty:
        box_data = [terra_data[terra_data['method'] == m]['mean_peff'].values 
                   for m in PEFF_METHODS]
        bp = ax.boxplot(box_data, labels=[method_info[m]['name'] for m in PEFF_METHODS],
                       patch_artist=True)
        for patch, method in zip(bp['boxes'], PEFF_METHODS):
            patch.set_facecolor(method_info[method]['color'])
            patch.set_alpha(0.7)
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('TerraClimate - Method Distribution')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Bar chart comparing methods
    ax = axes[1, 0]
    method_means = df.groupby(['dataset', 'method'])['mean_peff'].mean().unstack()
    x = np.arange(len(PEFF_METHODS))
    width = 0.35
    
    if 'ERA5Land' in method_means.index:
        bars1 = ax.bar(x - width/2, 
                      [method_means.loc['ERA5Land', m] for m in PEFF_METHODS],
                      width, label='ERA5-Land', color='steelblue', alpha=0.8)
    if 'TerraClimate' in method_means.index:
        bars2 = ax.bar(x + width/2,
                      [method_means.loc['TerraClimate', m] for m in PEFF_METHODS],
                      width, label='TerraClimate', color='coral', alpha=0.8)
    
    ax.set_ylabel('Mean Effective Precipitation [mm]')
    ax.set_title('Dataset × Method Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels([method_info[m]['name'] for m in PEFF_METHODS], rotation=15)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Efficiency by method
    ax = axes[1, 1]
    df['efficiency'] = df['mean_peff'] / df['mean_precip']
    eff_means = df.groupby(['dataset', 'method'])['efficiency'].mean().unstack()
    
    if 'ERA5Land' in eff_means.index:
        ax.bar(x - width/2,
              [eff_means.loc['ERA5Land', m] for m in PEFF_METHODS],
              width, label='ERA5-Land', color='steelblue', alpha=0.8)
    if 'TerraClimate' in eff_means.index:
        ax.bar(x + width/2,
              [eff_means.loc['TerraClimate', m] for m in PEFF_METHODS],
              width, label='TerraClimate', color='coral', alpha=0.8)
    
    ax.set_ylabel('Mean Efficiency (Peff/P)')
    ax.set_title('Rainfall Effectiveness by Method')
    ax.set_xticks(x)
    ax.set_xticklabels([method_info[m]['name'] for m in PEFF_METHODS], rotation=15)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1)
    
    plt.suptitle('Effective Precipitation Method Comparison', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = METHOD_COMPARISON_DIR / 'method_comparison_summary.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved method comparison summary: {output_path}")
    
    # Save statistics to CSV
    summary_stats = df.groupby(['dataset', 'method']).agg({
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
        Name of the dataset ('ERA5Land' or 'TerraClimate').
    """
    config = DATASETS[dataset_name]
    input_dir = Path(config['output_dir'])
    output_file = ANALYSIS_DIR / f'{dataset_name}_effective_precip_{START_YEAR}_{END_YEAR}.nc'
    
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
    Run the complete pyCropWat workflow.
    
    Parameters
    ----------
    skip_processing : bool
        If True, skip GEE processing if data already exists.
        Set to False to force reprocessing.
    n_workers : int
        Number of parallel workers for GEE processing.
    """
    logger.info("=" * 60)
    logger.info("pyCropWat Complete Workflow")
    logger.info("Study Area: Rio de la Plata Basin")
    logger.info(f"Period: {START_YEAR}-{END_YEAR}")
    logger.info("=" * 60)
    
    # Create output directories
    create_output_directories()
    
    # Create sample zones for zonal statistics
    create_sample_zones()
    
    # Process each dataset
    for dataset_name in ['ERA5Land', 'TerraClimate']:
        logger.info(f"\n{'=' * 40}")
        logger.info(f"Processing {dataset_name}")
        logger.info("=" * 40)
        
        # Step 1: Process effective precipitation
        process_effective_precipitation(dataset_name, skip_if_exists=skip_processing, n_workers=n_workers)
        
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
        
        # Step 6: Export to NetCDF
        export_netcdf(dataset_name)
    
    # Step 5: Dataset comparison (includes zonal comparison)
    logger.info(f"\n{'=' * 40}")
    logger.info("Dataset Comparison")
    logger.info("=" * 40)
    compare_datasets()
    
    # Step 5b: Method comparison
    logger.info(f"\n{'=' * 40}")
    logger.info("Effective Precipitation Method Comparison")
    logger.info("=" * 40)
    compare_peff_methods()
    
    logger.info("\n" + "=" * 60)
    logger.info("Workflow Complete!")
    logger.info("=" * 60)
    logger.info(f"Analysis outputs saved to: {ANALYSIS_DIR.absolute()}")


def run_analysis_only():
    """
    Run analysis workflow only (without GEE processing).
    
    Use this if you already have processed effective precipitation data
    in the RDP_ERA5Land and RDP_TerraClimate directories.
    """
    logger.info("Running analysis-only workflow...")
    
    create_output_directories()
    
    # Create sample zones for zonal statistics
    create_sample_zones()
    
    for dataset_name in ['ERA5Land', 'TerraClimate']:
        config = DATASETS[dataset_name]
        input_dir = Path(config['output_dir'])
        
        if not input_dir.exists():
            logger.warning(f"Input directory not found: {input_dir}")
            logger.warning(f"Run full workflow or process {dataset_name} data first")
            continue
        
        logger.info(f"\nAnalyzing {dataset_name}...")
        
        agg = run_temporal_aggregation(dataset_name)
        stats = run_statistical_analysis(dataset_name)
        viz = run_visualization(dataset_name)
        
        # Anomaly, climatology, and trend maps
        plot_anomaly_maps(dataset_name)
        plot_climatology_maps(dataset_name)
        plot_trend_maps(dataset_name)
        
        export_netcdf(dataset_name)
    
    compare_datasets()
    compare_peff_methods()
    
    logger.info("\nAnalysis workflow complete!")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='pyCropWat Complete Workflow Example',
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

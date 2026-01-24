"""
pyCropWat Western U.S. PCML Example
====================================

This script demonstrates the Physics-Constrained Machine Learning (PCML) 
effective precipitation method for the Western United States.

The PCML method provides pre-computed effective precipitation values from a 
trained ML model that combines satellite data with physics-based constraints.
This example downloads PCML data and calculates water year totals.

PCML Coverage:
- **Region**: Western U.S. (17 states: AZ, CA, CO, ID, KS, MT, NE, NV, NM, 
              ND, OK, OR, SD, TX, UT, WA, WY)
- **Temporal**: January 2000 - September 2024 (monthly)
- **Resolution**: ~2 km (native scale retrieved dynamically from GEE asset)
- **GEE Asset**: projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked

PCML Geometry Options:
- **No geometry provided**: Downloads the entire PCML asset (full Western U.S. - 17 states)
- **User provides geometry**: PCML data is clipped/subsetted to that geometry. 
  Note: Only Western U.S. vectors overlapping the 17-state extent can be used
  (e.g., AZ.geojson, pacific_northwest.geojson)

Water Year Definition:
- Western U.S. water year runs from October to September
- Water Year 2001 = October 2000 - September 2001
- Water Year 2024 = October 2023 - September 2024
- PCML annual fractions are water year (Oct-Sep) fractions

Output:
- Monthly effective precipitation rasters (PCML Peff)
- Annual (water year, Oct-Sep) effective precipitation fraction rasters (from GEE asset, WY 2000-2024)
- Water year total effective precipitation (WY 2001-2024, since Peff starts Jan 2000)
- Mean annual fraction maps (WY 2000-2024)

Reference:
Hasan, M. F., Smith, R. G., Majumdar, S., Huntington, J. L., Alves Meira Neto, A., 
& Minor, B. A. (2025). Satellite data and physics-constrained machine learning for 
estimating effective precipitation in the Western United States and application for 
monitoring groundwater irrigation. Agricultural Water Management, 319, 109821.
https://doi.org/10.1016/j.agwat.2025.109821

Requirements:
- Google Earth Engine account with access to the Western U.S. geometry asset
- pycropwat package installed

Usage:
    python western_us_pcml_example.py
    python western_us_pcml_example.py -w 8  # Use 8 workers
    python western_us_pcml_example.py --gee-project your-project-id
    python western_us_pcml_example.py --analysis-only  # Skip download, run analysis only
"""

import argparse
import logging
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import rioxarray

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import pyCropWat modules
from pycropwat import EffectivePrecipitation


# =============================================================================
# Configuration
# =============================================================================

# GEE Configuration
GEE_PROJECT = None  # Set to your GEE project ID, or None to use default

# Time period (PCML Peff available: Jan 2000 - Sep 2024)
# Note: Water year totals can be calculated for WY 2001-2024
START_YEAR = 2000
END_YEAR = 2024

# Output directories
OUTPUT_DIR = Path('./WesternUS_PCML')
FIGURES_DIR = OUTPUT_DIR / 'figures'
WATER_YEAR_DIR = OUTPUT_DIR / 'water_year'
FRACTION_DIR = OUTPUT_DIR / 'fractions'


def create_output_directories():
    """Create all necessary output directories."""
    directories = [OUTPUT_DIR, FIGURES_DIR, WATER_YEAR_DIR, FRACTION_DIR]
    for d in directories:
        d.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directories created in {OUTPUT_DIR}")


# =============================================================================
# Step 1: Download PCML Data
# =============================================================================

def download_pcml_data(skip_if_exists: bool = True, n_workers: int = 4):
    """
    Download PCML effective precipitation data for Western U.S.
    
    PCML method automatically:
    - Uses the default PCML asset (covers entire Western U.S.)
    - Retrieves native scale (~2km) from the asset using nominalScale()
    - Uses the PCML asset's geometry (no external geometry needed)
    - Downloads annual fractions directly from GEE asset (WY 2000-2024)
    
    Parameters
    ----------
    skip_if_exists : bool
        Skip download if files already exist
    n_workers : int
        Number of parallel workers
    """
    # Check if data already exists
    if skip_if_exists and OUTPUT_DIR.exists():
        existing_files = list(OUTPUT_DIR.glob('effective_precip_[0-9]*.tif'))
        # PCML Peff: Jan 2000 - Sep 2024 = 24*12 + 9 = 297 months
        expected_months = (END_YEAR - START_YEAR) * 12 + 9
        if len(existing_files) >= expected_months * 0.9:
            logger.info(f"Skipping download - data already exists ({len(existing_files)} files)")
            return
    
    logger.info("=" * 60)
    logger.info("Downloading PCML Effective Precipitation Data")
    logger.info("=" * 60)
    logger.info(f"Period: {START_YEAR} - {END_YEAR}")
    logger.info("Method: PCML (Physics-Constrained Machine Learning)")
    logger.info("Region: Western United States (17 states)")
    logger.info("Note: PCML asset covers entire Western U.S. - no geometry needed")
    
    # Initialize PCML processor
    # For PCML method:
    # - Asset and band are automatically set to PCML defaults when method='pcml'
    # - Geometry is derived from the PCML asset itself (no external geometry needed)
    # Note: PCML Peff asset ends at Sep 2024; download will skip non-existent bands
    ep = EffectivePrecipitation(
        asset_id=None,  # Auto-set to PCML default asset
        precip_band=None,  # Auto-set for PCML
        geometry_path=None,  # PCML uses its own geometry (entire Western U.S.)
        start_year=START_YEAR,
        end_year=END_YEAR,
        gee_project=GEE_PROJECT,
        method='pcml'  # This triggers PCML-specific handling
    )
    
    # Process and download
    logger.info(f"Processing with {n_workers} workers...")
    results = ep.process(
        output_dir=str(OUTPUT_DIR),
        n_workers=n_workers,
        save_inputs=False,  # PCML doesn't need gridMET precip (Peff comes from asset)
        input_dir=str(OUTPUT_DIR / 'inputs')
    )
    
    successful = sum(1 for r in results if r[0] is not None)
    logger.info(f"Download complete: {successful}/{len(results)} months processed")


# =============================================================================
# Step 2: Calculate Water Year Totals
# =============================================================================

def calculate_water_year_totals():
    """
    Calculate water year total effective precipitation.
    
    Water year in Western U.S. runs from October to September.
    - Water Year 2001 = Oct 2000 - Sep 2001
    - Water Year 2024 = Oct 2023 - Sep 2024
    
    This function aggregates monthly PCML Peff to water year totals.
    """
    logger.info("=" * 60)
    logger.info("Calculating Water Year Totals (Oct-Sep)")
    logger.info("=" * 60)
    
    WATER_YEAR_DIR.mkdir(parents=True, exist_ok=True)
    
    # Water years available: WY 2001 (Oct 2000 - Sep 2001) to WY 2024 (Oct 2023 - Sep 2024)
    water_years = range(2001, 2025)  # WY 2001 to WY 2024
    
    wy_totals = {}
    
    for wy in water_years:
        logger.info(f"Processing Water Year {wy} (Oct {wy-1} - Sep {wy})...")
        
        # Collect monthly files for this water year
        monthly_data = []
        template_da = None
        
        # October, November, December of previous year
        for month in [10, 11, 12]:
            month_file = OUTPUT_DIR / f'effective_precip_{wy-1}_{month:02d}.tif'
            if month_file.exists():
                da = rioxarray.open_rasterio(month_file).squeeze('band', drop=True)
                monthly_data.append(da.values)
                if template_da is None:
                    template_da = da
        
        # January to September of water year
        for month in range(1, 10):  # Jan-Sep
            month_file = OUTPUT_DIR / f'effective_precip_{wy}_{month:02d}.tif'
            if month_file.exists():
                da = rioxarray.open_rasterio(month_file).squeeze('band', drop=True)
                monthly_data.append(da.values)
                if template_da is None:
                    template_da = da
        
        if len(monthly_data) < 10:  # Need at least 10 months
            logger.warning(f"  Insufficient data for WY {wy} ({len(monthly_data)} months)")
            continue
        
        # Sum to water year total
        wy_total = np.nansum(monthly_data, axis=0)
        wy_totals[wy] = np.nanmean(wy_total)
        
        # Save water year total raster
        output_file = WATER_YEAR_DIR / f'pcml_peff_wy{wy}.tif'
        
        da_wy = template_da.copy(data=wy_total)
        da_wy.attrs = {
            'units': 'mm',
            'long_name': f'PCML effective precipitation water year {wy}',
            'water_year': wy,
            'period': f'Oct {wy-1} - Sep {wy}',
            'method': 'PCML'
        }
        da_wy = da_wy.rio.write_crs("EPSG:4326")
        da_wy.rio.to_raster(output_file)
        
        logger.info(f"  WY {wy}: Mean Peff = {wy_totals[wy]:.1f} mm")
    
    logger.info(f"Saved {len(wy_totals)} water year totals to {WATER_YEAR_DIR}")
    return wy_totals


# =============================================================================
# Step 3: Calculate Mean Fractions
# =============================================================================

def calculate_mean_fractions():
    """
    Calculate mean effective precipitation fractions.
    
    Creates:
    - Long-term mean annual (water year, Oct-Sep) fraction (average of WY 2000-2024)
    - Annual fraction time series
    
    Note: PCML annual fractions are water year (Oct-Sep) fractions, loaded
    directly from the PCML annual fraction GEE asset. Band format: bYYYY,
    where bYYYY represents WY YYYY (Oct YYYY-1 to Sep YYYY).
    """
    logger.info("=" * 60)
    logger.info("Calculating Mean Effective Precipitation Fractions (Water Year Oct-Sep)")
    logger.info("=" * 60)
    
    FRACTION_DIR.mkdir(parents=True, exist_ok=True)
    
    # Collect PCML annual fraction files (format: effective_precip_fraction_YYYY.tif)
    # Each file contains the water year fraction (band bYYYY = WY YYYY)
    fraction_files = sorted(OUTPUT_DIR.glob('effective_precip_fraction_[0-9][0-9][0-9][0-9].tif'))
    
    if not fraction_files:
        logger.warning("No fraction files found. Run download first.")
        return None
    
    logger.info(f"Found {len(fraction_files)} annual fraction files")
    
    # Calculate overall mean fraction and collect annual values
    all_fractions = []
    annual_fractions = {}  # water year -> mean fraction value
    template_da = None
    
    for frac_file in fraction_files:
        da = rioxarray.open_rasterio(frac_file).squeeze('band', drop=True)
        values = da.values
        
        # Extract year from filename (format: effective_precip_fraction_YYYY.tif)
        # The year in the filename corresponds to the water year (band bYYYY = WY YYYY)
        parts = frac_file.stem.split('_')
        water_year = int(parts[-1])  # Last part is the year (no month suffix for PCML)
        
        annual_fractions[water_year] = np.nanmean(values[values > 0])  # Exclude nodata
        all_fractions.append(values)
        
        if template_da is None:
            template_da = da
    
    # Determine water year range
    min_wy = min(annual_fractions.keys())
    max_wy = max(annual_fractions.keys())
    
    # Calculate and save overall mean fraction
    mean_fraction = np.nanmean(all_fractions, axis=0)
    
    output_file = FRACTION_DIR / f'pcml_mean_fraction_wy{min_wy}_{max_wy}.tif'
    da_mean = template_da.copy(data=mean_fraction)
    da_mean.attrs = {
        'units': 'fraction',
        'long_name': 'Mean PCML effective precipitation fraction (water year, Oct-Sep)',
        'period': f'WY {min_wy}-{max_wy}',
        'source': 'PCML annual fraction GEE asset (water year, Oct-Sep)'
    }
    da_mean = da_mean.rio.write_crs("EPSG:4326")
    da_mean = da_mean.rio.write_nodata(0)
    da_mean.rio.to_raster(output_file)
    
    logger.info(f"Overall mean fraction (WY {min_wy}-{max_wy}): {np.nanmean(mean_fraction[mean_fraction > 0]):.3f}")
    logger.info(f"Saved: {output_file}")
    
    # Log annual fraction statistics
    logger.info(f"\nAnnual (Water Year) Fraction Statistics (WY {min_wy}-{max_wy}):")
    for year in sorted(annual_fractions.keys()):
        logger.info(f"  WY {year}: {annual_fractions[year]:.3f}")
    
    logger.info(f"\nMean of annual fractions: {np.mean(list(annual_fractions.values())):.3f}")
    
    return mean_fraction, annual_fractions


# =============================================================================
# Step 4: Visualization
# =============================================================================

def plot_water_year_time_series(wy_totals: dict):
    """
    Plot water year total effective precipitation time series.
    
    Parameters
    ----------
    wy_totals : dict
        Dictionary with water year as key and mean Peff as value
    """
    logger.info("Creating water year time series plot...")
    
    years = sorted(wy_totals.keys())
    values = [wy_totals[y] for y in years]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot bars with color gradient based on value
    colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(years)))
    bars = ax.bar(years, values, color=colors, edgecolor='navy', linewidth=0.5)
    
    # Add trend line
    z = np.polyfit(years, values, 1)
    p = np.poly1d(z)
    ax.plot(years, p(years), 'r--', linewidth=2, label=f'Trend: {z[0]:.2f} mm/year')
    
    # Add mean line
    mean_val = np.mean(values)
    ax.axhline(y=mean_val, color='green', linestyle=':', linewidth=2, 
               label=f'Mean: {mean_val:.1f} mm')
    
    ax.set_xlabel('Water Year', fontsize=12)
    ax.set_ylabel('Total Effective Precipitation [mm]', fontsize=12)
    ax.set_title('PCML Water Year Total Effective Precipitation (Oct-Sep)\nWestern United States WY 2001-2024', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Rotate x-axis labels
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    output_path = FIGURES_DIR / 'water_year_time_series.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved: {output_path}")


def plot_annual_fraction_time_series(annual_fractions: dict):
    """
    Plot annual (water year) effective precipitation fraction time series.
    
    Parameters
    ----------
    annual_fractions : dict
        Dictionary with water year as key and mean fraction as value
    """
    logger.info("Creating annual fraction time series plot...")
    
    years = sorted(annual_fractions.keys())
    values = [annual_fractions[y] for y in years]
    min_wy = min(years)
    max_wy = max(years)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot bars with color gradient based on value
    colors = plt.cm.YlGnBu(np.linspace(0.3, 0.9, len(years)))
    bars = ax.bar(years, values, color=colors, edgecolor='darkblue', linewidth=0.5)
    
    # Add trend line
    z = np.polyfit(years, values, 1)
    p = np.poly1d(z)
    ax.plot(years, p(years), 'r--', linewidth=2, label=f'Trend: {z[0]*100:.4f}% /year')
    
    # Add mean line
    mean_val = np.mean(values)
    ax.axhline(y=mean_val, color='green', linestyle=':', linewidth=2, 
               label=f'Mean: {mean_val:.3f}')
    
    ax.set_xlabel('Water Year', fontsize=12)
    ax.set_ylabel('Effective Precipitation Fraction', fontsize=12)
    ax.set_title(f'PCML Annual Peff Fraction Time Series (Water Year Oct-Sep)\nWestern United States WY {min_wy}-{max_wy}', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Rotate x-axis labels
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    output_path = FIGURES_DIR / 'annual_fraction_time_series.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved: {output_path}")


def plot_mean_fraction_map():
    """
    Plot spatial map of mean effective precipitation fraction.
    """
    logger.info("Creating mean fraction map...")
    
    # Find the mean fraction file (dynamically named based on water year range)
    frac_files = list(FRACTION_DIR.glob('pcml_mean_fraction_wy*.tif'))
    if not frac_files:
        logger.warning("Mean fraction file not found")
        return
    frac_file = frac_files[0]  # Use first match
    
    da = rioxarray.open_rasterio(frac_file).squeeze('band', drop=True)
    # Mask out nodata (0) values so they don't appear in the plot
    da = da.where(da != 0)
    
    # Calculate 2-98% percentile for color scale
    vmin = np.nanpercentile(da.values, 2)
    vmax = np.nanpercentile(da.values, 98)
    
    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': None})
    
    # Plot the fraction map
    im = da.plot(ax=ax, cmap='YlGnBu', vmin=vmin, vmax=vmax,
                 cbar_kwargs={'label': 'Effective Precipitation Fraction',
                             'shrink': 0.8})
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(f'PCML Mean Effective Precipitation Fraction\nWestern U.S. ({START_YEAR}-{END_YEAR})',
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    output_path = FIGURES_DIR / 'mean_fraction_map.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved: {output_path}")


def plot_mean_annual_water_year_map():
    """
    Plot spatial map of mean annual water year total across all years.
    """
    logger.info("Creating mean annual water year map...")
    
    # Find all water year files
    wy_files = sorted(WATER_YEAR_DIR.glob('pcml_peff_wy*.tif'))
    
    if not wy_files:
        logger.warning("No water year files found")
        return None
    
    # Stack all water years and compute mean
    arrays = []
    years = []
    for wy_file in wy_files:
        da = rioxarray.open_rasterio(wy_file).squeeze('band', drop=True)
        arrays.append(da.values)
        # Extract year from filename
        year = int(wy_file.stem.replace('pcml_peff_wy', ''))
        years.append(year)
    
    # Compute mean across all years
    stacked = np.stack(arrays, axis=0)
    mean_arr = np.nanmean(stacked, axis=0)
    
    # Create DataArray with same coordinates as input
    da_mean = da.copy(data=mean_arr)
    
    # Save the mean annual map
    mean_output = WATER_YEAR_DIR / f'pcml_mean_annual_peff_wy{min(years)}-{max(years)}.tif'
    da_mean.rio.to_raster(mean_output)
    logger.info(f"Saved mean annual Peff: {mean_output}")
    
    # Calculate 2-98% percentile for color scale
    vmin = np.nanpercentile(mean_arr, 2)
    vmax = np.nanpercentile(mean_arr, 98)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot the mean Peff map
    im = da_mean.plot(ax=ax, cmap='Blues', vmin=vmin, vmax=vmax,
                      cbar_kwargs={'label': 'Effective Precipitation [mm]',
                                  'shrink': 0.8})
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(f'PCML Mean Annual Effective Precipitation\n'
                 f'Water Years {min(years)}-{max(years)} (n={len(years)} years)',
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    output_path = FIGURES_DIR / f'mean_annual_water_year_map_wy{min(years)}-{max(years)}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved: {output_path}")
    
    return da_mean, min(years), max(years)


def create_summary_figure(wy_totals: dict, annual_fractions: dict):
    """
    Create a summary figure with multiple panels.
    
    Parameters
    ----------
    wy_totals : dict
        Water year totals
    annual_fractions : dict
        Annual (water year) mean fractions
    """
    logger.info("Creating summary figure...")
    
    fig = plt.figure(figsize=(16, 12))
    
    # Panel 1: Water year time series (top left)
    ax1 = fig.add_subplot(2, 2, 1)
    years = sorted(wy_totals.keys())
    values = [wy_totals[y] for y in years]
    ax1.bar(years, values, color='steelblue', edgecolor='navy', linewidth=0.5)
    ax1.axhline(y=np.mean(values), color='red', linestyle='--', 
                label=f'Mean: {np.mean(values):.1f} mm')
    ax1.set_xlabel('Water Year')
    ax1.set_ylabel('Total Peff [mm]')
    ax1.set_title('(a) Water Year Total Effective Precipitation (Oct-Sep)')
    ax1.legend()
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Panel 2: Annual fraction time series (top right)
    ax2 = fig.add_subplot(2, 2, 2)
    frac_years = sorted(annual_fractions.keys())
    frac_values = [annual_fractions[y] for y in frac_years]
    ax2.bar(frac_years, frac_values, color='teal', edgecolor='darkslategray', linewidth=0.5)
    ax2.axhline(y=np.mean(frac_values), color='red', linestyle='--',
                label=f'Mean: {np.mean(frac_values):.3f}')
    ax2.set_xlabel('Water Year')
    ax2.set_ylabel('Mean Fraction')
    ax2.set_title('(b) Annual Mean Peff Fraction (Water Year Oct-Sep)')
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_ylim(0, 1)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Panel 3: Mean fraction map (bottom left)
    ax3 = fig.add_subplot(2, 2, 3)
    frac_files = list(FRACTION_DIR.glob('pcml_mean_fraction_wy*.tif'))
    if frac_files:
        frac_file = frac_files[0]
        da_frac = rioxarray.open_rasterio(frac_file).squeeze('band', drop=True)
        # Mask out nodata (0) values so they don't appear in the plot
        da_frac = da_frac.where(da_frac > 0)
        # Calculate 2-98% percentile for color scale
        vmin3 = np.nanpercentile(da_frac.values, 2)
        vmax3 = np.nanpercentile(da_frac.values, 98)
        da_frac.plot(ax=ax3, cmap='YlGnBu', vmin=vmin3, vmax=vmax3,
                     cbar_kwargs={'label': 'Fraction', 'shrink': 0.8})
        ax3.set_xlabel('Longitude (°)')
        ax3.set_ylabel('Latitude (°)')
        # Extract WY range from filename
        wy_range = frac_file.stem.replace('pcml_mean_fraction_', '').upper()
        ax3.set_title(f'(c) Mean Effective Precipitation Fraction ({wy_range})')
    
    # Panel 4: Mean annual water year map (bottom right)
    ax4 = fig.add_subplot(2, 2, 4)
    # Find all water year files and compute mean
    wy_files = sorted(WATER_YEAR_DIR.glob('pcml_peff_wy*.tif'))
    if wy_files:
        arrays = []
        years = []
        for wy_file in wy_files:
            # Skip mean annual file if it exists
            if 'mean_annual' in wy_file.stem:
                continue
            da_wy = rioxarray.open_rasterio(wy_file).squeeze('band', drop=True)
            arrays.append(da_wy.values)
            year = int(wy_file.stem.replace('pcml_peff_wy', ''))
            years.append(year)
        
        if arrays:
            mean_arr = np.nanmean(np.stack(arrays, axis=0), axis=0)
            da_mean = da_wy.copy(data=mean_arr)
            # Calculate 2-98% percentile for color scale
            vmin4 = np.nanpercentile(mean_arr, 2)
            vmax4 = np.nanpercentile(mean_arr, 98)
            da_mean.plot(ax=ax4, cmap='Blues', vmin=vmin4, vmax=vmax4,
                         cbar_kwargs={'label': 'Peff [mm]', 'shrink': 0.8})
            ax4.set_xlabel('Longitude (°)')
            ax4.set_ylabel('Latitude (°)')
            ax4.set_title(f'(d) Mean Annual Peff (WY {min(years)}-{max(years)})')
    
    plt.suptitle('PCML Effective Precipitation Analysis - Western United States\n'
                 f'Period: {START_YEAR}-{END_YEAR}', 
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    output_path = FIGURES_DIR / 'pcml_summary.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved: {output_path}")


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Western U.S. PCML Effective Precipitation Example'
    )
    parser.add_argument('-w', '--workers', type=int, default=4,
                        help='Number of parallel workers (default: 4)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Force reprocessing even if files exist')
    parser.add_argument('--analysis-only', action='store_true',
                        help='Skip download, run analysis only')
    parser.add_argument('--gee-project', type=str, default=None,
                        help='Google Earth Engine project ID')
    
    args = parser.parse_args()
    
    # Update global config if GEE project provided
    global GEE_PROJECT
    if args.gee_project:
        GEE_PROJECT = args.gee_project
    
    logger.info("=" * 70)
    logger.info("PCML Effective Precipitation Analysis - Western United States")
    logger.info("=" * 70)
    logger.info(f"Period: {START_YEAR} - {END_YEAR}")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    
    # Create directories
    create_output_directories()
    
    # Step 1: Download PCML data
    if not args.analysis_only:
        download_pcml_data(
            skip_if_exists=not args.force,
            n_workers=args.workers
        )
    
    # Step 2: Calculate water year totals
    wy_totals = calculate_water_year_totals()
    
    # Step 3: Calculate mean fractions
    result = calculate_mean_fractions()
    if result:
        mean_fraction, annual_fractions = result
    else:
        annual_fractions = {}
    
    # Step 4: Create visualizations
    logger.info("=" * 60)
    logger.info("Creating Visualizations")
    logger.info("=" * 60)
    
    if wy_totals:
        plot_water_year_time_series(wy_totals)
        plot_mean_annual_water_year_map()
    
    if annual_fractions:
        plot_annual_fraction_time_series(annual_fractions)
    
    plot_mean_fraction_map()
    
    if wy_totals and annual_fractions:
        create_summary_figure(wy_totals, annual_fractions)
    
    # Print summary
    logger.info("=" * 70)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Water year totals: {WATER_YEAR_DIR}")
    logger.info(f"Mean fractions: {FRACTION_DIR}")
    logger.info(f"Figures: {FIGURES_DIR}")
    
    if wy_totals:
        logger.info(f"\nWater Year Statistics (Oct-Sep):")
        logger.info(f"  Years processed: WY {min(wy_totals.keys())} - WY {max(wy_totals.keys())} ({len(wy_totals)} years)")
        logger.info(f"  Mean annual Peff: {np.mean(list(wy_totals.values())):.1f} mm")
        logger.info(f"  Min: WY {min(wy_totals, key=wy_totals.get)} = {min(wy_totals.values()):.1f} mm")
        logger.info(f"  Max: WY {max(wy_totals, key=wy_totals.get)} = {max(wy_totals.values()):.1f} mm")
    
    if annual_fractions:
        logger.info(f"\nAnnual Fraction Statistics (Water Year Oct-Sep):")
        logger.info(f"  Years: WY {min(annual_fractions.keys())} - WY {max(annual_fractions.keys())}")
        logger.info(f"  Mean annual fraction: {np.mean(list(annual_fractions.values())):.3f}")
        logger.info(f"  Min: WY {min(annual_fractions, key=annual_fractions.get)} = {min(annual_fractions.values()):.3f}")
        logger.info(f"  Max: WY {max(annual_fractions, key=annual_fractions.get)} = {max(annual_fractions.values()):.3f}")


if __name__ == '__main__':
    main()

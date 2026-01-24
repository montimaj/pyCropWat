"""
Command-line interface for pyCropWat.

Provides subcommands for:
- process: Calculate effective precipitation from GEE data
- aggregate: Temporal aggregation (seasonal, annual, growing season)
- analyze: Statistical analysis (anomaly, trend, zonal statistics)
- export: Export to NetCDF or Cloud-Optimized GeoTIFF
- plot: Visualization (time series, climatology, maps)
"""

import argparse
import logging
import sys
from pathlib import Path

from .core import EffectivePrecipitation
from .methods import list_available_methods

# Default PCML asset for Western U.S. effective precipitation
# Note: Only Western U.S. vectors overlapping the 17-state extent can be used with PCML
# (AZ, CA, CO, ID, KS, MT, NE, NV, NM, ND, OK, OR, SD, TX, UT, WA, WY)
PCML_DEFAULT_ASSET = 'projects/ee-peff-westus-unmasked/assets/effective_precip_monthly_unmasked'
PCML_DEFAULT_BAND = 'pcml'  # Special marker - actual bands are bYYYY_M format (e.g., b2015_9, b2016_10)
PCML_DEFAULT_SCALE = None  # Retrieved dynamically from asset using nominalScale()

# PCML annual fraction asset (available WY 2000-2024)
# Note: Only annual (water year, Oct-Sep) fractions are available for PCML, not monthly. Band format: bYYYY
PCML_FRACTION_ASSET = 'projects/ee-peff-westus-unmasked/assets/effective_precip_fraction_unmasked'


def get_pcml_band_name(year: int, month: int) -> str:
    """Get PCML band name for a specific year and month.
    
    PCML bands are formatted as bYYYY_M where months 1-9 do not have a preceding zero.
    Examples: b2015_9, b2016_10
    """
    return f"b{year}_{month}"


def setup_logging(verbose: bool = False) -> None:
    """Configure logging based on verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add common arguments to a parser."""
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )


def cmd_process(args):
    """Handle the 'process' subcommand."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    # Validate arguments
    # Note: PCML method doesn't require asset/band/geometry - it uses the PCML asset defaults
    if args.method != 'pcml':
        if args.asset is None or args.band is None:
            logger.error("--asset and --band are required (except for PCML method)")
            sys.exit(1)
        if args.geometry is None and args.gee_geometry is None:
            logger.error("Either --geometry or --gee-geometry must be provided (except for PCML method)")
            sys.exit(1)
    
    if args.geometry is not None and not Path(args.geometry).exists():
        from .utils import is_gee_asset
        if not is_gee_asset(args.geometry):
            logger.error(f"Geometry file not found: {args.geometry}")
            sys.exit(1)
    
    if args.start_year > args.end_year:
        logger.error("Start year must be less than or equal to end year")
        sys.exit(1)
    
    if args.months:
        invalid_months = [m for m in args.months if m < 1 or m > 12]
        if invalid_months:
            logger.error(f"Invalid months: {invalid_months}. Must be between 1 and 12.")
            sys.exit(1)
    
    # Validate USDA-SCS method requirements
    if args.method == 'usda_scs':
        if not args.awc_asset:
            logger.error("USDA-SCS method requires --awc-asset. "
                        "U.S.: projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite, "
                        "Global: projects/sat-io/open-datasets/FAO/HWSD_V2_SMU")
            sys.exit(1)
        if not args.eto_asset:
            logger.error("USDA-SCS method requires --eto-asset. "
                        "U.S.: projects/openet/assets/reference_et/conus/gridmet/monthly/v1, "
                        "Global: projects/climate-engine-pro/assets/ce-ag-era5-v2/daily (set --eto-is-daily)")
            sys.exit(1)
    
    # Validate ensemble method requirements (needs AWC and ETo for USDA-SCS component)
    if args.method == 'ensemble':
        if not args.awc_asset:
            logger.error("Ensemble method requires --awc-asset for USDA-SCS component. "
                        "U.S.: projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite, "
                        "Global: projects/sat-io/open-datasets/FAO/HWSD_V2_SMU")
            sys.exit(1)
        if not args.eto_asset:
            logger.error("Ensemble method requires --eto-asset for USDA-SCS/SuET components. "
                        "U.S.: projects/openet/assets/reference_et/conus/gridmet/monthly/v1, "
                        "Global: projects/climate-engine-pro/assets/ce-ag-era5-v2/daily (set --eto-is-daily)")
            sys.exit(1)
    
    # Validate SuET method requirements
    if args.method == 'suet':
        if not args.eto_asset:
            logger.error("SuET method requires --eto-asset. "
                        "U.S.: projects/openet/assets/reference_et/conus/gridmet/monthly/v1, "
                        "Global: projects/climate-engine-pro/assets/ce-ag-era5-v2/daily (set --eto-is-daily)")
            sys.exit(1)
    
    # Handle PCML method - always use default asset (overrides any user-provided values)
    if args.method == 'pcml':
        args.asset = PCML_DEFAULT_ASSET
        args.band = PCML_DEFAULT_BAND
        # Scale will be retrieved dynamically from PCML asset using nominalScale()
        args.scale = None  # Let core.py get it from asset
        logger.info(f"Using PCML default asset: {PCML_DEFAULT_ASSET}")
        logger.info(f"Using PCML default band: {PCML_DEFAULT_BAND}")
        logger.info(f"PCML native scale will be retrieved from asset using nominalScale()")
        logger.info(f"PCML annual fractions will be loaded from GEE asset: {PCML_FRACTION_ASSET}")
        if args.geometry is None and args.gee_geometry is None:
            logger.info(f"Using PCML asset geometry (entire Western U.S.)")
    
    try:
        logger.info(f"Initializing effective precipitation processor...")
        logger.info(f"Asset: {args.asset}")
        logger.info(f"Band: {args.band}")
        logger.info(f"Method: {args.method}")
        if args.gee_geometry:
            logger.info(f"GEE Geometry Asset: {args.gee_geometry}")
        elif args.geometry:
            logger.info(f"Geometry: {args.geometry}")
        elif args.method == 'pcml':
            logger.info(f"Geometry: Using PCML asset's built-in geometry (Western U.S.)")
        logger.info(f"Date range: {args.start_year} - {args.end_year}")
        
        # Build method parameters
        method_params = {}
        if args.method == 'fixed_percentage':
            method_params['percentage'] = args.percentage
        elif args.method == 'dependable_rainfall':
            method_params['probability'] = args.probability
        elif args.method == 'usda_scs':
            method_params['awc_asset'] = args.awc_asset
            method_params['awc_band'] = args.awc_band
            method_params['eto_asset'] = args.eto_asset
            method_params['eto_band'] = args.eto_band
            method_params['eto_is_daily'] = args.eto_is_daily
            method_params['rooting_depth'] = args.rooting_depth
            method_params['mad_factor'] = args.mad_factor
            band_info = f"band: {args.awc_band}" if args.awc_band else "single-band"
            logger.info(f"AWC Asset: {args.awc_asset} ({band_info})")
            logger.info(f"ETo Asset: {args.eto_asset} (band: {args.eto_band})")
            logger.info(f"Rooting Depth: {args.rooting_depth} m")
            logger.info(f"MAD Factor: {args.mad_factor}")
        elif args.method == 'suet':
            method_params['eto_asset'] = args.eto_asset
            method_params['eto_band'] = args.eto_band
            method_params['eto_is_daily'] = args.eto_is_daily
            logger.info(f"ETo Asset: {args.eto_asset} (band: {args.eto_band})")
        elif args.method == 'ensemble':
            method_params['awc_asset'] = args.awc_asset
            method_params['awc_band'] = args.awc_band
            method_params['eto_asset'] = args.eto_asset
            method_params['eto_band'] = args.eto_band
            method_params['eto_is_daily'] = args.eto_is_daily
            method_params['rooting_depth'] = args.rooting_depth
            method_params['mad_factor'] = args.mad_factor
            method_params['percentage'] = args.percentage
            method_params['probability'] = args.probability
            band_info = f"band: {args.awc_band}" if args.awc_band else "single-band"
            logger.info(f"AWC Asset: {args.awc_asset} ({band_info})")
            logger.info(f"ETo Asset: {args.eto_asset} (band: {args.eto_band})")
            logger.info(f"Rooting Depth: {args.rooting_depth} m")
            logger.info(f"MAD Factor: {args.mad_factor}")
        
        ep = EffectivePrecipitation(
            asset_id=args.asset,
            precip_band=args.band,
            geometry_path=args.geometry,
            start_year=args.start_year,
            end_year=args.end_year,
            scale=args.scale,
            precip_scale_factor=args.scale_factor,
            gee_project=args.project,
            gee_geometry_asset=args.gee_geometry,
            method=args.method,
            method_params=method_params
        )
        
        if args.sequential:
            logger.info("Processing sequentially...")
            results = ep.process_sequential(args.output, months=args.months)
        else:
            logger.info(f"Processing with {args.workers} workers...")
            results = ep.process(args.output, n_workers=args.workers, months=args.months)
        
        successful = sum(1 for r in results if r[0] is not None)
        logger.info(f"Processing complete. {successful}/{len(results)} months processed successfully.")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def cmd_aggregate(args):
    """Handle the 'aggregate' subcommand."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    try:
        from .analysis import TemporalAggregator
        
        agg = TemporalAggregator(args.input, pattern=args.pattern)
        
        if args.type == 'annual':
            logger.info(f"Calculating annual {args.method} for {args.year}...")
            result = agg.annual_aggregate(
                args.year,
                method=args.method,
                output_path=args.output
            )
            
        elif args.type == 'seasonal':
            if not args.season:
                logger.error("--season is required for seasonal aggregation")
                sys.exit(1)
            logger.info(f"Calculating {args.season} {args.method} for {args.year}...")
            result = agg.seasonal_aggregate(
                args.year,
                args.season,
                method=args.method,
                output_path=args.output
            )
            
        elif args.type == 'growing-season':
            logger.info(f"Calculating growing season ({args.start_month}-{args.end_month}) {args.method} for {args.year}...")
            result = agg.growing_season_aggregate(
                args.year,
                start_month=args.start_month,
                end_month=args.end_month,
                method=args.method,
                output_path=args.output
            )
            
        elif args.type == 'custom':
            if not args.months:
                logger.error("--months is required for custom aggregation")
                sys.exit(1)
            logger.info(f"Calculating custom aggregation for {args.year}, months {args.months}...")
            result = agg.custom_aggregate(
                args.year,
                months=args.months,
                method=args.method,
                output_path=args.output
            )
            
        elif args.type == 'climatology':
            if not args.start_year or not args.end_year:
                logger.error("--start-year and --end-year are required for climatology")
                sys.exit(1)
            logger.info(f"Calculating climatology for {args.start_year}-{args.end_year}...")
            result = agg.multi_year_climatology(
                args.start_year,
                args.end_year,
                months=args.months,
                output_dir=args.output
            )
        
        if result is not None:
            logger.info("Aggregation complete.")
        else:
            logger.warning("No data available for aggregation.")
            
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def cmd_analyze(args):
    """Handle the 'analyze' subcommand."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    try:
        from .analysis import StatisticalAnalyzer
        
        stats = StatisticalAnalyzer(args.input, pattern=args.pattern)
        
        if args.analysis_type == 'anomaly':
            if not all([args.year, args.month, args.clim_start, args.clim_end]):
                logger.error("--year, --month, --clim-start, and --clim-end are required for anomaly")
                sys.exit(1)
            logger.info(f"Calculating {args.anomaly_type} anomaly for {args.year}-{args.month:02d}...")
            result = stats.calculate_anomaly(
                args.year,
                args.month,
                args.clim_start,
                args.clim_end,
                anomaly_type=args.anomaly_type,
                output_path=args.output
            )
            
        elif args.analysis_type == 'trend':
            if not args.start_year or not args.end_year:
                logger.error("--start-year and --end-year are required for trend analysis")
                sys.exit(1)
            logger.info(f"Calculating {args.trend_method} trend for {args.start_year}-{args.end_year}...")
            slope, pvalue = stats.calculate_trend(
                args.start_year,
                args.end_year,
                month=args.month,
                method=args.trend_method,
                output_dir=args.output
            )
            if slope is not None:
                logger.info("Trend analysis complete. Saved slope and p-value rasters.")
            
        elif args.analysis_type == 'zonal':
            if not args.zones:
                logger.error("--zones is required for zonal statistics")
                sys.exit(1)
            if not args.start_year or not args.end_year:
                logger.error("--start-year and --end-year are required for zonal statistics")
                sys.exit(1)
            logger.info(f"Calculating zonal statistics for {args.start_year}-{args.end_year}...")
            df = stats.zonal_statistics(
                args.zones,
                args.start_year,
                args.end_year,
                months=args.months,
                stats=args.stats.split(',') if args.stats else ['mean', 'sum', 'min', 'max', 'std'],
                output_path=args.output
            )
            logger.info(f"Zonal statistics complete. {len(df)} records processed.")
            
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def cmd_export(args):
    """Handle the 'export' subcommand."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    try:
        if args.format == 'netcdf':
            from .analysis import export_to_netcdf
            logger.info(f"Exporting to NetCDF: {args.output}")
            export_to_netcdf(
                args.input,
                args.output,
                pattern=args.pattern,
                variable_name=args.variable or 'effective_precipitation',
                compression=not args.no_compression
            )
            logger.info("NetCDF export complete.")
            
        elif args.format == 'cog':
            from .analysis import export_to_cog
            logger.info(f"Converting to COG: {args.output}")
            export_to_cog(
                args.input,
                args.output,
                overview_levels=[2, 4, 8, 16]
            )
            logger.info("COG conversion complete.")
            
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def cmd_plot(args):
    """Handle the 'plot' subcommand."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    try:
        from .analysis import Visualizer
        
        viz = Visualizer(args.input, pattern=args.pattern)
        
        if args.plot_type == 'timeseries':
            if not args.start_year or not args.end_year:
                logger.error("--start-year and --end-year are required for time series plot")
                sys.exit(1)
            logger.info(f"Creating time series plot for {args.start_year}-{args.end_year}...")
            viz.plot_time_series(
                args.start_year,
                args.end_year,
                months=args.months,
                stat=args.stat,
                title=args.title,
                output_path=args.output,
                figsize=(args.width, args.height)
            )
            
        elif args.plot_type == 'climatology':
            if not args.start_year or not args.end_year:
                logger.error("--start-year and --end-year are required for climatology plot")
                sys.exit(1)
            logger.info(f"Creating climatology plot for {args.start_year}-{args.end_year}...")
            viz.plot_monthly_climatology(
                args.start_year,
                args.end_year,
                stat=args.stat,
                title=args.title,
                output_path=args.output,
                figsize=(args.width, args.height)
            )
            
        elif args.plot_type == 'map':
            if not args.year or not args.month:
                logger.error("--year and --month are required for map plot")
                sys.exit(1)
            logger.info(f"Creating map for {args.year}-{args.month:02d}...")
            viz.plot_raster(
                args.year,
                args.month,
                cmap=args.cmap,
                title=args.title,
                output_path=args.output,
                figsize=(args.width, args.height),
                vmin=args.vmin,
                vmax=args.vmax
            )
        
        elif args.plot_type == 'interactive':
            if not args.year or not args.month:
                logger.error("--year and --month are required for interactive map")
                sys.exit(1)
            logger.info(f"Creating interactive map for {args.year}-{args.month:02d}...")
            viz.plot_interactive_map(
                args.year,
                args.month,
                cmap=args.cmap,
                title=args.title,
                output_path=args.output,
                zoom_start=args.zoom,
                opacity=args.opacity
            )
        
        elif args.plot_type == 'compare':
            if not args.year or not args.month:
                logger.error("--year and --month are required for comparison plot")
                sys.exit(1)
            if not args.other_input:
                logger.error("--other-input is required for comparison plot")
                sys.exit(1)
            logger.info(f"Creating comparison plot for {args.year}-{args.month:02d}...")
            viz.plot_comparison(
                args.year,
                args.month,
                other_dir=args.other_input,
                other_pattern=args.other_pattern,
                labels=(args.label1, args.label2),
                cmap=args.cmap,
                diff_cmap=args.diff_cmap,
                title=args.title,
                output_path=args.output,
                figsize=(args.width, args.height)
            )
        
        elif args.plot_type == 'scatter':
            if not args.start_year or not args.end_year:
                logger.error("--start-year and --end-year are required for scatter plot")
                sys.exit(1)
            if not args.other_input:
                logger.error("--other-input is required for scatter plot")
                sys.exit(1)
            logger.info(f"Creating scatter comparison for {args.start_year}-{args.end_year}...")
            viz.plot_scatter_comparison(
                args.start_year,
                args.end_year,
                other_dir=args.other_input,
                other_pattern=args.other_pattern,
                labels=(args.label1, args.label2),
                months=args.months,
                sample_size=args.sample_size,
                title=args.title,
                output_path=args.output,
                figsize=(args.width, args.height)
            )
        
        elif args.plot_type == 'annual-compare':
            if not args.start_year or not args.end_year:
                logger.error("--start-year and --end-year are required for annual comparison")
                sys.exit(1)
            if not args.other_input:
                logger.error("--other-input is required for annual comparison")
                sys.exit(1)
            logger.info(f"Creating annual comparison for {args.start_year}-{args.end_year}...")
            viz.plot_annual_comparison(
                args.start_year,
                args.end_year,
                other_dir=args.other_input,
                other_pattern=args.other_pattern,
                labels=(args.label1, args.label2),
                stat=args.stat,
                title=args.title,
                output_path=args.output,
                figsize=(args.width, args.height)
            )
        
        logger.info(f"Plot saved to: {args.output}")
            
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def create_parser():
    """Create the argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog='pycropwat',
        description='pyCropWat - Calculate effective precipitation from GEE climate data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.2'
    )
    
    parser.add_argument(
        '--list-methods',
        action='store_true',
        help='List available effective precipitation methods and exit'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # =========================================================================
    # PROCESS subcommand
    # =========================================================================
    process_parser = subparsers.add_parser(
        'process',
        help='Calculate effective precipitation from GEE climate data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pycropwat process --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \\
                    --band total_precipitation_sum \\
                    --gee-geometry projects/my-project/assets/roi \\
                    --start-year 2020 --end-year 2023 \\
                    --scale-factor 1000 --output ./output

  pycropwat process --asset IDAHO_EPSCOR/TERRACLIMATE --band pr \\
                    --geometry roi.geojson \\
                    --start-year 2015 --end-year 2020 \\
                    --method fao_aglw --output ./output

  # PCML method (Western U.S. only - no asset/band/geometry required)
  pycropwat process --method pcml --start-year 2020 --end-year 2024 \\
                    --output ./WesternUS_PCML --workers 8
        """
    )
    process_parser.add_argument('--asset', '-a', help='GEE ImageCollection asset ID (not required for PCML method)')
    process_parser.add_argument('--band', '-b', help='Precipitation band name (not required for PCML method)')
    process_parser.add_argument('--geometry', '-g', type=str, help='Path to shapefile or GeoJSON')
    process_parser.add_argument('--gee-geometry', '-G', type=str, help='GEE FeatureCollection asset ID')
    process_parser.add_argument('--start-year', '-s', required=True, type=int, help='Start year')
    process_parser.add_argument('--end-year', '-e', required=True, type=int, help='End year')
    process_parser.add_argument('--output', '-o', required=True, type=Path, help='Output directory')
    process_parser.add_argument('--scale-factor', '-f', type=float, default=1.0, help='Precipitation scale factor')
    process_parser.add_argument('--scale', '-r', type=float, help='Output resolution in meters')
    process_parser.add_argument('--workers', '-w', type=int, default=4, help='Number of parallel workers')
    process_parser.add_argument('--months', '-m', type=int, nargs='+', help='Specific months to process (1-12)')
    process_parser.add_argument('--project', '-p', type=str, help='GEE project ID')
    process_parser.add_argument('--method', type=str, default='ensemble',
                               choices=['cropwat', 'fao_aglw', 'fixed_percentage', 'dependable_rainfall', 'farmwest', 'usda_scs', 'suet', 'pcml', 'ensemble'],
                               help='Effective precipitation method (pcml uses default PCML asset for Western U.S.)')
    process_parser.add_argument('--percentage', type=float, default=0.7, help='Percentage for fixed_percentage method')
    process_parser.add_argument('--probability', type=float, default=0.75, help='Probability for dependable_rainfall method')
    # USDA-SCS method parameters
    process_parser.add_argument('--awc-asset', type=str, 
                               help='GEE AWC asset for usda_scs/ensemble method. U.S.: projects/openet/soil/ssurgo_AWC_WTA_0to152cm_composite, Global: projects/sat-io/open-datasets/FAO/HWSD_V2_SMU')
    process_parser.add_argument('--awc-band', type=str, default=None,
                               help='AWC band name. Omit for SSURGO (single-band), use "AWC" for HWSD')
    process_parser.add_argument('--eto-asset', type=str,
                               help='GEE ETo asset for usda_scs/suet methods. U.S.: projects/openet/assets/reference_et/conus/gridmet/monthly/v1, Global: projects/climate-engine-pro/assets/ce-ag-era5-v2/daily')
    process_parser.add_argument('--eto-band', type=str, default='eto',
                               help='ETo band name. GridMET: "eto", ERA5: "ReferenceET_PenmanMonteith_FAO56"')
    process_parser.add_argument('--eto-is-daily', action='store_true',
                               help='Set if ETo asset is daily (will aggregate to monthly)')
    process_parser.add_argument('--rooting-depth', type=float, default=1.0,
                               help='Crop rooting depth in meters for usda_scs method (default: 1.0)')
    process_parser.add_argument('--mad-factor', type=float, default=0.5,
                               help='Management Allowed Depletion factor (0-1) for usda_scs method (default: 0.5)')
    process_parser.add_argument('--sequential', action='store_true', help='Process sequentially')
    add_common_args(process_parser)
    process_parser.set_defaults(func=cmd_process)
    
    # =========================================================================
    # AGGREGATE subcommand
    # =========================================================================
    agg_parser = subparsers.add_parser(
        'aggregate',
        help='Temporal aggregation of effective precipitation rasters',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Annual total
  pycropwat aggregate --input ./output --type annual --year 2020 --output ./annual_2020.tif

  # Seasonal (summer)
  pycropwat aggregate --input ./output --type seasonal --year 2020 --season JJA --output ./summer_2020.tif

  # Growing season (April-October)
  pycropwat aggregate --input ./output --type growing-season --year 2020 \\
                      --start-month 4 --end-month 10 --output ./growing_2020.tif

  # Multi-year climatology
  pycropwat aggregate --input ./output --type climatology \\
                      --start-year 2000 --end-year 2020 --output ./climatology/
        """
    )
    agg_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory with monthly rasters')
    agg_parser.add_argument('--type', '-t', required=True,
                           choices=['annual', 'seasonal', 'growing-season', 'custom', 'climatology'],
                           help='Aggregation type')
    agg_parser.add_argument('--year', '-y', type=int, help='Year to aggregate')
    agg_parser.add_argument('--start-year', type=int, help='Start year (for climatology)')
    agg_parser.add_argument('--end-year', type=int, help='End year (for climatology)')
    agg_parser.add_argument('--season', choices=['DJF', 'MAM', 'JJA', 'SON'], help='Season code')
    agg_parser.add_argument('--start-month', type=int, default=4, help='Growing season start month')
    agg_parser.add_argument('--end-month', type=int, default=10, help='Growing season end month')
    agg_parser.add_argument('--months', '-m', type=int, nargs='+', help='Specific months for custom aggregation')
    agg_parser.add_argument('--method', default='sum', choices=['sum', 'mean', 'min', 'max', 'std'],
                           help='Aggregation method')
    agg_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    agg_parser.add_argument('--output', '-o', required=True, help='Output file or directory')
    add_common_args(agg_parser)
    agg_parser.set_defaults(func=cmd_aggregate)
    
    # =========================================================================
    # ANALYZE subcommand
    # =========================================================================
    analyze_parser = subparsers.add_parser(
        'analyze',
        help='Statistical analysis of effective precipitation data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate anomaly
  pycropwat analyze anomaly --input ./output --year 2020 --month 6 \\
                            --clim-start 1990 --clim-end 2020 --output ./anomaly_2020_06.tif

  # Calculate trend
  pycropwat analyze trend --input ./output --start-year 2000 --end-year 2020 \\
                          --method sen --output ./trend/

  # Zonal statistics
  pycropwat analyze zonal --input ./output --zones ./regions.shp \\
                          --start-year 2000 --end-year 2020 --output ./zonal_stats.csv
        """
    )
    analyze_subparsers = analyze_parser.add_subparsers(dest='analysis_type', help='Analysis type')
    
    # Anomaly
    anomaly_parser = analyze_subparsers.add_parser('anomaly', help='Calculate anomaly relative to climatology')
    anomaly_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    anomaly_parser.add_argument('--year', '-y', required=True, type=int, help='Target year')
    anomaly_parser.add_argument('--month', '-m', required=True, type=int, help='Target month (1-12)')
    anomaly_parser.add_argument('--clim-start', required=True, type=int, help='Climatology start year')
    anomaly_parser.add_argument('--clim-end', required=True, type=int, help='Climatology end year')
    anomaly_parser.add_argument('--anomaly-type', default='absolute',
                               choices=['absolute', 'percent', 'standardized'],
                               help='Type of anomaly')
    anomaly_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    anomaly_parser.add_argument('--output', '-o', required=True, help='Output file')
    add_common_args(anomaly_parser)
    anomaly_parser.set_defaults(func=cmd_analyze)
    
    # Trend
    trend_parser = analyze_subparsers.add_parser('trend', help='Calculate temporal trend')
    trend_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    trend_parser.add_argument('--start-year', required=True, type=int, help='Start year')
    trend_parser.add_argument('--end-year', required=True, type=int, help='End year')
    trend_parser.add_argument('--month', '-m', type=int, help='Specific month (or annual if omitted)')
    trend_parser.add_argument('--trend-method', default='linear', choices=['linear', 'sen'],
                             help='Trend method (linear or Theil-Sen)')
    trend_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    trend_parser.add_argument('--output', '-o', required=True, help='Output directory')
    add_common_args(trend_parser)
    trend_parser.set_defaults(func=cmd_analyze)
    
    # Zonal
    zonal_parser = analyze_subparsers.add_parser('zonal', help='Calculate zonal statistics')
    zonal_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    zonal_parser.add_argument('--zones', '-z', required=True, help='Path to zones shapefile/GeoJSON')
    zonal_parser.add_argument('--start-year', required=True, type=int, help='Start year')
    zonal_parser.add_argument('--end-year', required=True, type=int, help='End year')
    zonal_parser.add_argument('--months', '-m', type=int, nargs='+', help='Specific months')
    zonal_parser.add_argument('--stats', default='mean,sum,min,max,std',
                             help='Statistics to compute (comma-separated)')
    zonal_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    zonal_parser.add_argument('--output', '-o', required=True, help='Output CSV file')
    add_common_args(zonal_parser)
    zonal_parser.set_defaults(func=cmd_analyze)
    
    # =========================================================================
    # EXPORT subcommand
    # =========================================================================
    export_parser = subparsers.add_parser(
        'export',
        help='Export data to different formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Export to NetCDF
  pycropwat export netcdf --input ./output --output ./data.nc

  # Convert to Cloud-Optimized GeoTIFF
  pycropwat export cog --input ./output/effective_precip_2020_06.tif --output ./cog_2020_06.tif
        """
    )
    export_subparsers = export_parser.add_subparsers(dest='format', help='Export format')
    
    # NetCDF
    netcdf_parser = export_subparsers.add_parser('netcdf', help='Export to NetCDF')
    netcdf_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    netcdf_parser.add_argument('--output', '-o', required=True, help='Output NetCDF file')
    netcdf_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    netcdf_parser.add_argument('--variable', help='Variable name in NetCDF')
    netcdf_parser.add_argument('--no-compression', action='store_true', help='Disable compression')
    add_common_args(netcdf_parser)
    netcdf_parser.set_defaults(func=cmd_export)
    
    # COG
    cog_parser = export_subparsers.add_parser('cog', help='Convert to Cloud-Optimized GeoTIFF')
    cog_parser.add_argument('--input', '-i', required=True, help='Input GeoTIFF file')
    cog_parser.add_argument('--output', '-o', required=True, help='Output COG file')
    add_common_args(cog_parser)
    cog_parser.set_defaults(func=cmd_export)
    
    # =========================================================================
    # PLOT subcommand
    # =========================================================================
    plot_parser = subparsers.add_parser(
        'plot',
        help='Create visualizations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Time series plot
  pycropwat plot timeseries --input ./output --start-year 2000 --end-year 2020 \\
                            --output ./timeseries.png

  # Monthly climatology bar chart
  pycropwat plot climatology --input ./output --start-year 2000 --end-year 2020 \\
                             --output ./climatology.png

  # Single month map
  pycropwat plot map --input ./output --year 2020 --month 6 --output ./map_2020_06.png

  # Interactive map (requires leafmap or folium)
  pycropwat plot interactive --input ./output --year 2020 --month 6 --output ./map.html

  # Compare two datasets side-by-side
  pycropwat plot compare --input ./era5_output --other-input ./terraclimate_output \\
                         --year 2020 --month 6 --label1 ERA5 --label2 TerraClimate \\
                         --output ./comparison.png

  # Scatter plot comparison
  pycropwat plot scatter --input ./era5_output --other-input ./terraclimate_output \\
                         --start-year 2000 --end-year 2020 --output ./scatter.png

  # Annual comparison bar chart
  pycropwat plot annual-compare --input ./era5_output --other-input ./terraclimate_output \\
                                --start-year 2000 --end-year 2020 --output ./annual.png
        """
    )
    plot_subparsers = plot_parser.add_subparsers(dest='plot_type', help='Plot type')
    
    # Time series
    ts_parser = plot_subparsers.add_parser('timeseries', help='Time series plot')
    ts_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    ts_parser.add_argument('--start-year', required=True, type=int, help='Start year')
    ts_parser.add_argument('--end-year', required=True, type=int, help='End year')
    ts_parser.add_argument('--months', '-m', type=int, nargs='+', help='Specific months')
    ts_parser.add_argument('--stat', default='mean', choices=['mean', 'sum', 'min', 'max'],
                          help='Spatial aggregation statistic')
    ts_parser.add_argument('--title', help='Plot title')
    ts_parser.add_argument('--width', type=int, default=12, help='Figure width in inches')
    ts_parser.add_argument('--height', type=int, default=6, help='Figure height in inches')
    ts_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    ts_parser.add_argument('--output', '-o', required=True, help='Output image file')
    add_common_args(ts_parser)
    ts_parser.set_defaults(func=cmd_plot)
    
    # Climatology
    clim_parser = plot_subparsers.add_parser('climatology', help='Monthly climatology bar chart')
    clim_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    clim_parser.add_argument('--start-year', required=True, type=int, help='Start year')
    clim_parser.add_argument('--end-year', required=True, type=int, help='End year')
    clim_parser.add_argument('--stat', default='mean', choices=['mean', 'sum'],
                            help='Spatial aggregation statistic')
    clim_parser.add_argument('--title', help='Plot title')
    clim_parser.add_argument('--width', type=int, default=10, help='Figure width in inches')
    clim_parser.add_argument('--height', type=int, default=6, help='Figure height in inches')
    clim_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    clim_parser.add_argument('--output', '-o', required=True, help='Output image file')
    add_common_args(clim_parser)
    clim_parser.set_defaults(func=cmd_plot)
    
    # Map
    map_parser = plot_subparsers.add_parser('map', help='Single raster map')
    map_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    map_parser.add_argument('--year', '-y', required=True, type=int, help='Year')
    map_parser.add_argument('--month', '-m', required=True, type=int, help='Month (1-12)')
    map_parser.add_argument('--cmap', default='YlGnBu', help='Colormap name')
    map_parser.add_argument('--title', help='Plot title')
    map_parser.add_argument('--vmin', type=float, help='Minimum color scale value')
    map_parser.add_argument('--vmax', type=float, help='Maximum color scale value')
    map_parser.add_argument('--width', type=int, default=10, help='Figure width in inches')
    map_parser.add_argument('--height', type=int, default=8, help='Figure height in inches')
    map_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    map_parser.add_argument('--output', '-o', required=True, help='Output image file')
    add_common_args(map_parser)
    map_parser.set_defaults(func=cmd_plot)
    
    # Interactive map
    interactive_parser = plot_subparsers.add_parser('interactive', help='Interactive map (HTML)')
    interactive_parser.add_argument('--input', '-i', required=True, type=Path, help='Input directory')
    interactive_parser.add_argument('--year', '-y', required=True, type=int, help='Year')
    interactive_parser.add_argument('--month', '-m', required=True, type=int, help='Month (1-12)')
    interactive_parser.add_argument('--cmap', default='YlGnBu', help='Colormap name')
    interactive_parser.add_argument('--title', help='Map title')
    interactive_parser.add_argument('--zoom', type=int, default=6, help='Initial zoom level')
    interactive_parser.add_argument('--opacity', type=float, default=0.7, help='Layer opacity (0-1)')
    interactive_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    interactive_parser.add_argument('--output', '-o', required=True, help='Output HTML file')
    add_common_args(interactive_parser)
    interactive_parser.set_defaults(func=cmd_plot)
    
    # Comparison plot
    compare_parser = plot_subparsers.add_parser('compare', help='Side-by-side dataset comparison')
    compare_parser.add_argument('--input', '-i', required=True, type=Path, help='Primary input directory')
    compare_parser.add_argument('--other-input', required=True, type=Path, help='Secondary input directory')
    compare_parser.add_argument('--year', '-y', required=True, type=int, help='Year')
    compare_parser.add_argument('--month', '-m', required=True, type=int, help='Month (1-12)')
    compare_parser.add_argument('--label1', default='Dataset 1', help='Label for primary dataset')
    compare_parser.add_argument('--label2', default='Dataset 2', help='Label for secondary dataset')
    compare_parser.add_argument('--cmap', default='YlGnBu', help='Colormap for data')
    compare_parser.add_argument('--diff-cmap', default='RdBu', help='Colormap for difference')
    compare_parser.add_argument('--title', help='Plot title')
    compare_parser.add_argument('--width', type=int, default=16, help='Figure width in inches')
    compare_parser.add_argument('--height', type=int, default=5, help='Figure height in inches')
    compare_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    compare_parser.add_argument('--other-pattern', default='effective_precip_[0-9]*.tif', help='Pattern for secondary dataset')
    compare_parser.add_argument('--output', '-o', required=True, help='Output image file')
    add_common_args(compare_parser)
    compare_parser.set_defaults(func=cmd_plot)
    
    # Scatter plot comparison
    scatter_parser = plot_subparsers.add_parser('scatter', help='Scatter plot comparison')
    scatter_parser.add_argument('--input', '-i', required=True, type=Path, help='Primary input directory')
    scatter_parser.add_argument('--other-input', required=True, type=Path, help='Secondary input directory')
    scatter_parser.add_argument('--start-year', required=True, type=int, help='Start year')
    scatter_parser.add_argument('--end-year', required=True, type=int, help='End year')
    scatter_parser.add_argument('--months', '-m', type=int, nargs='+', help='Specific months')
    scatter_parser.add_argument('--label1', default='Dataset 1', help='Label for primary dataset')
    scatter_parser.add_argument('--label2', default='Dataset 2', help='Label for secondary dataset')
    scatter_parser.add_argument('--sample-size', type=int, default=10000, help='Max samples for plot')
    scatter_parser.add_argument('--title', help='Plot title')
    scatter_parser.add_argument('--width', type=int, default=8, help='Figure width in inches')
    scatter_parser.add_argument('--height', type=int, default=8, help='Figure height in inches')
    scatter_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    scatter_parser.add_argument('--other-pattern', default='effective_precip_[0-9]*.tif', help='Pattern for secondary dataset')
    scatter_parser.add_argument('--output', '-o', required=True, help='Output image file')
    add_common_args(scatter_parser)
    scatter_parser.set_defaults(func=cmd_plot)
    
    # Annual comparison
    annual_compare_parser = plot_subparsers.add_parser('annual-compare', help='Annual totals comparison')
    annual_compare_parser.add_argument('--input', '-i', required=True, type=Path, help='Primary input directory')
    annual_compare_parser.add_argument('--other-input', required=True, type=Path, help='Secondary input directory')
    annual_compare_parser.add_argument('--start-year', required=True, type=int, help='Start year')
    annual_compare_parser.add_argument('--end-year', required=True, type=int, help='End year')
    annual_compare_parser.add_argument('--label1', default='Dataset 1', help='Label for primary dataset')
    annual_compare_parser.add_argument('--label2', default='Dataset 2', help='Label for secondary dataset')
    annual_compare_parser.add_argument('--stat', default='mean', choices=['mean', 'sum'], help='Spatial statistic')
    annual_compare_parser.add_argument('--title', help='Plot title')
    annual_compare_parser.add_argument('--width', type=int, default=12, help='Figure width in inches')
    annual_compare_parser.add_argument('--height', type=int, default=6, help='Figure height in inches')
    annual_compare_parser.add_argument('--pattern', default='effective_precip_[0-9]*.tif', help='File glob pattern (use [0-9]* to exclude fraction files)')
    annual_compare_parser.add_argument('--other-pattern', default='effective_precip_[0-9]*.tif', help='Pattern for secondary dataset')
    annual_compare_parser.add_argument('--output', '-o', required=True, help='Output image file')
    add_common_args(annual_compare_parser)
    annual_compare_parser.set_defaults(func=cmd_plot)
    
    return parser


def main():
    """Main entry point for the CLI."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Handle --list-methods at top level
    if args.list_methods:
        print("\nAvailable effective precipitation methods:\n")
        for name, desc in list_available_methods().items():
            print(f"  {name:25s} {desc}")
        print()
        sys.exit(0)
    
    # If no command provided, check for legacy usage or show help
    if args.command is None:
        # Check if legacy arguments are provided
        if len(sys.argv) > 1 and sys.argv[1].startswith('--'):
            # Legacy mode - redirect to process command
            print("Note: Running in legacy mode. Consider using 'pycropwat process' subcommand.")
            print("Run 'pycropwat --help' for full command reference.\n")
            # Re-parse with process as default
            sys.argv.insert(1, 'process')
            args = parser.parse_args()
        else:
            parser.print_help()
            sys.exit(0)
    
    # Execute the subcommand
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(0)


if __name__ == '__main__':
    main()

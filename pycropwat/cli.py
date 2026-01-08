"""
Command-line interface for pyCropWat.
"""

import argparse
import logging
import sys
from pathlib import Path

from .core import EffectivePrecipitation


def setup_logging(verbose: bool = False) -> None:
    """Configure logging based on verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description='Calculate CROPWAT effective precipitation from GEE climate data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process ERA5-Land data for a region
  pycropwat --asset ECMWF/ERA5_LAND/MONTHLY_AGGR \\
            --band total_precipitation_sum \\
            --geometry roi.geojson \\
            --start-year 2015 \\
            --end-year 2020 \\
            --output ./output \\
            --scale-factor 1000

  # Process CHIRPS data
  pycropwat --asset UCSB-CHG/CHIRPS/DAILY \\
            --band precipitation \\
            --geometry watershed.shp \\
            --start-year 2010 \\
            --end-year 2015 \\
            --output ./chirps_output \\
            --scale-factor 1
        """
    )
    
    parser.add_argument(
        '--asset', '-a',
        required=True,
        help='GEE ImageCollection asset ID (e.g., ECMWF/ERA5_LAND/MONTHLY_AGGR)'
    )
    
    parser.add_argument(
        '--band', '-b',
        required=True,
        help='Name of the precipitation band in the asset'
    )
    
    parser.add_argument(
        '--geometry', '-g',
        type=str,
        default=None,
        help='Path to shapefile or GeoJSON file for region of interest'
    )
    
    parser.add_argument(
        '--gee-geometry', '-G',
        type=str,
        default=None,
        help='GEE FeatureCollection asset ID for region of interest (e.g., projects/my-project/assets/boundary)'
    )
    
    parser.add_argument(
        '--start-year', '-s',
        required=True,
        type=int,
        help='Start year (inclusive)'
    )
    
    parser.add_argument(
        '--end-year', '-e',
        required=True,
        type=int,
        help='End year (inclusive)'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        type=Path,
        help='Output directory for rasters'
    )
    
    parser.add_argument(
        '--scale-factor', '-f',
        type=float,
        default=1.0,
        help='Factor to convert precipitation to mm (default: 1.0)'
    )
    
    parser.add_argument(
        '--scale', '-r',
        type=float,
        default=None,
        help='Output resolution in meters. If not specified, uses native resolution of input data.'
    )
    
    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=4,
        help='Number of parallel workers (default: 4)'
    )
    
    parser.add_argument(
        '--months', '-m',
        type=int,
        nargs='+',
        help='Specific months to process (1-12). If not specified, all months are processed.'
    )
    
    parser.add_argument(
        '--project', '-p',
        type=str,
        default=None,
        help='GEE project ID for authentication'
    )
    
    parser.add_argument(
        '--sequential',
        action='store_true',
        help='Process sequentially instead of in parallel (useful for debugging)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    # Validate arguments
    if args.geometry is None and args.gee_geometry is None:
        logger.error("Either --geometry or --gee-geometry must be provided")
        sys.exit(1)
    
    if args.geometry is not None and not Path(args.geometry).exists():
        # Check if it might be a GEE asset path
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
    
    try:
        # Initialize processor
        logger.info(f"Initializing effective precipitation processor...")
        logger.info(f"Asset: {args.asset}")
        logger.info(f"Band: {args.band}")
        if args.gee_geometry:
            logger.info(f"GEE Geometry Asset: {args.gee_geometry}")
        else:
            logger.info(f"Geometry: {args.geometry}")
        logger.info(f"Date range: {args.start_year} - {args.end_year}")
        
        ep = EffectivePrecipitation(
            asset_id=args.asset,
            precip_band=args.band,
            geometry_path=args.geometry,
            start_year=args.start_year,
            end_year=args.end_year,
            scale=args.scale,
            precip_scale_factor=args.scale_factor,
            gee_project=args.project,
            gee_geometry_asset=args.gee_geometry
        )
        
        # Process
        if args.sequential:
            logger.info("Processing sequentially...")
            results = ep.process_sequential(args.output, months=args.months)
        else:
            logger.info(f"Processing with {args.workers} workers...")
            results = ep.process(args.output, n_workers=args.workers, months=args.months)
        
        # Summary
        successful = sum(1 for r in results if r[0] is not None)
        logger.info(f"Processing complete. {successful}/{len(results)} months processed successfully.")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

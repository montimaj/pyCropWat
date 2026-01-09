"""
Core module for effective precipitation calculations using Google Earth Engine.
"""

import logging
import tempfile
from pathlib import Path
from typing import Union, Optional, List, Tuple
import ee
import numpy as np
import xarray as xr
import rioxarray
from rioxarray.merge import merge_arrays
import dask
from dask import delayed, compute
from dask.diagnostics import ProgressBar

from .utils import load_geometry, get_date_range, get_monthly_dates, initialize_gee
from .methods import (
    cropwat_effective_precip,
    fao_aglw_effective_precip,
    fixed_percentage_effective_precip,
    dependable_rainfall_effective_precip,
    get_method_function,
    PeffMethod
)

logger = logging.getLogger(__name__)

# Maximum pixels per tile for GEE sampleRectangle (conservative limit)
MAX_PIXELS_PER_TILE = 65536  # 256 x 256


class EffectivePrecipitation:
    """
    Calculate CROPWAT effective precipitation from GEE climate data.
    
    The CROPWAT effective precipitation formula:
    - If pr <= 250mm: ep = pr * (125 - 0.2 * pr) / 125
    - If pr > 250mm: ep = 0.1 * pr + 125
    
    Parameters
    ----------
    asset_id : str
        GEE ImageCollection asset ID (e.g., 'ECMWF/ERA5_LAND/MONTHLY_AGGR').
    precip_band : str
        Name of the precipitation band in the asset.
    geometry_path : str, Path, or None
        Path to shapefile or GeoJSON file defining the region of interest.
        Can also be a GEE FeatureCollection asset ID. Set to None if using
        gee_geometry_asset instead.
    start_year : int
        Start year for processing (inclusive).
    end_year : int
        End year for processing (inclusive).
    scale : float, optional
        Output resolution in meters. Default is native resolution (None).
    precip_scale_factor : float, optional
        Factor to convert precipitation to mm. Default is 1000 (for ERA5 m to mm).
    gee_project : str, optional
        GEE project ID for authentication.
    gee_geometry_asset : str, optional
        GEE FeatureCollection asset ID for the region of interest.
        Takes precedence over geometry_path if both are provided.
        
    Examples
    --------
    >>> from pycropwat import EffectivePrecipitation
    >>> # Using local file
    >>> ep = EffectivePrecipitation(
    ...     asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    ...     precip_band='total_precipitation_sum',
    ...     geometry_path='roi.geojson',
    ...     start_year=2015,
    ...     end_year=2020,
    ...     precip_scale_factor=1000
    ... )
    >>> 
    >>> # Using GEE FeatureCollection asset
    >>> ep = EffectivePrecipitation(
    ...     asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR',
    ...     precip_band='total_precipitation_sum',
    ...     gee_geometry_asset='projects/my-project/assets/study_area',
    ...     start_year=2015,
    ...     end_year=2020,
    ...     precip_scale_factor=1000
    ... )
    >>> ep.process(output_dir='./output', n_workers=4)
    """
    
    def __init__(
        self,
        asset_id: str,
        precip_band: str,
        geometry_path: Optional[Union[str, Path]] = None,
        start_year: int = None,
        end_year: int = None,
        scale: Optional[float] = None,
        precip_scale_factor: float = 1.0,
        gee_project: Optional[str] = None,
        gee_geometry_asset: Optional[str] = None,
        method: PeffMethod = 'cropwat',
        method_params: Optional[dict] = None,
    ):
        self.asset_id = asset_id
        self.precip_band = precip_band
        self.geometry_path = geometry_path
        self.gee_geometry_asset = gee_geometry_asset
        self.start_year = start_year
        self.end_year = end_year
        self.scale = scale  # None means use native resolution
        self.precip_scale_factor = precip_scale_factor
        self.gee_project = gee_project
        self.method = method
        self.method_params = method_params or {}
        
        # Get the effective precipitation function
        self._peff_function = get_method_function(method)
        
        # USDA-SCS specific: cache for AWC data (loaded once)
        self._awc_cache = None
        
        # Input directory for saving downloaded data (set during process())
        self._input_dir = None
        
        # Validate that at least one geometry source is provided
        if geometry_path is None and gee_geometry_asset is None:
            raise ValueError("Either geometry_path or gee_geometry_asset must be provided")
        
        # Initialize GEE
        initialize_gee(self.gee_project)
        
        # Load geometry from GEE asset or local file
        self.geometry = load_geometry(geometry_path, gee_asset=gee_geometry_asset)
        self.bounds = self.geometry.bounds().getInfo()['coordinates'][0]
        
        # Get date range
        self.start_date, self.end_date = get_date_range(start_year, end_year)
        
        # Load and filter image collection
        self._load_collection()
        
    def _load_collection(self) -> None:
        """Load and prepare the GEE ImageCollection."""
        self.collection = (
            ee.ImageCollection(self.asset_id)
            .select(self.precip_band)
            .filterDate(self.start_date, self.end_date)
            .filterBounds(self.geometry)
        )
        
        # Apply scale factor and rename band to 'pr'
        def scale_and_rename(img):
            return (
                img.multiply(self.precip_scale_factor)
                .rename('pr')
                .copyProperties(img, ['system:time_start', 'system:time_end'])
            )
        
        self.collection = self.collection.map(scale_and_rename)
        
    @staticmethod
    def cropwat_effective_precip(pr: np.ndarray) -> np.ndarray:
        """
        Calculate CROPWAT effective precipitation.
        
        Parameters
        ----------
        pr : np.ndarray
            Precipitation in mm.
            
        Returns
        -------
        np.ndarray
            Effective precipitation in mm.
        """
        ep = np.where(
            pr <= 250,
            pr * (125 - 0.2 * pr) / 125,
            0.1 * pr + 125
        )
        return ep
    
    def _get_native_scale(self) -> float:
        """
        Get the native scale (resolution) of the image collection in meters.
        
        Returns
        -------
        float
            Native scale in meters.
        """
        try:
            # Get the first image from the collection to determine native scale
            first_img = self.collection.first()
            projection = first_img.projection()
            native_scale = projection.nominalScale().getInfo()
            logger.info(f"Native scale: {native_scale} meters")
            return native_scale
        except Exception as e:
            logger.warning(f"Could not determine native scale, defaulting to 10000m: {e}")
            return 10000.0
    
    def _estimate_pixel_count(self, bounds_coords: List, scale_meters: float) -> int:
        """
        Estimate the number of pixels for a given bounds and scale.
        
        Parameters
        ----------
        bounds_coords : list
            Bounding box coordinates [[min_lon, min_lat], [max_lon, min_lat], ...]
        scale_meters : float
            Resolution in meters.
            
        Returns
        -------
        int
            Estimated number of pixels.
        """
        min_lon = min(c[0] for c in bounds_coords)
        max_lon = max(c[0] for c in bounds_coords)
        min_lat = min(c[1] for c in bounds_coords)
        max_lat = max(c[1] for c in bounds_coords)
        
        # Approximate width and height in meters (at mid-latitude)
        mid_lat = (min_lat + max_lat) / 2
        lat_meters_per_degree = 111320  # meters per degree latitude
        lon_meters_per_degree = 111320 * np.cos(np.radians(mid_lat))
        
        width_meters = (max_lon - min_lon) * lon_meters_per_degree
        height_meters = (max_lat - min_lat) * lat_meters_per_degree
        
        n_cols = int(np.ceil(width_meters / scale_meters))
        n_rows = int(np.ceil(height_meters / scale_meters))
        
        return n_cols * n_rows
    
    def _create_tile_grid(self, bounds_coords: List, scale_meters: float) -> List[ee.Geometry]:
        """
        Create a grid of tiles that cover the bounding box.
        
        Parameters
        ----------
        bounds_coords : list
            Bounding box coordinates.
        scale_meters : float
            Resolution in meters.
            
        Returns
        -------
        list
            List of ee.Geometry.Rectangle tiles.
        """
        min_lon = min(c[0] for c in bounds_coords)
        max_lon = max(c[0] for c in bounds_coords)
        min_lat = min(c[1] for c in bounds_coords)
        max_lat = max(c[1] for c in bounds_coords)
        
        # Calculate tile size in degrees based on MAX_PIXELS_PER_TILE
        tile_pixels = int(np.sqrt(MAX_PIXELS_PER_TILE))  # e.g., 256 pixels per side
        
        mid_lat = (min_lat + max_lat) / 2
        lat_meters_per_degree = 111320
        lon_meters_per_degree = 111320 * np.cos(np.radians(mid_lat))
        
        # Tile size in degrees
        tile_height_deg = (tile_pixels * scale_meters) / lat_meters_per_degree
        tile_width_deg = (tile_pixels * scale_meters) / lon_meters_per_degree
        
        tiles = []
        lat = min_lat
        while lat < max_lat:
            lon = min_lon
            while lon < max_lon:
                tile_max_lat = min(lat + tile_height_deg, max_lat)
                tile_max_lon = min(lon + tile_width_deg, max_lon)
                
                tile = ee.Geometry.Rectangle([lon, lat, tile_max_lon, tile_max_lat])
                tiles.append(tile)
                
                lon += tile_width_deg
            lat += tile_height_deg
        
        logger.info(f"Created {len(tiles)} tiles for download")
        return tiles
    
    def _download_tile(
        self,
        img: ee.Image,
        tile: ee.Geometry,
        default_value: float = 0
    ) -> Optional[Tuple[np.ndarray, List]]:
        """
        Download a single tile from GEE.
        
        Parameters
        ----------
        img : ee.Image
            Image to download.
        tile : ee.Geometry
            Tile geometry.
        default_value : float
            Default value for missing data.
            
        Returns
        -------
        tuple or None
            Tuple of (array, coordinates) or None if download fails.
        """
        try:
            arr = img.sampleRectangle(
                region=tile,
                defaultValue=default_value
            ).get('pr').getInfo()
            
            if arr is None:
                return None
            
            arr = np.array(arr, dtype=np.float32)
            coords = tile.getInfo()['coordinates'][0]
            
            return arr, coords
            
        except Exception as e:
            logger.warning(f"Failed to download tile: {e}")
            return None
    
    def _mosaic_tiles(
        self,
        tile_files: List[Path],
        output_path: Path,
        cleanup: bool = True
    ) -> bool:
        """
        Mosaic tile files using rioxarray.
        
        Parameters
        ----------
        tile_files : list
            List of tile file paths.
        output_path : Path
            Output mosaic file path.
        cleanup : bool
            Whether to delete tile files after mosaicking.
            
        Returns
        -------
        bool
            True if successful, False otherwise.
        """
        if len(tile_files) == 0:
            return False
        
        if len(tile_files) == 1:
            # Just copy the single tile
            import shutil
            shutil.copy(tile_files[0], output_path)
            if cleanup:
                tile_files[0].unlink()
            return True
        
        try:
            # Read all tiles as DataArrays
            tile_arrays = []
            for f in tile_files:
                da = rioxarray.open_rasterio(f).squeeze('band', drop=True)
                tile_arrays.append(da)
            
            # Merge tiles using rioxarray
            merged = merge_arrays(tile_arrays)
            
            # Save merged result
            merged.rio.to_raster(output_path, compress='LZW')
            
            # Close arrays and cleanup
            for da in tile_arrays:
                da.close()
            
            if cleanup:
                for f in tile_files:
                    f.unlink()
            
            return True
            
        except Exception as e:
            logger.error(f"Mosaicking failed: {e}")
            return False
    
    def _get_monthly_image(self, year: int, month: int) -> ee.Image:
        """
        Get a single monthly image from the collection.
        
        Parameters
        ----------
        year : int
            Year.
        month : int
            Month (1-12).
            
        Returns
        -------
        ee.Image
            Monthly precipitation image (sum of all images in that month).
        """
        monthly_img = (
            self.collection
            .filter(ee.Filter.calendarRange(year, year, 'year'))
            .filter(ee.Filter.calendarRange(month, month, 'month'))
            .sum()  # Sum all images to get monthly total precipitation
        )
        return monthly_img.clip(self.geometry)
    
    def _download_monthly_precip(self, year: int, month: int, temp_dir: Optional[Path] = None) -> Optional[xr.DataArray]:
        """
        Download monthly precipitation data from GEE, using chunked download if needed.
        
        Parameters
        ----------
        year : int
            Year.
        month : int
            Month (1-12).
        temp_dir : Path, optional
            Temporary directory for tile files. If None, uses system temp.
            
        Returns
        -------
        xr.DataArray or None
            Precipitation data array, or None if download fails.
        """
        try:
            img = self._get_monthly_image(year, month)
            region = self.geometry.bounds()
            bounds_coords = region.getInfo()['coordinates'][0]
            
            # Determine the scale to use (native or specified)
            if self.scale is not None:
                scale_meters = self.scale
            else:
                # Use native scale from the dataset
                scale_meters = self._get_native_scale()
            
            # Always reproject to ensure consistent resolution
            img = img.reproject(
                crs='EPSG:4326',
                scale=scale_meters
            )
            
            # Estimate pixel count
            estimated_pixels = self._estimate_pixel_count(bounds_coords, scale_meters)
            logger.debug(f"Estimated pixels: {estimated_pixels}, max allowed: {MAX_PIXELS_PER_TILE}")
            
            # Check if we need chunked download
            if estimated_pixels <= MAX_PIXELS_PER_TILE:
                # Direct download (small region)
                return self._download_single_tile(img, region, year, month)
            else:
                # Chunked download (large region)
                logger.info(f"Large region detected ({estimated_pixels} pixels), using chunked download...")
                return self._download_chunked(img, bounds_coords, scale_meters, year, month, temp_dir)
            
        except Exception as e:
            logger.error(f"Error downloading data for {year}-{month:02d}: {e}")
            return None
    
    def _download_single_tile(
        self,
        img: ee.Image,
        region: ee.Geometry,
        year: int,
        month: int
    ) -> Optional[xr.DataArray]:
        """
        Download a single tile directly (for small regions).
        """
        arr = img.sampleRectangle(
            region=region,
            defaultValue=0
        ).get('pr').getInfo()
        
        if arr is None:
            logger.warning(f"No data for {year}-{month:02d}")
            return None
        
        arr = np.array(arr, dtype=np.float32)
        
        # Get coordinates
        coords = region.getInfo()['coordinates'][0]
        min_lon = min(c[0] for c in coords)
        max_lon = max(c[0] for c in coords)
        min_lat = min(c[1] for c in coords)
        max_lat = max(c[1] for c in coords)
        
        # Create coordinate arrays
        lats = np.linspace(max_lat, min_lat, arr.shape[0])
        lons = np.linspace(min_lon, max_lon, arr.shape[1])
        
        # Create xarray DataArray
        da = xr.DataArray(
            arr,
            dims=['y', 'x'],
            coords={
                'y': lats,
                'x': lons
            },
            attrs={
                'units': 'mm',
                'long_name': 'precipitation',
                'year': year,
                'month': month
            }
        )
        da = da.rio.write_crs("EPSG:4326")
        
        # Save precipitation input if input_dir is set
        if self._input_dir is not None:
            pr_path = self._input_dir / f"precip_{year}_{month:02d}.tif"
            if not pr_path.exists():
                da.rio.to_raster(pr_path)
                logger.info(f"Saved input precipitation: {pr_path.name}")
        
        return da
    
    def _load_awc_data(self, template_da: xr.DataArray) -> np.ndarray:
        """
        Load Available Water Capacity (AWC) data for USDA-SCS method.
        
        Parameters
        ----------
        template_da : xr.DataArray
            Template DataArray to match spatial extent and resolution.
            
        Returns
        -------
        np.ndarray
            AWC values resampled to match template grid.
        """
        if self._awc_cache is not None:
            return self._awc_cache
        
        awc_asset = self.method_params.get('awc_asset')
        awc_band = self.method_params.get('awc_band', 'AWC')
        
        logger.info(f"Loading AWC data from {awc_asset}")
        
        try:
            # Load AWC image (static, not a time series)
            awc_img = ee.Image(awc_asset)
            
            # Select band if specified
            if awc_band:
                awc_img = awc_img.select(awc_band)
            
            # Get the scale to use
            if self.scale is not None:
                scale_meters = self.scale
            else:
                scale_meters = self._get_native_scale()
            
            # Reproject to match template
            awc_img = awc_img.reproject(
                crs='EPSG:4326',
                scale=scale_meters
            ).rename('awc')
            
            # Download AWC data
            region = self.geometry.bounds()
            arr = awc_img.sampleRectangle(
                region=region,
                defaultValue=0.15  # Default AWC if missing
            ).get('awc').getInfo()
            
            if arr is None:
                logger.warning("No AWC data available, using default value of 0.15")
                awc_arr = np.full(template_da.shape, 0.15, dtype=np.float32)
            else:
                awc_arr = np.array(arr, dtype=np.float32)
                
                # Resample to match template if shapes don't match
                if awc_arr.shape != template_da.shape:
                    from scipy.ndimage import zoom
                    zoom_factors = (template_da.shape[0] / awc_arr.shape[0],
                                   template_da.shape[1] / awc_arr.shape[1])
                    awc_arr = zoom(awc_arr, zoom_factors, order=1)
            
            # Cache the result
            self._awc_cache = awc_arr
            
            # Save AWC input if input_dir is set
            if self._input_dir is not None:
                awc_path = self._input_dir / "awc.tif"
                if not awc_path.exists():
                    # Create DataArray for saving - use template coords but ensure shape matches
                    # Resize array if needed to exactly match template coordinates
                    if awc_arr.shape != (len(template_da.coords['y']), len(template_da.coords['x'])):
                        from scipy.ndimage import zoom
                        target_shape = (len(template_da.coords['y']), len(template_da.coords['x']))
                        zoom_factors = (target_shape[0] / awc_arr.shape[0],
                                       target_shape[1] / awc_arr.shape[1])
                        awc_arr = zoom(awc_arr, zoom_factors, order=1)
                    
                    awc_da = xr.DataArray(
                        awc_arr,
                        dims=template_da.dims,
                        coords=template_da.coords,
                        attrs={
                            'units': 'fraction',
                            'long_name': 'available_water_capacity',
                            'source': awc_asset
                        }
                    )
                    awc_da = awc_da.rio.write_crs("EPSG:4326")
                    awc_da.rio.to_raster(awc_path)
                    logger.info(f"Saved input AWC: {awc_path.name}")
            
            return awc_arr
            
        except Exception as e:
            logger.warning(f"Error loading AWC data: {e}. Using default value of 0.15")
            awc_arr = np.full(template_da.shape, 0.15, dtype=np.float32)
            self._awc_cache = awc_arr
            return awc_arr
    
    def _load_monthly_eto(self, year: int, month: int, template_da: xr.DataArray) -> np.ndarray:
        """
        Load Reference Evapotranspiration (ETo) data for USDA-SCS method.
        
        Parameters
        ----------
        year : int
            Year.
        month : int
            Month (1-12).
        template_da : xr.DataArray
            Template DataArray to match spatial extent and resolution.
            
        Returns
        -------
        np.ndarray
            Monthly ETo values in mm.
        """
        eto_asset = self.method_params.get('eto_asset')
        eto_band = self.method_params.get('eto_band', 'eto')
        eto_is_daily = self.method_params.get('eto_is_daily', False)
        eto_scale_factor = self.method_params.get('eto_scale_factor', 1.0)
        
        logger.info(f"Loading ETo data from {eto_asset} for {year}-{month:02d}")
        
        try:
            # Get date range for this month
            import calendar
            _, days_in_month = calendar.monthrange(year, month)
            start_date = f"{year}-{month:02d}-01"
            end_date = f"{year}-{month:02d}-{days_in_month}"
            
            # Load ETo collection
            eto_coll = (
                ee.ImageCollection(eto_asset)
                .select(eto_band)
                .filterDate(start_date, end_date)
                .filterBounds(self.geometry)
            )
            
            # Sum to monthly (for daily data) or just get the image (for monthly)
            if eto_is_daily:
                eto_img = eto_coll.sum().rename('eto')
                logger.debug(f"Aggregated {days_in_month} daily ETo images to monthly")
            else:
                eto_img = eto_coll.first().rename('eto')
            
            # Get the scale to use
            if self.scale is not None:
                scale_meters = self.scale
            else:
                scale_meters = self._get_native_scale()
            
            # Reproject to match template
            eto_img = eto_img.reproject(
                crs='EPSG:4326',
                scale=scale_meters
            )
            
            # Download ETo data
            region = self.geometry.bounds()
            arr = eto_img.sampleRectangle(
                region=region,
                defaultValue=0
            ).get('eto').getInfo()
            
            if arr is None:
                logger.warning(f"No ETo data for {year}-{month:02d}, using default 100 mm")
                eto_arr = np.full(template_da.shape, 100.0, dtype=np.float32)
            else:
                eto_arr = np.array(arr, dtype=np.float32)
                
                # Apply scale factor (e.g., convert units or handle negative values)
                if eto_scale_factor != 1.0:
                    eto_arr = eto_arr * eto_scale_factor
                    logger.debug(f"Applied ETo scale factor: {eto_scale_factor}")
                
                # Resample to match template if shapes don't match
                if eto_arr.shape != template_da.shape:
                    from scipy.ndimage import zoom
                    zoom_factors = (template_da.shape[0] / eto_arr.shape[0],
                                   template_da.shape[1] / eto_arr.shape[1])
                    eto_arr = zoom(eto_arr, zoom_factors, order=1)
            
            # Save ETo input if input_dir is set
            if self._input_dir is not None:
                eto_path = self._input_dir / f"eto_{year}_{month:02d}.tif"
                if not eto_path.exists():
                    # Create DataArray for saving - use template coords but ensure shape matches
                    # Resize array if needed to exactly match template coordinates
                    eto_arr_save = eto_arr
                    if eto_arr_save.shape != (len(template_da.coords['y']), len(template_da.coords['x'])):
                        from scipy.ndimage import zoom
                        target_shape = (len(template_da.coords['y']), len(template_da.coords['x']))
                        zoom_factors = (target_shape[0] / eto_arr_save.shape[0],
                                       target_shape[1] / eto_arr_save.shape[1])
                        eto_arr_save = zoom(eto_arr_save, zoom_factors, order=1)
                    
                    eto_da = xr.DataArray(
                        eto_arr_save,
                        dims=template_da.dims,
                        coords=template_da.coords,
                        attrs={
                            'units': 'mm',
                            'long_name': 'reference_evapotranspiration',
                            'year': year,
                            'month': month,
                            'source': eto_asset
                        }
                    )
                    eto_da = eto_da.rio.write_crs("EPSG:4326")
                    eto_da.rio.to_raster(eto_path)
                    logger.info(f"Saved input ETo: {eto_path.name}")
            
            return eto_arr
            
        except Exception as e:
            logger.warning(f"Error loading ETo data: {e}. Using default value of 100 mm")
            return np.full(template_da.shape, 100.0, dtype=np.float32)

    def _download_chunked(
        self,
        img: ee.Image,
        bounds_coords: List,
        scale_meters: float,
        year: int,
        month: int,
        temp_dir: Optional[Path] = None
    ) -> Optional[xr.DataArray]:
        """
        Download image in chunks and mosaic together.
        
        Parameters
        ----------
        img : ee.Image
            Image to download.
        bounds_coords : list
            Bounding box coordinates.
        scale_meters : float
            Resolution in meters.
        year : int
            Year.
        month : int
            Month.
        temp_dir : Path, optional
            Temporary directory for tiles.
            
        Returns
        -------
        xr.DataArray or None
            Mosaicked data array.
        """
        # Create tiles
        tiles = self._create_tile_grid(bounds_coords, scale_meters)
        
        # Create temp directory for tiles
        if temp_dir is None:
            temp_dir = Path(tempfile.mkdtemp(prefix='pycropwat_'))
        else:
            temp_dir.mkdir(parents=True, exist_ok=True)
        
        tile_files = []
        
        try:
            # Download each tile
            for i, tile in enumerate(tiles):
                logger.debug(f"Downloading tile {i+1}/{len(tiles)}")
                result = self._download_tile(img, tile)
                
                if result is None:
                    continue
                
                arr, coords = result
                
                # Get tile bounds
                min_lon = min(c[0] for c in coords)
                max_lon = max(c[0] for c in coords)
                min_lat = min(c[1] for c in coords)
                max_lat = max(c[1] for c in coords)
                
                # Create coordinate arrays for this tile
                lats = np.linspace(max_lat, min_lat, arr.shape[0])
                lons = np.linspace(min_lon, max_lon, arr.shape[1])
                
                # Create DataArray for this tile
                tile_da = xr.DataArray(
                    arr,
                    dims=['y', 'x'],
                    coords={'y': lats, 'x': lons},
                    attrs={'units': 'mm', 'long_name': 'precipitation'}
                )
                tile_da = tile_da.rio.write_crs("EPSG:4326")
                
                # Save tile to temp file
                tile_path = temp_dir / f"tile_{i:04d}.tif"
                tile_da.rio.to_raster(tile_path)
                tile_files.append(tile_path)
            
            if not tile_files:
                logger.warning(f"No tiles downloaded for {year}-{month:02d}")
                return None
            
            # Mosaic tiles
            mosaic_path = temp_dir / f"mosaic_{year}_{month:02d}.tif"
            success = self._mosaic_tiles(tile_files, mosaic_path, cleanup=True)
            
            if not success:
                logger.error("Failed to mosaic tiles")
                return None
            
            # Read mosaicked file
            da = rioxarray.open_rasterio(mosaic_path).squeeze('band', drop=True)
            da.attrs = {
                'units': 'mm',
                'long_name': 'precipitation',
                'year': year,
                'month': month
            }
            
            # Load into memory and cleanup
            da = da.load()
            mosaic_path.unlink()
            
            # Save precipitation input if input_dir is set
            if self._input_dir is not None:
                pr_path = self._input_dir / f"precip_{year}_{month:02d}.tif"
                if not pr_path.exists():
                    da.rio.to_raster(pr_path)
                    logger.info(f"Saved input precipitation: {pr_path.name}")
            
            return da
            
        finally:
            # Cleanup temp directory if empty
            try:
                if temp_dir.exists() and not any(temp_dir.iterdir()):
                    temp_dir.rmdir()
            except Exception:
                pass
    
    def _process_single_month(
        self,
        year: int,
        month: int,
        output_dir: Path
    ) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Process a single month: download, calculate effective precipitation, save.
        
        Parameters
        ----------
        year : int
            Year.
        month : int
            Month (1-12).
        output_dir : Path
            Output directory.
            
        Returns
        -------
        tuple
            Tuple of (ep_path, epf_path) or (None, None) if processing fails.
        """
        logger.info(f"Processing {year}-{month:02d}")
        
        # Download precipitation data
        pr_da = self._download_monthly_precip(year, month)
        
        if pr_da is None:
            return None, None
        
        # Calculate effective precipitation using the configured method
        if self.method == 'usda_scs':
            # USDA-SCS method requires AWC and ETo data
            awc_arr = self._load_awc_data(pr_da)
            eto_arr = self._load_monthly_eto(year, month, pr_da)
            rooting_depth = self.method_params.get('rooting_depth', 1.0)
            
            ep_arr = self._peff_function(
                pr_da.values,
                eto_arr,
                awc_arr,
                rooting_depth
            )
        elif self.method_params:
            # Other methods with parameters (fixed_percentage, dependable_rainfall)
            # Filter to only pass valid parameters for the method
            valid_params = {}
            if self.method == 'fixed_percentage':
                valid_params['percentage'] = self.method_params.get('percentage', 0.7)
            elif self.method == 'dependable_rainfall':
                valid_params['probability'] = self.method_params.get('probability', 0.75)
            ep_arr = self._peff_function(pr_da.values, **valid_params)
        else:
            ep_arr = self._peff_function(pr_da.values)
        
        # Calculate effective precipitation fraction
        with np.errstate(divide='ignore', invalid='ignore'):
            epf_arr = np.where(pr_da.values > 0, ep_arr / pr_da.values, 0)
        
        # Create effective precipitation DataArray
        ep_da = xr.DataArray(
            ep_arr,
            dims=pr_da.dims,
            coords=pr_da.coords,
            attrs={
                'units': 'mm',
                'long_name': 'effective_precipitation',
                'year': year,
                'month': month,
                'method': self.method.upper()
            }
        )
        ep_da = ep_da.rio.write_crs("EPSG:4326")
        
        # Create effective precipitation fraction DataArray
        epf_da = xr.DataArray(
            epf_arr.astype(np.float32),
            dims=pr_da.dims,
            coords=pr_da.coords,
            attrs={
                'units': 'fraction',
                'long_name': 'effective_precipitation_fraction',
                'year': year,
                'month': month,
                'method': self.method.upper()
            }
        )
        epf_da = epf_da.rio.write_crs("EPSG:4326")
        
        # Save to GeoTIFF
        ep_path = output_dir / f"effective_precip_{year}_{month:02d}.tif"
        epf_path = output_dir / f"effective_precip_fraction_{year}_{month:02d}.tif"
        
        ep_da.rio.to_raster(ep_path)
        epf_da.rio.to_raster(epf_path)
        
        logger.info(f"Saved: {ep_path.name}, {epf_path.name}")
        
        return ep_path, epf_path
    
    def process(
        self,
        output_dir: Union[str, Path],
        n_workers: int = 4,
        months: Optional[List[int]] = None,
        input_dir: Optional[Union[str, Path]] = None,
        save_inputs: bool = False
    ) -> List[Tuple[Optional[Path], Optional[Path]]]:
        """
        Process all months and save effective precipitation rasters.
        
        Uses dask for parallel processing of multiple months.
        
        Parameters
        ----------
        output_dir : str or Path
            Directory to save output rasters.
        n_workers : int, optional
            Number of parallel workers for dask. Default is 4.
        months : list of int, optional
            List of months to process (1-12). If None, processes all months.
        input_dir : str or Path, optional
            Directory to save downloaded input data (precipitation, AWC, ETo).
            If None and save_inputs is True, uses output_dir/../analysis_inputs.
        save_inputs : bool, optional
            Whether to save downloaded input data. Default is False.
            
        Returns
        -------
        list
            List of tuples containing paths to saved files (ep_path, epf_path).
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up input directory for saving downloaded data
        if save_inputs:
            if input_dir is not None:
                self._input_dir = Path(input_dir)
            else:
                # Default: parallel to output_dir in analysis_inputs
                self._input_dir = output_dir.parent / 'analysis_inputs' / output_dir.name
            self._input_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Input data will be saved to: {self._input_dir}")
        else:
            self._input_dir = None
        
        # Generate list of (year, month) to process
        all_dates = get_monthly_dates(self.start_year, self.end_year)
        
        if months is not None:
            all_dates = [(y, m) for y, m in all_dates if m in months]
        
        logger.info(f"Processing {len(all_dates)} months with {n_workers} workers")
        
        # Create delayed tasks
        tasks = [
            delayed(self._process_single_month)(year, month, output_dir)
            for year, month in all_dates
        ]
        
        # Execute in parallel with progress bar
        with ProgressBar():
            results = compute(*tasks, num_workers=n_workers)
        
        return list(results)
    
    def process_sequential(
        self,
        output_dir: Union[str, Path],
        months: Optional[List[int]] = None,
        input_dir: Optional[Union[str, Path]] = None,
        save_inputs: bool = False
    ) -> List[Tuple[Optional[Path], Optional[Path]]]:
        """
        Process all months sequentially (useful for debugging).
        
        Parameters
        ----------
        output_dir : str or Path
            Directory to save output rasters.
        months : list of int, optional
            List of months to process (1-12). If None, processes all months.
        input_dir : str or Path, optional
            Directory to save downloaded input data (precipitation, AWC, ETo).
            If None and save_inputs is True, uses output_dir/../analysis_inputs.
        save_inputs : bool, optional
            Whether to save downloaded input data. Default is False.
            
        Returns
        -------
        list
            List of tuples containing paths to saved files (ep_path, epf_path).
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up input directory for saving downloaded data
        if save_inputs:
            if input_dir is not None:
                self._input_dir = Path(input_dir)
            else:
                # Default: parallel to output_dir in analysis_inputs
                self._input_dir = output_dir.parent / 'analysis_inputs' / output_dir.name
            self._input_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Input data will be saved to: {self._input_dir}")
        else:
            self._input_dir = None
        
        all_dates = get_monthly_dates(self.start_year, self.end_year)
        
        if months is not None:
            all_dates = [(y, m) for y, m in all_dates if m in months]
        
        results = []
        for year, month in all_dates:
            result = self._process_single_month(year, month, output_dir)
            results.append(result)
        
        return results

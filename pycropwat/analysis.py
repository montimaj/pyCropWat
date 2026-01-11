"""
Analysis module for temporal aggregation, statistics, and visualization.

This module provides comprehensive tools for post-processing effective precipitation
data, including temporal aggregation, statistical analysis, and visualization.

Classes
-------
TemporalAggregator
    Aggregate effective precipitation rasters over time periods (seasonal,
    annual, growing season, custom date ranges).
    
StatisticalAnalyzer
    Calculate statistics, trends, anomalies, and spatial summaries from
    effective precipitation time series.
    
Visualizer
    Create publication-quality maps, time series plots, and comparison
    visualizations.

Functions
---------
export_to_netcdf
    Export raster data to NetCDF format with CF-compliant metadata.
    
export_to_cog
    Export raster data to Cloud-Optimized GeoTIFF format.

Example
-------
```python
from pycropwat.analysis import (
    TemporalAggregator,
    StatisticalAnalyzer,
    Visualizer
)

# Temporal aggregation
agg = TemporalAggregator('./output')
annual = agg.annual_aggregate(2020, method='sum')
seasonal = agg.seasonal_aggregate(2020, 'JJA')

# Statistical analysis
stats = StatisticalAnalyzer('./output')
trend = stats.calculate_trend(2010, 2020)
anomaly = stats.calculate_anomaly(2020)

# Visualization
viz = Visualizer()
viz.plot_map(annual, title='Annual Effective Precipitation 2020')
viz.plot_time_series(stats.get_time_series(2010, 2020))
```

Notes
-----
All raster operations preserve geospatial metadata (CRS, transform) and
support both in-memory and file-based workflows.

See Also
--------
pycropwat.core : Core effective precipitation calculations.
pycropwat.methods : Effective precipitation calculation methods.
"""

import logging
import warnings
from pathlib import Path
from typing import Union, Optional, List, Tuple, Literal, Dict
import numpy as np
import xarray as xr
import rioxarray
import pandas as pd

logger = logging.getLogger(__name__)

# Suppress rioxarray nodata warning globally
warnings.filterwarnings('ignore', message='.*nodata.*', category=UserWarning)

# Season definitions (Northern Hemisphere by default)
SEASONS = {
    'DJF': [12, 1, 2],    # Winter
    'MAM': [3, 4, 5],     # Spring
    'JJA': [6, 7, 8],     # Summer
    'SON': [9, 10, 11],   # Fall/Autumn
}

# Aggregation types
AggregationType = Literal['sum', 'mean', 'min', 'max', 'std']


class TemporalAggregator:
    """
    Aggregate effective precipitation rasters over time.
    
    Supports seasonal, annual, growing season, and custom date range
    aggregations.
    
    Parameters
    ----------
    
    input_dir : str or Path
        Directory containing monthly effective precipitation rasters.
    
    pattern : str, optional
        Glob pattern for finding input files. Default is 'effective_precip_[0-9]*.tif'
        which excludes fraction files. Use 'effective_precip_fraction_*.tif' to work
        with fraction files instead.
        
    Examples
    --------
    Basic usage:
    
    ```python
    from pycropwat.analysis import TemporalAggregator
    agg = TemporalAggregator('./output')
    # Annual totals
    annual = agg.annual_aggregate(2020, 'sum')
    # Seasonal aggregation
    seasonal = agg.seasonal_aggregate(2020, 'JJA', 'sum')
    # Growing season (April-October)
    growing = agg.custom_aggregate(2020, months=[4, 5, 6, 7, 8, 9, 10])
    ```
    """
    
    def __init__(
        self,
        input_dir: Union[str, Path],
        pattern: str = 'effective_precip_[0-9]*.tif'
    ):
        self.input_dir = Path(input_dir)
        self.pattern = pattern
        self._files = None
        self._index = None
        
        if not self.input_dir.exists():
            raise FileNotFoundError(f"Input directory not found: {input_dir}")
        
        self._build_file_index()
    
    def _build_file_index(self) -> None:
        """Build an index of available files by year and month."""
        self._files = sorted(self.input_dir.glob(self.pattern))
        self._index = {}
        
        for f in self._files:
            # Parse year and month from filename
            # Expected format: effective_precip_YYYY_MM.tif
            parts = f.stem.split('_')
            try:
                year = int(parts[-2])
                month = int(parts[-1])
                if year not in self._index:
                    self._index[year] = {}
                self._index[year][month] = f
            except (ValueError, IndexError):
                logger.warning(f"Could not parse date from filename: {f.name}")
        
        logger.info(f"Found {len(self._files)} files spanning years {list(self._index.keys())}")
    
    def _load_rasters(
        self,
        year: int,
        months: List[int],
        year_for_month: Optional[Dict[int, int]] = None
    ) -> Optional[xr.DataArray]:
        """
        Load and stack rasters for specified months.
        
        Parameters
        ----------
        
        year : int
            Primary year to load (used for months not in year_for_month).
        
        months : list of int
            Months to load (1-12).
        
        year_for_month : dict, optional
            Mapping of month -> year for cross-year aggregations.
            E.g., {10: 2020, 11: 2020, 12: 2020, 1: 2021, 2: 2021, 3: 2021}
            If None, all months are loaded from the primary year.
            
        Returns
        -------
        xr.DataArray or None
            Stacked DataArray with time dimension, or None if no data.
        """
        arrays = []
        valid_months = []
        time_labels = []
        
        for month in months:
            # Determine which year to use for this month
            if year_for_month and month in year_for_month:
                load_year = year_for_month[month]
            else:
                load_year = year
            
            if load_year not in self._index:
                logger.warning(f"No data available for year {load_year}")
                continue
                
            if month in self._index[load_year]:
                da = rioxarray.open_rasterio(self._index[load_year][month])
                da = da.squeeze('band', drop=True)
                arrays.append(da)
                valid_months.append(month)
                time_labels.append(f"{load_year}-{month:02d}")
            else:
                logger.warning(f"No data for {load_year}-{month:02d}")
        
        if not arrays:
            return None
        
        # Stack along new time dimension
        # Use join='override' to handle slight coordinate differences between files
        stacked = xr.concat(arrays, dim='time', join='override')
        stacked = stacked.assign_coords(time=time_labels)
        
        return stacked
    
    def annual_aggregate(
        self,
        year: int,
        method: AggregationType = 'sum',
        output_path: Optional[Union[str, Path]] = None
    ) -> Optional[xr.DataArray]:
        """
        Calculate annual aggregate of effective precipitation.
        
        Parameters
        ----------
        
        year : int
            Year to aggregate.
        
        method : str, optional
            Aggregation method: 'sum', 'mean', 'min', 'max', 'std'.
            Default is 'sum'.
        
        output_path : str or Path, optional
            Path to save output raster. If None, returns DataArray only.
            
        Returns
        -------
        xr.DataArray or None
            Aggregated data, or None if insufficient data.
        """
        return self.custom_aggregate(
            year,
            months=list(range(1, 13)),
            method=method,
            output_path=output_path,
            output_name=f'annual_{method}'
        )
    
    def seasonal_aggregate(
        self,
        year: int,
        season: str,
        method: AggregationType = 'sum',
        output_path: Optional[Union[str, Path]] = None
    ) -> Optional[xr.DataArray]:
        """
        Calculate seasonal aggregate of effective precipitation.
        
        Parameters
        ----------
        
        year : int
            Year to aggregate. For DJF, December is from the previous year.
        
        season : str
            Season code: 'DJF', 'MAM', 'JJA', or 'SON'.
        
        method : str, optional
            Aggregation method: 'sum', 'mean', 'min', 'max', 'std'.
            Default is 'sum'.
        
        output_path : str or Path, optional
            Path to save output raster.
            
        Returns
        -------
        xr.DataArray or None
            Aggregated data, or None if insufficient data.
        """
        if season not in SEASONS:
            raise ValueError(f"Unknown season '{season}'. Use: {list(SEASONS.keys())}")
        
        months = SEASONS[season]
        
        # For DJF, December is from the previous year
        if season == 'DJF':
            arrays = []
            # December from previous year
            if year - 1 in self._index and 12 in self._index[year - 1]:
                da = rioxarray.open_rasterio(self._index[year - 1][12]).squeeze('band', drop=True)
                arrays.append(da)
            # January and February from current year
            for m in [1, 2]:
                if year in self._index and m in self._index[year]:
                    da = rioxarray.open_rasterio(self._index[year][m]).squeeze('band', drop=True)
                    arrays.append(da)
            
            if not arrays:
                return None
            
            stacked = xr.concat(arrays, dim='time', join='override')
        else:
            stacked = self._load_rasters(year, months)
            if stacked is None:
                return None
        
        # Apply aggregation
        result = self._apply_aggregation(stacked, method)
        result.attrs['season'] = season
        result.attrs['year'] = year
        result.attrs['aggregation'] = method
        
        if output_path:
            result.rio.to_raster(output_path, compress='LZW')
            logger.info(f"Saved: {output_path}")
        
        return result
    
    def custom_aggregate(
        self,
        year: int,
        months: List[int],
        method: AggregationType = 'sum',
        output_path: Optional[Union[str, Path]] = None,
        output_name: Optional[str] = None,
        cross_year: bool = False
    ) -> Optional[xr.DataArray]:
        """
        Calculate custom temporal aggregate of effective precipitation.
        
        Parameters
        ----------
        
        year : int
            Starting year to aggregate. For cross-year aggregations (e.g.,
            Southern Hemisphere growing season Oct-Mar), this is the year
            containing the first months.
        
        months : list of int
            List of months to include (1-12). For cross-year aggregations,
            list months in chronological order (e.g., [10, 11, 12, 1, 2, 3]).
        
        method : str, optional
            Aggregation method: 'sum', 'mean', 'min', 'max', 'std'.
            Default is 'sum'.
        
        output_path : str or Path, optional
            Path to save output raster.
        
        output_name : str, optional
            Name for output attributes.
        
        cross_year : bool, optional
            If True, handles seasons that span two calendar years.
            Months after the "wrap point" (where month number decreases)
            are loaded from year+1. Default is False.
            
        Returns
        -------
        xr.DataArray or None
            Aggregated data, or None if insufficient data.
            
        Examples
        --------
        Northern Hemisphere growing season (same year):
        >>> agg.custom_aggregate(2020, months=[4, 5, 6, 7, 8, 9])
        
        Southern Hemisphere growing season (cross-year, Oct 2020 - Mar 2021):
        >>> agg.custom_aggregate(2020, months=[10, 11, 12, 1, 2, 3], cross_year=True)
        """
        year_for_month = None
        
        if cross_year:
            # Build year mapping for cross-year aggregations
            # Detect where months wrap around (e.g., 12 -> 1)
            year_for_month = {}
            current_year = year
            prev_month = 0
            
            for month in months:
                if month < prev_month:
                    # Month wrapped around to next year
                    current_year = year + 1
                year_for_month[month] = current_year
                prev_month = month
        
        stacked = self._load_rasters(year, months, year_for_month=year_for_month)
        if stacked is None:
            return None
        
        result = self._apply_aggregation(stacked, method)
        result.attrs['year'] = year
        if cross_year:
            result.attrs['end_year'] = year + 1
        result.attrs['months'] = months
        result.attrs['aggregation'] = method
        result.attrs['cross_year'] = cross_year
        if output_name:
            result.attrs['name'] = output_name
        
        if output_path:
            result.rio.to_raster(output_path, compress='LZW')
            logger.info(f"Saved: {output_path}")
        
        return result
    
    def growing_season_aggregate(
        self,
        year: int,
        start_month: int = 4,
        end_month: int = 10,
        method: AggregationType = 'sum',
        output_path: Optional[Union[str, Path]] = None
    ) -> Optional[xr.DataArray]:
        """
        Calculate growing season aggregate.
        
        Automatically handles cross-year seasons (e.g., Southern Hemisphere
        Oct-Mar) when start_month > end_month.
        
        Parameters
        ----------
        
        year : int
            Starting year to aggregate. For cross-year seasons, this is the
            year containing start_month (e.g., 2020 for Oct 2020 - Mar 2021).
        
        start_month : int, optional
            Growing season start month (1-12). Default is 4 (April).
        
        end_month : int, optional
            Growing season end month (1-12). Default is 10 (October).
            If end_month < start_month, assumes cross-year season
            (e.g., start=10, end=3 means Oct-Mar spanning two years).
        
        method : str, optional
            Aggregation method. Default is 'sum'.
        
        output_path : str or Path, optional
            Path to save output raster.
            
        Returns
        -------
        xr.DataArray or None
            Aggregated data.
            
        Examples
        --------
        Northern Hemisphere (Apr-Oct, same year):
        >>> agg.growing_season_aggregate(2020, start_month=4, end_month=10)
        
        Southern Hemisphere (Oct-Mar, cross-year):
        >>> agg.growing_season_aggregate(2020, start_month=10, end_month=3)
        # This aggregates Oct 2020 - Mar 2021
        """
        # Detect cross-year season (start_month > end_month)
        if start_month > end_month:
            # Cross-year season (e.g., Oct-Mar for Southern Hemisphere)
            # Build month list: [10, 11, 12, 1, 2, 3]
            months = list(range(start_month, 13)) + list(range(1, end_month + 1))
            cross_year = True
        else:
            # Same-year season (e.g., Apr-Oct for Northern Hemisphere)
            months = list(range(start_month, end_month + 1))
            cross_year = False
        
        return self.custom_aggregate(
            year,
            months=months,
            method=method,
            output_path=output_path,
            output_name=f'growing_season_{start_month:02d}_{end_month:02d}',
            cross_year=cross_year
        )
    
    def _apply_aggregation(
        self,
        data: xr.DataArray,
        method: AggregationType
    ) -> xr.DataArray:
        """Apply aggregation method along time dimension."""
        if method == 'sum':
            return data.sum(dim='time')
        elif method == 'mean':
            return data.mean(dim='time')
        elif method == 'min':
            return data.min(dim='time')
        elif method == 'max':
            return data.max(dim='time')
        elif method == 'std':
            return data.std(dim='time')
        else:
            raise ValueError(f"Unknown aggregation method: {method}")
    
    def multi_year_climatology(
        self,
        start_year: int,
        end_year: int,
        months: Optional[List[int]] = None,
        output_dir: Optional[Union[str, Path]] = None
    ) -> dict:
        """
        Calculate multi-year climatology (long-term averages).
        
        Parameters
        ----------
        
        start_year : int
            Start year (inclusive).
        
        end_year : int
            End year (inclusive).
        
        months : list of int, optional
            Months to include. If None, calculates for each month.
        
        output_dir : str or Path, optional
            Directory to save output rasters.
            
        Returns
        -------
        dict
            Dictionary mapping month to climatology DataArray.
        """
        if months is None:
            months = list(range(1, 13))
        
        climatology = {}
        
        for month in months:
            arrays = []
            for year in range(start_year, end_year + 1):
                if year in self._index and month in self._index[year]:
                    da = rioxarray.open_rasterio(self._index[year][month])
                    da = da.squeeze('band', drop=True)
                    arrays.append(da)
            
            if arrays:
                stacked = xr.concat(arrays, dim='year', join='override')
                clim = stacked.mean(dim='year')
                clim.attrs = {
                    'long_name': f'climatology_month_{month:02d}',
                    'start_year': start_year,
                    'end_year': end_year,
                    'n_years': len(arrays)
                }
                climatology[month] = clim
                
                if output_dir:
                    output_dir = Path(output_dir)
                    output_dir.mkdir(parents=True, exist_ok=True)
                    output_path = output_dir / f'climatology_{start_year}_{end_year}_month_{month:02d}.tif'
                    clim.rio.to_raster(output_path, compress='LZW')
                    logger.info(f"Saved climatology: {output_path}")
        
        return climatology


class StatisticalAnalyzer:
    """
    Statistical analysis tools for effective precipitation data.
    
    Provides methods for calculating anomalies, trends, and other
    statistical metrics from time series of raster data.
    
    Parameters
    ----------
    
    input_dir : str or Path
        Directory containing monthly effective precipitation rasters.
    
    pattern : str, optional
        Glob pattern for finding input files.
        
    Examples
    --------
    Basic usage:
    
    ```python
    from pycropwat.analysis import StatisticalAnalyzer
    stats = StatisticalAnalyzer('./output')
    # Calculate anomaly
    anomaly = stats.calculate_anomaly(2020, 6, clim_start=1990, clim_end=2020)
    # Calculate trend
    trend, pvalue = stats.calculate_trend(start_year=2000, end_year=2020, month=6)
    ```
    """
    
    def __init__(
        self,
        input_dir: Union[str, Path],
        pattern: str = 'effective_precip_[0-9]*.tif'
    ):
        self.aggregator = TemporalAggregator(input_dir, pattern)
        self.input_dir = Path(input_dir)
        self._index = self.aggregator._index
    
    def calculate_anomaly(
        self,
        year: int,
        month: int,
        clim_start: int,
        clim_end: int,
        anomaly_type: Literal['absolute', 'percent', 'standardized'] = 'absolute',
        output_path: Optional[Union[str, Path]] = None
    ) -> Optional[xr.DataArray]:
        """
        Calculate anomaly relative to climatology.
        
        Parameters
        ----------
        
        year : int
            Year of interest.
        
        month : int
            Month of interest (1-12).
        
        clim_start : int
            Climatology start year.
        
        clim_end : int
            Climatology end year.
        
        anomaly_type : str, optional
            Type of anomaly: 'absolute', 'percent', or 'standardized'.
            Default is 'absolute'.
        
        output_path : str or Path, optional
            Path to save output raster.
            
        Returns
        -------
        xr.DataArray or None
            Anomaly data, or None if insufficient data.
        """
        # Load target data
        if year not in self._index or month not in self._index[year]:
            logger.warning(f"No data for {year}-{month:02d}")
            return None
        
        target = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
        
        # Calculate climatology
        clim_arrays = []
        for y in range(clim_start, clim_end + 1):
            if y in self._index and month in self._index[y]:
                da = rioxarray.open_rasterio(self._index[y][month]).squeeze('band', drop=True)
                clim_arrays.append(da)
        
        if len(clim_arrays) < 3:
            logger.warning(f"Insufficient data for climatology (need at least 3 years)")
            return None
        
        clim_stacked = xr.concat(clim_arrays, dim='year', join='override')
        clim_mean = clim_stacked.mean(dim='year')
        clim_std = clim_stacked.std(dim='year')
        
        # Calculate anomaly
        if anomaly_type == 'absolute':
            anomaly = target - clim_mean
            anomaly.attrs['units'] = 'mm'
        elif anomaly_type == 'percent':
            with np.errstate(divide='ignore', invalid='ignore'):
                anomaly = ((target - clim_mean) / clim_mean) * 100
                anomaly = xr.where(clim_mean > 0, anomaly, np.nan)
            anomaly.attrs['units'] = '%'
        elif anomaly_type == 'standardized':
            with np.errstate(divide='ignore', invalid='ignore'):
                anomaly = (target - clim_mean) / clim_std
                anomaly = xr.where(clim_std > 0, anomaly, np.nan)
            anomaly.attrs['units'] = 'standard deviations'
        else:
            raise ValueError(f"Unknown anomaly type: {anomaly_type}")
        
        anomaly.attrs['long_name'] = f'{anomaly_type}_anomaly'
        anomaly.attrs['year'] = year
        anomaly.attrs['month'] = month
        anomaly.attrs['climatology_period'] = f'{clim_start}-{clim_end}'
        
        if output_path:
            anomaly.rio.to_raster(output_path, compress='LZW')
            logger.info(f"Saved anomaly: {output_path}")
        
        return anomaly
    
    def calculate_trend(
        self,
        start_year: int,
        end_year: int,
        month: Optional[int] = None,
        method: Literal['linear', 'sen'] = 'linear',
        output_dir: Optional[Union[str, Path]] = None
    ) -> Tuple[Optional[xr.DataArray], Optional[xr.DataArray]]:
        """
        Calculate temporal trend in effective precipitation.
        
        Parameters
        ----------
        
        start_year : int
            Start year for trend analysis.
        
        end_year : int
            End year for trend analysis.
        
        month : int, optional
            Specific month to analyze. If None, uses annual totals.
        
        method : str, optional
            Trend method: 'linear' (OLS) or 'sen' (Theil-Sen).
            Default is 'linear'.
        
        output_dir : str or Path, optional
            Directory to save output rasters.
            
        Returns
        -------
        tuple
            (slope, pvalue) DataArrays, or (None, None) if insufficient data.
            Slope units are mm/year.
        """
        # Collect time series
        arrays = []
        years = []
        
        for year in range(start_year, end_year + 1):
            if month is not None:
                # Single month
                if year in self._index and month in self._index[year]:
                    da = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
                    arrays.append(da)
                    years.append(year)
            else:
                # Annual total
                annual = self.aggregator.annual_aggregate(year, 'sum')
                if annual is not None:
                    arrays.append(annual)
                    years.append(year)
        
        if len(arrays) < 5:
            logger.warning(f"Insufficient data for trend analysis (need at least 5 years)")
            return None, None
        
        # Stack into time series
        stacked = xr.concat(arrays, dim='year', join='override')
        stacked = stacked.assign_coords(year=years)
        
        # Calculate trend using vectorized approach
        if method == 'linear':
            slope, pvalue = self._linear_trend(stacked, years)
        elif method == 'sen':
            slope, pvalue = self._sen_trend(stacked, years)
        else:
            raise ValueError(f"Unknown trend method: {method}")
        
        slope.attrs = {
            'long_name': 'trend_slope',
            'units': 'mm/year',
            'method': method,
            'period': f'{start_year}-{end_year}'
        }
        pvalue.attrs = {
            'long_name': 'trend_significance',
            'units': 'p-value',
            'method': method
        }
        
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            suffix = f'month_{month:02d}' if month else 'annual'
            slope.rio.to_raster(
                output_dir / f'trend_slope_{start_year}_{end_year}_{suffix}.tif',
                compress='LZW'
            )
            pvalue.rio.to_raster(
                output_dir / f'trend_pvalue_{start_year}_{end_year}_{suffix}.tif',
                compress='LZW'
            )
            logger.info(f"Saved trend analysis to {output_dir}")
        
        return slope, pvalue
    
    def _linear_trend(
        self,
        data: xr.DataArray,
        years: List[int]
    ) -> Tuple[xr.DataArray, xr.DataArray]:
        """Calculate linear trend using ordinary least squares."""
        from scipy import stats as scipy_stats
        
        years_arr = np.array(years, dtype=np.float64)
        n = len(years)
        
        # Reshape for vectorized computation
        shape = data.shape[1:]  # Spatial shape
        flat_data = data.values.reshape(n, -1)
        
        # Calculate slope and p-value for each pixel
        slopes = np.zeros(flat_data.shape[1])
        pvalues = np.ones(flat_data.shape[1])
        
        for i in range(flat_data.shape[1]):
            y = flat_data[:, i]
            if np.all(np.isfinite(y)):
                result = scipy_stats.linregress(years_arr, y)
                slopes[i] = result.slope
                pvalues[i] = result.pvalue
            else:
                slopes[i] = np.nan
                pvalues[i] = np.nan
        
        # Reshape back to spatial
        slopes = slopes.reshape(shape)
        pvalues = pvalues.reshape(shape)
        
        slope_da = xr.DataArray(
            slopes,
            dims=data.dims[1:],
            coords={k: v for k, v in data.coords.items() if k != 'year'}
        )
        slope_da = slope_da.rio.write_crs(data.rio.crs)
        
        pvalue_da = xr.DataArray(
            pvalues,
            dims=data.dims[1:],
            coords={k: v for k, v in data.coords.items() if k != 'year'}
        )
        pvalue_da = pvalue_da.rio.write_crs(data.rio.crs)
        
        return slope_da, pvalue_da
    
    def _sen_trend(
        self,
        data: xr.DataArray,
        years: List[int]
    ) -> Tuple[xr.DataArray, xr.DataArray]:
        """Calculate Theil-Sen slope with Mann-Kendall significance."""
        from scipy import stats as scipy_stats
        
        years_arr = np.array(years, dtype=np.float64)
        n = len(years)
        
        shape = data.shape[1:]
        flat_data = data.values.reshape(n, -1)
        
        slopes = np.zeros(flat_data.shape[1])
        pvalues = np.ones(flat_data.shape[1])
        
        for i in range(flat_data.shape[1]):
            y = flat_data[:, i]
            if np.all(np.isfinite(y)):
                # Theil-Sen slope
                result = scipy_stats.theilslopes(y, years_arr)
                slopes[i] = result.slope
                
                # Mann-Kendall test for significance
                try:
                    from scipy.stats import kendalltau
                    tau, p = kendalltau(years_arr, y)
                    pvalues[i] = p
                except Exception:
                    pvalues[i] = np.nan
            else:
                slopes[i] = np.nan
                pvalues[i] = np.nan
        
        slopes = slopes.reshape(shape)
        pvalues = pvalues.reshape(shape)
        
        slope_da = xr.DataArray(
            slopes,
            dims=data.dims[1:],
            coords={k: v for k, v in data.coords.items() if k != 'year'}
        )
        slope_da = slope_da.rio.write_crs(data.rio.crs)
        
        pvalue_da = xr.DataArray(
            pvalues,
            dims=data.dims[1:],
            coords={k: v for k, v in data.coords.items() if k != 'year'}
        )
        pvalue_da = pvalue_da.rio.write_crs(data.rio.crs)
        
        return slope_da, pvalue_da
    
    def zonal_statistics(
        self,
        geometry_path: Union[str, Path],
        start_year: int,
        end_year: int,
        months: Optional[List[int]] = None,
        stats: List[str] = ['mean', 'sum', 'min', 'max', 'std'],
        output_path: Optional[Union[str, Path]] = None
    ) -> pd.DataFrame:
        """
        Calculate zonal statistics for polygons.
        
        Parameters
        ----------
        
        geometry_path : str or Path
            Path to shapefile or GeoJSON with zones.
        
        start_year : int
            Start year.
        
        end_year : int
            End year.
        
        months : list of int, optional
            Months to include. If None, includes all available.
        
        stats : list of str, optional
            Statistics to calculate: 'mean', 'sum', 'min', 'max', 'std', 'count'.
        
        output_path : str or Path, optional
            Path to save CSV output.
            
        Returns
        -------
        pd.DataFrame
            DataFrame with zonal statistics.
        """
        import geopandas as gpd
        from rasterstats import zonal_stats
        
        gdf = gpd.read_file(geometry_path)
        if gdf.crs is not None and gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(epsg=4326)
        
        results = []
        
        if months is None:
            months = list(range(1, 13))
        
        for year in range(start_year, end_year + 1):
            for month in months:
                if year not in self._index or month not in self._index[year]:
                    continue
                
                raster_path = self._index[year][month]
                
                # Calculate zonal statistics
                zone_stats = zonal_stats(
                    gdf,
                    str(raster_path),
                    stats=stats
                )
                
                for i, zs in enumerate(zone_stats):
                    # Get zone_id from GeoJSON properties if available, otherwise use index
                    zone_id = gdf.iloc[i].get('zone_id', i) if 'zone_id' in gdf.columns else i
                    row = {
                        'zone_id': zone_id,
                        'year': year,
                        'month': month
                    }
                    row.update(zs)
                    results.append(row)
        
        df = pd.DataFrame(results)
        
        if output_path:
            df.to_csv(output_path, index=False)
            logger.info(f"Saved zonal statistics: {output_path}")
        
        return df


class Visualizer:
    """
    Visualization tools for effective precipitation data.
    
    Provides methods for creating time series plots, maps, and
    comparison visualizations.
    
    Parameters
    ----------
    
    input_dir : str or Path
        Directory containing effective precipitation rasters.
    
    pattern : str, optional
        Glob pattern for finding input files.
    """
    
    def __init__(
        self,
        input_dir: Union[str, Path],
        pattern: str = 'effective_precip_[0-9]*.tif'
    ):
        self.input_dir = Path(input_dir)
        self.aggregator = TemporalAggregator(input_dir, pattern)
        self._index = self.aggregator._index
    
    def plot_time_series(
        self,
        start_year: int,
        end_year: int,
        months: Optional[List[int]] = None,
        geometry_path: Optional[Union[str, Path]] = None,
        stat: str = 'mean',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (12, 6)
    ):
        """
        Plot time series of effective precipitation.
        
        Parameters
        ----------
        
        start_year : int
            Start year.
        
        end_year : int
            End year.
        
        months : list of int, optional
            Months to include. If None, includes all.
        
        geometry_path : str or Path, optional
            Path to geometry for spatial averaging. If None, uses entire raster.
        
        stat : str, optional
            Statistic for spatial aggregation: 'mean', 'sum', 'min', 'max'.
            Default is 'mean'.
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size (width, height) in inches.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        from datetime import datetime
        
        if months is None:
            months = list(range(1, 13))
        
        dates = []
        values = []
        
        for year in range(start_year, end_year + 1):
            for month in months:
                if year not in self._index or month not in self._index[year]:
                    continue
                
                da = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
                
                # Spatial aggregation
                if stat == 'mean':
                    val = float(da.mean())
                elif stat == 'sum':
                    val = float(da.sum())
                elif stat == 'min':
                    val = float(da.min())
                elif stat == 'max':
                    val = float(da.max())
                else:
                    raise ValueError(f"Unknown stat: {stat}")
                
                dates.append(datetime(year, month, 15))
                values.append(val)
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(dates, values, 'b-o', linewidth=1, markersize=3)
        ax.fill_between(dates, values, alpha=0.3)
        
        ax.set_xlabel('Date')
        ax.set_ylabel(f'Effective Precipitation ({stat}) [mm]')
        ax.set_title(title or f'Effective Precipitation Time Series ({start_year}-{end_year})')
        
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved figure: {output_path}")
        
        return fig
    
    def plot_monthly_climatology(
        self,
        start_year: int,
        end_year: int,
        stat: str = 'mean',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (10, 6)
    ):
        """
        Plot monthly climatology bar chart.
        
        Parameters
        ----------
        
        start_year : int
            Climatology start year.
        
        end_year : int
            Climatology end year.
        
        stat : str, optional
            Statistic: 'mean', 'sum'. Default is 'mean'.
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        monthly_means = []
        monthly_stds = []
        
        for month in range(1, 13):
            month_values = []
            for year in range(start_year, end_year + 1):
                if year in self._index and month in self._index[year]:
                    da = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
                    if stat == 'mean':
                        month_values.append(float(da.mean()))
                    elif stat == 'sum':
                        month_values.append(float(da.sum()))
            
            if month_values:
                monthly_means.append(np.mean(month_values))
                monthly_stds.append(np.std(month_values))
            else:
                monthly_means.append(0)
                monthly_stds.append(0)
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        x = np.arange(12)
        
        bars = ax.bar(x, monthly_means, yerr=monthly_stds, capsize=3,
                      color='steelblue', edgecolor='navy', alpha=0.8)
        
        ax.set_xlabel('Month')
        ax.set_ylabel(f'Effective Precipitation ({stat}) [mm]')
        ax.set_title(title or f'Monthly Climatology ({start_year}-{end_year})')
        ax.set_xticks(x)
        ax.set_xticklabels(months)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved figure: {output_path}")
        
        return fig
    
    def plot_raster(
        self,
        year: int,
        month: int,
        cmap: str = 'YlGnBu',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (10, 8),
        vmin: Optional[float] = None,
        vmax: Optional[float] = None
    ):
        """
        Plot a single month's effective precipitation raster.
        
        Parameters
        ----------
        
        year : int
            Year.
        
        month : int
            Month (1-12).
        
        cmap : str, optional
            Colormap name. Default is 'YlGnBu'.
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size.
        
        vmin, vmax : float, optional
            Color scale limits.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        if year not in self._index or month not in self._index[year]:
            raise ValueError(f"No data for {year}-{month:02d}")
        
        da = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
        
        # Mask 0 values as nodata
        da = da.where(da != 0)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        im = da.plot(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
                     cbar_kwargs={'label': 'Effective Precipitation [mm]'})
        
        ax.set_title(title or f'Effective Precipitation - {year}/{month:02d}')
        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved figure: {output_path}")
        
        return fig

    def plot_interactive_map(
        self,
        year: int,
        month: int,
        cmap: str = 'YlGnBu',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        zoom_start: int = 6,
        opacity: float = 0.7
    ):
        """
        Create an interactive map using folium/leafmap.
        
        Parameters
        ----------
        
        year : int
            Year.
        
        month : int
            Month (1-12).
        
        cmap : str, optional
            Colormap name. Default is 'YlGnBu'.
        
        title : str, optional
            Map title.
        
        output_path : str or Path, optional
            Path to save HTML file.
        
        zoom_start : int, optional
            Initial zoom level. Default is 6.
        
        opacity : float, optional
            Raster layer opacity (0-1). Default is 0.7.
            
        Returns
        -------
        folium.Map or leafmap.Map
            The interactive map object.
        """
        try:
            import leafmap.foliumap as leafmap
        except ImportError:
            try:
                import folium
                from folium import raster_layers
                USE_FOLIUM = True
            except ImportError:
                raise ImportError(
                    "Interactive maps require 'leafmap' or 'folium'. "
                    "Install with: pip install leafmap folium"
                )
            USE_FOLIUM = True
        else:
            USE_FOLIUM = False
        
        if year not in self._index or month not in self._index[year]:
            raise ValueError(f"No data for {year}-{month:02d}")
        
        raster_path = self._index[year][month]
        da = rioxarray.open_rasterio(raster_path).squeeze('band', drop=True)
        
        # Mask 0 values as nodata
        da = da.where(da != 0)
        
        # Get bounds for centering map
        bounds = da.rio.bounds()
        center_lat = (bounds[1] + bounds[3]) / 2
        center_lon = (bounds[0] + bounds[2]) / 2
        
        map_title = title or f'Effective Precipitation - {year}/{month:02d}'
        
        if not USE_FOLIUM:
            # Use leafmap for better raster support
            m = leafmap.Map(center=[center_lat, center_lon], zoom=zoom_start)
            m.add_raster(
                str(raster_path),
                colormap=cmap,
                layer_name=map_title,
                opacity=opacity
            )
            # Add colorbar with required colors parameter
            try:
                m.add_colorbar(
                    colors=None,
                    cmap=cmap,
                    vmin=float(da.min()),
                    vmax=float(da.max()),
                    label='Effective Precipitation [mm]'
                )
            except Exception:
                # Some versions of leafmap have different API
                pass
        else:
            # Fallback to folium
            import folium
            m = folium.Map(location=[center_lat, center_lon], zoom_start=zoom_start)
            
            # Add a simple tile layer - full raster overlay requires additional setup
            folium.TileLayer('CartoDB positron').add_to(m)
            
            # Add bounds rectangle to show coverage
            folium.Rectangle(
                bounds=[[bounds[1], bounds[0]], [bounds[3], bounds[2]]],
                color='blue',
                fill=True,
                fill_opacity=0.1,
                popup=f'{map_title}<br>Min: {float(da.min()):.1f} mm<br>Max: {float(da.max()):.1f} mm'
            ).add_to(m)
            
            logger.warning(
                "Full raster visualization requires leafmap. "
                "Install with: pip install leafmap"
            )
        
        if output_path:
            m.save(str(output_path))
            logger.info(f"Saved interactive map: {output_path}")
        
        return m

    def plot_comparison(
        self,
        year: int,
        month: int,
        other_dir: Union[str, Path],
        other_pattern: str = 'effective_precip_*.tif',
        labels: Tuple[str, str] = ('Dataset 1', 'Dataset 2'),
        cmap: str = 'YlGnBu',
        diff_cmap: str = 'RdBu',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (16, 5)
    ):
        """
        Create a comparison plot between two datasets.
        
        Parameters
        ----------
        
        year : int
            Year.
        
        month : int
            Month (1-12).
        
        other_dir : str or Path
            Directory containing the second dataset.
        
        other_pattern : str, optional
            Glob pattern for the second dataset.
        
        labels : tuple of str, optional
            Labels for the two datasets.
        
        cmap : str, optional
            Colormap for the datasets. Default is 'YlGnBu'.
        
        diff_cmap : str, optional
            Colormap for the difference. Default is 'RdBu'.
        
        title : str, optional
            Overall plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        if year not in self._index or month not in self._index[year]:
            raise ValueError(f"No data for {year}-{month:02d} in primary dataset")
        
        # Load primary dataset
        da1 = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
        # Mask 0 values as nodata
        da1 = da1.where(da1 != 0)
        
        # Load secondary dataset
        other_dir = Path(other_dir)
        other_agg = TemporalAggregator(other_dir, other_pattern)
        
        if year not in other_agg._index or month not in other_agg._index[year]:
            raise ValueError(f"No data for {year}-{month:02d} in secondary dataset")
        
        da2 = rioxarray.open_rasterio(other_agg._index[year][month]).squeeze('band', drop=True)
        # Mask 0 values as nodata
        da2 = da2.where(da2 != 0)
        
        # Align grids if necessary
        if da1.shape != da2.shape:
            logger.warning("Datasets have different shapes. Resampling to match...")
            da2 = da2.rio.reproject_match(da1)
        
        # Calculate difference
        diff = da1 - da2
        
        # Determine common color scale
        vmin = min(float(da1.min()), float(da2.min()))
        vmax = max(float(da1.max()), float(da2.max()))
        
        # Symmetric difference scale
        diff_abs_max = max(abs(float(diff.min())), abs(float(diff.max())))
        
        # Create figure
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        # Plot dataset 1
        im1 = da1.plot(ax=axes[0], cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)
        axes[0].set_title(f'{labels[0]}\n{year}/{month:02d}')
        axes[0].set_xlabel('Longitude [°]')
        axes[0].set_ylabel('Latitude [°]')
        plt.colorbar(im1, ax=axes[0], label='Peff [mm]', shrink=0.8)
        
        # Plot dataset 2
        im2 = da2.plot(ax=axes[1], cmap=cmap, vmin=vmin, vmax=vmax, add_colorbar=False)
        axes[1].set_title(f'{labels[1]}\n{year}/{month:02d}')
        axes[1].set_xlabel('Longitude [°]')
        axes[1].set_ylabel('Latitude [°]')
        plt.colorbar(im2, ax=axes[1], label='Peff [mm]', shrink=0.8)
        
        # Plot difference
        im3 = diff.plot(ax=axes[2], cmap=diff_cmap, vmin=-diff_abs_max, vmax=diff_abs_max, add_colorbar=False)
        axes[2].set_title(f'Difference\n({labels[0]} - {labels[1]})')
        axes[2].set_xlabel('Longitude [°]')
        axes[2].set_ylabel('Latitude [°]')
        plt.colorbar(im3, ax=axes[2], label='Δ Peff [mm]', shrink=0.8)
        
        if title:
            fig.suptitle(title, y=1.02, fontsize=14)
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved comparison figure: {output_path}")
        
        return fig

    def plot_scatter_comparison(
        self,
        start_year: int,
        end_year: int,
        other_dir: Union[str, Path],
        other_pattern: str = 'effective_precip_*.tif',
        labels: Tuple[str, str] = ('Dataset 1', 'Dataset 2'),
        months: Optional[List[int]] = None,
        sample_size: int = 10000,
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (8, 8)
    ):
        """
        Create a scatter plot comparing two datasets.
        
        Parameters
        ----------
        
        start_year : int
            Start year.
        
        end_year : int
            End year.
        
        other_dir : str or Path
            Directory containing the second dataset.
        
        other_pattern : str, optional
            Glob pattern for the second dataset.
        
        labels : tuple of str, optional
            Labels for the two datasets.
        
        months : list of int, optional
            Months to include. If None, includes all.
        
        sample_size : int, optional
            Number of random pixels to sample for scatter plot.
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure with scatter plot and statistics.
        """
        import matplotlib.pyplot as plt
        from scipy import stats
        
        if months is None:
            months = list(range(1, 13))
        
        other_dir = Path(other_dir)
        other_agg = TemporalAggregator(other_dir, other_pattern)
        
        all_vals1 = []
        all_vals2 = []
        
        for year in range(start_year, end_year + 1):
            for month in months:
                if (year not in self._index or month not in self._index[year] or
                    year not in other_agg._index or month not in other_agg._index[year]):
                    continue
                
                da1 = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
                da2 = rioxarray.open_rasterio(other_agg._index[year][month]).squeeze('band', drop=True)
                
                # Align if necessary
                if da1.shape != da2.shape:
                    da2 = da2.rio.reproject_match(da1)
                
                # Flatten and remove NaN
                v1 = da1.values.flatten()
                v2 = da2.values.flatten()
                
                mask = ~(np.isnan(v1) | np.isnan(v2))
                all_vals1.extend(v1[mask])
                all_vals2.extend(v2[mask])
        
        all_vals1 = np.array(all_vals1)
        all_vals2 = np.array(all_vals2)
        
        # Random sampling if too many points
        if len(all_vals1) > sample_size:
            idx = np.random.choice(len(all_vals1), sample_size, replace=False)
            plot_vals1 = all_vals1[idx]
            plot_vals2 = all_vals2[idx]
        else:
            plot_vals1 = all_vals1
            plot_vals2 = all_vals2
        
        # Calculate statistics
        r, p_value = stats.pearsonr(all_vals1, all_vals2)
        rmse = np.sqrt(np.mean((all_vals1 - all_vals2) ** 2))
        bias = np.mean(all_vals1 - all_vals2)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        ax.scatter(plot_vals1, plot_vals2, alpha=0.3, s=5, c='steelblue')
        
        # 1:1 line
        max_val = max(all_vals1.max(), all_vals2.max())
        min_val = min(all_vals1.min(), all_vals2.min())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=2, label='1:1 line')
        
        # Regression line
        slope, intercept, _, _, _ = stats.linregress(all_vals1, all_vals2)
        ax.plot([min_val, max_val], 
                [slope * min_val + intercept, slope * max_val + intercept],
                'r-', linewidth=2, label=f'Fit: y={slope:.2f}x+{intercept:.1f}')
        
        ax.set_xlabel(f'{labels[0]} Effective Precipitation [mm]')
        ax.set_ylabel(f'{labels[1]} Effective Precipitation [mm]')
        ax.set_title(title or f'Dataset Comparison ({start_year}-{end_year})')
        
        # Statistics text box
        stats_text = f'R² = {r**2:.3f}\nRMSE = {rmse:.2f} mm\nBias = {bias:.2f} mm\nn = {len(all_vals1):,}'
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        ax.legend(loc='lower right')
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved scatter comparison: {output_path}")
        
        return fig

    def plot_annual_comparison(
        self,
        start_year: int,
        end_year: int,
        other_dir: Union[str, Path],
        other_pattern: str = 'effective_precip_*.tif',
        labels: Tuple[str, str] = ('Dataset 1', 'Dataset 2'),
        stat: str = 'mean',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (12, 6)
    ):
        """
        Compare annual totals between two datasets.
        
        Parameters
        ----------
        
        start_year : int
            Start year.
        
        end_year : int
            End year.
        
        other_dir : str or Path
            Directory containing the second dataset.
        
        other_pattern : str, optional
            Glob pattern for the second dataset.
        
        labels : tuple of str, optional
            Labels for the two datasets.
        
        stat : str, optional
            Spatial statistic: 'mean', 'sum'. Default is 'mean'.
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        other_dir = Path(other_dir)
        other_agg = TemporalAggregator(other_dir, other_pattern)
        
        years = list(range(start_year, end_year + 1))
        vals1 = []
        vals2 = []
        
        for year in years:
            annual1 = 0
            annual2 = 0
            count1 = 0
            count2 = 0
            
            for month in range(1, 13):
                if year in self._index and month in self._index[year]:
                    da = rioxarray.open_rasterio(self._index[year][month]).squeeze('band', drop=True)
                    if stat == 'mean':
                        annual1 += float(da.mean())
                    else:
                        annual1 += float(da.sum())
                    count1 += 1
                
                if year in other_agg._index and month in other_agg._index[year]:
                    da = rioxarray.open_rasterio(other_agg._index[year][month]).squeeze('band', drop=True)
                    if stat == 'mean':
                        annual2 += float(da.mean())
                    else:
                        annual2 += float(da.sum())
                    count2 += 1
            
            vals1.append(annual1 if count1 == 12 else np.nan)
            vals2.append(annual2 if count2 == 12 else np.nan)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        x = np.arange(len(years))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, vals1, width, label=labels[0], color='steelblue', alpha=0.8)
        bars2 = ax.bar(x + width/2, vals2, width, label=labels[1], color='coral', alpha=0.8)
        
        ax.set_xlabel('Year')
        ax.set_ylabel(f'Annual Effective Precipitation ({stat}) [mm]')
        ax.set_title(title or f'Annual Comparison ({start_year}-{end_year})')
        ax.set_xticks(x)
        ax.set_xticklabels(years, rotation=45)
        ax.legend()
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved annual comparison: {output_path}")
        
        return fig

    def plot_anomaly_map(
        self,
        anomaly_path: Union[str, Path],
        cmap: str = 'RdBu',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (10, 8),
        center_value: float = 100.0
    ):
        """
        Plot an anomaly map from a GeoTIFF file.
        
        Uses a diverging colormap centered at the specified value
        (default 100 for percent anomalies).
        
        Parameters
        ----------
        
        anomaly_path : str or Path
            Path to the anomaly GeoTIFF file.
        
        cmap : str, optional
            Colormap name. Default is 'RdBu' (red=dry, blue=wet).
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size (width, height) in inches.
        
        center_value : float, optional
            Center value for the diverging colormap. Default is 100 (for percent).
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        anomaly_path = Path(anomaly_path)
        if not anomaly_path.exists():
            raise FileNotFoundError(f"Anomaly file not found: {anomaly_path}")
        
        da = rioxarray.open_rasterio(anomaly_path).squeeze('band', drop=True)
        da = da.where(da != 0)  # Mask nodata
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Calculate symmetric limits around center value
        vmax = max(
            abs(float(da.min()) - center_value), 
            abs(float(da.max()) - center_value), 
            50
        )
        
        im = da.plot(
            ax=ax,
            cmap=cmap,
            vmin=center_value - vmax,
            vmax=center_value + vmax,
            cbar_kwargs={'label': 'Anomaly [% of Climatology]'}
        )
        
        ax.set_title(title or f'Effective Precipitation Anomaly')
        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved anomaly map: {output_path}")
        
        return fig

    def plot_climatology_map(
        self,
        climatology_path: Union[str, Path],
        cmap: str = 'YlGnBu',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (10, 8),
        vmin: Optional[float] = None,
        vmax: Optional[float] = None
    ):
        """
        Plot a climatology map from a GeoTIFF file.
        
        Parameters
        ----------
        
        climatology_path : str or Path
            Path to the climatology GeoTIFF file.
        
        cmap : str, optional
            Colormap name. Default is 'YlGnBu'.
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size (width, height) in inches.
        
        vmin, vmax : float, optional
            Color scale limits.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        climatology_path = Path(climatology_path)
        if not climatology_path.exists():
            raise FileNotFoundError(f"Climatology file not found: {climatology_path}")
        
        da = rioxarray.open_rasterio(climatology_path).squeeze('band', drop=True)
        da = da.where(da != 0)  # Mask nodata
        
        fig, ax = plt.subplots(figsize=figsize)
        
        im = da.plot(
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            cbar_kwargs={'label': 'Climatology [mm]'}
        )
        
        ax.set_title(title or f'Effective Precipitation Climatology')
        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved climatology map: {output_path}")
        
        return fig

    def plot_trend_map(
        self,
        slope_path: Union[str, Path],
        pvalue_path: Optional[Union[str, Path]] = None,
        cmap: str = 'RdBu',
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (10, 8),
        show_significance: bool = True,
        significance_level: float = 0.05
    ):
        """
        Plot a trend map from slope GeoTIFF, optionally with significance overlay.
        
        Parameters
        ----------
        
        slope_path : str or Path
            Path to the slope (trend) GeoTIFF file.
        
        pvalue_path : str or Path, optional
            Path to the p-value GeoTIFF file for significance overlay.
        
        cmap : str, optional
            Colormap name. Default is 'RdBu' (red=decreasing, blue=increasing).
        
        title : str, optional
            Plot title.
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size (width, height) in inches.
        
        show_significance : bool, optional
            Whether to overlay significance stippling. Default is True.
        
        significance_level : float, optional
            P-value threshold for significance. Default is 0.05.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        slope_path = Path(slope_path)
        if not slope_path.exists():
            raise FileNotFoundError(f"Slope file not found: {slope_path}")
        
        da_slope = rioxarray.open_rasterio(slope_path).squeeze('band', drop=True)
        da_slope = da_slope.where(da_slope != 0)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Symmetric color limits centered at 0
        vmax_slope = max(abs(float(da_slope.min())), abs(float(da_slope.max())), 0.5)
        
        im = da_slope.plot(
            ax=ax,
            cmap=cmap,
            vmin=-vmax_slope,
            vmax=vmax_slope,
            cbar_kwargs={'label': 'Trend Slope [mm/month/year]'}
        )
        
        # Add significance stippling if p-value file provided
        if show_significance and pvalue_path:
            pvalue_path = Path(pvalue_path)
            if pvalue_path.exists():
                da_pvalue = rioxarray.open_rasterio(pvalue_path).squeeze('band', drop=True)
                significant = (da_pvalue < significance_level) & (da_pvalue != 0)
                
                if significant.sum() > 0:
                    # Subsample for stippling
                    step = max(1, min(da_pvalue.shape) // 50)
                    sig_y, sig_x = np.where(significant.values[::step, ::step])
                    
                    if len(sig_x) > 0:
                        x_coords = da_pvalue.x.values[::step][sig_x]
                        y_coords = da_pvalue.y.values[::step][sig_y]
                        
                        ax.scatter(x_coords, y_coords, marker='.', s=1, c='black', 
                                  alpha=0.3, label=f'p < {significance_level}')
                        ax.legend(loc='lower right')
        
        ax.set_title(title or 'Effective Precipitation Trend')
        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved trend map: {output_path}")
        
        return fig

    def plot_trend_panel(
        self,
        slope_path: Union[str, Path],
        pvalue_path: Union[str, Path],
        title: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: Tuple[int, int] = (16, 6)
    ):
        """
        Plot a two-panel figure with slope and p-value maps side by side.
        
        Parameters
        ----------
        
        slope_path : str or Path
            Path to the slope (trend) GeoTIFF file.
        
        pvalue_path : str or Path
            Path to the p-value GeoTIFF file.
        
        title : str, optional
            Overall figure title (used for slope panel).
        
        output_path : str or Path, optional
            Path to save figure.
        
        figsize : tuple, optional
            Figure size (width, height) in inches.
            
        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        import matplotlib.pyplot as plt
        
        slope_path = Path(slope_path)
        pvalue_path = Path(pvalue_path)
        
        if not slope_path.exists():
            raise FileNotFoundError(f"Slope file not found: {slope_path}")
        if not pvalue_path.exists():
            raise FileNotFoundError(f"P-value file not found: {pvalue_path}")
        
        da_slope = rioxarray.open_rasterio(slope_path).squeeze('band', drop=True)
        da_slope = da_slope.where(da_slope != 0)
        
        da_pvalue = rioxarray.open_rasterio(pvalue_path).squeeze('band', drop=True)
        da_pvalue = da_pvalue.where(da_pvalue != 0)
        
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        # Plot 1: Slope (trend magnitude)
        vmax_slope = max(abs(float(da_slope.min())), abs(float(da_slope.max())), 0.5)
        
        da_slope.plot(
            ax=axes[0],
            cmap='RdBu',
            vmin=-vmax_slope,
            vmax=vmax_slope,
            cbar_kwargs={'label': 'Trend Slope [mm/month/year]'}
        )
        axes[0].set_title(title or 'Effective Precipitation Trend')
        axes[0].set_xlabel('Longitude [°]')
        axes[0].set_ylabel('Latitude [°]')
        
        # Plot 2: P-value (significance)
        da_pvalue.plot(
            ax=axes[1],
            cmap='RdYlGn_r',  # Green for low p-values (significant)
            vmin=0,
            vmax=0.1,
            cbar_kwargs={'label': 'P-value'}
        )
        axes[1].set_title('Trend Significance (p < 0.05 = significant)')
        axes[1].set_xlabel('Longitude [°]')
        axes[1].set_ylabel('Latitude [°]')
        
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved trend panel: {output_path}")
        
        return fig


def export_to_netcdf(
    input_dir: Union[str, Path],
    output_path: Union[str, Path],
    pattern: str = 'effective_precip_*.tif',
    variable_name: str = 'effective_precipitation',
    compression: bool = True
) -> None:
    """
    Export monthly rasters to a single NetCDF file.
    
    Parameters
    ----------
    
    input_dir : str or Path
        Directory containing monthly rasters.
    
    output_path : str or Path
        Output NetCDF file path.
    
    pattern : str, optional
        Glob pattern for finding input files.
    
    variable_name : str, optional
        Name for the variable in NetCDF.
    
    compression : bool, optional
        Whether to compress the output. Default is True.
    """
    import pandas as pd
    
    input_dir = Path(input_dir)
    files = sorted(input_dir.glob(pattern))
    
    if not files:
        raise FileNotFoundError(f"No files found matching {pattern} in {input_dir}")
    
    arrays = []
    times = []
    
    for f in files:
        # Parse date from filename
        parts = f.stem.split('_')
        try:
            year = int(parts[-2])
            month = int(parts[-1])
            times.append(pd.Timestamp(year=year, month=month, day=1))
            
            da = rioxarray.open_rasterio(f).squeeze('band', drop=True)
            # Set nodata explicitly to avoid warning
            if da.rio.nodata is None:
                da = da.rio.write_nodata(0)
            arrays.append(da)
        except (ValueError, IndexError):
            logger.warning(f"Could not parse date from: {f.name}")
    
    if not arrays:
        raise ValueError("No valid files found")
    
    # Stack into time series
    stacked = xr.concat(arrays, dim='time', join='override')
    stacked = stacked.assign_coords(time=times)
    stacked.name = variable_name
    stacked.attrs['units'] = 'mm'
    stacked.attrs['long_name'] = 'Monthly effective precipitation'
    
    # Create dataset
    ds = stacked.to_dataset()
    
    # Set encoding - try netCDF4 engine with compression first
    encoding = {}
    engine = None
    if compression:
        try:
            import netCDF4  # noqa: F401
            encoding[variable_name] = {
                'zlib': True,
                'complevel': 4
            }
            engine = 'netcdf4'
        except ImportError:
            # Fall back to scipy without compression
            logger.warning("netCDF4 not installed, saving without compression")
            engine = 'scipy'
    
    ds.to_netcdf(output_path, encoding=encoding if engine == 'netcdf4' else None, engine=engine)
    logger.info(f"Exported {len(arrays)} time steps to {output_path}")


def export_to_cog(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    overview_levels: List[int] = [2, 4, 8, 16]
) -> None:
    """
    Convert a GeoTIFF to Cloud-Optimized GeoTIFF (COG).
    
    Parameters
    ----------
    
    input_path : str or Path
        Input GeoTIFF path.
    
    output_path : str or Path
        Output COG path.
    
    overview_levels : list of int, optional
        Overview levels to create.
    """
    import rasterio
    from rasterio.enums import Resampling
    
    with rasterio.open(input_path) as src:
        profile = src.profile.copy()
        
        # Update profile for COG
        profile.update(
            driver='GTiff',
            tiled=True,
            blockxsize=512,
            blockysize=512,
            compress='LZW',
            interleave='band'
        )
        
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(src.read())
            
            # Build overviews
            dst.build_overviews(overview_levels, Resampling.average)
            dst.update_tags(ns='rio_overview', resampling='average')
    
    logger.info(f"Created COG: {output_path}")

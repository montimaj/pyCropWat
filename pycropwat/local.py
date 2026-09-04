"""
Local precipitation data sources for pyCropWat.

This module lets you drive ``pycropwat.EffectivePrecipitation`` with your **own**
precipitation rasters instead of downloading precipitation from Google Earth Engine.
Two on-disk layouts are supported:

* **Raster mode** - a directory (or glob) of one GeoTIFF per month, where the file name
  carries the year and month, e.g. ``Precip_2005_07.tif``.
* **NetCDF mode** - one or more NetCDF/HDF5 files holding a precipitation variable with a
  ``time`` coordinate (or a single month per file, dated by file name).

Precipitation is returned as a 2-D ``xarray.DataArray`` with dims ``('y', 'x')``,
``float32`` values in millimetres, nodata as ``NaN``, a north-up ``y`` coordinate and a CRS
attached via ``rioxarray``. That is exactly what the rest of pyCropWat expects, so the
precipitation-only methods (``cropwat``, ``fao_aglw``, ``fixed_percentage``,
``dependable_rainfall``, ``farmwest``) run with no Earth Engine involvement at all.

Classes
-------
LocalPrecipitationSource
    Indexed, lazily-opened reader for local monthly precipitation data.

Functions
---------
parse_year_month
    Extract ``(year, month)`` from a file name or arbitrary string.

open_local_precipitation
    Convenience factory returning a ``LocalPrecipitationSource``.

Example
-------
```python
from pycropwat.local import LocalPrecipitationSource

# A directory of monthly WRF GeoTIFFs named Precip_YYYY_MM.tif
src = LocalPrecipitationSource(
    '../pyCropWat_Data/Precip',
    pattern='Precip_*.tif',
    nodata=-9999,
)
print(src.kind, len(src), src.year_range, src.crs, src.shape)

da = src.get_month(2005, 7)     # xr.DataArray, dims ('y', 'x'), mm
print(float(da.mean()))
```

NetCDF, with an explicit variable and a unit conversion (m -> mm):

```python
from pycropwat.local import open_local_precipitation

with open_local_precipitation('precip_2000_2020.nc',
                              variable='tp',
                              scale_factor=1000.0) as src:
    for year, month in src.available_months():
        da = src.get_month(year, month)
```

See Also
--------
pycropwat.core : Effective precipitation workflow that consumes these arrays.
pycropwat.methods : Precipitation-only effective precipitation formulas.
"""

import glob as _glob
import logging
import re
from pathlib import Path
from typing import Union, Optional, List, Tuple, Dict, Any, Sequence

import numpy as np
import pandas as pd
import xarray as xr
import rioxarray
import rasterio
from rasterio.crs import CRS
from rasterio.transform import array_bounds
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)


#: File suffixes treated as single-month rasters.
RASTER_SUFFIXES = (
    '.tif', '.tiff', '.vrt', '.img', '.bil', '.bsq', '.bip',
    '.asc', '.jp2', '.grd', '.dat', '.hgt'
)

#: File suffixes treated as NetCDF/HDF stacks.
NETCDF_SUFFIXES = ('.nc', '.nc4', '.cdf', '.netcdf', '.h5', '.hdf5')

#: Variable names searched (case-insensitively) when auto-detecting NetCDF precipitation.
PRECIP_VARIABLE_NAMES = (
    'precip', 'precipitation', 'pr', 'prcp', 'ppt', 'tp',
    'total_precipitation', 'rain', 'rainfall', 'RAINNC', 'PRECIP'
)

#: Candidate names for the X (longitude / easting) dimension.
X_DIM_NAMES = ('x', 'lon', 'longitude', 'west_east', 'XLONG')

#: Candidate names for the Y (latitude / northing) dimension.
Y_DIM_NAMES = ('y', 'lat', 'latitude', 'south_north', 'XLAT')

#: Candidate names for the time coordinate.
TIME_DIM_NAMES = ('time', 'Time', 'valid_time', 't')

# Variables that are never the precipitation variable.
_NON_DATA_VARS = frozenset({
    'crs', 'spatial_ref', 'lat', 'lon', 'latitude', 'longitude', 'xlat', 'xlong',
    'time_bnds', 'time_bounds', 'lat_bnds', 'lon_bnds', 'x_bnds', 'y_bnds',
    'transverse_mercator', 'lambert_conformal_conic', 'polar_stereographic',
    'albers_conical_equal_area', 'latitude_longitude'
})

# Plausible calendar-year window used to validate parsed dates.
_MIN_YEAR = 1700
_MAX_YEAR = 2200

# Ordered most-specific-first. Every pattern exposes named groups 'year' and 'month'.
_DATE_PATTERNS = (
    re.compile(r'(?<!\d)(?P<year>\d{4})_(?P<month>\d{1,2})(?!\d)'),    # YYYY_MM
    re.compile(r'(?<!\d)(?P<year>\d{4})-(?P<month>\d{1,2})(?!\d)'),    # YYYY-MM
    re.compile(r'(?<!\d)(?P<year>\d{4})\.(?P<month>\d{1,2})(?!\d)'),   # YYYY.MM
    re.compile(r'(?<!\d)(?P<year>\d{4})(?P<month>\d{2})(?!\d)'),       # YYYYMM
)

# Used only to detect the "one multi-band raster per year" layout we deliberately reject.
_YEAR_ONLY_PATTERN = re.compile(r'(?<!\d)(?P<year>\d{4})(?!\d)')

_GLOB_CHARS = ('*', '?', '[')

# Affine coefficients of two files count as the same grid when they agree to within
# this fraction of a pixel.
_TRANSFORM_PIXEL_TOL = 1e-6

# Collapsing a 2-D curvilinear coordinate to 1-D is reported at WARNING once the
# worst-placed pixel moves by more than this many cells ...
_CURVILINEAR_WARN_CELLS = 0.5

# ... and rejected outright once it moves by more than this fraction of the axis span.
_CURVILINEAR_MAX_REL_DEV = 0.2


def parse_year_month(
    name: Union[str, Path],
    date_regex: Optional[str] = None
) -> Optional[Tuple[int, int]]:
    """
    Extract a ``(year, month)`` pair from a file name or arbitrary string.

    Layouts are tried most-specific-first and the **last** match in the string wins, so
    directory noise and prefixes such as ``Precip_`` or ``effective_precip_`` never
    confuse the parser. Recognised layouts:

    | Layout    | Example                    |
    |-----------|----------------------------|
    | ``YYYY_MM`` | ``Precip_2005_07.tif``   |
    | ``YYYY-MM`` | ``precip-2005-07.nc``    |
    | ``YYYYMM``  | ``pr200507.tif``         |
    | ``YYYY.MM`` | ``pr.2005.07.tif``       |

    Parameters
    ----------

    name : str or Path
        File name, file stem, full path or any other string to parse.

    date_regex : str, optional
        Custom regular expression exposing named groups ``year`` and ``month``.
        When given, the built-in layouts are not used.

    Returns
    -------
    tuple or None
        ``(year, month)`` with ``1 <= month <= 12`` and a plausible year
        (1700-2200), or ``None`` when nothing matched.

    Raises
    ------
    ValueError
        If ``date_regex`` is not a valid regular expression or lacks the
        ``year`` / ``month`` named groups.

    Examples
    --------
    ```python
    from pycropwat.local import parse_year_month

    parse_year_month('Precip_2005_07.tif')          # (2005, 7)
    parse_year_month('pr200507')                    # (2005, 7)
    parse_year_month('/data/2019/x-2005-07.nc')     # (2005, 7)
    parse_year_month('no_date_here.tif')            # None

    # Custom layout: MM_YYYY
    parse_year_month('rain_07_2005.tif',
                     date_regex=r'(?P<month>\\d{2})_(?P<year>\\d{4})')  # (2005, 7)
    ```
    """
    text = str(name)

    if date_regex is not None:
        try:
            pattern = re.compile(date_regex)
        except re.error as exc:
            raise ValueError(f"Invalid date_regex {date_regex!r}: {exc}") from exc
        group_names = pattern.groupindex
        if 'year' not in group_names or 'month' not in group_names:
            raise ValueError(
                f"date_regex must define named groups 'year' and 'month'. "
                f"Got groups: {sorted(group_names)}. "
                r"Example: r'(?P<year>\d{4})_(?P<month>\d{2})'"
            )
        patterns = (pattern,)
    else:
        patterns = _DATE_PATTERNS

    for pattern in patterns:
        matches = list(pattern.finditer(text))
        # Anchor on the LAST match, but fall back to earlier ones if it is implausible.
        for match in reversed(matches):
            try:
                year = int(match.group('year'))
                month = int(match.group('month'))
            except (TypeError, ValueError):
                continue
            if 1 <= month <= 12 and _MIN_YEAR <= year <= _MAX_YEAR:
                return year, month

    return None


def _is_glob(text: str) -> bool:
    """Return True when ``text`` contains shell glob wildcards."""
    return any(char in text for char in _GLOB_CHARS)


def _classify_suffix(path: Path) -> Optional[str]:
    """Return ``'raster'``, ``'netcdf'`` or None for a path based on its suffix."""
    suffix = path.suffix.lower()
    if suffix in NETCDF_SUFFIXES:
        return 'netcdf'
    if suffix in RASTER_SUFFIXES:
        return 'raster'
    return None


def _times_to_year_month(values: Any) -> List[Tuple[int, int]]:
    """
    Convert an array of datetimes to a list of ``(year, month)`` tuples.

    Handles ``numpy.datetime64``, ``datetime.datetime`` and ``cftime`` objects.

    Parameters
    ----------

    values : array-like
        Decoded time values.

    Returns
    -------
    list
        List of ``(year, month)`` integer tuples.
    """
    array = np.atleast_1d(np.asarray(values))

    try:
        index = pd.to_datetime(array)
        return [(int(stamp.year), int(stamp.month)) for stamp in np.atleast_1d(index)]
    except Exception:
        # cftime (non-standard calendars) and other duck-typed datetimes.
        pairs = []
        for value in array.ravel().tolist():
            if hasattr(value, 'year') and hasattr(value, 'month'):
                pairs.append((int(value.year), int(value.month)))
            else:
                stamp = pd.Timestamp(value)
                pairs.append((int(stamp.year), int(stamp.month)))
        return pairs


class LocalPrecipitationSource:
    """
    Reader for user-supplied local monthly precipitation data.

    Wraps a directory of single-month rasters, or one/many NetCDF files, and exposes
    them as analysis-ready monthly ``xarray.DataArray`` objects. Files are indexed
    up front (cheap) and pixel data is only read when ``get_month`` is called.

    Parameters
    ----------

    path : str or Path
        A directory of rasters, a single raster/NetCDF file, or a glob string
        (e.g. ``'/data/precip/Precip_*.tif'``).

    pattern : str, default '*.tif'
        Glob applied when ``path`` is a directory. Files whose suffix is not a
        recognised raster/NetCDF suffix are ignored, so stray companions such as
        ``.aux.xml``, ``metadata.csv`` or ``.geeup-state.json`` never break indexing.

    variable : str, optional
        NetCDF variable holding precipitation. Auto-detected when None.

    scale_factor : float, default 1.0
        Constant multiplier applied to every returned array, for converting the source
        units to millimetres per month: ``1000`` for metres, ``10`` for centimetres,
        ``25.4`` for inches, ``0.1`` for tenths of a millimetre. A single multiplier
        cannot convert a *rate* such as kg m-2 s-1 to a monthly total, because the
        factor depends on the number of days in each month - accumulate such data to
        monthly totals before handing it to pyCropWat.

    nodata : float or sequence of float, optional
        Extra nodata sentinel(s) to mask, **in addition** to the nodata declared in the
        file metadata. Compared with ``numpy.isclose`` and NaN-safe, e.g. ``-9999``.

    crs : str or CRS, optional
        CRS to attach to the returned arrays (e.g. ``'EPSG:4326'``). It **overrides**
        whatever the file declares, and logs at INFO when it replaces a declared CRS.
        Overriding *relabels* the grid only: the pixels and their coordinates are used
        exactly as they are on disk and nothing is reprojected, so this is the fix for
        missing or wrong CRS metadata, not a way to change projection. When None, the
        file's own CRS is used, falling back to EPSG:4326 with a warning. Passing it
        does **not** switch off the mixed-CRS check: raster mode still compares the
        files' declared CRSs with each other, and reports a directory whose files
        disagree - see the grid note below.

    date_regex : str, optional
        Custom regex with named groups ``year`` and ``month`` used to date file names.
        See ``parse_year_month``.

    time_dim : str, default 'time'
        Name of the NetCDF time coordinate. Falls back to ``TIME_DIM_NAMES`` when the
        given name is absent.

    x_dim : str, optional
        Name of the NetCDF X dimension. Auto-detected from ``X_DIM_NAMES`` when None.

    y_dim : str, optional
        Name of the NetCDF Y dimension. Auto-detected from ``Y_DIM_NAMES`` when None.

    Attributes
    ----------
    path : Path
        The path/glob the source was constructed from.

    scale_factor : float
        Multiplier applied to returned data.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist and matches nothing.

    ValueError
        If no usable files are found, raster and NetCDF files are mixed, no file name
        could be dated, or the NetCDF precipitation variable is ambiguous.

    Notes
    -----
    **Raster layout limitation.** Raster mode expects *one file per month*, dated by
    file name. A single multi-band raster per year (12 bands = 12 months) is **not**
    supported and raises a ``ValueError``: split it into monthly files, or supply
    a ``date_regex`` matching a one-file-per-month naming scheme. A raster that *is*
    dated to a single month but happens to carry extra bands is accepted - the first
    band is used and a warning is logged.

    Every array returned by ``get_month`` has dims exactly ``('y', 'x')``, dtype
    ``float32``, a descending (north-up) ``y`` coordinate, nodata as ``NaN``, the source
    CRS written via ``rioxarray``, and attrs
    ``{'units', 'long_name', 'year', 'month', 'source'}``.

    **One grid for every month.** All months are stacked into a single time series
    downstream, so every file must share one CRS, shape and affine transform. Raster
    mode checks that up front from file metadata only (no pixels are read) and raises a
    ``ValueError`` naming the first offending file and how it differs; NetCDF mode
    catches a mismatched grid when the month is read. That check always compares the
    CRSs the files *declare* with each other, never with the ``crs`` argument, so a
    mixed-CRS directory is caught whether or not an override was supplied: without
    ``crs`` it is one of the grid differences that raise, and with ``crs`` the override
    still wins (that is what it is for) but the disagreement is logged as a WARNING,
    since relabelling never reprojects. Hence ``shape``, ``resolution``, ``bounds`` and
    ``native_bounds`` describe *every* month, and ``bounds`` / ``native_bounds`` are
    always ordered ``(minx, miny, maxx, maxy)``, including for south-up files.

    **Curvilinear NetCDF coordinates.** A file with 2-D coordinates (the WRF
    ``XLAT`` / ``XLONG`` layout on ``south_north`` / ``west_east`` dims) is supported,
    whether xarray has promoted them to coordinates (the file gives the precipitation
    variable a CF ``coordinates = "XLONG XLAT"`` attribute) or they sit in the file as
    plain 2-D *data variables* on the same dims, as hand-cut and post-processed subsets
    often leave them. Either way the 2-D coordinates are collapsed to 1-D along the
    middle row/column. Map-projected
    model output is always mildly curvilinear in lon/lat, so the resulting displacement
    is measured - beyond half a cell it is logged as a WARNING quantifying the error, and
    beyond 20% of the axis span (or if the collapsed coordinate is not monotonic) it
    raises rather than silently returning a wrong grid. A dimension with *no* usable
    coordinate at all falls back to pixel indices with a WARNING: the shape survives but
    the grid is not georeferenced.

    Examples
    --------
    Directory of monthly GeoTIFFs:

    ```python
    from pycropwat.local import LocalPrecipitationSource

    src = LocalPrecipitationSource('../pyCropWat_Data/Precip',
                                   pattern='Precip_*.tif',
                                   nodata=-9999)
    print(src.kind, len(src), src.year_range)     # raster 264 (2000, 2021)
    da = src.get_month(2005, 7)
    print(da.dims, da.shape)                      # ('y', 'x') (689, 799)
    src.close()
    ```

    NetCDF stack used as a context manager:

    ```python
    with LocalPrecipitationSource('precip.nc', variable='precip') as src:
        print(src.available_months()[:3])
        da = src.get_month(2005, 7)
    ```
    """

    def __init__(
        self,
        path: Union[str, Path],
        pattern: str = '*.tif',
        variable: Optional[str] = None,
        scale_factor: float = 1.0,
        nodata: Optional[Union[float, Sequence[float]]] = None,
        crs: Optional[Union[str, CRS]] = None,
        date_regex: Optional[str] = None,
        time_dim: str = 'time',
        x_dim: Optional[str] = None,
        y_dim: Optional[str] = None
    ):
        self.path = Path(str(path))
        self.pattern = pattern
        self.scale_factor = float(scale_factor)

        self._variable = variable
        self._crs_arg = crs
        self._date_regex = date_regex
        self._time_dim = time_dim
        self._x_dim = x_dim
        self._y_dim = y_dim

        if nodata is None:
            self._nodata_values = ()
        elif isinstance(nodata, (list, tuple, set, np.ndarray)):
            self._nodata_values = tuple(float(value) for value in nodata)
        else:
            self._nodata_values = (float(nodata),)

        # Populated by the mode-specific setup below.
        self._files = self._resolve_files(path, pattern)
        self._kind = self._resolve_kind(self._files)

        self._index: Dict[Tuple[int, int], Any] = {}
        self._dataset: Optional[xr.Dataset] = None
        self._data_array: Optional[xr.DataArray] = None
        self._file_cache: Dict[Path, xr.Dataset] = {}
        self._time_name: Optional[str] = None
        self._var_name: Optional[str] = None
        self._crs: Optional[CRS] = None
        self._shape: Optional[Tuple[int, int]] = None
        self._resolution: Optional[Tuple[float, float]] = None
        self._native_bounds: Optional[Tuple[float, float, float, float]] = None
        self._bounds: Optional[Tuple[float, float, float, float]] = None

        if self._kind == 'raster':
            self._setup_raster()
        else:
            self._setup_netcdf()

        if not self._index:
            raise ValueError(
                f"No dated months could be resolved from {self.path}. "
                f"Expected file names such as 'Precip_2005_07.tif' or a NetCDF time "
                f"coordinate. Pass date_regex=r'(?P<year>\\d{{4}})_(?P<month>\\d{{2}})' "
                f"(adapted to your naming) to override the built-in parsing."
            )

        self._compute_bounds()

        years = self.year_range
        logger.info(
            "Local precipitation source: kind=%s, files=%d, months=%d, years=%d-%d, "
            "crs=%s, shape=%s",
            self._kind, len(self._files), len(self._index), years[0], years[1],
            self._crs, self._shape
        )

    # ------------------------------------------------------------------
    # File resolution
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_files(path: Union[str, Path], pattern: str) -> List[Path]:
        """
        Expand ``path`` into a sorted list of usable data files.

        Parameters
        ----------

        path : str or Path
            Directory, single file or glob string.

        pattern : str
            Glob applied when ``path`` is a directory.

        Returns
        -------
        list
            Sorted list of ``pathlib.Path`` objects.
        """
        path_str = str(path)
        candidate = Path(path_str)

        if candidate.is_dir():
            matches = [item for item in candidate.glob(pattern) if item.is_file()]
            if not matches:
                raise ValueError(
                    f"No files matching pattern '{pattern}' in directory {candidate}. "
                    f"Adjust local_precip_pattern (e.g. '*.tif', 'Precip_*.tif', '*.nc')."
                )
        elif candidate.is_file():
            matches = [candidate]
        elif _is_glob(path_str):
            matches = [Path(item) for item in _glob.glob(path_str, recursive=True)]
            matches = [item for item in matches if item.is_file()]
            if not matches:
                raise ValueError(
                    f"Glob '{path_str}' matched no files. Check the path and wildcards."
                )
        else:
            raise FileNotFoundError(
                f"Local precipitation path not found: {path_str}. "
                f"Provide a directory of monthly rasters, a NetCDF file, or a glob "
                f"such as '/data/precip/Precip_*.tif'."
            )

        usable = [item for item in matches if _classify_suffix(item) is not None]
        skipped = len(matches) - len(usable)
        if skipped:
            logger.debug("Ignored %d non-raster/non-NetCDF file(s) under %s", skipped, path_str)

        if not usable:
            raise ValueError(
                f"No usable raster or NetCDF files found at {path_str}. "
                f"Supported suffixes: {', '.join(RASTER_SUFFIXES + NETCDF_SUFFIXES)}."
            )

        return sorted(usable, key=lambda item: str(item))

    @staticmethod
    def _resolve_kind(files: List[Path]) -> str:
        """
        Determine whether the resolved files are rasters or NetCDF stacks.

        Parameters
        ----------

        files : list
            Resolved data files.

        Returns
        -------
        str
            Either ``'raster'`` or ``'netcdf'``.
        """
        kinds = {_classify_suffix(item) for item in files}
        if len(kinds) > 1:
            raster_examples = [item.name for item in files if _classify_suffix(item) == 'raster']
            netcdf_examples = [item.name for item in files if _classify_suffix(item) == 'netcdf']
            raise ValueError(
                "Mixed raster and NetCDF inputs are not supported. Found "
                f"{len(raster_examples)} raster file(s) (e.g. {raster_examples[0]}) and "
                f"{len(netcdf_examples)} NetCDF file(s) (e.g. {netcdf_examples[0]}). "
                "Narrow the selection with local_precip_pattern (e.g. '*.tif' or '*.nc')."
            )
        return kinds.pop()

    # ------------------------------------------------------------------
    # Raster mode
    # ------------------------------------------------------------------

    def _setup_raster(self) -> None:
        """Index one-file-per-month rasters and read their spatial metadata."""
        undated: List[Path] = []

        for item in self._files:
            key = parse_year_month(item.stem, self._date_regex)
            if key is None:
                undated.append(item)
                continue
            if key in self._index:
                logger.warning(
                    "Duplicate month %d-%02d: keeping %s, ignoring %s",
                    key[0], key[1], Path(self._index[key]).name, item.name
                )
                continue
            self._index[key] = item

        for item in undated:
            self._check_multiband_year_raster(item)
        if undated:
            logger.warning(
                "Skipped %d file(s) with no parsable year/month (e.g. %s)",
                len(undated), undated[0].name
            )

        if not self._index:
            return

        self._validate_raster_grid()

    def _validate_raster_grid(self) -> None:
        """
        Check that every indexed raster shares one grid and record that common grid.

        Only file *metadata* is read (``rasterio.open`` + ``.profile``), never pixel data,
        so indexing a directory of a few hundred monthly files stays in the millisecond
        range. The earliest month is the reference; every other month must agree with it
        on CRS, ``(n_rows, n_cols)`` and affine transform.

        The CRS comparison is always between the CRSs the *files themselves declare*,
        never against the ``crs`` argument, so a mixed-CRS directory is still detected
        when an override is supplied. What happens next depends on that argument:

        * no ``crs`` argument - a disagreement is part of the grid mismatch and raises;
        * ``crs`` supplied - the override still wins (it is documented to replace
          whatever a file declares, and relabelling files whose CRS metadata is wrong or
          inconsistent is exactly what it is for), but the disagreement is logged as a
          loud WARNING naming the files and both CRSs, because relabelling does not
          reproject: the result is only correct if the pixels really do share one grid.

        Under an override, two *declarations* must differ to count as a conflict, and
        they are compared against the first file that declares one (not necessarily the
        earliest month): a file declaring no CRS at all is not reported, because
        supplying the missing label is the other thing the override exists for.

        Raises
        ------
        ValueError
            If any file differs from the reference grid. The message names the first
            offending file and exactly how it differs.
        """
        reference = Path(self._index[min(self._index)])
        with rasterio.open(reference) as src:
            profile = src.profile
        ref_crs = profile.get('crs')
        ref_shape = (int(profile['height']), int(profile['width']))
        ref_transform = profile['transform']

        # First declared CRS seen, used only to compare declarations with each other.
        declared_crs = ref_crs
        declared_by = reference.name
        crs_conflicts: List[str] = []

        for key in sorted(self._index):
            path = Path(self._index[key])
            if path == reference:
                continue
            with rasterio.open(path) as src:
                profile = src.profile          # metadata only - no pixels are read
            shape = (int(profile['height']), int(profile['width']))
            transform = profile['transform']
            crs = profile.get('crs')

            differences = []
            if shape != ref_shape:
                differences.append(
                    f"shape (n_rows, n_cols) is {shape}, expected {ref_shape}"
                )
            if not self._transforms_match(transform, ref_transform):
                differences.append(
                    f"affine transform is {tuple(round(v, 12) for v in transform[:6])}, "
                    f"expected {tuple(round(v, 12) for v in ref_transform[:6])}"
                )
            # Declared CRSs are always compared with each other, never with the
            # override, so a genuinely mixed directory is detected either way.
            if self._crs_arg is None:
                if crs != ref_crs:
                    differences.append(f"CRS is {crs}, expected {ref_crs}")
            elif crs is not None:
                if declared_crs is None:
                    declared_crs, declared_by = crs, path.name
                elif crs != declared_crs:
                    crs_conflicts.append(
                        f"{path.name} ({key[0]}-{key[1]:02d}) declares {crs}"
                    )

            if differences:
                raise ValueError(
                    f"Local precipitation files do not share a single grid. "
                    f"{path.name} ({key[0]}-{key[1]:02d}) disagrees with "
                    f"{reference.name}: " + '; '.join(differences) + ". "
                    f"Every monthly file must sit on the same grid, because the monthly "
                    f"outputs are stacked into one time series. Regrid the odd file(s) "
                    f"(e.g. gdalwarp -t_srs ... -te ... -ts ...) or narrow "
                    f"local_precip_pattern to a consistent subset."
                )

        if crs_conflicts:
            shown = '; '.join(crs_conflicts[:3])
            if len(crs_conflicts) > 3:
                shown += f"; ... ({len(crs_conflicts)} file(s) in total)"
            logger.warning(
                "Local precipitation files declare more than one CRS: %s declares %s, "
                "but %s. The requested crs=%s is being applied to all of them, which "
                "relabels the grids without reprojecting anything, so the output is "
                "only correct if these files really do sit on the same grid. If they do "
                "not, regrid the odd file(s) (e.g. gdalwarp -t_srs ... -te ... -ts ...) "
                "or narrow local_precip_pattern to a consistent subset.",
                declared_by, declared_crs, shown, self._crs_arg
            )

        self._shape = ref_shape
        self._resolution = (abs(float(ref_transform.a)), abs(float(ref_transform.e)))
        self._native_bounds = self._normalize_bounds(
            array_bounds(ref_shape[0], ref_shape[1], ref_transform)
        )
        self._crs = self._resolve_crs(ref_crs, reference)

    @staticmethod
    def _transforms_match(transform: Any, reference: Any) -> bool:
        """
        Compare two affine transforms to within a fraction of a pixel.

        Parameters
        ----------

        transform : affine.Affine
            Transform under test.

        reference : affine.Affine
            Transform of the reference file.

        Returns
        -------
        bool
            True when the two describe the same grid.
        """
        scale = max(abs(float(reference.a)), abs(float(reference.e)), 1e-12)
        tolerance = _TRANSFORM_PIXEL_TOL * scale
        return all(
            abs(float(left) - float(right)) <= tolerance
            for left, right in zip(tuple(transform)[:6], tuple(reference)[:6])
        )

    @staticmethod
    def _normalize_bounds(values: Sequence[float]) -> Tuple[float, float, float, float]:
        """
        Order a bounds 4-tuple as ``(minx, miny, maxx, maxy)``.

        South-up transforms (a positive north-south pixel size) yield ``bottom > top``,
        which downstream code reads as an inverted box; normalising here keeps
        ``bounds`` and ``native_bounds`` well-ordered for every grid.

        Parameters
        ----------

        values : sequence of float
            ``(left, bottom, right, top)`` in any order.

        Returns
        -------
        tuple
            ``(minx, miny, maxx, maxy)``.
        """
        left, bottom, right, top = (float(value) for value in values)
        return (min(left, right), min(bottom, top), max(left, right), max(bottom, top))

    def _check_multiband_year_raster(self, path: Path) -> None:
        """
        Reject the unsupported "one 12-band raster per year" layout with a clear error.

        Parameters
        ----------

        path : Path
            A raster whose file name yielded no ``(year, month)`` pair.
        """
        match = None
        for candidate in _YEAR_ONLY_PATTERN.finditer(path.stem):
            year = int(candidate.group('year'))
            if _MIN_YEAR <= year <= _MAX_YEAR:
                match = year
        if match is None:
            return

        try:
            with rasterio.open(path) as src:
                band_count = int(src.count)
        except Exception:  # pragma: no cover - unreadable file, handled by caller
            return

        if band_count >= 12:
            raise ValueError(
                f"{path.name} looks like a single multi-band raster for year {match} "
                f"({band_count} bands). One multi-band raster per year is not supported. "
                f"Split the data into one file per month (e.g. 'Precip_{match}_01.tif' ... "
                f"'Precip_{match}_12.tif'), or supply a date_regex matching a "
                f"one-file-per-month naming scheme."
            )

    def _read_raster_month(self, year: int, month: int) -> xr.DataArray:
        """
        Read one monthly raster into memory.

        Parameters
        ----------

        year : int
            Calendar year.

        month : int
            Calendar month (1-12).

        Returns
        -------
        xr.DataArray
            The raw (unfinalised) 2-D array for that month.
        """
        path = Path(self._index[(year, month)])
        raw = rioxarray.open_rasterio(path, masked=True)
        try:
            if 'band' in raw.dims:
                if raw.sizes['band'] == 1:
                    raw = raw.squeeze('band', drop=True)
                else:
                    logger.warning(
                        "%s has %d bands; using the first band for %d-%02d",
                        path.name, raw.sizes['band'], year, month
                    )
                    raw = raw.isel(band=0, drop=True)
            data = raw.load()
        finally:
            raw.close()

        file_crs = data.rio.crs
        if (self._crs_arg is None and file_crs is not None
                and self._crs is not None and file_crs != self._crs):
            logger.warning(
                "%s CRS (%s) differs from the source CRS (%s)", path.name, file_crs, self._crs
            )
        return self._finalize(data, year, month, path)

    # ------------------------------------------------------------------
    # NetCDF mode
    # ------------------------------------------------------------------

    def _setup_netcdf(self) -> None:
        """Open the NetCDF file(s), pick the variable, and index the months."""
        probe = xr.open_dataset(self._files[0], decode_coords='all')
        try:
            self._var_name = self._detect_variable(probe)
            time_name = self._detect_time(probe)
        except Exception:
            probe.close()
            raise

        if time_name is None:
            # One month per file, dated by file name.
            probe.close()
            logger.info(
                "No time coordinate found in %s; dating NetCDF files by name",
                self._files[0].name
            )
            undated = []
            for item in self._files:
                key = parse_year_month(item.stem, self._date_regex)
                if key is None:
                    undated.append(item)
                    continue
                if key in self._index:
                    logger.warning(
                        "Duplicate month %d-%02d: keeping %s, ignoring %s",
                        key[0], key[1], Path(self._index[key]).name, item.name
                    )
                    continue
                self._index[key] = item
            if undated:
                logger.warning(
                    "Skipped %d NetCDF file(s) with no parsable year/month (e.g. %s)",
                    len(undated), undated[0].name
                )
            if not self._index:
                return
            template_source = Path(self._index[min(self._index)])
            template_ds = self._get_file_dataset(template_source)
            template = self._prepare_spatial(template_ds[self._var_name], template_ds)
        else:
            try:
                if len(self._files) == 1:
                    dataset = xr.open_dataset(self._files[0], decode_coords='all')
                else:
                    dataset = xr.open_mfdataset(
                        [str(item) for item in self._files],
                        combine='by_coords',
                        decode_coords='all'
                    )
                self._dataset = dataset
                self._time_name = time_name
                # 2-D coordinate variables are looked for in ``probe`` - the first file,
                # unconcatenated - rather than in ``dataset``: open_mfdataset stacks plain
                # data variables along the time axis, so a static 2-D XLAT/XLONG that the
                # file left in data_vars is 3-D by the time it reaches here. Every file
                # must sit on one grid anyway, so the first file describes them all.
                # ``_promote_2d_coords`` materialises what it takes, so probe can close.
                data_array = self._prepare_spatial(dataset[self._var_name], probe)
                self._data_array = data_array
            finally:
                probe.close()

            times = data_array[time_name].values
            for position, key in enumerate(_times_to_year_month(times)):
                if key in self._index:
                    logger.warning(
                        "Duplicate month %d-%02d in the NetCDF time axis; keeping the first",
                        key[0], key[1]
                    )
                    continue
                self._index[key] = position
            if not self._index:
                return
            template = data_array.isel({time_name: 0}, drop=True)

        template = self._drop_singleton_dims(template)
        self._crs = self._resolve_crs(self._netcdf_crs(template), self._files[0])
        template = template.rio.write_crs(self._crs)
        self._shape = (int(template.sizes['y']), int(template.sizes['x']))
        try:
            x_res, y_res = template.rio.resolution()
            native_bounds = template.rio.bounds()
        except Exception as exc:
            raise ValueError(
                f"Could not derive a regular grid for variable '{self._var_name}' in "
                f"{self._files[0].name} (dims {tuple(str(dim) for dim in template.dims)}): "
                f"{exc}. Pass x_dim=/y_dim= naming the dimensions (or 1-D coordinate "
                f"variables) that describe the grid."
            ) from exc
        self._resolution = (abs(float(x_res)), abs(float(y_res)))
        self._native_bounds = self._normalize_bounds(native_bounds)

    def _detect_variable(self, dataset: xr.Dataset) -> str:
        """
        Choose the precipitation variable in a NetCDF dataset.

        Parameters
        ----------

        dataset : xr.Dataset
            Open dataset to inspect.

        Returns
        -------
        str
            Name of the precipitation variable.
        """
        available = [str(name) for name in dataset.data_vars]

        if self._variable is not None:
            if self._variable in dataset.data_vars:
                return self._variable
            raise ValueError(
                f"Variable '{self._variable}' not found in {self._files[0].name}. "
                f"Available variables: {available}"
            )

        lookup = {name.lower(): name for name in available}
        for candidate in PRECIP_VARIABLE_NAMES:
            if candidate.lower() in lookup:
                return lookup[candidate.lower()]

        candidates = [
            name for name in available
            if dataset[name].ndim >= 2 and name.lower() not in _NON_DATA_VARS
        ]
        if len(candidates) == 1:
            logger.info("Auto-detected NetCDF precipitation variable '%s'", candidates[0])
            return candidates[0]
        if not candidates:
            raise ValueError(
                f"No 2-D (or higher) data variable found in {self._files[0].name}. "
                f"Available variables: {available}"
            )
        raise ValueError(
            f"Could not determine which variable holds precipitation in "
            f"{self._files[0].name}. Candidates: {candidates}. "
            f"Pass local_precip_variable='<name>' (or variable='<name>')."
        )

    def _detect_time(self, obj: Union[xr.Dataset, xr.DataArray]) -> Optional[str]:
        """
        Find the time coordinate name, or None when the data has no time axis.

        Parameters
        ----------

        obj : xr.Dataset or xr.DataArray
            Object to inspect.

        Returns
        -------
        str or None
            Name of the time coordinate.
        """
        names: List[str] = []
        if self._time_dim:
            names.append(self._time_dim)
        names.extend(name for name in TIME_DIM_NAMES if name not in names)

        for name in names:
            if name in obj.coords:
                return name
        for name in names:
            if name in obj.dims:
                logger.warning(
                    "Dimension '%s' has no coordinate values; dating by file name instead",
                    name
                )
                return None
        return None

    def _get_file_dataset(self, path: Path) -> xr.Dataset:
        """
        Open (and cache) a single NetCDF file in per-file mode.

        Parameters
        ----------

        path : Path
            NetCDF file to open.

        Returns
        -------
        xr.Dataset
            The cached open dataset.
        """
        cached = self._file_cache.get(path)
        if cached is None:
            cached = xr.open_dataset(path, decode_coords='all')
            self._file_cache[path] = cached
        return cached

    def _prepare_spatial(
        self,
        data_array: xr.DataArray,
        dataset: Optional[xr.Dataset] = None
    ) -> xr.DataArray:
        """
        Rename the spatial dimensions of a NetCDF array to ``'x'`` / ``'y'``.

        Parameters
        ----------

        data_array : xr.DataArray
            Array straight out of the dataset.

        dataset : xr.Dataset, optional
            Dataset to search for 2-D coordinate variables - normally the one
            ``data_array`` came from. When given, 2-D latitude/longitude *data variables*
            that the file never advertised as coordinates (no CF ``coordinates``
            attribute) are promoted to coordinates first, so they are found by
            ``_ensure_1d_coords`` exactly like promoted ones.

        Returns
        -------
        xr.DataArray
            Array whose horizontal dims are named ``'x'`` and ``'y'``.
        """
        dims = [str(dim) for dim in data_array.dims]

        x_name = self._x_dim
        if x_name is None:
            for candidate in X_DIM_NAMES:
                x_name = next((dim for dim in dims if dim.lower() == candidate.lower()), None)
                if x_name is not None:
                    break
        y_name = self._y_dim
        if y_name is None:
            for candidate in Y_DIM_NAMES:
                y_name = next((dim for dim in dims if dim.lower() == candidate.lower()), None)
                if y_name is not None:
                    break

        if x_name is None or y_name is None:
            raise ValueError(
                f"Could not identify the spatial dimensions of variable "
                f"'{self._var_name}' (dims: {dims}). "
                f"Pass x_dim=/y_dim= explicitly."
            )

        if dataset is not None:
            data_array = self._promote_2d_coords(data_array, dataset, x_name, y_name)

        renames = {}
        if x_name != 'x':
            renames[x_name] = 'x'
        if y_name != 'y':
            renames[y_name] = 'y'
        if renames:
            data_array = data_array.rename(renames)

        return self._ensure_1d_coords(data_array)

    def _promote_2d_coords(
        self,
        data_array: xr.DataArray,
        dataset: xr.Dataset,
        x_name: str,
        y_name: str
    ) -> xr.DataArray:
        """
        Attach 2-D lat/lon *data variables* to the array as coordinates.

        xarray only promotes ``XLAT`` / ``XLONG`` to coordinates when the file says so,
        via a CF ``coordinates = "XLONG XLAT"`` attribute on the precipitation variable.
        Hand-cut and post-processed WRF subsets frequently drop that attribute, leaving
        the 2-D coordinates behind as ordinary data variables; without this step the grid
        would silently degrade to pixel indices even though the real coordinates are
        sitting in the file.

        A variable is only promoted when **both** its dims and its name fit: its dims
        must be exactly the array's own ``(y, x)`` dims (in either order), and its name
        must appear in ``X_DIM_NAMES`` / ``Y_DIM_NAMES``, matched case-insensitively.
        Requiring the dims to match is what keeps an unrelated 2-D field from ever being
        mistaken for a coordinate. An axis that already carries a 1-D coordinate is left
        untouched, and so is the precipitation variable itself. Whatever is promoted is
        read into memory, so ``dataset`` need not stay open afterwards.

        Parameters
        ----------

        data_array : xr.DataArray
            Array whose spatial dims are still named as they are in the file.

        dataset : xr.Dataset
            Dataset to search. Normally the one ``data_array`` came from; in multi-file
            mode it is the first file on its own, because ``open_mfdataset`` stacks plain
            data variables along the time axis and a static 2-D XLAT/XLONG left in
            data_vars is no longer 2-D once that has happened.

        x_name : str
            Name of the X dimension in the file.

        y_name : str
            Name of the Y dimension in the file.

        Returns
        -------
        xr.DataArray
            The array, with any such variables assigned as 2-D coordinates.
        """
        grid_dims = {str(x_name), str(y_name)}
        present = {str(name).lower() for name in data_array.coords}
        promoted: Dict[str, Any] = {}

        for dim_name, candidates in ((x_name, X_DIM_NAMES), (y_name, Y_DIM_NAMES)):
            # An axis that already has a real 1-D coordinate needs nothing.
            if dim_name in data_array.coords and data_array.coords[dim_name].ndim == 1:
                continue
            wanted = {str(candidate).lower() for candidate in candidates}
            for name in dataset.data_vars:
                text = str(name)
                if text == self._var_name or text.lower() in present:
                    continue
                if text.lower() not in wanted or text in promoted:
                    continue
                variable = dataset[name]
                if variable.ndim != 2:
                    continue
                if {str(dim) for dim in variable.dims} != grid_dims:
                    continue
                # Read now: the dataset searched may be closed straight afterwards.
                promoted[text] = xr.Variable(
                    variable.dims, np.asarray(variable.values), attrs=dict(variable.attrs)
                )

        if promoted:
            logger.info(
                "Promoted 2-D data variable(s) %s in %s to coordinates; the file "
                "declares no CF 'coordinates' attribute for them",
                ', '.join(f"'{name}'" for name in sorted(promoted)), self._files[0].name
            )
            data_array = data_array.assign_coords(promoted)
        return data_array

    def _ensure_1d_coords(self, data_array: xr.DataArray) -> xr.DataArray:
        """
        Give the ``x`` / ``y`` dims real 1-D coordinates, deriving them if necessary.

        Three cases are handled, in order:

        1. A real 1-D ``x`` / ``y`` coordinate already exists - nothing to do.
        2. A 2-D coordinate covers the axis (the WRF ``XLAT`` / ``XLONG`` layout). It is
           collapsed to 1-D along the middle row/column, and how much that collapse
           displaces the worst pixel is measured - see ``_collapse_2d_coord``. This
           covers files that never advertised those 2-D fields as coordinates too,
           because ``_prepare_spatial`` has already promoted them by this point.
        3. The dimension has no usable coordinate at all. Pixel indices
           ``0 .. n - 1`` are assigned and a WARNING is logged: the grid keeps its shape
           but carries no meaningful georeferencing.

        Membership is tested with ``axis in data_array.coords`` rather than
        ``coords.get(axis)``, because xarray synthesises a *virtual* 1-D integer index
        for a dimension that has no coordinate variable - reading it as an existing
        coordinate would skip cases 2 and 3 entirely.

        Parameters
        ----------

        data_array : xr.DataArray
            Array with dims named ``'x'`` and ``'y'``.

        Returns
        -------
        xr.DataArray
            Array carrying real 1-D ``x`` and ``y`` coordinates.

        Raises
        ------
        ValueError
            If a 2-D coordinate is too curvilinear to be represented by 1-D coordinates
            (see ``_collapse_2d_coord``).
        """
        for axis, candidates in (('x', X_DIM_NAMES), ('y', Y_DIM_NAMES)):
            if axis in data_array.coords and data_array.coords[axis].ndim == 1:
                continue

            source = self._find_2d_coord(data_array, axis, candidates)
            if source is not None:
                values = self._collapse_2d_coord(source, axis)
                data_array = data_array.assign_coords({axis: values})
                logger.info(
                    "Derived a 1-D '%s' coordinate from the 2-D coordinate '%s'",
                    axis, source.name
                )
                continue

            size = int(data_array.sizes[axis])
            logger.warning(
                "Dimension '%s' of variable '%s' in %s has no 1-D coordinate values, and "
                "no 2-D coordinate for it was found either (searched the variable's "
                "coordinates and the file's 2-D data variables on this grid for the "
                "names %s); falling back to pixel indices 0-%d. The grid keeps its shape "
                "but is not georeferenced. To fix it: pass %s_dim='<name>' to name the "
                "coordinate variable that describes this axis; or, if the file does hold "
                "2-D coordinates under some other name, give '%s' a CF 'coordinates' "
                "attribute pointing at them (e.g. coordinates = \"XLONG XLAT\"); or pass "
                "crs=... if the pixel grid is genuinely all the file has and you only "
                "need it labelled.",
                axis, self._var_name, self._files[0].name,
                ', '.join(str(name) for name in candidates), size - 1,
                axis, self._var_name
            )
            data_array = data_array.assign_coords(
                {axis: np.arange(size, dtype=np.float64)}
            )
        return data_array

    @staticmethod
    def _find_2d_coord(
        data_array: xr.DataArray,
        axis: str,
        candidates: Sequence[str]
    ) -> Optional[xr.DataArray]:
        """
        Find a 2-D coordinate spanning ``('y', 'x')`` that describes ``axis``.

        Only ``data_array.coords`` is searched. 2-D lat/lon variables that the file left
        in ``data_vars`` are handled earlier, by ``_promote_2d_coords``, which assigns
        them as coordinates so they arrive here like any other.

        Parameters
        ----------

        data_array : xr.DataArray
            Array whose horizontal dims are already named ``'x'`` / ``'y'``.

        axis : str
            Either ``'x'`` or ``'y'``.

        candidates : sequence of str
            Coordinate names to try, matched case-insensitively (e.g. ``X_DIM_NAMES``).

        Returns
        -------
        xr.DataArray or None
            The 2-D coordinate, or None when there is none.
        """
        present = {str(name).lower(): str(name) for name in data_array.coords}
        for candidate in (axis, *candidates):
            name = present.get(str(candidate).lower())
            if name is None:
                continue
            coord = data_array.coords[name]
            if coord.ndim == 2 and {str(dim) for dim in coord.dims} == {'x', 'y'}:
                return coord
        return None

    def _collapse_2d_coord(self, source: xr.DataArray, axis: str) -> np.ndarray:
        """
        Collapse a 2-D curvilinear coordinate to the 1-D coordinate of one axis.

        The middle row (for ``x``) or middle column (for ``y``) is taken as the
        representative line, which minimises the worst-case displacement over the domain.
        The displacement of every other row/column from that line is then measured:

        * more than ``_CURVILINEAR_WARN_CELLS`` (half a cell) - accepted, with a WARNING
          quantifying the distortion, since map-projected model output such as a WRF
          Lambert-conformal domain is always mildly curvilinear in lon/lat;
        * more than ``_CURVILINEAR_MAX_REL_DEV`` (20%) of the axis span, or a
          non-monotonic / non-finite line - rejected with a ``ValueError``, because
          no 1-D coordinate can describe such a grid.

        Parameters
        ----------

        source : xr.DataArray
            2-D coordinate with dims ``('y', 'x')`` in either order.

        axis : str
            Either ``'x'`` or ``'y'``.

        Returns
        -------
        np.ndarray
            The derived 1-D coordinate values.

        Raises
        ------
        ValueError
            If the grid is too curvilinear to be represented by 1-D coordinates.
        """
        grid = np.asarray(source.transpose('y', 'x').values, dtype=np.float64)
        if axis == 'x':
            line = grid[grid.shape[0] // 2, :]
            deviation = np.abs(grid - line[np.newaxis, :])
        else:
            line = grid[:, grid.shape[1] // 2]
            deviation = np.abs(grid - line[:, np.newaxis])

        def _reject(reason: str) -> None:
            raise ValueError(
                f"The 2-D coordinate '{source.name}' of variable '{self._var_name}' in "
                f"{self._files[0].name} cannot be represented by a 1-D '{axis}' "
                f"coordinate: {reason}. pyCropWat needs a regular grid. Regrid the file "
                f"first (e.g. with gdalwarp, cdo remapbil or xESMF), or pass "
                f"x_dim=/y_dim= naming 1-D coordinate variables that already describe "
                f"one."
            )

        if not np.all(np.isfinite(line)):
            _reject('the derived coordinate contains non-finite values')

        steps = np.diff(line)
        if line.size > 1 and not (np.all(steps > 0) or np.all(steps < 0)):
            _reject('the derived coordinate is not monotonic')

        span = float(abs(line[-1] - line[0])) if line.size > 1 else 0.0
        cell = float(np.mean(np.abs(steps))) if steps.size else 0.0
        worst = float(np.nanmax(deviation)) if deviation.size else 0.0

        if span > 0.0 and worst / span > _CURVILINEAR_MAX_REL_DEV:
            _reject(
                f'the grid is genuinely curvilinear - collapsing it would misplace '
                f'pixels by up to {worst:.6g} coordinate units, '
                f'{100.0 * worst / span:.1f}% of the {span:.6g} span'
            )

        if cell > 0.0 and worst > _CURVILINEAR_WARN_CELLS * cell:
            logger.warning(
                "'%s' is curvilinear: the 1-D '%s' coordinate derived from it misplaces "
                "pixels by up to %.6g units (%.1f cells, %.2f%% of the %.6g span). "
                "Regrid to a regular grid if that displacement matters.",
                source.name, axis, worst, worst / cell, 100.0 * worst / span, span
            )

        if cell > 0.0 and steps.size > 1:
            spread = float(np.max(np.abs(steps)) - np.min(np.abs(steps)))
            if spread > 0.01 * cell:
                logger.warning(
                    "The derived 1-D '%s' coordinate is not evenly spaced (cell size "
                    "varies by %.2f%%); pyCropWat treats it as a regular grid.",
                    axis, 100.0 * spread / cell
                )

        return line

    @staticmethod
    def _netcdf_crs(data_array: xr.DataArray) -> Optional[CRS]:
        """
        Extract a CRS from a NetCDF array via rioxarray or a grid-mapping variable.

        Parameters
        ----------

        data_array : xr.DataArray
            Spatially-prepared array.

        Returns
        -------
        CRS or None
            The CRS found in the file, or None.
        """
        try:
            found = data_array.rio.crs
        except Exception:
            found = None
        if found is not None:
            return found

        for name in ('spatial_ref', 'crs'):
            coord = data_array.coords.get(name)
            if coord is None:
                continue
            attrs = coord.attrs
            for key in ('crs_wkt', 'spatial_ref', 'esri_pe_string', 'proj4', 'proj4_params'):
                if attrs.get(key):
                    try:
                        return CRS.from_user_input(attrs[key])
                    except Exception:
                        continue
            if attrs.get('epsg_code'):
                try:
                    return CRS.from_user_input(attrs['epsg_code'])
                except Exception:
                    pass
        return None

    def _read_netcdf_month(self, year: int, month: int) -> xr.DataArray:
        """
        Read one month out of the NetCDF stack.

        Parameters
        ----------

        year : int
            Calendar year.

        month : int
            Calendar month (1-12).

        Returns
        -------
        xr.DataArray
            The finalised 2-D array for that month.
        """
        entry = self._index[(year, month)]

        if self._data_array is not None:
            selected = self._data_array.isel({self._time_name: int(entry)}, drop=True)
            source = self._files[0] if len(self._files) == 1 else self.path
        else:
            path = Path(entry)
            dataset = self._get_file_dataset(path)
            selected = self._prepare_spatial(dataset[self._var_name], dataset)
            source = path

        return self._finalize(selected.load(), year, month, source)

    # ------------------------------------------------------------------
    # Shared finalisation
    # ------------------------------------------------------------------

    @staticmethod
    def _drop_singleton_dims(data_array: xr.DataArray) -> xr.DataArray:
        """
        Drop every length-1 dimension other than ``'x'`` and ``'y'``.

        Parameters
        ----------

        data_array : xr.DataArray
            Array possibly carrying singleton band/time/level dims.

        Returns
        -------
        xr.DataArray
            Array without extraneous singleton dimensions.
        """
        extra = [
            str(dim) for dim in data_array.dims
            if str(dim) not in ('x', 'y') and data_array.sizes[dim] == 1
        ]
        if extra:
            data_array = data_array.squeeze(extra, drop=True)
        return data_array

    def _resolve_crs(self, found: Optional[CRS], source: Path) -> CRS:
        """
        Pick the CRS to attach to returned arrays.

        Parameters
        ----------

        found : CRS or None
            CRS read from the file, if any.

        source : Path
            File the CRS was read from (used for logging).

        Returns
        -------
        CRS
            The ``crs`` constructor argument when one was given, otherwise the CRS of
            the source, otherwise EPSG:4326.

        Notes
        -----
        The ``crs`` argument **overrides** whatever the file declares, as documented for
        ``local_precip_crs``. Overriding only relabels the grid - the pixels and their
        coordinates are left exactly as they are on disk, nothing is reprojected - so it
        is the right tool for a file whose CRS metadata is missing or wrong, and the
        wrong tool for changing the projection of correctly-labelled data.
        """
        if self._crs_arg is not None:
            override = CRS.from_user_input(self._crs_arg)
            if found is None:
                logger.info("%s declares no CRS; assigning %s", source.name, override)
            elif CRS.from_user_input(found) != override:
                logger.info(
                    "%s declares CRS %s; overriding it with the requested %s "
                    "(relabelled only - no reprojection)", source.name, found, override
                )
            return override
        if found is not None:
            return CRS.from_user_input(found)
        logger.warning(
            "%s declares no CRS and none was supplied; assuming EPSG:4326. "
            "Pass local_precip_crs=... if this is wrong.", source.name
        )
        return CRS.from_epsg(4326)

    def _mask_nodata(self, values: np.ndarray, extra: Sequence[float]) -> np.ndarray:
        """
        Replace nodata sentinels with NaN.

        Parameters
        ----------

        values : np.ndarray
            Float array to mask (modified in place).

        extra : sequence of float
            Sentinel values to mask.

        Returns
        -------
        np.ndarray
            The masked array.
        """
        for sentinel in extra:
            if sentinel is None:
                continue
            sentinel = float(sentinel)
            if np.isnan(sentinel):
                continue
            mask = np.isclose(values, sentinel, equal_nan=False)
            if mask.any():
                values[mask] = np.nan
        return values

    def _finalize(
        self,
        data_array: xr.DataArray,
        year: int,
        month: int,
        source: Union[str, Path]
    ) -> xr.DataArray:
        """
        Bring a raw monthly array into the canonical pyCropWat form.

        Squeezes singleton dims, orders dims as ``('y', 'x')``, casts to ``float32``,
        masks nodata, applies ``scale_factor``, flips to north-up, writes the CRS and
        sets the standard attributes.

        Parameters
        ----------

        data_array : xr.DataArray
            Raw array read from disk.

        year : int
            Calendar year.

        month : int
            Calendar month (1-12).

        source : str or Path
            File the array came from, recorded in ``attrs['source']``.

        Returns
        -------
        xr.DataArray
            Analysis-ready 2-D precipitation array in mm.
        """
        data_array = self._drop_singleton_dims(data_array)

        dims = tuple(str(dim) for dim in data_array.dims)
        if set(dims) != {'y', 'x'}:
            raise ValueError(
                f"Expected a 2-D ('y', 'x') array for {year}-{month:02d} but got dims "
                f"{dims}. Subset the extra dimension(s) or pass x_dim=/y_dim= "
                f"(and time_dim=) so the layout can be resolved."
            )
        if dims != ('y', 'x'):
            data_array = data_array.transpose('y', 'x')

        shape = (int(data_array.sizes['y']), int(data_array.sizes['x']))
        if self._shape is not None and shape != self._shape:
            raise ValueError(
                f"{Path(source).name} holds a {shape} grid for {year}-{month:02d} but "
                f"the source grid is {self._shape}. Every month must sit on the same "
                f"grid, because the monthly outputs are stacked into one time series. "
                f"Regrid the odd file(s) or narrow the selection to a consistent subset."
            )

        values = np.asarray(data_array.values, dtype=np.float32)

        # Honour _FillValue / missing_value left in attrs (xarray did not mask them).
        sentinels = list(self._nodata_values)
        for key in ('_FillValue', 'missing_value'):
            attr_value = data_array.attrs.get(key)
            if attr_value is None:
                continue
            for item in np.atleast_1d(attr_value).tolist():
                try:
                    sentinels.append(float(item))
                except (TypeError, ValueError):
                    continue
        try:
            file_nodata = data_array.rio.nodata
        except Exception:
            file_nodata = None
        if file_nodata is not None:
            sentinels.append(float(file_nodata))

        values = self._mask_nodata(values, sentinels)

        if self.scale_factor != 1.0:
            values = values * np.float32(self.scale_factor)

        result = xr.DataArray(
            values.astype(np.float32),
            dims=('y', 'x'),
            coords={
                'y': np.asarray(data_array['y'].values),
                'x': np.asarray(data_array['x'].values)
            },
            name='precipitation'
        )

        # GeoTIFF writers downstream expect a north-up (descending y) array.
        y_values = result['y'].values
        if y_values.size > 1 and y_values[0] < y_values[-1]:
            result = result.isel(y=slice(None, None, -1))

        result = result.rio.write_crs(self._crs)
        result.attrs = {
            'units': 'mm',
            'long_name': 'precipitation',
            'year': int(year),
            'month': int(month),
            'source': str(source)
        }
        return result

    def _compute_bounds(self) -> None:
        """Reproject the native bounds to EPSG:4326 for the ``bounds`` property."""
        if self._native_bounds is None:
            return
        target = CRS.from_epsg(4326)
        if self._crs is not None and self._crs != target:
            try:
                self._bounds = tuple(
                    float(value) for value in
                    transform_bounds(self._crs, target, *self._native_bounds)
                )
                return
            except Exception as exc:  # pragma: no cover - exotic CRS
                logger.warning("Could not reproject bounds to EPSG:4326: %s", exc)
        self._bounds = tuple(float(value) for value in self._native_bounds)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def kind(self) -> str:
        """str : ``'raster'`` for one-file-per-month rasters, ``'netcdf'`` otherwise."""
        return self._kind

    @property
    def files(self) -> List[Path]:
        """list of Path : The data files backing this source, sorted."""
        return list(self._files)

    @property
    def crs(self) -> CRS:
        """CRS : Coordinate reference system of the source data."""
        return self._crs

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """tuple : ``(minx, miny, maxx, maxy)`` in EPSG:4326 lon/lat degrees."""
        return self._bounds

    @property
    def native_bounds(self) -> Tuple[float, float, float, float]:
        """tuple : ``(minx, miny, maxx, maxy)`` in the source CRS."""
        return self._native_bounds

    @property
    def year_range(self) -> Tuple[int, int]:
        """tuple : ``(first_year, last_year)`` present in the source."""
        years = [year for year, _ in self._index]
        return min(years), max(years)

    @property
    def resolution(self) -> Tuple[float, float]:
        """tuple : ``(x_resolution, y_resolution)`` as positive values in CRS units."""
        return self._resolution

    @property
    def shape(self) -> Tuple[int, int]:
        """
        tuple : ``(n_rows, n_cols)`` of every monthly grid.

        In raster mode this is validated across all indexed files when the source is
        built; a file on a different grid raises rather than being read.
        """
        return self._shape

    def available_months(self) -> List[Tuple[int, int]]:
        """
        List every month present in the source.

        Returns
        -------
        list
            Sorted list of ``(year, month)`` tuples.

        Examples
        --------
        ```python
        src.available_months()[:2]   # [(2000, 1), (2000, 2)]
        ```
        """
        return sorted(self._index.keys())

    def has_month(self, year: int, month: int) -> bool:
        """
        Check whether a given month is available.

        Parameters
        ----------

        year : int
            Calendar year.

        month : int
            Calendar month (1-12).

        Returns
        -------
        bool
            True when the month can be read.
        """
        return (int(year), int(month)) in self._index

    def get_month(self, year: int, month: int) -> Optional[xr.DataArray]:
        """
        Read the precipitation grid for one month.

        Parameters
        ----------

        year : int
            Calendar year.

        month : int
            Calendar month (1-12).

        Returns
        -------
        xr.DataArray or None
            2-D array with dims ``('y', 'x')``, dtype ``float32``, values in mm
            (``scale_factor`` applied), nodata as NaN, a descending ``y`` coordinate
            and the source CRS attached. Returns None when the month is missing.

        Examples
        --------
        ```python
        da = src.get_month(2005, 7)
        if da is not None:
            print(da.dims, da.shape, float(da.mean()))
        ```
        """
        year = int(year)
        month = int(month)
        if not 1 <= month <= 12:
            raise ValueError(f"month must be between 1 and 12, got {month}")

        if (year, month) not in self._index:
            logger.warning("No local precipitation data for %d-%02d", year, month)
            return None

        if self._kind == 'raster':
            return self._read_raster_month(year, month)
        return self._read_netcdf_month(year, month)

    def close(self) -> None:
        """
        Release any open NetCDF file handles.

        Safe to call more than once. Rasters are opened and closed per call to
        ``get_month``, so this is a no-op in raster mode.
        """
        if self._dataset is not None:
            try:
                self._dataset.close()
            except Exception:  # pragma: no cover - already closed
                pass
            self._dataset = None
        self._data_array = None
        for dataset in self._file_cache.values():
            try:
                dataset.close()
            except Exception:  # pragma: no cover - already closed
                pass
        self._file_cache = {}

    def __enter__(self) -> "LocalPrecipitationSource":
        """Return self so the source can be used in a ``with`` block."""
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        """Close open file handles on exit; never suppresses exceptions."""
        self.close()
        return False

    def __len__(self) -> int:
        """int : Number of months available in this source."""
        return len(self._index)

    def __repr__(self) -> str:
        """Return a concise, informative representation."""
        years = self.year_range
        return (
            f"LocalPrecipitationSource(kind='{self._kind}', files={len(self._files)}, "
            f"months={len(self._index)}, years={years[0]}-{years[1]}, "
            f"crs='{self._crs}', shape={self._shape})"
        )


def open_local_precipitation(path: Union[str, Path], **kwargs) -> LocalPrecipitationSource:
    """
    Open a local precipitation dataset.

    Thin convenience factory around ``LocalPrecipitationSource``; every keyword
    argument is forwarded unchanged.

    Parameters
    ----------

    path : str or Path
        A directory of monthly rasters, a single raster/NetCDF file, or a glob string.

    **kwargs
        Keyword arguments forwarded to ``LocalPrecipitationSource``
        (``pattern``, ``variable``, ``scale_factor``, ``nodata``, ``crs``,
        ``date_regex``, ``time_dim``, ``x_dim``, ``y_dim``).

    Returns
    -------
    LocalPrecipitationSource
        An indexed, ready-to-read precipitation source.

    Examples
    --------
    ```python
    from pycropwat.local import open_local_precipitation

    src = open_local_precipitation('../pyCropWat_Data/Precip',
                                   pattern='Precip_*.tif',
                                   nodata=-9999)
    da = src.get_month(2005, 7)
    src.close()
    ```
    """
    return LocalPrecipitationSource(path, **kwargs)

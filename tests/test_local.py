"""
Tests for local precipitation input (``pycropwat.local`` and the local mode of
:class:`pycropwat.core.EffectivePrecipitation`).

Every fixture is synthesised on disk in ``tmp_path``: the suite never touches the
network, never initialises Earth Engine and does not depend on any external data.
"""

import json
import logging
import re
from pathlib import Path

import numpy as np
import pytest
import rasterio
import rioxarray  # noqa: F401  (registers the .rio accessor used below)
import xarray as xr
from rasterio.transform import Affine, from_origin

import pycropwat
import pycropwat.core as core
import pycropwat.utils as utils
from pycropwat.local import (
    LocalPrecipitationSource,
    open_local_precipitation,
    parse_year_month,
)


# ---------------------------------------------------------------------------
# Synthetic data helpers
# ---------------------------------------------------------------------------

#: Sentinel written into the GeoTIFF/NetCDF fixtures for the nodata pixel.
NODATA = -9999.0

#: 3x3 base grid in mm. NaN marks the pixel that is stored as nodata on disk.
BASE_PRECIP = np.array([
    [100.0, 300.0, 0.0],
    [50.0, np.nan, 250.0],
    [10.0, 20.0, 30.0],
], dtype=np.float32)

#: North-up EPSG:4326 transform: 1 degree cells with the origin at (-100, 40).
BASE_TRANSFORM = from_origin(-100.0, 40.0, 1.0, 1.0)

#: (minx, miny, maxx, maxy) of :data:`BASE_TRANSFORM` for a 3x3 grid.
BASE_BOUNDS = (-100.0, 37.0, -97.0, 40.0)

#: Descending (north-up) cell-centre latitudes of the 3x3 fixture grid.
BASE_Y = np.array([39.5, 38.5, 37.5])

#: Cell-centre longitudes of the 3x3 fixture grid.
BASE_X = np.array([-99.5, -98.5, -97.5])

#: Months written by the shared raster/NetCDF fixtures.
FIXTURE_MONTHS = [(2001, 1), (2001, 2), (2001, 3), (2002, 1), (2002, 2), (2002, 3)]

#: First NetCDF engine available in this environment, or None.
NETCDF_ENGINE = next(
    (name for name in ('netcdf4', 'h5netcdf', 'scipy') if name in xr.backends.list_engines()),
    None
)

#: Skip marker applied to every NetCDF test when no backend is installed.
requires_netcdf = pytest.mark.skipif(
    NETCDF_ENGINE is None, reason='No NetCDF backend (netcdf4/h5netcdf/scipy) available'
)


def month_array(year, month):
    """Build the deterministic precipitation grid for one month (NaN = nodata)."""
    offset = np.float32((year - 2001) * 12 + (month - 1))
    return (BASE_PRECIP + offset).astype(np.float32)


def expected_cropwat(pr):
    """Re-implement the CROPWAT formula independently of ``pycropwat.methods``."""
    return np.where(pr <= 250.0, pr * (125.0 - 0.2 * pr) / 125.0, 0.1 * pr + 125.0)


def write_raster(path, values, transform=BASE_TRANSFORM, crs='EPSG:4326', nodata=NODATA,
                 declare_nodata=True):
    """Write a single-band GeoTIFF, replacing NaN with ``nodata`` on disk."""
    values = np.asarray(values, dtype=np.float32)
    if nodata is not None:
        values = np.where(np.isnan(values), np.float32(nodata), values)
    profile = {
        'driver': 'GTiff',
        'height': values.shape[0],
        'width': values.shape[1],
        'count': 1,
        'dtype': 'float32',
        'transform': transform,
    }
    if crs is not None:
        profile['crs'] = crs
    if nodata is not None and declare_nodata:
        profile['nodata'] = nodata
    with rasterio.open(path, 'w', **profile) as dataset:
        dataset.write(values, 1)
    return Path(path)


def write_netcdf(dataset, path):
    """Write a Dataset with whichever NetCDF backend this environment provides."""
    dataset.to_netcdf(path, engine=NETCDF_ENGINE)
    return Path(path)


def build_stack_dataset(months, y_name='y', x_name='x', var_name='precip', ascending_y=False):
    """Build a (time, y, x) Dataset holding the fixture months."""
    y_values = BASE_Y[::-1] if ascending_y else BASE_Y
    stack = []
    for year, month in months:
        values = month_array(year, month)
        stack.append(values[::-1] if ascending_y else values)
    data = np.stack(stack).astype(np.float32)
    data = np.where(np.isnan(data), np.float32(NODATA), data)
    times = np.array(
        [f'{year:04d}-{month:02d}-01' for year, month in months], dtype='datetime64[ns]'
    )
    dataset = xr.Dataset(
        {var_name: (('time', y_name, x_name), data)},
        coords={'time': times, y_name: y_values, x_name: BASE_X}
    )
    dataset[var_name].attrs['_FillValue'] = NODATA
    return dataset


def write_geojson(path, minx, miny, maxx, maxy):
    """Write a one-feature GeoJSON polygon covering the given WGS84 box."""
    ring = [[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy], [minx, miny]]
    payload = {
        'type': 'FeatureCollection',
        'features': [{
            'type': 'Feature',
            'properties': {},
            'geometry': {'type': 'Polygon', 'coordinates': [ring]}
        }]
    }
    Path(path).write_text(json.dumps(payload))
    return Path(path)


def read_band(path):
    """Read a single-band GeoTIFF as a masked 2-D DataArray."""
    return rioxarray.open_rasterio(path, masked=True).squeeze('band', drop=True)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def raster_dir(tmp_path_factory):
    """Directory of monthly GeoTIFFs plus a non-raster decoy file."""
    directory = tmp_path_factory.mktemp('precip_rasters')
    for year, month in FIXTURE_MONTHS:
        write_raster(directory / f'Precip_{year}_{month:02d}.tif', month_array(year, month))
    # Mirrors the stray companion file that sits in the real WRF data directory.
    (directory / '.geeup-state.json').write_text('{}')
    return directory


@pytest.fixture(scope='module')
def netcdf_stack(tmp_path_factory):
    """Single NetCDF file holding every fixture month on a time axis."""
    if NETCDF_ENGINE is None:
        pytest.skip('No NetCDF backend available')
    directory = tmp_path_factory.mktemp('precip_netcdf')
    return write_netcdf(build_stack_dataset(FIXTURE_MONTHS), directory / 'precip.nc')


@pytest.fixture
def no_earth_engine(monkeypatch):
    """Turn any Earth Engine initialisation into an immediate test failure."""
    def _boom(*args, **kwargs):
        raise AssertionError(
            'Earth Engine was initialised during a local, precipitation-only run'
        )

    monkeypatch.setattr(core, 'initialize_gee', _boom)
    monkeypatch.setattr(utils, 'initialize_gee', _boom)
    monkeypatch.setattr(core.ee, 'Initialize', _boom)
    monkeypatch.setattr(core.ee, 'Authenticate', _boom, raising=False)
    return _boom


# ---------------------------------------------------------------------------
# parse_year_month
# ---------------------------------------------------------------------------

class TestParseYearMonth:
    """Tests for the file-name date parser."""

    @pytest.mark.parametrize('name, expected', [
        ('Precip_2005_07.tif', (2005, 7)),
        ('effective_precip_2005_07.tif', (2005, 7)),
        ('precip-2005-07.nc', (2005, 7)),
        ('pr200507.tif', (2005, 7)),
        ('pr.2005.07.tif', (2005, 7)),
        ('Precip_2005_7.tif', (2005, 7)),
        ('2005_07', (2005, 7)),
    ])
    def test_supported_layouts(self, name, expected):
        """Test that every documented naming layout is recognised."""
        assert parse_year_month(name) == expected

    def test_accepts_path_objects(self):
        """Test that a Path is parsed the same way as a string."""
        assert parse_year_month(Path('/data/Precip_2005_07.tif')) == (2005, 7)

    def test_directory_noise_is_ignored(self):
        """Test that a year-like directory component does not confuse the parser."""
        assert parse_year_month('/data/2019/x-2005-07.nc') == (2005, 7)

    def test_last_match_wins(self):
        """Test that the last date in the string is the one that is used."""
        assert parse_year_month('/precip/2005_07/data_2010_03.tif') == (2010, 3)

    @pytest.mark.parametrize('name', ['Precip_2005_00.tif', 'Precip_2005_13.tif'])
    def test_month_out_of_range_rejected(self, name):
        """Test that month 0 and month 13 are not accepted as dates."""
        assert parse_year_month(name) is None

    @pytest.mark.parametrize('name', ['Precip_1500_07.tif', 'Precip_2500_07.tif'])
    def test_implausible_years_rejected(self, name):
        """Test that years outside the 1700-2200 window are rejected."""
        assert parse_year_month(name) is None

    @pytest.mark.parametrize('name', ['no_date_here.tif', 'precipitation.tif', ''])
    def test_no_date_returns_none(self, name):
        """Test that strings carrying no date return None."""
        assert parse_year_month(name) is None

    def test_custom_date_regex(self):
        """Test a custom MM_YYYY layout supplied through date_regex."""
        parsed = parse_year_month(
            'rain_07_2005.tif', date_regex=r'(?P<month>\d{2})_(?P<year>\d{4})'
        )
        assert parsed == (2005, 7)

    def test_custom_date_regex_overrides_builtins(self):
        """Test that a custom regex that matches nothing does not fall back."""
        assert parse_year_month('Precip_2005_07.tif', date_regex=r'(?P<year>y)(?P<month>m)') is None

    def test_date_regex_requires_named_groups(self):
        """Test that a regex without 'year'/'month' groups raises ValueError."""
        with pytest.raises(ValueError, match='named groups'):
            parse_year_month('x', date_regex=r'(?P<yr>\d{4})')

    def test_invalid_date_regex(self):
        """Test that an uncompilable regex raises ValueError."""
        with pytest.raises(ValueError, match='Invalid date_regex'):
            parse_year_month('x', date_regex=r'(?P<year>\d{4}')


# ---------------------------------------------------------------------------
# LocalPrecipitationSource - raster mode
# ---------------------------------------------------------------------------

class TestRasterSource:
    """Tests for one-file-per-month raster sources."""

    def test_from_directory(self, raster_dir):
        """Test construction from a directory of monthly GeoTIFFs."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        assert src.kind == 'raster'
        assert len(src.files) == len(FIXTURE_MONTHS)
        assert len(src) == len(FIXTURE_MONTHS)
        src.close()

    def test_from_glob(self, raster_dir):
        """Test construction from a glob string instead of a directory."""
        src = LocalPrecipitationSource(str(raster_dir / 'Precip_2001_*.tif'))
        assert [path.name for path in src.files] == [
            'Precip_2001_01.tif', 'Precip_2001_02.tif', 'Precip_2001_03.tif'
        ]
        assert src.year_range == (2001, 2001)
        src.close()

    def test_from_single_file(self, raster_dir):
        """Test construction from a single raster file."""
        src = LocalPrecipitationSource(raster_dir / 'Precip_2002_03.tif')
        assert len(src) == 1
        assert src.available_months() == [(2002, 3)]
        src.close()

    def test_files_are_sorted(self, raster_dir):
        """Test that the resolved file list comes back sorted."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        names = [path.name for path in src.files]
        assert names == sorted(names)
        src.close()

    def test_index_metadata(self, raster_dir):
        """Test year_range, available_months, has_month, shape and resolution."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        assert src.year_range == (2001, 2002)
        assert src.available_months() == FIXTURE_MONTHS
        assert src.has_month(2001, 2) is True
        assert src.has_month(2001, 12) is False
        assert src.shape == (3, 3)
        assert src.resolution == pytest.approx((1.0, 1.0))
        src.close()

    def test_crs_and_bounds(self, raster_dir):
        """Test that the CRS and both bounds properties reflect the files."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        assert src.crs.to_epsg() == 4326
        assert src.native_bounds == pytest.approx(BASE_BOUNDS)
        # Source is already WGS84, so bounds and native_bounds agree.
        assert src.bounds == pytest.approx(BASE_BOUNDS)
        src.close()

    def test_bounds_reprojected_to_wgs84(self, tmp_path):
        """Test that a projected source reports lon/lat bounds via bounds."""
        transform = from_origin(-11131949.0, 4865942.0, 100000.0, 100000.0)
        write_raster(
            tmp_path / 'Precip_2001_01.tif',
            np.zeros((3, 4), dtype=np.float32),
            transform=transform,
            crs='EPSG:3857'
        )
        src = LocalPrecipitationSource(tmp_path)
        assert src.crs.to_epsg() == 3857
        assert src.native_bounds == pytest.approx(
            (-11131949.0, 4565942.0, -10731949.0, 4865942.0)
        )
        minx, miny, maxx, maxy = src.bounds
        assert minx == pytest.approx(-100.0, abs=0.01)
        assert maxy == pytest.approx(40.0, abs=0.01)
        assert -180.0 < minx < maxx < 180.0
        assert -90.0 < miny < maxy < 90.0
        src.close()

    def test_get_month_shape_and_metadata(self, raster_dir):
        """Test the canonical form of an array returned by get_month."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        data = src.get_month(2001, 1)

        assert data.dims == ('y', 'x')
        assert data.dtype == np.float32
        assert data.shape == (3, 3)
        assert data.rio.crs.to_epsg() == 4326
        # y must be descending (north-up) for the GeoTIFF writers downstream.
        assert np.all(np.diff(data['y'].values) < 0)
        assert data.attrs['units'] == 'mm'
        assert data.attrs['long_name'] == 'precipitation'
        assert data.attrs['year'] == 2001
        assert data.attrs['month'] == 1
        assert Path(data.attrs['source']).name == 'Precip_2001_01.tif'
        src.close()

    @pytest.mark.parametrize('year, month', FIXTURE_MONTHS)
    def test_get_month_values(self, raster_dir, year, month):
        """Test that every month reads back the values that were written."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        data = src.get_month(year, month)
        np.testing.assert_allclose(data.values, month_array(year, month), equal_nan=True)
        src.close()

    def test_declared_nodata_becomes_nan(self, raster_dir):
        """Test that the nodata value stored in the file metadata is masked."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        data = src.get_month(2001, 1)
        assert np.isnan(data.values[1, 1])
        assert int(np.isnan(data.values).sum()) == 1
        src.close()

    def test_explicit_nodata_sentinel(self, tmp_path):
        """Test that a nodata sentinel passed in is masked even when undeclared."""
        write_raster(
            tmp_path / 'Precip_2001_01.tif', month_array(2001, 1), declare_nodata=False
        )

        without = LocalPrecipitationSource(tmp_path)
        assert not np.isnan(without.get_month(2001, 1).values).any()
        assert without.get_month(2001, 1).values[1, 1] == pytest.approx(NODATA)
        without.close()

        with_sentinel = LocalPrecipitationSource(tmp_path, nodata=NODATA)
        assert np.isnan(with_sentinel.get_month(2001, 1).values[1, 1])
        with_sentinel.close()

    def test_scale_factor(self, raster_dir):
        """Test that scale_factor multiplies the returned values."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif', scale_factor=10.0)
        data = src.get_month(2001, 1)
        np.testing.assert_allclose(
            data.values, month_array(2001, 1) * 10.0, equal_nan=True, rtol=1e-6
        )
        assert np.isnan(data.values[1, 1])
        src.close()

    def test_missing_month_returns_none(self, raster_dir):
        """Test that a month with no file returns None instead of raising."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        assert src.get_month(2001, 12) is None
        src.close()

    @pytest.mark.parametrize('month', [0, 13])
    def test_get_month_rejects_invalid_month(self, raster_dir, month):
        """Test that get_month validates the month number."""
        src = LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif')
        with pytest.raises(ValueError, match='month must be between 1 and 12'):
            src.get_month(2001, month)
        src.close()

    def test_south_up_raster_is_flipped(self, tmp_path):
        """Test that a south-up source comes back north-up."""
        transform = Affine(1.0, 0.0, -100.0, 0.0, 1.0, 37.0)
        values = np.arange(9, dtype=np.float32).reshape(3, 3)
        write_raster(tmp_path / 'Precip_2001_01.tif', values, transform=transform, nodata=None)

        src = LocalPrecipitationSource(tmp_path)
        data = src.get_month(2001, 1)
        assert np.all(np.diff(data['y'].values) < 0)
        # Row 0 of the returned array is the northernmost row, i.e. the last row on disk.
        np.testing.assert_array_equal(data.values[0], values[2])
        src.close()

    def test_non_raster_files_are_ignored(self, raster_dir):
        """Test that a '.geeup-state.json' decoy never reaches the reader."""
        src = LocalPrecipitationSource(raster_dir, pattern='*')
        names = [path.name for path in src.files]
        assert '.geeup-state.json' not in names
        assert len(names) == len(FIXTURE_MONTHS)
        src.close()

    def test_default_pattern_skips_decoy(self, raster_dir):
        """Test that the default '*.tif' pattern picks up only the rasters."""
        src = LocalPrecipitationSource(raster_dir)
        assert len(src) == len(FIXTURE_MONTHS)
        src.close()

    def test_missing_path_raises(self, tmp_path):
        """Test that a non-existent path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match='not found'):
            LocalPrecipitationSource(tmp_path / 'does_not_exist')

    def test_empty_pattern_raises(self, raster_dir):
        """Test that a pattern matching nothing raises ValueError."""
        with pytest.raises(ValueError, match='No files matching pattern'):
            LocalPrecipitationSource(raster_dir, pattern='*.nc')

    def test_glob_matching_nothing_raises(self, raster_dir):
        """Test that a glob string matching no file raises ValueError."""
        with pytest.raises(ValueError, match='matched no files'):
            LocalPrecipitationSource(str(raster_dir / 'Nope_*.tif'))

    def test_only_unusable_files_raises(self, tmp_path):
        """Test that a directory holding only non-raster files raises ValueError."""
        (tmp_path / '.geeup-state.json').write_text('{}')
        with pytest.raises(ValueError, match='No usable raster or NetCDF files'):
            LocalPrecipitationSource(tmp_path, pattern='*')

    def test_undated_files_raise(self, tmp_path):
        """Test that files with no parsable date leave the index empty."""
        write_raster(tmp_path / 'precipitation.tif', month_array(2001, 1))
        with pytest.raises(ValueError, match='No dated months could be resolved'):
            LocalPrecipitationSource(tmp_path)

    def test_custom_date_regex_indexes_files(self, tmp_path):
        """Test that date_regex drives the raster index for an unusual layout."""
        write_raster(tmp_path / 'rain_07_2005.tif', month_array(2001, 1))
        src = LocalPrecipitationSource(
            tmp_path, date_regex=r'(?P<month>\d{2})_(?P<year>\d{4})'
        )
        assert src.available_months() == [(2005, 7)]
        src.close()

    def test_crs_argument_assigned_when_file_has_none(self, tmp_path):
        """Test that the crs argument is used when the file declares none."""
        write_raster(
            tmp_path / 'Precip_2001_01.tif', month_array(2001, 1), crs=None, nodata=None
        )
        src = LocalPrecipitationSource(tmp_path, crs='EPSG:3857')
        assert src.crs.to_epsg() == 3857
        assert src.get_month(2001, 1).rio.crs.to_epsg() == 3857
        src.close()

    def test_context_manager_and_repr(self, raster_dir):
        """Test the context-manager protocol and the repr summary."""
        with LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif') as src:
            text = repr(src)
        assert "kind='raster'" in text
        assert 'months=6' in text
        assert 'years=2001-2002' in text

    def test_open_local_precipitation_factory(self, raster_dir):
        """Test that the convenience factory forwards its keyword arguments."""
        src = open_local_precipitation(raster_dir, pattern='Precip_*.tif', scale_factor=2.0)
        assert isinstance(src, LocalPrecipitationSource)
        assert src.scale_factor == 2.0
        assert src.kind == 'raster'
        src.close()

    def test_exports_are_public(self):
        """Test that the local API is re-exported from the package root."""
        assert pycropwat.LocalPrecipitationSource is LocalPrecipitationSource
        assert pycropwat.open_local_precipitation is open_local_precipitation
        assert pycropwat.parse_year_month is parse_year_month
        for name in ('LocalPrecipitationSource', 'open_local_precipitation', 'parse_year_month'):
            assert name in pycropwat.__all__


# ---------------------------------------------------------------------------
# LocalPrecipitationSource - NetCDF mode
# ---------------------------------------------------------------------------

@requires_netcdf
class TestNetCDFSource:
    """Tests for NetCDF precipitation sources."""

    def test_single_file_with_time(self, netcdf_stack):
        """Test indexing a single NetCDF file through its time coordinate."""
        src = LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326')
        assert src.kind == 'netcdf'
        assert len(src.files) == 1
        assert src.available_months() == FIXTURE_MONTHS
        assert src.year_range == (2001, 2002)
        assert src.shape == (3, 3)
        assert src.resolution == pytest.approx((1.0, 1.0))
        assert src.native_bounds == pytest.approx(BASE_BOUNDS)
        src.close()

    def test_get_month_canonical_form(self, netcdf_stack):
        """Test that NetCDF months come back in the canonical ('y', 'x') form."""
        src = LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326')
        data = src.get_month(2002, 3)
        assert data.dims == ('y', 'x')
        assert data.dtype == np.float32
        assert data.rio.crs.to_epsg() == 4326
        assert np.all(np.diff(data['y'].values) < 0)
        assert data.attrs['year'] == 2002
        assert data.attrs['month'] == 3
        np.testing.assert_allclose(data.values, month_array(2002, 3), equal_nan=True)
        src.close()

    def test_fill_value_becomes_nan(self, netcdf_stack):
        """Test that the NetCDF _FillValue is masked to NaN."""
        src = LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326')
        data = src.get_month(2001, 1)
        assert np.isnan(data.values[1, 1])
        assert int(np.isnan(data.values).sum()) == 1
        src.close()

    def test_missing_month_returns_none(self, netcdf_stack):
        """Test that a month outside the time axis returns None."""
        src = LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326')
        assert src.get_month(2003, 1) is None
        src.close()

    def test_multi_file(self, tmp_path):
        """Test combining several NetCDF files into one source."""
        write_netcdf(build_stack_dataset(FIXTURE_MONTHS[:3]), tmp_path / 'precip_2001.nc')
        write_netcdf(build_stack_dataset(FIXTURE_MONTHS[3:]), tmp_path / 'precip_2002.nc')
        src = LocalPrecipitationSource(tmp_path, pattern='*.nc', crs='EPSG:4326')
        assert len(src.files) == 2
        assert src.available_months() == FIXTURE_MONTHS
        np.testing.assert_allclose(
            src.get_month(2002, 2).values, month_array(2002, 2), equal_nan=True
        )
        src.close()

    def test_explicit_variable(self, tmp_path):
        """Test selecting the precipitation variable explicitly."""
        dataset = build_stack_dataset(FIXTURE_MONTHS[:2], var_name='wanted')
        dataset['ignored'] = dataset['wanted'] * 0.0
        write_netcdf(dataset, tmp_path / 'multi_var.nc')
        src = LocalPrecipitationSource(
            tmp_path / 'multi_var.nc', variable='wanted', crs='EPSG:4326'
        )
        np.testing.assert_allclose(
            src.get_month(2001, 1).values, month_array(2001, 1), equal_nan=True
        )
        src.close()

    def test_unknown_variable_raises(self, netcdf_stack):
        """Test that naming a variable that is absent raises ValueError."""
        with pytest.raises(ValueError, match="Variable 'nope' not found"):
            LocalPrecipitationSource(netcdf_stack, variable='nope', crs='EPSG:4326')

    @pytest.mark.parametrize('var_name', ['precip', 'pr', 'tp', 'RAINNC'])
    def test_variable_auto_detection_by_name(self, tmp_path, var_name):
        """Test that well-known precipitation variable names are auto-detected."""
        write_netcdf(
            build_stack_dataset(FIXTURE_MONTHS[:1], var_name=var_name),
            tmp_path / f'{var_name}.nc'
        )
        src = LocalPrecipitationSource(tmp_path / f'{var_name}.nc', crs='EPSG:4326')
        np.testing.assert_allclose(
            src.get_month(2001, 1).values, month_array(2001, 1), equal_nan=True
        )
        src.close()

    def test_single_unknown_variable_auto_detected(self, tmp_path):
        """Test that a lone unrecognised 2-D variable is used as precipitation."""
        write_netcdf(
            build_stack_dataset(FIXTURE_MONTHS[:1], var_name='wrf_output'),
            tmp_path / 'single.nc'
        )
        src = LocalPrecipitationSource(tmp_path / 'single.nc', crs='EPSG:4326')
        np.testing.assert_allclose(
            src.get_month(2001, 1).values, month_array(2001, 1), equal_nan=True
        )
        src.close()

    def test_ambiguous_variables_raise(self, tmp_path):
        """Test that two unrecognised candidates raise a helpful ValueError."""
        dataset = build_stack_dataset(FIXTURE_MONTHS[:1], var_name='foo')
        dataset['bar'] = dataset['foo'] * 2.0
        write_netcdf(dataset, tmp_path / 'ambiguous.nc')
        with pytest.raises(ValueError, match='Could not determine which variable'):
            LocalPrecipitationSource(tmp_path / 'ambiguous.nc', crs='EPSG:4326')

    def test_lon_lat_dimensions_renamed(self, tmp_path):
        """Test that lon/lat dimensions are auto-renamed to x/y."""
        dataset = build_stack_dataset(FIXTURE_MONTHS[:2], y_name='lat', x_name='lon')
        write_netcdf(dataset, tmp_path / 'lonlat.nc')
        src = LocalPrecipitationSource(tmp_path / 'lonlat.nc', crs='EPSG:4326')
        data = src.get_month(2001, 2)
        assert data.dims == ('y', 'x')
        np.testing.assert_allclose(data['x'].values, BASE_X)
        np.testing.assert_allclose(data.values, month_array(2001, 2), equal_nan=True)
        src.close()

    def test_wrf_dimensions_with_explicit_names(self, tmp_path):
        """Test west_east/south_north dims renamed via x_dim and y_dim."""
        dataset = build_stack_dataset(
            FIXTURE_MONTHS[:2], y_name='south_north', x_name='west_east', ascending_y=True
        )
        write_netcdf(dataset, tmp_path / 'wrf.nc')
        src = LocalPrecipitationSource(
            tmp_path / 'wrf.nc', x_dim='west_east', y_dim='south_north', crs='EPSG:4326'
        )
        data = src.get_month(2001, 1)
        assert data.dims == ('y', 'x')
        # The file stores latitudes ascending; the reader must flip them north-up.
        assert np.all(np.diff(data['y'].values) < 0)
        np.testing.assert_allclose(data.values, month_array(2001, 1), equal_nan=True)
        src.close()

    def test_no_time_coordinate_dates_by_filename(self, tmp_path):
        """Test a one-month-per-file NetCDF layout with no time coordinate."""
        for year, month in FIXTURE_MONTHS[:2]:
            values = np.where(
                np.isnan(month_array(year, month)), np.float32(NODATA), month_array(year, month)
            )
            dataset = xr.Dataset(
                {'precip': (('y', 'x'), values)},
                coords={'y': BASE_Y, 'x': BASE_X}
            )
            dataset['precip'].attrs['_FillValue'] = NODATA
            write_netcdf(dataset, tmp_path / f'Precip_{year}_{month:02d}.nc')

        src = LocalPrecipitationSource(tmp_path, pattern='*.nc', crs='EPSG:4326')
        assert src.kind == 'netcdf'
        assert src.available_months() == FIXTURE_MONTHS[:2]
        np.testing.assert_allclose(
            src.get_month(2001, 2).values, month_array(2001, 2), equal_nan=True
        )
        src.close()

    def test_default_crs_when_absent(self, netcdf_stack):
        """Test that a NetCDF with no CRS falls back to EPSG:4326."""
        src = LocalPrecipitationSource(netcdf_stack)
        assert src.crs.to_epsg() == 4326
        src.close()

    def test_close_is_idempotent(self, netcdf_stack):
        """Test that close() can safely be called more than once."""
        src = LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326')
        src.get_month(2001, 1)
        src.close()
        src.close()

    def test_context_manager_closes(self, netcdf_stack):
        """Test that the context manager yields the source and closes it."""
        with LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326') as src:
            assert isinstance(src, LocalPrecipitationSource)
            assert src.get_month(2001, 1) is not None
        assert src._dataset is None

    def test_context_manager_does_not_swallow_errors(self, netcdf_stack):
        """Test that __exit__ never suppresses an exception."""
        with pytest.raises(RuntimeError, match='boom'):
            with LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326'):
                raise RuntimeError('boom')

    def test_mixed_inputs_raise(self, tmp_path):
        """Test that mixing rasters and NetCDF files raises ValueError."""
        write_raster(tmp_path / 'Precip_2001_01.tif', month_array(2001, 1))
        write_netcdf(build_stack_dataset(FIXTURE_MONTHS[:1]), tmp_path / 'Precip_2001_02.nc')
        with pytest.raises(ValueError, match='Mixed raster and NetCDF inputs'):
            LocalPrecipitationSource(tmp_path, pattern='*')


def build_curvilinear_dataset(promote, decoy=False, coord_names=('XLAT', 'XLONG')):
    """
    Build a WRF-style Dataset with 2-D latitude/longitude arrays.

    Parameters
    ----------

    promote : bool
        If True, name the 2-D arrays in the variable's CF ``coordinates``
        attribute so xarray promotes them to coordinates. If False, leave
        them as plain data variables, which is the layout that used to
        degrade silently to pixel indices.

    decoy : bool
        If True, add an unrelated 2-D field on the same dims, which must
        never be mistaken for a coordinate.

    coord_names : tuple of str
        Names for the (latitude, longitude) arrays.
    """
    lat_name, lon_name = coord_names
    lon_2d, lat_2d = np.meshgrid(BASE_X, BASE_Y)
    values = np.stack([month_array(*ym) for ym in FIXTURE_MONTHS[:2]]).astype(np.float32)
    values = np.where(np.isnan(values), np.float32(NODATA), values)
    data_vars = {
        'RAINNC': (('time', 'south_north', 'west_east'), values),
        lat_name: (('south_north', 'west_east'), lat_2d.astype(np.float32)),
        lon_name: (('south_north', 'west_east'), lon_2d.astype(np.float32)),
    }
    if decoy:
        data_vars['HGT'] = (
            ('south_north', 'west_east'), np.zeros_like(lat_2d, dtype=np.float32)
        )
    times = [np.datetime64(f'{y:04d}-{m:02d}-01') for y, m in FIXTURE_MONTHS[:2]]
    dataset = xr.Dataset(data_vars, coords={'time': times})
    dataset['RAINNC'].attrs['_FillValue'] = NODATA
    if promote:
        dataset['RAINNC'].attrs['coordinates'] = f'{lon_name} {lat_name}'
    return dataset


@requires_netcdf
class TestCurvilinearCoordinates:
    """
    Tests for 2-D XLAT/XLONG coordinate handling.

    A WRF file may name its 2-D coordinate arrays in a CF ``coordinates``
    attribute, in which case xarray promotes them, or leave them as plain
    data variables. Both layouts must produce the same georeferenced grid.
    """

    def _open(self, tmp_path, name, dataset):
        write_netcdf(dataset, tmp_path / name)
        return LocalPrecipitationSource(
            tmp_path / name, variable='RAINNC', crs='EPSG:4326', nodata=NODATA
        )

    def test_2d_coords_as_data_variables_are_georeferenced(self, tmp_path):
        """Test 2-D lat/lon data variables with no CF attribute are still found."""
        src = self._open(tmp_path, 'plain.nc', build_curvilinear_dataset(promote=False))
        # Before the fix this silently returned pixel indices (-0.5 .. n-0.5).
        assert src.native_bounds[0] < -90.0
        assert src.native_bounds[3] > 30.0
        np.testing.assert_allclose(src.get_month(2001, 1)['x'].values, BASE_X, atol=1e-4)
        src.close()

    def test_promoted_and_plain_layouts_agree(self, tmp_path):
        """Test the CF-promoted and plain-data-variable layouts give identical grids."""
        plain = self._open(tmp_path, 'plain.nc', build_curvilinear_dataset(promote=False))
        promoted = self._open(tmp_path, 'promoted.nc', build_curvilinear_dataset(promote=True))
        assert plain.shape == promoted.shape
        assert plain.native_bounds == promoted.native_bounds
        np.testing.assert_allclose(
            plain.get_month(2001, 1).values, promoted.get_month(2001, 1).values,
            equal_nan=True
        )
        plain.close()
        promoted.close()

    def test_unrelated_2d_field_is_not_used_as_coordinate(self, tmp_path):
        """Test a 2-D field on the same dims is not mistaken for a coordinate."""
        src = self._open(
            tmp_path, 'decoy.nc', build_curvilinear_dataset(promote=False, decoy=True)
        )
        np.testing.assert_allclose(src.get_month(2001, 1)['x'].values, BASE_X, atol=1e-4)
        src.close()

    def test_lowercase_lon_lat_data_variables_are_found(self, tmp_path):
        """Test the lower-case lon/lat naming convention is also detected."""
        src = self._open(
            tmp_path, 'lower.nc',
            build_curvilinear_dataset(promote=False, coord_names=('lat', 'lon'))
        )
        np.testing.assert_allclose(src.get_month(2001, 1)['x'].values, BASE_X, atol=1e-4)
        src.close()

    def test_no_coordinates_falls_back_to_pixel_indices(self, tmp_path):
        """Test a file with no coordinates of any kind still opens, ungeoreferenced."""
        dataset = build_curvilinear_dataset(promote=False)
        dataset = dataset.drop_vars(['XLAT', 'XLONG'])
        src = self._open(tmp_path, 'bare.nc', dataset)
        # No coordinates to derive: the grid keeps its shape but is pixel-indexed.
        assert src.native_bounds == (-0.5, -0.5, len(BASE_X) - 0.5, len(BASE_Y) - 0.5)
        src.close()


@requires_netcdf
class TestRasterNetCDFEquivalence:
    """Tests that both on-disk layouts yield identical arrays."""

    @pytest.mark.parametrize('year, month', [(2001, 1), (2002, 3)])
    def test_identical_output(self, raster_dir, netcdf_stack, year, month):
        """Test that the same data read as GeoTIFF and NetCDF matches exactly."""
        with LocalPrecipitationSource(raster_dir, pattern='Precip_*.tif') as raster_src, \
                LocalPrecipitationSource(netcdf_stack, crs='EPSG:4326') as netcdf_src:
            from_raster = raster_src.get_month(year, month)
            from_netcdf = netcdf_src.get_month(year, month)

        assert from_raster.dims == from_netcdf.dims
        assert from_raster.shape == from_netcdf.shape
        assert from_raster.dtype == from_netcdf.dtype
        np.testing.assert_array_equal(from_raster.values, from_netcdf.values)
        np.testing.assert_allclose(from_raster['y'].values, from_netcdf['y'].values)
        np.testing.assert_allclose(from_raster['x'].values, from_netcdf['x'].values)


# ---------------------------------------------------------------------------
# EffectivePrecipitation with local precipitation and no Earth Engine
# ---------------------------------------------------------------------------

class TestEffectivePrecipitationLocal:
    """Tests for the local precipitation workflow in EffectivePrecipitation."""

    def test_construction_never_touches_earth_engine(self, raster_dir, no_earth_engine):
        """Test that a precipitation-only local run skips GEE initialisation."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir),
            local_precip_pattern='Precip_*.tif',
            method='cropwat'
        )
        assert ep._is_local is True
        assert ep._needs_gee is False
        assert ep.geometry is None
        assert ep.collection is None

    def test_awc_driven_method_still_needs_earth_engine(self, raster_dir, no_earth_engine):
        """Test the guard itself: 'ensemble' needs GEE for AWC/ETo even locally."""
        with pytest.raises(AssertionError, match='Earth Engine was initialised'):
            core.EffectivePrecipitation(
                local_precip=str(raster_dir),
                local_precip_pattern='Precip_*.tif',
                method='ensemble'
            )

    def test_years_inferred_from_files(self, raster_dir, no_earth_engine):
        """Test that start_year/end_year default to the range on disk."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif', method='cropwat'
        )
        assert (ep.start_year, ep.end_year) == (2001, 2002)
        # get_date_range() returns an exclusive end date one year past end_year.
        assert (ep.start_date, ep.end_date) == ('2001-01-01', '2003-01-01')

    def test_requested_years_are_clamped(self, raster_dir, no_earth_engine, caplog):
        """Test that a wider requested range is clamped to the available years."""
        with caplog.at_level(logging.WARNING, logger='pycropwat.core'):
            ep = core.EffectivePrecipitation(
                local_precip=str(raster_dir),
                local_precip_pattern='Precip_*.tif',
                start_year=1995,
                end_year=2030,
                method='cropwat'
            )
        assert (ep.start_year, ep.end_year) == (2001, 2002)
        assert 'clamped' in caplog.text

    def test_non_overlapping_years_raise(self, raster_dir, no_earth_engine):
        """Test that a year range with no data at all raises ValueError."""
        with pytest.raises(ValueError, match='do not overlap'):
            core.EffectivePrecipitation(
                local_precip=str(raster_dir),
                local_precip_pattern='Precip_*.tif',
                start_year=1990,
                end_year=1995,
                method='cropwat'
            )

    def test_bounds_from_local_extent(self, raster_dir, no_earth_engine):
        """Test that the ROI falls back to the extent of the local files."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif', method='cropwat'
        )
        ring = np.asarray(ep.bounds, dtype=float)
        assert ring.shape == (5, 2)
        assert ring[:, 0].min() == pytest.approx(BASE_BOUNDS[0])
        assert ring[:, 1].min() == pytest.approx(BASE_BOUNDS[1])
        assert ring[:, 0].max() == pytest.approx(BASE_BOUNDS[2])
        assert ring[:, 1].max() == pytest.approx(BASE_BOUNDS[3])

    def test_process_sequential_outputs(self, raster_dir, tmp_path, no_earth_engine):
        """Test the file names and hand-computed values of a CROPWAT run."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif', method='cropwat'
        )
        output_dir = tmp_path / 'out'
        results = ep.process_sequential(output_dir=output_dir, months=[1])

        assert len(results) == 2
        expected_names = {
            'effective_precip_2001_01.tif', 'effective_precip_fraction_2001_01.tif',
            'effective_precip_2002_01.tif', 'effective_precip_fraction_2002_01.tif'
        }
        assert {path.name for path in output_dir.iterdir()} == expected_names

        precip = month_array(2001, 1)
        peff = read_band(output_dir / 'effective_precip_2001_01.tif').values
        frac = read_band(output_dir / 'effective_precip_fraction_2001_01.tif').values

        np.testing.assert_allclose(peff, expected_cropwat(precip), equal_nan=True, rtol=1e-6)
        # Hand-computed spot checks: 100 mm -> 100 * (125 - 20) / 125 = 84 mm;
        # 300 mm -> 0.1 * 300 + 125 = 155 mm.
        assert peff[0, 0] == pytest.approx(84.0)
        assert peff[0, 1] == pytest.approx(155.0)
        assert frac[0, 0] == pytest.approx(0.84)
        assert frac[0, 1] == pytest.approx(155.0 / 300.0, rel=1e-6)
        # Zero precipitation yields zero effective precipitation and a zero fraction.
        assert peff[0, 2] == pytest.approx(0.0)
        assert frac[0, 2] == pytest.approx(0.0)

    def test_nodata_propagates_as_nan(self, raster_dir, tmp_path, no_earth_engine):
        """Test that nodata pixels become NaN in both outputs, not zero."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif',
            start_year=2001, end_year=2001, method='cropwat'
        )
        output_dir = tmp_path / 'out'
        ep.process_sequential(output_dir=output_dir, months=[1])

        peff_da = read_band(output_dir / 'effective_precip_2001_01.tif')
        frac_da = read_band(output_dir / 'effective_precip_fraction_2001_01.tif')

        assert np.isnan(peff_da.values[1, 1])
        assert np.isnan(frac_da.values[1, 1])
        assert int(np.isnan(peff_da.values).sum()) == 1
        assert int(np.isnan(frac_da.values).sum()) == 1
        assert np.isnan(peff_da.rio.nodata)
        assert np.isnan(frac_da.rio.nodata)
        assert peff_da.rio.crs.to_epsg() == 4326

    def test_save_inputs_writes_precipitation(self, raster_dir, tmp_path, no_earth_engine):
        """Test that save_inputs=True writes precip_YYYY_MM.tif."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif',
            start_year=2001, end_year=2001, method='cropwat'
        )
        input_dir = tmp_path / 'inputs'
        ep.process_sequential(
            output_dir=tmp_path / 'out', months=[1, 2], save_inputs=True, input_dir=input_dir
        )

        assert sorted(path.name for path in input_dir.iterdir()) == [
            'precip_2001_01.tif', 'precip_2001_02.tif'
        ]
        saved = read_band(input_dir / 'precip_2001_01.tif')
        np.testing.assert_allclose(saved.values, month_array(2001, 1), equal_nan=True, rtol=1e-6)

    def test_parallel_matches_sequential(self, raster_dir, tmp_path, no_earth_engine):
        """Test that process() with several workers reproduces process_sequential()."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif', method='cropwat'
        )
        sequential_dir = tmp_path / 'sequential'
        parallel_dir = tmp_path / 'parallel'

        sequential = ep.process_sequential(output_dir=sequential_dir, months=[1, 2])
        parallel = ep.process(output_dir=parallel_dir, n_workers=3, months=[1, 2])

        assert len(sequential) == len(parallel) == 4
        names = sorted(path.name for path in sequential_dir.iterdir())
        assert names == sorted(path.name for path in parallel_dir.iterdir())
        for name in names:
            np.testing.assert_array_equal(
                read_band(sequential_dir / name).values,
                read_band(parallel_dir / name).values
            )

    @requires_netcdf
    def test_netcdf_source_end_to_end(self, netcdf_stack, tmp_path, no_earth_engine):
        """Test a full local run backed by a NetCDF stack, through Dask workers."""
        ep = core.EffectivePrecipitation(
            local_precip=str(netcdf_stack), local_precip_crs='EPSG:4326', method='farmwest'
        )
        assert ep._needs_gee is False
        results = ep.process(output_dir=tmp_path / 'out', n_workers=2, months=[3])
        assert len(results) == 2
        assert all(paths[0] is not None and paths[1] is not None for paths in results)

    def test_missing_months_are_skipped(self, raster_dir, tmp_path, no_earth_engine, caplog):
        """Test that requesting months with no file degrades gracefully."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif', method='cropwat'
        )
        with caplog.at_level(logging.WARNING, logger='pycropwat.core'):
            results = ep.process_sequential(output_dir=tmp_path / 'out', months=[7])

        assert results == []
        assert 'No local precipitation file' in caplog.text
        assert not any((tmp_path / 'out').iterdir())

    def test_partially_available_months(self, raster_dir, tmp_path, no_earth_engine):
        """Test that available months are still processed when others are missing."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif', method='cropwat'
        )
        results = ep.process_sequential(output_dir=tmp_path / 'out', months=[3, 9])
        assert len(results) == 2
        assert {path.name for path in (tmp_path / 'out').iterdir()} == {
            'effective_precip_2001_03.tif', 'effective_precip_fraction_2001_03.tif',
            'effective_precip_2002_03.tif', 'effective_precip_fraction_2002_03.tif'
        }

    def test_clip_to_geometry_shrinks_extent(self, raster_dir, tmp_path, no_earth_engine):
        """Test that a local vector geometry clips the local rasters."""
        geometry_path = write_geojson(tmp_path / 'roi.geojson', -99.9, 38.1, -98.1, 39.9)
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir),
            local_precip_pattern='Precip_*.tif',
            geometry_path=str(geometry_path),
            start_year=2001,
            end_year=2001,
            method='cropwat'
        )
        assert ep._needs_gee is False
        assert ep._clip_gdf is not None

        ep.process_sequential(output_dir=tmp_path / 'out', months=[1])
        clipped = read_band(tmp_path / 'out' / 'effective_precip_2001_01.tif')
        assert clipped.shape == (2, 2)
        np.testing.assert_allclose(
            clipped.values, expected_cropwat(month_array(2001, 1))[:2, :2],
            equal_nan=True, rtol=1e-6
        )

    def test_clip_can_be_disabled(self, raster_dir, tmp_path, no_earth_engine):
        """Test that clip_to_geometry=False keeps the full local extent."""
        geometry_path = write_geojson(tmp_path / 'roi.geojson', -99.9, 38.1, -98.1, 39.9)
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir),
            local_precip_pattern='Precip_*.tif',
            geometry_path=str(geometry_path),
            start_year=2001,
            end_year=2001,
            method='cropwat',
            clip_to_geometry=False
        )
        ep.process_sequential(output_dir=tmp_path / 'out', months=[1])
        assert read_band(tmp_path / 'out' / 'effective_precip_2001_01.tif').shape == (3, 3)

    def test_precip_scale_factor_applied(self, tmp_path, no_earth_engine):
        """Test that precip_scale_factor converts local units before the formula."""
        source_dir = tmp_path / 'metres'
        source_dir.mkdir()
        # 0.1 m == 100 mm, the value the CROPWAT spot check above is built on.
        write_raster(
            source_dir / 'Precip_2001_01.tif',
            np.full((3, 3), 0.1, dtype=np.float32),
            nodata=None
        )
        ep = core.EffectivePrecipitation(
            local_precip=str(source_dir), precip_scale_factor=1000.0, method='cropwat'
        )
        ep.process_sequential(output_dir=tmp_path / 'out', months=[1])
        peff = read_band(tmp_path / 'out' / 'effective_precip_2001_01.tif').values
        np.testing.assert_allclose(peff, np.full((3, 3), 84.0), rtol=1e-5)


class TestLocalValidation:
    """Tests for the argument validation added by the local precipitation feature."""

    def test_pcml_rejects_local_precip(self, raster_dir, no_earth_engine):
        """Test that method='pcml' cannot be combined with local files."""
        with pytest.raises(ValueError, match="local_precip cannot be combined with method='pcml'"):
            core.EffectivePrecipitation(
                local_precip=str(raster_dir),
                local_precip_pattern='Precip_*.tif',
                method='pcml'
            )

    def test_no_precipitation_source_raises(self, no_earth_engine):
        """Test that omitting every source names both ways of supplying one."""
        with pytest.raises(ValueError, match='No precipitation source provided') as excinfo:
            core.EffectivePrecipitation(method='cropwat')
        message = str(excinfo.value)
        assert 'asset_id and precip_band' in message
        assert 'local_precip' in message

    def test_asset_without_band_raises(self, no_earth_engine):
        """Test that asset_id without precip_band blames the band, not the source."""
        with pytest.raises(
            ValueError,
            match=re.escape('precip_band is required when asset_id is given')
        ) as excinfo:
            core.EffectivePrecipitation(asset_id='IDAHO_EPSCOR/TERRACLIMATE', method='cropwat')
        message = str(excinfo.value)
        # The asset that *was* supplied is named, with an example band for it.
        assert "asset 'IDAHO_EPSCOR/TERRACLIMATE' was provided" in message
        assert "precip_band='pr'" in message
        # The old wording claimed no source at all, which sent users asset hunting.
        assert 'No precipitation source provided' not in message

    def test_asset_without_band_unknown_asset_raises(self, no_earth_engine):
        """Test that an unrecognised asset still gets the band-specific message."""
        with pytest.raises(
            ValueError,
            match=re.escape('precip_band is required when asset_id is given')
        ) as excinfo:
            core.EffectivePrecipitation(
                asset_id='projects/my-project/assets/my_precip', method='cropwat'
            )
        message = str(excinfo.value)
        assert "asset 'projects/my-project/assets/my_precip' was provided" in message
        # No band is known for it, so the user is pointed at the catalog instead.
        assert 'Earth Engine Data Catalog' in message

    def test_band_without_asset_raises(self, no_earth_engine):
        """Test that precip_band without asset_id blames the asset, not the source."""
        with pytest.raises(
            ValueError,
            match=re.escape('asset_id is required when precip_band is given')
        ) as excinfo:
            core.EffectivePrecipitation(precip_band='total_precipitation_sum', method='cropwat')
        message = str(excinfo.value)
        assert "band 'total_precipitation_sum' was provided" in message
        # An asset documented as carrying that band is offered as the example.
        assert "asset_id='ECMWF/ERA5_LAND/MONTHLY_AGGR'" in message
        assert 'No precipitation source provided' not in message

    def test_bad_local_path_raises(self, tmp_path, no_earth_engine):
        """Test that a non-existent local_precip path fails fast."""
        with pytest.raises(FileNotFoundError):
            core.EffectivePrecipitation(
                local_precip=str(tmp_path / 'missing_dir'), method='cropwat'
            )

    def test_local_source_is_exposed(self, raster_dir, no_earth_engine):
        """Test that the opened source and its native scale are available."""
        ep = core.EffectivePrecipitation(
            local_precip=str(raster_dir), local_precip_pattern='Precip_*.tif', method='cropwat'
        )
        assert isinstance(ep._local_source, LocalPrecipitationSource)
        # 1 degree cells in a geographic CRS -> roughly 111 km.
        assert ep._get_native_scale() == pytest.approx(111320.0)

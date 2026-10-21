"""Regression checks for GeodePy #35, not a claim of physical GPS accuracy."""

import math

import pytest
from geographiclib.geodesic import Geodesic

from geodepy import geodesy
from geodepy.angles import dec2dms, dec2ddm
from geodepy.constants import grs80, wgs84, ans, intl24
from geodepy.convert import geo2grid, grid2geo


@pytest.mark.parametrize("latitude", [-75.0, -35.0, 0.0, 35.0, 75.0])
@pytest.mark.parametrize("longitude,bearing", [(179.9, 90.0), (-179.9, 270.0)])
def test_crossing_in_both_directions(latitude, longitude, bearing):
    actual = geodesy.vincdir(latitude, longitude, bearing, 100000.0)
    assert -180.0 <= actual[1] <= 180.0
    ref = Geodesic(grs80.semimaj, grs80.f)
    expected = ref.Direct(latitude, longitude, bearing, 100000.0)
    assert ref.Inverse(*actual[:2], expected["lat2"], expected["lon2"])["s12"] <= 0.001


@pytest.mark.parametrize("longitude", [-180.0, 180.0])
@pytest.mark.parametrize("latitude", [-80.0, -35.0, 0.0, 35.0, 84.0])
def test_sub_rounding_boundary_motion_keeps_existing_endpoint(latitude, longitude):
    # 1 mm northward with a tiny west component. Raw longitude can cross -180
    # by less than the existing eleven-decimal-degree return precision.
    actual = geodesy.vincdir(latitude, longitude, 359.999, 0.001)
    assert actual[1] == longitude


@pytest.mark.parametrize(
    "longitude",
    [-180.0, math.nextafter(-180.0, 0.0), math.nextafter(180.0, 0.0), 180.0],
)
@pytest.mark.parametrize("bearing", [90.0, 270.0])
def test_neighbouring_floats(longitude, bearing):
    for distance in [0.0, 0.001, 1.0]:
        actual = geodesy.vincdir(35.0, longitude, bearing, distance)
        assert -180.0 <= actual[1] <= 180.0


@pytest.mark.parametrize("adapter", [float, dec2dms, dec2ddm])
@pytest.mark.parametrize("longitude,bearing", [(179.9, 90.0), (-179.9, 270.0)])
def test_angle_objects(adapter, longitude, bearing):
    expected = geodesy.vincdir(35.0, longitude, bearing, 100000.0)
    actual = geodesy.vincdir(
        adapter(35.0), adapter(longitude), adapter(bearing), 100000.0
    )
    assert actual == expected
    assert -180.0 <= actual[1] <= 180.0


@pytest.mark.parametrize("ellipsoid", [grs80, wgs84, ans, intl24])
def test_other_ellipsoids(ellipsoid):
    ref = Geodesic(ellipsoid.semimaj, ellipsoid.f)
    for lon, az in [(179.9, 90.0), (-179.9, 270.0)]:
        actual = geodesy.vincdir(-35.0, lon, az, 100000.0, ellipsoid=ellipsoid)
        expected = ref.Direct(-35.0, lon, az, 100000.0)
        assert -180.0 <= actual[1] <= 180.0
        assert (
            ref.Inverse(*actual[:2], expected["lat2"], expected["lon2"])["s12"] <= 0.001
        )


@pytest.mark.parametrize("latitude", [-75.0, -35.0, 0.0, 35.0, 75.0])
@pytest.mark.parametrize("longitude,bearing", [(179.9, 90.0), (-179.9, 270.0)])
def test_actual_auto_grid_conversion_chain(latitude, longitude, bearing):
    lat, lon, _ = geodesy.vincdir(latitude, longitude, bearing, 100000.0)
    hemi, zone, east, north, *_ = geo2grid(lat, lon)
    assert 1 <= zone <= 60
    round_lat, round_lon, *_ = grid2geo(zone, east, north, hemisphere=hemi)
    ref = Geodesic(grs80.semimaj, grs80.f)
    assert ref.Inverse(lat, lon, round_lat, round_lon)["s12"] <= 0.001

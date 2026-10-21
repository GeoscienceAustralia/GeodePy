"""Regression tests for the UTM zone limit at positive 180 degrees."""

import math

import pytest

from geodepy.convert import geo2grid, grid2geo


@pytest.mark.parametrize("latitude", [-75.0, -35.0, 0.0, 35.0, 75.0])
@pytest.mark.parametrize("longitude", [math.nextafter(180.0, -math.inf), 180.0])
def test_positive_antimeridian_selects_zone_60(latitude, longitude):
    # The same location using an explicit valid zone already works.
    actual = geo2grid(latitude, longitude)
    assert actual == geo2grid(latitude, longitude, zone=60)
    assert actual[1] == 60


@pytest.mark.parametrize("latitude", [-75.0, -35.0, 0.0, 35.0, 75.0])
@pytest.mark.parametrize("longitude", [-180.0, math.nextafter(-180.0, math.inf)])
def test_negative_antimeridian_remains_zone_1(latitude, longitude):
    actual = geo2grid(latitude, longitude)
    assert actual == geo2grid(latitude, longitude, zone=1)
    assert actual[1] == 1


@pytest.mark.parametrize("latitude", [-35.0, 0.0, 35.0])
@pytest.mark.parametrize("longitude", [179.9, -179.9, 174.5, -174.5, 0.0])
def test_unaffected_locations_keep_their_zone(latitude, longitude):
    zone = int((longitude + 180.0) // 6.0) + 1
    assert geo2grid(latitude, longitude) == geo2grid(
        latitude, longitude, zone=zone
    )


@pytest.mark.parametrize("latitude", [-35.0, 0.0, 35.0])
@pytest.mark.parametrize(
    "longitude", [180.0, math.nextafter(180.0, -math.inf), -180.0]
)
def test_antimeridian_grid_round_trip(latitude, longitude):
    hemisphere, zone, east, north, *_ = geo2grid(latitude, longitude)
    lat2, lon2, *_ = grid2geo(zone, east, north, hemisphere=hemisphere)
    assert abs(lat2 - latitude) < 1e-8
    longitude_error = ((lon2 - longitude + 180.0) % 360.0) - 180.0
    assert abs(longitude_error) < 1e-8


def test_out_of_range_longitudes_remain_invalid():
    with pytest.raises(ValueError, match="Invalid Longitude"):
        geo2grid(0.0, math.nextafter(180.0, math.inf))
    with pytest.raises(ValueError, match="Invalid Longitude"):
        geo2grid(0.0, math.nextafter(-180.0, -math.inf))

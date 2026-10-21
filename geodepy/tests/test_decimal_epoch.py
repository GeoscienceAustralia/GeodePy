"""Regression tests for decimal-year epoch support (GeodePy #104)."""

from datetime import date

import pytest

from geodepy.constants import _epoch_datetime, itrf2014_to_gda2020
from geodepy.transform import (
    conform14,
    plate_motion_transformation,
    transform_atrf2014_to_gda2020,
)

COORD = (-4052052.7373, 4212835.9835, -2545104.5867)
MODEL = itrf2014_to_gda2020


@pytest.mark.parametrize("year", [1994, 2000, 2018, 2020, 2024, 2025, 2026])
def test_integer_decimal_year_matches_january_first(year):
    assert conform14(*COORD, float(year), MODEL) == conform14(
        *COORD, date(year, 1, 1), MODEL
    )


@pytest.mark.parametrize(
    "from_year,to_year", [(2018, 2025), (2025, 2018), (2020, 2020), (2024, 2026)]
)
def test_plate_motion_integer_years_match_dates(from_year, to_year):
    expected = plate_motion_transformation(
        *COORD, date(from_year, 1, 1), date(to_year, 1, 1), MODEL
    )
    assert plate_motion_transformation(
        *COORD, float(from_year), float(to_year), MODEL
    ) == expected


@pytest.mark.parametrize(
    "from_epoch,to_epoch",
    [(2018.0, date(2025, 1, 1)), (date(2018, 1, 1), 2025.0)],
)
def test_mixed_input_types(from_epoch, to_epoch):
    expected = plate_motion_transformation(
        *COORD, date(2018, 1, 1), date(2025, 1, 1), MODEL
    )
    assert plate_motion_transformation(*COORD, from_epoch, to_epoch, MODEL) == expected


@pytest.mark.parametrize(
    "year,fraction,expected",
    [
        (2024, 0.0, (2024, 1, 1, 0)),
        (2024, 0.5, (2024, 7, 2, 0)),
        (2025, 0.5, (2025, 7, 2, 12)),
        (2026, 0.0, (2026, 1, 1, 0)),
        (2026, 0.5, (2026, 7, 2, 12)),
    ],
)
def test_decimal_calendar_fraction(year, fraction, expected):
    observed = _epoch_datetime(year + fraction)
    assert (observed.year, observed.month, observed.day, observed.hour) == expected


def test_leap_halfway_matches_same_calendar_date():
    assert conform14(*COORD, 2024.5, MODEL) == conform14(
        *COORD, date(2024, 7, 2), MODEL
    )


def test_fractional_year_affects_transform():
    jan = conform14(*COORD, 2025.0, MODEL)
    mid = conform14(*COORD, 2025.5, MODEL)
    next_year = conform14(*COORD, 2026.0, MODEL)
    for a, b, c in zip(jan[:3], mid[:3], next_year[:3]):
        assert abs(a - c) > 0
        assert min(a, c) <= b <= max(a, c)


def test_fractional_plate_motion_reverse_round_trip():
    endpoint = plate_motion_transformation(*COORD, 2018.25, 2025.75, MODEL)
    reverse = plate_motion_transformation(*endpoint[:3], 2025.75, 2018.25, MODEL)
    for actual, original in zip(reverse[:3], COORD):
        assert abs(actual - original) < 0.002


def test_wrapper_supports_decimal_year():
    assert transform_atrf2014_to_gda2020(*COORD, 2018.0) == (
        transform_atrf2014_to_gda2020(*COORD, date(2018, 1, 1))
    )


@pytest.mark.parametrize(
    "bad", [float("nan"), float("inf"), float("-inf"), 0.5, -2026.1, 10000.0]
)
def test_invalid_year_fails(bad):
    with pytest.raises(ValueError):
        conform14(*COORD, bad, MODEL)


@pytest.mark.parametrize("bad", [True, False, "2025.5", None])
def test_unsupported_year_type_fails(bad):
    with pytest.raises(ValueError):
        conform14(*COORD, bad, MODEL)


@pytest.mark.parametrize(
    "invalid",
    [float("nan"), float("inf"), float("-inf"), True, "2025.5", None, 0.1, 10000.0],
)
@pytest.mark.parametrize("slot", [0, 1])
def test_plate_motion_invalid_decimal_years(invalid, slot):
    epochs = [2018.25, 2025.75]
    epochs[slot] = invalid
    with pytest.raises(ValueError):
        plate_motion_transformation(*COORD, *epochs, MODEL)


def test_invalid_transform_addition_fails_loudly():
    with pytest.raises(ValueError):
        MODEL + float("nan")
    with pytest.raises(ValueError):
        MODEL + "2018.0"


def test_ref_epoch_preserved_on_date_path():
    result = MODEL + date(2026, 5, 11)
    assert result.ref_epoch == date(2026, 5, 11)
    assert type(result.ref_epoch) is date


def test_ref_epoch_is_numeric_on_numeric_path():
    result = MODEL + 2026.25
    assert result.ref_epoch == 2026.25

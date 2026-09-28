import struct
import tempfile
import unittest
from pathlib import Path

from geodepy.ntv2reader import read_ntv2_file
from geodepy.transform import ntv2_2d

# NTv2 longitudes are positive west and grid coordinates are in arcseconds.
SOUTH = -32.0 * 3600
EAST = -142.0 * 3600
INTERVAL = 360.0
SIZE = 7
METHODS = ("bilinear", "bicubic")
# Interpolated shifts are rounded to six decimal places in arcseconds.
TOLERANCE = 1e-6 / 3600 + 1e-13


def shifts_at(y, x, sign=1):
    """Affine fields with exactly representable float32 values at grid nodes."""
    return sign * (0.5 + y / 16 + x / 32), sign * (-0.25 + y / 64 - x / 32)


def coordinate_at(y, x):
    return (SOUTH + y * INTERVAL) / 3600, -(EAST + x * INTERVAL) / 3600


def write_grid(path, sign=1):
    """Write a small NTv2 file so tests exercise the actual interpolators."""

    def record(name, value):
        if isinstance(value, str):
            data = value.encode("ascii").ljust(8)
        elif isinstance(value, int):
            data = struct.pack("<i4x", value)
        else:
            data = struct.pack("<d", value)
        return name.encode("ascii").ljust(8) + data

    headers = (
        ("NUM_OREC", 11),
        ("NUM_SREC", 11),
        ("NUM_FILE", 1),
        ("GS_TYPE", "SECONDS"),
        ("VERSION", "TEST"),
        ("SYSTEM_F", "SOURCE"),
        ("SYSTEM_T", "TARGET"),
        ("MAJOR_F", 6378137.0),
        ("MINOR_F", 6356752.314),
        ("MAJOR_T", 6378137.0),
        ("MINOR_T", 6356752.314),
        ("SUB_NAME", "AFFINE"),
        ("PARENT", "NONE"),
        ("CREATED", "01012026"),
        ("UPDATED", "01012026"),
        ("S_LAT", SOUTH),
        ("N_LAT", SOUTH + (SIZE - 1) * INTERVAL),
        ("E_LONG", EAST),
        ("W_LONG", EAST + (SIZE - 1) * INTERVAL),
        ("LAT_INC", INTERVAL),
        ("LONG_INC", INTERVAL),
        ("GS_COUNT", SIZE * SIZE),
    )
    with path.open("wb") as stream:
        for name, value in headers:
            stream.write(record(name, value))
        for y in range(SIZE):
            for x in range(SIZE):
                stream.write(struct.pack("<4f", *shifts_at(y, x, sign), 0.0, 0.0))


class TestNTv2Transform(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.grids = {}
        for sign in (-1, 1):
            path = Path(directory.name) / f"affine_{sign}.gsb"
            write_grid(path, sign)
            self.grids[sign] = read_ntv2_file(str(path))

    def assertCoordinateAlmostEqual(self, actual, expected):
        for value, reference in zip(actual, expected):
            self.assertAlmostEqual(value, reference, delta=TOLERANCE)

    def test_forward_matches_analytic_field(self):
        for method in METHODS:
            for sign, grid in self.grids.items():
                with self.subTest(method=method, sign=sign):
                    lat, lon = coordinate_at(2.25, 2.75)
                    dlat, dlon = shifts_at(2.25, 2.75, sign)
                    expected = lat + dlat / 3600, lon - dlon / 3600
                    self.assertCoordinateAlmostEqual(
                        ntv2_2d(grid, lat, lon, method=method), expected
                    )

    def test_reverse_matches_known_source(self):
        # Construct targets analytically, independently of ntv2_2d's forward path.
        # The second point crosses an interpolation cell boundary when shifted.
        for method in METHODS:
            for sign, grid in self.grids.items():
                for y, x in ((2.25, 2.75), (2.9999, 3.0001)):
                    with self.subTest(method=method, sign=sign, y=y, x=x):
                        source = coordinate_at(y, x)
                        dlat, dlon = shifts_at(y, x, sign)
                        target = source[0] + dlat / 3600, source[1] - dlon / 3600
                        self.assertCoordinateAlmostEqual(
                            ntv2_2d(grid, *target, forward_tf=False, method=method),
                            source,
                        )

    def test_repeated_round_trips(self):
        for method in METHODS:
            for sign, grid in self.grids.items():
                with self.subTest(method=method, sign=sign):
                    source = coordinate_at(2.25, 2.75)
                    result = source
                    for _ in range(10):
                        target = ntv2_2d(grid, *result, method=method)
                        result = ntv2_2d(grid, *target, forward_tf=False, method=method)
                    self.assertCoordinateAlmostEqual(result, source)

    def test_reverse_estimate_outside_grid(self):
        # The input is in the grid, but subtracting the shift leaves its south edge.
        # Bilinear interpolation is supported here without a bicubic halo.
        with self.assertRaisesRegex(ValueError, "Coordinate outside of grid extents"):
            ntv2_2d(
                self.grids[1],
                *coordinate_at(0.0001, 2.75),
                forward_tf=False,
                method="bilinear",
            )


if __name__ == "__main__":
    unittest.main()

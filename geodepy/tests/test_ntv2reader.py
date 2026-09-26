"""NTv2 interpolation regression tests using small, generated binary grids."""

import struct
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from geodepy.ntv2reader import interpolate_ntv2, read_ntv2_file
from geodepy.transform import ntv2_2d


def write_grid(path, sizes=(6,), nonlinear=False):
    """Write valid little-endian NTv2 headers and independent planar fields."""

    def record(name, value):
        if isinstance(value, str):
            data = value.encode("ascii").ljust(8, b" ")
        elif isinstance(value, int):
            data = struct.pack("<i4x", value)
        else:
            data = struct.pack("<d", value)
        return name.encode("ascii").ljust(8, b" ") + data

    with path.open("wb") as stream:
        for name, value in (
            ("NUM_OREC", 11),
            ("NUM_SREC", 11),
            ("NUM_FILE", len(sizes)),
            ("GS_TYPE", "SECONDS"),
            ("VERSION", "TEST"),
            ("SYSTEM_F", "FROM"),
            ("SYSTEM_T", "TO"),
            ("MAJOR_F", 6378137.0),
            ("MINOR_F", 6356752.314),
            ("MAJOR_T", 6378137.0),
            ("MINOR_T", 6356752.314),
        ):
            stream.write(record(name, value))
        for index, size in enumerate(sizes):
            south = float(index * 10)
            for name, value in (
                ("SUB_NAME", f"GRID{index}"),
                ("PARENT", "NONE"),
                ("CREATED", "01012026"),
                ("UPDATED", "01012026"),
                ("S_LAT", south),
                ("N_LAT", south + size - 1),
                ("E_LONG", 0.0),
                ("W_LONG", float(size - 1)),
                ("LAT_INC", 1.0),
                ("LONG_INC", 1.0),
                ("GS_COUNT", size * size),
            ):
                stream.write(record(name, value))
            for y in range(size):
                for x in range(size):
                    first = x * x + y * y if nonlinear else x + 10 * y
                    stream.write(
                        struct.pack(
                            "<4f", first + index * 1000, 2 * x + 3 * y, x - y, 1.0
                        )
                    )


class TestNTv2Boundaries(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "test.gsb"

    def grid(self, sizes=(6,), nonlinear=False):
        write_grid(self.path, sizes, nonlinear)
        return read_ntv2_file(self.path)

    def test_bicubic_rejects_all_edge_bands_and_corners(self):
        grid = self.grid((6, 6))
        # Exercise first and last subgrids: neither headers nor neighbours
        # may be used to fill in the missing interpolation stencil.
        for index in range(2):
            for x, y in (
                (2.5, 0.5),
                (2.5, 4.5),
                (0.5, 2.5),
                (4.5, 2.5),
                (0.5, 0.5),
                (0.5, 4.5),
                (4.5, 0.5),
                (4.5, 4.5),
                (0.0, 2.5),
                (2.5, 0.0),
            ):
                with self.subTest(grid=index, x=x, y=y):
                    lat, lon = (index * 10 + y) / 3600, -x / 3600
                    with self.assertRaisesRegex(ValueError, "bilinear"):
                        interpolate_ntv2(grid, lat, lon)
                    self.assertEqual(
                        interpolate_ntv2(grid, lat, lon, method="bilinear"),
                        (x + 10 * y + index * 1000, 2 * x + 3 * y, x - y, 1.0),
                    )

    def test_bicubic_guard_runs_before_file_access(self):
        subgrid = self.grid().subgrids["GRID0"]
        for row, col in ((0, 2), (4, 2), (2, 0), (2, 4)):
            with self.subTest(row=row, col=col):
                stream = Mock()
                with self.assertRaisesRegex(ValueError, "bilinear"):
                    subgrid.ntv2_bicubic(row + 0.5, col + 0.5, 6, row, col, stream, 352)
                self.assertEqual(stream.mock_calls, [])

    def test_small_grids_support_bilinear_only(self):
        for size in (2, 3):
            with self.subTest(size=size):
                grid = self.grid((size,))
                with self.assertRaisesRegex(ValueError, "bilinear"):
                    interpolate_ntv2(grid, 0.5 / 3600, -0.5 / 3600)
                self.assertEqual(
                    interpolate_ntv2(grid, 0.5 / 3600, -0.5 / 3600, "bilinear"),
                    (5.5, 2.5, 0.0, 1.0),
                )

    def test_interior_retains_bicubic_interpolation(self):
        grid = self.grid(nonlinear=True)
        for x, y in ((1.0, 1.0), (1.5, 1.5), (2.5, 2.5), (3.5, 3.5)):
            with self.subTest(x=x, y=y):
                self.assertEqual(
                    interpolate_ntv2(grid, y / 3600, -x / 3600),
                    (x * x + y * y, 2 * x + 3 * y, x - y, 1.0),
                )
        self.assertNotEqual(
            interpolate_ntv2(grid, 2.5 / 3600, -2.5 / 3600),
            interpolate_ntv2(grid, 2.5 / 3600, -2.5 / 3600, "bilinear"),
        )

    def test_minimum_bicubic_grid(self):
        grid = self.grid((4,), nonlinear=True)
        self.assertEqual(
            interpolate_ntv2(grid, 1.5 / 3600, -1.5 / 3600), (4.5, 7.5, 0.0, 1.0)
        )

    def test_north_and_west_extent_exclusions_are_unchanged(self):
        grid = self.grid()
        for lat, lon in ((5 / 3600, -2.5 / 3600), (2.5 / 3600, -5 / 3600)):
            for method in ("bicubic", "bilinear"):
                with self.subTest(lat=lat, lon=lon, method=method):
                    self.assertEqual(
                        interpolate_ntv2(grid, lat, lon, method),
                        (None, None, None, None),
                    )

    def test_transform_exposes_boundary_error_and_bilinear_alternative(self):
        grid = self.grid()
        lat, lon = 0.5 / 3600, -2.5 / 3600
        with self.assertRaisesRegex(ValueError, "bilinear"):
            ntv2_2d(grid, lat, lon)
        self.assertEqual(
            ntv2_2d(grid, lat, lon, method="bilinear"),
            (lat + 7.5 / 3600, lon - 6.5 / 3600),
        )

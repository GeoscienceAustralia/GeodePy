import unittest

from geodepy.transform import (geo2grid,
                               grid2geo,
                               llh2xyz,
                               xyz2llh,
                               conform7,
                               conform14,
                               atrftogda2020,
                               gda2020toatrf,
                               mga94to2020,
                               mga2020to94)
from geodepy.convert import DMSAngle, hp2dec_v
from geodepy.constants import itrf14togda20, gda94to20
from geodepy.fileio import read_dnacoord
from datetime import date
import numpy as np
import os.path


class TestTransforms(unittest.TestCase):
    # Tests equality between values produced using grid2geo and geo2grid
    # TODO: Change Source Data input to natadjust_rvs_example.dat
    def test_geo_grid_transform_interoperability(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        test_geo_coords = np.genfromtxt(os.path.join(abs_path, 'resources/Test_Conversion_Geo.csv'),
                                        delimiter=',',
                                        dtype='S4,f8,f8',
                                        names=['site', 'lat', 'lon'])

        test_grid_coords = np.genfromtxt(os.path.join(abs_path, 'resources/Test_Conversion_Grid.csv'),
                                         delimiter=',',
                                         dtype='S4,i4,f8,f8',
                                         names=['site', 'zone', 'east', 'north'])

        geoed_grid = np.array(list(grid2geo(*x) for x in test_grid_coords[['zone', 'east', 'north']]))
        np.testing.assert_almost_equal(geoed_grid[:, :2], hp2dec_v(np.array(test_geo_coords[['lat', 'lon']].tolist())),
                                       decimal=8)

        gridded_geo = np.stack(geo2grid(*x) for x in hp2dec_v(np.array(test_geo_coords[['lat', 'lon']].tolist())))
        np.testing.assert_almost_equal(gridded_geo[:, 2:4].astype(float),
                                       np.array(test_grid_coords[['east', 'north']].tolist()),
                                       decimal=3)

    # Tests comparing results from DynAdjust Adjusted Coordinate File including 109 Points Across Australia
    def test_llh2xyz(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        testdata = read_dnacoord(os.path.join(abs_path, 'resources/natadjust_rvs_example.dat'))
        for coord in testdata:
            coord.converthptodd()
            xcomp, ycomp, zcomp = llh2xyz(coord.lat, coord.long, coord.ell_ht)
            assert (abs(coord.x - xcomp) < 2e-4)
            assert (abs(coord.y - ycomp) < 2e-4)
            assert (abs(coord.z - zcomp) < 2e-4)

    def test_xyz2llh(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        testdata = read_dnacoord(os.path.join(abs_path, 'resources/natadjust_rvs_example.dat'))
        for coord in testdata:
            coord.converthptodd()
            latcomp, longcomp, ell_htcomp = xyz2llh(coord.x, coord.y, coord.z)
            assert (abs(latcomp - coord.lat) < 2e-9)
            assert (abs(longcomp - coord.long) < 2e-9)
            assert (abs(ell_htcomp - coord.ell_ht) < 1e-4)

    def test_geo2grid(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        # Test various coordinates in Australia
        testdata = read_dnacoord(os.path.join(abs_path, 'resources/natadjust_rvs_example.dat'))
        for coord in testdata:
            coord.converthptodd()
            hem, zonecomp, eastcomp, northcomp, psf, grid_conv = geo2grid(coord.lat, coord.long)
            self.assertEqual(zonecomp, coord.zone)
            self.assertLess(abs(eastcomp - coord.easting), 4e-4)
            self.assertLess((northcomp - coord.northing), 4e-4)

        # Test North and South Hemisphere Output
        north_ex = (DMSAngle(34, 57, 00.79653).dec(), DMSAngle(117, 48, 36.68783).dec())
        south_ex = (DMSAngle(-34, 57, 00.79653).dec(), DMSAngle(117, 48, 36.68783).dec())
        north_grid = geo2grid(north_ex[0], north_ex[1])
        south_grid = geo2grid(south_ex[0], south_ex[1])
        self.assertEqual(north_grid[0], 'North')
        self.assertEqual(south_grid[0], 'South')
        self.assertEqual(north_grid[1], south_grid[1])  # Zone
        self.assertEqual(north_grid[2], south_grid[2])  # Easting
        self.assertEqual(north_grid[3], 10000000 - south_grid[3])  # Northing
        self.assertEqual(north_grid[4], south_grid[4])  # PSF
        self.assertEqual(north_grid[5], -south_grid[5])  # Grid Convergence

        # Test Input Validation
        with self.assertRaises(ValueError):
            geo2grid(0, 45, -1)
        with self.assertRaises(ValueError):
            geo2grid(0, 45, 61)
        with self.assertRaises(ValueError):
            geo2grid(-81, 45, 0)
        with self.assertRaises(ValueError):
            geo2grid(85, 45, 0)
        with self.assertRaises(ValueError):
            geo2grid(0, -181, 0)
        with self.assertRaises(ValueError):
            geo2grid(0, 181, 0)

    def test_grid2geo(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        testdata = read_dnacoord(os.path.join(abs_path, 'resources/natadjust_rvs_example.dat'))
        for coord in testdata:
            coord.converthptodd()
            latcomp, longcomp, psf, grid_conv = grid2geo(coord.zone, coord.easting, coord.northing)
            self.assertLess(abs(latcomp - coord.lat), 5e-9)
            self.assertLess(abs(longcomp - coord.long), 5e-9)

        # Test North and South Hemisphere Output
        north_ex = (50, 573976.8747, 3867822.4539, 'North')
        south_ex = (50, 573976.8747, 6132177.5461, 'South')
        north_geo = grid2geo(*north_ex)
        south_geo = grid2geo(*south_ex)
        self.assertEqual(north_geo[0], -south_geo[0])
        self.assertEqual(north_geo[1], south_geo[1])
        self.assertEqual(north_geo[2], south_geo[2])
        self.assertEqual(north_geo[3], -south_geo[3])

        # Test Input Validation
        with self.assertRaises(ValueError):
            grid2geo(-1, 0, 500000)
        with self.assertRaises(ValueError):
            grid2geo(61, 0, 500000)
        with self.assertRaises(ValueError):
            grid2geo(0, -2830001, 500000)
        with self.assertRaises(ValueError):
            grid2geo(0, 3830001, 500000)
        with self.assertRaises(ValueError):
            grid2geo(0, 0, -1)
        with self.assertRaises(ValueError):
            grid2geo(0, 0, 10000001)
        with self.assertRaises(ValueError):
            grid2geo(0, 0, 500000, 'fail')

    def test_conform7(self):
        # Replication of tests in GDA2020 Tech Manual v1.2 - Sect 3.1.1
        alic_gda1994 = (-4052051.7643, 4212836.2017, -2545106.0245)
        alic_gda2020 = (-4052052.7379, 4212835.9897, -2545104.5898)
        alic_gda2020_comp = conform7(*alic_gda1994, gda94to20)
        alic_gda1994_comp = conform7(*alic_gda2020, -gda94to20)
        assert (abs(alic_gda2020_comp[0] - alic_gda2020[0]) < 5e-5)
        assert (abs(alic_gda2020_comp[1] - alic_gda2020[1]) < 5e-5)
        assert (abs(alic_gda2020_comp[2] - alic_gda2020[2]) < 5e-5)
        assert (abs(alic_gda1994_comp[0] - alic_gda1994[0]) < 5e-5)
        assert (abs(alic_gda1994_comp[1] - alic_gda1994[1]) < 5e-5)
        assert (abs(alic_gda1994_comp[2] - alic_gda1994[2]) < 5e-5)

    def test_conform14(self):
        # Replication of tests in GDA2020 Tech Manual v1.2 - Sect 3.3.1
        alic_gda2020 = (-4052052.7373, 4212835.9835, -2545104.5867)
        alic_itrf14at2018 = (-4052052.6588, 4212835.9938, -2545104.6946)
        alic_itrf14at2018_comp = conform14(*alic_gda2020, date(2018, 1, 1), -itrf14togda20)
        alic_gda2020_comp = conform14(*alic_itrf14at2018, date(2018, 1, 1), itrf14togda20)
        assert (abs(alic_itrf14at2018_comp[0] - alic_itrf14at2018[0]) < 5e-5)
        assert (abs(alic_itrf14at2018_comp[1] - alic_itrf14at2018[1]) < 5e-5)
        assert (abs(alic_itrf14at2018_comp[2] - alic_itrf14at2018[2]) < 5e-5)
        assert (abs(alic_gda2020_comp[0] - alic_gda2020[0]) < 5e-5)
        assert (abs(alic_gda2020_comp[1] - alic_gda2020[1]) < 5e-5)
        assert (abs(alic_gda2020_comp[2] - alic_gda2020[2]) < 5e-5)

    def test_atrftogda2020(self):
        alic_gda2020 = (-4052052.7373, 4212835.9835, -2545104.5867)
        alic_itrf14at2018 = (-4052052.6588, 4212835.9938, -2545104.6946)
        alic_itrf14at2018_comp = gda2020toatrf(*alic_gda2020, date(2018, 1, 1))
        alic_gda2020_comp = atrftogda2020(*alic_itrf14at2018, date(2018, 1, 1))
        assert (abs(alic_itrf14at2018_comp[0] - alic_itrf14at2018[0]) < 5e-5)
        assert (abs(alic_itrf14at2018_comp[1] - alic_itrf14at2018[1]) < 5e-5)
        assert (abs(alic_itrf14at2018_comp[2] - alic_itrf14at2018[2]) < 5e-5)
        assert (abs(alic_gda2020_comp[0] - alic_gda2020[0]) < 5e-5)
        assert (abs(alic_gda2020_comp[1] - alic_gda2020[1]) < 5e-5)
        assert (abs(alic_gda2020_comp[2] - alic_gda2020[2]) < 5e-5)

    def test_mga94to2020(self):
        alic_mga94 = (53, 386352.3979, 7381850.7689, 603.3466)
        alic_mga20 = (53, 386353.2343, 7381852.2986, 603.2489)
        # Test with no ellipsoid height supplied
        alic_mga20_noellht_comp = mga94to2020(alic_mga94[0], alic_mga94[1], alic_mga94[2])
        assert ((alic_mga20_noellht_comp[0] - alic_mga20[0]) == 0)
        assert (abs(alic_mga20_noellht_comp[1] - alic_mga20[1]) < 5e-5)
        assert (abs(alic_mga20_noellht_comp[2] - alic_mga20[2]) < 5e-5)
        assert (alic_mga20_noellht_comp[3] == 0)
        # Test with ellipsoid height supplied
        alic_mga20_ellht_comp = mga94to2020(alic_mga94[0], alic_mga94[1], alic_mga94[2], alic_mga94[3])
        assert ((alic_mga20_ellht_comp[0] - alic_mga20[0]) == 0)
        assert (abs(alic_mga20_ellht_comp[1] - alic_mga20[1]) < 5e-5)
        assert (abs(alic_mga20_ellht_comp[2] - alic_mga20[2]) < 5e-5)
        assert (abs(alic_mga20_ellht_comp[3] - alic_mga20[3]) < 5e-5)

    def test_mga2020to94(self):
        alic_mga94 = (53, 386352.3979, 7381850.7689, 603.3466)
        alic_mga20 = (53, 386353.2343, 7381852.2986, 603.2489)
        # Test with no ellipsoid height supplied
        alic_mga94_noellht_comp = mga2020to94(alic_mga20[0], alic_mga20[1], alic_mga20[2])
        assert ((alic_mga94_noellht_comp[0] - alic_mga94[0]) == 0)
        assert (abs(alic_mga94_noellht_comp[1] - alic_mga94[1]) < 5e-5)
        assert (abs(alic_mga94_noellht_comp[2] - alic_mga94[2]) < 5e-5)
        assert (alic_mga94_noellht_comp[3] == 0)
        # Test with ellipsoid height supplied
        alic_mga94_ellht_comp = mga2020to94(alic_mga20[0], alic_mga20[1], alic_mga20[2], alic_mga20[3])
        assert ((alic_mga94_ellht_comp[0] - alic_mga94[0]) == 0)
        assert (abs(alic_mga94_ellht_comp[1] - alic_mga94[1]) < 5e-5)
        assert (abs(alic_mga94_ellht_comp[2] - alic_mga94[2]) < 5e-5)
        assert (abs(alic_mga94_ellht_comp[3] - alic_mga94[3]) < 5e-5)


if __name__ == '__main__':
    unittest.main()

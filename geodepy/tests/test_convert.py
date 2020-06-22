import os
import unittest
import datetime
import numpy as np
from geodepy.fileio import read_dnacoord
from geodepy.convert import (dec2hp, hp2dec, DMSAngle, DDMAngle, dec2dms,
                             dec2ddm, hp2dms, hp2ddm, dd2sec,
                             yyyydoy_to_date, date_to_yyyydoy, grid2geo,
                             hp2dec_v, geo2grid, llh2xyz, xyz2llh)

dec_ex = 123.74875
dec_ex2 = 12.575
dec_ex3 = -12.575
dec_ex4 = 0.0525
dec_ex5 = 0.005

hp_ex = 123.44555
hp_ex2 = 12.3430
hp_ex3 = -12.3430
hp_ex4 = 0.0309
hp_ex5 = 0.0018

dms_ex = DMSAngle(123, 44, 55.5)
dms_ex2 = DMSAngle(12, 34, 30)
dms_ex3 = DMSAngle(-12, -34, -30)
dms_ex4 = DMSAngle(0, 3, 9)
dms_ex5 = DMSAngle(0, 0, 18)

ddm_ex = DDMAngle(123, 44.925)
ddm_ex2 = DDMAngle(12, 34.5)
ddm_ex3 = DDMAngle(-12, -34.5)
ddm_ex4 = DDMAngle(0, 3.15)
ddm_ex5 = DDMAngle(0, 0.3)


class TestConvert(unittest.TestCase):
    def test_dec2hp(self):
        self.assertAlmostEqual(hp_ex, dec2hp(dec_ex), 13)
        self.assertAlmostEqual(-hp_ex, dec2hp(-dec_ex), 13)

    def test_hp2dec(self):
        self.assertAlmostEqual(dec_ex, hp2dec(hp_ex), 13)
        self.assertAlmostEqual(-dec_ex, hp2dec(-hp_ex), 13)
        self.assertAlmostEqual(hp2dec(hp_ex) +
                               hp2dec(hp_ex2), dec_ex + dec_ex2, 13)

    def test_DMSAngle(self):
        # Test DMSAngle Methods
        self.assertEqual(dec_ex, dms_ex.dec())
        self.assertEqual(hp_ex, dms_ex.hp())
        self.assertEqual(hp_ex3, dms_ex3.hp())
        self.assertEqual(ddm_ex, dms_ex.ddm())
        self.assertEqual(-ddm_ex, -dms_ex.ddm())
        self.assertEqual(ddm_ex3, dms_ex3.ddm())

        # Test DMSAngle Sign Conventions
        self.assertEqual(-dec_ex, DMSAngle(-dms_ex.degree, dms_ex.minute,
                                           dms_ex.second).dec())
        self.assertEqual(dec_ex, DMSAngle(dms_ex.degree, -dms_ex.minute,
                                          -dms_ex.second).dec())
        self.assertAlmostEqual(-dec_ex4, DMSAngle(0, -dms_ex4.minute,
                                                  dms_ex4.second).dec(), 9)
        self.assertAlmostEqual(dec_ex4, DMSAngle(0, dms_ex4.minute,
                                                 dms_ex4.second).dec(), 9)
        self.assertEqual(-dec_ex5, DMSAngle(0, 0, -dms_ex5.second).dec())
        self.assertEqual(dec_ex5, DMSAngle(0, 0, dms_ex5.second).dec())
        self.assertEqual(-dms_ex3, DMSAngle(12, 34, -30))
        self.assertEqual(dms_ex.sign, 1)
        self.assertEqual(-dms_ex.sign, -1)
        self.assertEqual(dms_ex4.sign, 1)
        self.assertEqual(-dms_ex4.sign, -1)
        self.assertEqual(dms_ex5.sign, 1)
        self.assertEqual(-dms_ex5.sign, -1)
        self.assertEqual(DMSAngle(-1, 2, 3).sign, -1)
        self.assertEqual(DMSAngle(1, -2, 3).sign, 1)
        self.assertEqual(DMSAngle(1, 2, -3).sign, 1)
        self.assertEqual(DMSAngle(0, -1, 2).sign, -1)
        self.assertEqual(DMSAngle(0, 0, -3).sign, -1)
        self.assertEqual(DMSAngle(-0, 1, 2).sign, 1)
        self.assertEqual(DMSAngle(-0.0, 1, 2).sign, -1)
        self.assertEqual(repr(dms_ex), '{DMSAngle: +123d 44m 55.5s}')
        self.assertEqual(repr(dms_ex3), '{DMSAngle: -12d 34m 30s}')

        # Test DMSAngle Overloads
        self.assertEqual(dec_ex + dec_ex2, (dms_ex + dms_ex2).dec())
        self.assertEqual(dec_ex2 + dec_ex, (dms_ex2 + dms_ex).dec())
        self.assertEqual(dec_ex - dec_ex2, (dms_ex - dms_ex2).dec())
        self.assertEqual(dec_ex2 - dec_ex, (dms_ex2 - dms_ex).dec())
        self.assertEqual(dec_ex * 5, (dms_ex * 5).dec())
        self.assertEqual(5 * dec_ex, (5 * dms_ex).dec())
        self.assertEqual(dec_ex / 3, (dms_ex / 3).dec())
        self.assertEqual(abs(-dms_ex), dms_ex)
        self.assertEqual(-dms_ex2, dms_ex3)
        self.assertEqual(dms_ex2, abs(dms_ex3))
        self.assertEqual(dms_ex, ddm_ex)
        self.assertTrue(dms_ex == dms_ex)
        self.assertFalse(dms_ex == dms_ex2)
        self.assertTrue(dms_ex != dms_ex2)
        self.assertFalse(dms_ex != dms_ex)
        self.assertTrue(dms_ex > dms_ex2)
        self.assertFalse(dms_ex2 > dms_ex)
        self.assertTrue(dms_ex2 < dms_ex)
        self.assertFalse(dms_ex < dms_ex2)
        with self.assertRaises(TypeError):
            dms_ex * 'a'
        with self.assertRaises(TypeError):
            'a' * dms_ex
        with self.assertRaises(TypeError):
            dms_ex / 'a'
        with self.assertRaises(TypeError):
            dms_ex + 'a'
        with self.assertRaises(TypeError):
            'a' + dms_ex
        with self.assertRaises(TypeError):
            dms_ex - 'a'
        with self.assertRaises(TypeError):
            'a' - dms_ex

        # Test Class Interoperability
        self.assertEqual(DMSAngle(1, 2, 3) + DDMAngle(2, 3), DMSAngle(3, 5, 3))
        self.assertEqual(DMSAngle(3, 2, 0) - DDMAngle(2, 2.5),
                         DMSAngle(0, 59, 30))
        self.assertEqual(DDMAngle(2, 3) + DMSAngle(1, 2, 3), DDMAngle(3, 5.05))
        self.assertEqual(DDMAngle(3, 2) - DMSAngle(2, 2, 30),
                         DDMAngle(0, 59.5))

    def test_DDMAngle(self):
        # Test DDMAngle Methods
        self.assertEqual(dec_ex, ddm_ex.dec())
        self.assertEqual(hp_ex, ddm_ex.hp())
        self.assertEqual(dms_ex, ddm_ex.dms())
        self.assertEqual(hp_ex3, ddm_ex3.hp())

        # Test DMSAngle Sign Conventions
        self.assertEqual(-dec_ex, DDMAngle(-dms_ex.degree,
                                           ddm_ex.minute).dec())
        self.assertEqual(dec_ex, DDMAngle(dms_ex.degree, -ddm_ex.minute).dec())
        self.assertAlmostEqual(-dec_ex4, DDMAngle(0, -ddm_ex4.minute).dec(), 9)
        self.assertAlmostEqual(dec_ex4, DDMAngle(0, ddm_ex4.minute).dec(), 9)
        self.assertEqual(-ddm_ex3, DDMAngle(12, 34.5))
        self.assertEqual(ddm_ex.sign, 1)
        self.assertEqual(-ddm_ex.sign, -1)
        self.assertEqual(ddm_ex4.sign, 1)
        self.assertEqual(-ddm_ex4.sign, -1)
        self.assertEqual(ddm_ex5.sign, 1)
        self.assertEqual(-ddm_ex5.sign, -1)
        self.assertEqual(DDMAngle(-1, 2).sign, -1)
        self.assertEqual(DDMAngle(1, -2).sign, 1)
        self.assertEqual(DDMAngle(1, 2).sign, 1)
        self.assertEqual(DDMAngle(0, -1).sign, -1)
        self.assertEqual(DDMAngle(-0, 1).sign, 1)
        self.assertEqual(DDMAngle(-0.0, 1).sign, -1)
        self.assertEqual(repr(ddm_ex), '{DDMAngle: +123d 44.925m}')
        self.assertEqual(repr(ddm_ex3), '{DDMAngle: -12d 34.5m}')

        # Test DDMAngle Overloads
        self.assertEqual(dec_ex + dec_ex2, (ddm_ex + ddm_ex2).dec())
        self.assertEqual(dec_ex2 + dec_ex, (ddm_ex2 + ddm_ex).dec())
        self.assertEqual(dec_ex - dec_ex2, (ddm_ex - ddm_ex2).dec())
        self.assertEqual(dec_ex2 - dec_ex, (ddm_ex2 - ddm_ex).dec())
        self.assertEqual(dec_ex * 5, (ddm_ex * 5).dec())
        self.assertEqual(5 * dec_ex, (5 * ddm_ex).dec())
        self.assertEqual(dec_ex / 3, (ddm_ex / 3).dec())
        self.assertEqual(abs(-ddm_ex), ddm_ex)
        self.assertEqual(-ddm_ex2, ddm_ex3)
        self.assertEqual(ddm_ex2, abs(ddm_ex3))
        self.assertEqual(ddm_ex, dms_ex)
        self.assertTrue(ddm_ex == ddm_ex)
        self.assertFalse(ddm_ex == ddm_ex2)
        self.assertTrue(ddm_ex != ddm_ex2)
        self.assertFalse(ddm_ex != ddm_ex)
        self.assertTrue(ddm_ex > ddm_ex2)
        self.assertFalse(ddm_ex2 > ddm_ex)
        self.assertTrue(ddm_ex2 < ddm_ex)
        self.assertFalse(ddm_ex < ddm_ex2)
        with self.assertRaises(TypeError):
            ddm_ex * 'a'
        with self.assertRaises(TypeError):
            'a' * ddm_ex
        with self.assertRaises(TypeError):
            ddm_ex / 'a'
        with self.assertRaises(TypeError):
            ddm_ex + 'a'
        with self.assertRaises(TypeError):
            'a' + ddm_ex
        with self.assertRaises(TypeError):
            ddm_ex - 'a'
        with self.assertRaises(TypeError):
            'a' - ddm_ex

    def test_dec2dms(self):
        self.assertEqual(dms_ex, dec2dms(dec_ex))
        self.assertEqual(-dms_ex, dec2dms(-dec_ex))

    def test_dec2ddm(self):
        self.assertEqual(ddm_ex, dec2ddm(dec_ex))
        self.assertEqual(-ddm_ex, dec2ddm(-dec_ex))

    def test_hp2dms(self):
        self.assertEqual(dms_ex.degree, hp2dms(hp_ex).degree)
        self.assertEqual(dms_ex.minute, hp2dms(hp_ex).minute)
        self.assertAlmostEqual(dms_ex.second, hp2dms(hp_ex).second, 10)

        self.assertEqual(-dms_ex.sign, hp2dms(-hp_ex).sign)
        self.assertEqual(dms_ex.degree, hp2dms(-hp_ex).degree)
        self.assertEqual(dms_ex.minute, hp2dms(-hp_ex).minute)
        self.assertAlmostEqual(dms_ex.second, hp2dms(-hp_ex).second, 10)

    def test_hp2ddm(self):
        self.assertEqual(ddm_ex, hp2ddm(hp_ex))
        self.assertEqual(-ddm_ex, hp2ddm(-hp_ex))

    def test_dd2sec(self):
        self.assertEqual(dd2sec(1), 3600)
        self.assertEqual(dd2sec(-1), -3600)
        self.assertEqual(dd2sec(hp2dec(0.0001)), 1)
        self.assertEqual(dd2sec(hp2dec(-0.0001)), -1)
        self.assertEqual(dd2sec(hp2dec(0.00001)), 0.1)
        self.assertEqual(dd2sec(dec_ex4), 189)
        self.assertEqual(dd2sec(-dec_ex4), -189)
        self.assertEqual(dd2sec(dec_ex2), 45270)
        self.assertEqual(dd2sec(-dec_ex2), -45270)

    def test_date_to_yyyydoy(self):
        self.assertEqual(date_to_yyyydoy(datetime.date(2020, 1, 4)),
                         '2020.004')
        self.assertEqual(date_to_yyyydoy(datetime.date(2020, 10, 12)),
                         '2020.286')
        self.assertEqual(date_to_yyyydoy(datetime.date(1998, 4, 7)),
                         '1998.097')
        self.assertEqual(date_to_yyyydoy(datetime.date(2000, 11, 22)),
                         '2000.327')
        self.assertEqual(date_to_yyyydoy(datetime.date(2008, 2, 29)),
                         '2008.060')
        with self.assertRaises(AttributeError):
            date_to_yyyydoy('a')
        with self.assertRaises(AttributeError):
            date_to_yyyydoy('2020123')

    def test_yyyydoy_to_date(self):
        self.assertEqual(yyyydoy_to_date('2020.004'),
                         datetime.date(2020, 1, 4))
        self.assertEqual(yyyydoy_to_date('2020.286'),
                         datetime.date(2020, 10, 12))
        self.assertEqual(yyyydoy_to_date('1998.097'),
                         datetime.date(1998, 4, 7))
        self.assertEqual(yyyydoy_to_date('2000.327'),
                         datetime.date(2000, 11, 22))
        self.assertEqual(yyyydoy_to_date('2008.060'),
                         datetime.date(2008, 2, 29))
        self.assertEqual(yyyydoy_to_date('2020004'),
                         datetime.date(2020, 1, 4))
        self.assertEqual(yyyydoy_to_date('2020286'),
                         datetime.date(2020, 10, 12))
        self.assertEqual(yyyydoy_to_date('1998097'),
                         datetime.date(1998, 4, 7))
        self.assertEqual(yyyydoy_to_date('2000327'),
                         datetime.date(2000, 11, 22))
        self.assertEqual(yyyydoy_to_date('2008060'),
                         datetime.date(2008, 2, 29))
        with self.assertRaises(ValueError):
            yyyydoy_to_date('a')
        with self.assertRaises(ValueError):
            yyyydoy_to_date('20201234')
        with self.assertRaises(ValueError):
            yyyydoy_to_date('2020.1234')
        with self.assertRaises(ValueError):
            yyyydoy_to_date('202012')
        with self.assertRaises(ValueError):
            yyyydoy_to_date('2020.12')

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

    def test_llh2xyz(self):

        # Test of single point
        x, y, z = llh2xyz(hp2dec(-37.482667598), hp2dec(144.581644114), 39.6514)
        self.assertAlmostEqual(x, -4131654.2815, 3)
        self.assertAlmostEqual(y, 2896107.9738, 3)
        self.assertAlmostEqual(z, -3888601.3067, 3)

        # Test DMSAngle input
        x, y, z = llh2xyz(DMSAngle(-37, 48, 26.67598),
                          DMSAngle(144, 58, 16.44114),
                          39.6514)
        self.assertAlmostEqual(x, -4131654.2815, 3)
        self.assertAlmostEqual(y, 2896107.9738, 3)
        self.assertAlmostEqual(z, -3888601.3067, 3)

        # Test DDMAngle input
        x, y, z = llh2xyz(DDMAngle(-37, 48.4445996),
                          DDMAngle(144, 58.274019),
                          39.6514)
        self.assertAlmostEqual(x, -4131654.2815, 3)
        self.assertAlmostEqual(y, 2896107.9738, 3)
        self.assertAlmostEqual(z, -3888601.3067, 3)

        # Tests comparing results from DynAdjust Adjusted Coordinate File
        # including 109 Points Across Australia
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
        # Single Point Test
        hem, zone, east, north, psf, grid_conv = geo2grid(hp2dec(-37.482667598),
                                                          hp2dec(144.581644114))
        self.assertEqual(hem, 'South')
        self.assertEqual(zone, 55)
        self.assertAlmostEqual(east, 321405.5592, 3)
        self.assertAlmostEqual(north, 5813614.1613, 3)
        self.assertAlmostEqual(psf, 0.99999287, 8)
        self.assertAlmostEqual(grid_conv, -1.2439811331, 9)

        # Test DMSAngle Input
        (hem, zone, east,
         north, psf, grid_conv) = geo2grid(DMSAngle(-37, 48, 26.67598),
                                           DMSAngle(144, 58, 16.44114))
        self.assertEqual(hem, 'South')
        self.assertEqual(zone, 55)
        self.assertAlmostEqual(east, 321405.5592, 3)
        self.assertAlmostEqual(north, 5813614.1613, 3)
        self.assertAlmostEqual(psf, 0.99999287, 8)
        self.assertAlmostEqual(grid_conv, -1.2439811331, 9)

        # Test DDMAngle Input
        (hem, zone, east,
         north, psf, grid_conv) = geo2grid(DDMAngle(-37, 48.4445997),
                                           DDMAngle(144, 58.274019))
        self.assertEqual(hem, 'South')
        self.assertEqual(zone, 55)
        self.assertAlmostEqual(east, 321405.5592, 3)
        self.assertAlmostEqual(north, 5813614.1613, 3)
        self.assertAlmostEqual(psf, 0.99999287, 8)
        self.assertAlmostEqual(grid_conv, -1.2439811331, 9)

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


if __name__ == '__main__':
    unittest.main()

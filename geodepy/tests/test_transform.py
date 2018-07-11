import unittest

from geodepy.transform import geo2grid, grid2geo, llh2xyz, xyz2llh
from geodepy.convert import dms2dd_v, read_dnacoord
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
        np.testing.assert_almost_equal(geoed_grid[:, :2], dms2dd_v(np.array(test_geo_coords[['lat', 'lon']].tolist())),
                                       decimal=3)

        gridded_geo = np.stack(geo2grid(*x) for x in dms2dd_v(np.array(test_geo_coords[['lat', 'lon']].tolist())))
        np.testing.assert_almost_equal(gridded_geo[:, 2:4].astype(float), np.array(test_grid_coords[['east', 'north']].tolist()),
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

        testdata = read_dnacoord(os.path.join(abs_path, 'resources/natadjust_rvs_example.dat'))
        for coord in testdata:
            coord.converthptodd()
            hem, zonecomp, eastcomp, northcomp, psf, grid_conv = geo2grid(coord.lat, coord.long)
            assert (zonecomp == coord.zone)
            assert (abs(eastcomp - coord.easting) < 4e-4)
            assert (abs(northcomp - coord.northing) < 4e-4)

    def test_grid2geo(self):
        abs_path = os.path.abspath(os.path.dirname(__file__))

        testdata = read_dnacoord(os.path.join(abs_path, 'resources/natadjust_rvs_example.dat'))
        for coord in testdata:
            coord.converthptodd()
            latcomp, longcomp, psf, grid_conv = grid2geo(coord.zone, coord.easting, coord.northing)
            assert (abs(latcomp - coord.lat) < 5e-9)
            assert (abs(longcomp - coord.long) < 5e-9)


if __name__ == '__main__':
    unittest.main()

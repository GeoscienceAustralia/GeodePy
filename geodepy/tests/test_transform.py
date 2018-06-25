import unittest

from geodepy.transform import geo2grid, grid2geo
from geodepy.convert import dms2dd_v
import numpy as np
import os.path


class TestTransforms(unittest.TestCase):
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
        np.testing.assert_almost_equal(geoed_grid[:, :2], dms2dd_v(np.array(test_geo_coords[['lat', 'lon']].tolist())))

        gridded_geo = np.stack(geo2grid(*x) for x in dms2dd_v(np.array(test_geo_coords[['lat', 'lon']].tolist())))
        np.testing.assert_almost_equal(gridded_geo[:, 2:4].astype(float), np.array(test_grid_coords[['east', 'north']].tolist()))


if __name__ == '__main__':
    unittest.main()

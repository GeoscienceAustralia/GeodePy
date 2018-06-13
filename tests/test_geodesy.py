import unittest
import os.path
import numpy as np
import numpy.lib.recfunctions as rfn
from conversions import dms2dd, dd2dms
from geodesy import vincinv, vincdir

class TestGeodesy(unittest.TestCase):
    def test_vincinv(self):
        # Flinders Peak
        lat1 = dms2dd(-37.57037203)
        long1 = dms2dd(144.25295244)

        # Buninyong
        lat2 = dms2dd(-37.39101561)
        long2 = dms2dd(143.55353839)

        ell_dist, azimuth1to2, azimuth2to1 = vincinv(lat1, long1, lat2, long2)
        self.assertEqual(round(ell_dist, 3), 54972.271)
        self.assertEqual(round(dd2dms(azimuth1to2), 6), 306.520537)
        self.assertEqual(round(dd2dms(azimuth2to1), 6), 127.102507)

    def test_vincdir(self):
        # Flinders Peak
        lat1 = dms2dd(-37.57037203)
        long1 = dms2dd(144.25295244)

        # To Buninyong
        azimuth1to2 = dms2dd(306.520537)
        ell_dist = 54972.271

        lat2, long2, azimuth2to1 = vincdir(lat1, long1, azimuth1to2, ell_dist)
        self.assertEqual(round(dd2dms(lat2), 8), -37.39101562)
        self.assertEqual(round(dd2dms(long2), 8), 143.55353840)
        self.assertEqual(round(dd2dms(azimuth2to1), 6), 127.102507)

    # def test_vincentys(self):
        # Test multiple point-to-point vincinv calculations
        abs_path = os.path.abspath(os.path.dirname(__file__))

        test_geo_coords = np.genfromtxt(os.path.join(abs_path, 'resources/Test_Conversion_Geo.csv'),
                                        delimiter=',',
                                        dtype='S4,f8,f8',
                                        names=['site', 'lat1', 'long1'],
                                        usecols=('lat1', 'long1'))

        test_geo_coord2 = np.genfromtxt(os.path.join(abs_path, 'resources/Test_Conversion_Geo.csv'),
                                        delimiter=',',
                                        dtype='S4,f8,f8',
                                        names=['site', 'lat2', 'long2'],
                                        usecols=('lat2', 'long2'))

        # Form array with point pairs from test file
        np.roll(test_geo_coord2, 1)
        test_pairs = rfn.merge_arrays([test_geo_coords, test_geo_coord2], flatten=True)

        # Calculate Vincenty's Inverse Result using Lat, Long Pairs
        vincinv_result = np.array(list(vincinv(*x) for x in test_pairs[['lat1', 'long1', 'lat2', 'long2']]))

        # Then do a similar thing with vincdir, using the az1to2 and dist generated from vincinv to gen new latlong2
        # np.testing.assert_almost_equal(original lat2, long2 comp w vincdir lat2, long2)
        # np.testing.assert_almost_equal(vincinv az2to1 comp w vincdir az2to1)




if __name__ == '__main__':
    unittest.main()

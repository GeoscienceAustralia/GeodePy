import unittest
import os.path
import numpy as np
import numpy.lib.recfunctions as rfn
from geodepy.convert import hp2dec, dec2hp
from geodepy.geodesy import vincinv, vincdir


class TestGeodesy(unittest.TestCase):
    def test_vincinv(self):
        # Flinders Peak
        lat1 = hp2dec(-37.57037203)
        long1 = hp2dec(144.25295244)
        # Buninyong
        lat2 = hp2dec(-37.39101561)
        long2 = hp2dec(143.55353839)
        ell_dist, azimuth1to2, azimuth2to1 = vincinv(lat1, long1, lat2, long2)
        self.assertEqual(round(ell_dist, 3), 54972.271)
        self.assertEqual(round(dec2hp(azimuth1to2), 6), 306.520537)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

    def test_vincdir(self):
        # Flinders Peak
        lat1 = hp2dec(-37.57037203)
        long1 = hp2dec(144.25295244)
        # To Buninyong
        azimuth1to2 = hp2dec(306.520537)
        ell_dist = 54972.271
        lat2, long2, azimuth2to1 = vincdir(lat1, long1, azimuth1to2, ell_dist)
        self.assertEqual(round(dec2hp(lat2), 8), -37.39101561)
        self.assertEqual(round(dec2hp(long2), 8), 143.55353839)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

    def test_equality_vincentys(self):
        # Test multiple point-to-point vincinv calculations
        abs_path = os.path.abspath(os.path.dirname(__file__))

        test_geo_coords =\
            np.genfromtxt(os.path.join(abs_path,
                                       'resources/Test_Conversion_Geo.csv'),
                                        delimiter=',',
                                        dtype='S4,f8,f8',
                                        names=['site', 'lat1', 'long1'],
                                        usecols=('lat1', 'long1'))

        test_geo_coord2 = \
            np.genfromtxt(os.path.join(abs_path,
                                       'resources/Test_Conversion_Geo.csv'),
                                        delimiter=',',
                                        dtype='S4,f8,f8',
                                        names=['site', 'lat2', 'long2'],
                                        usecols=('lat2', 'long2'))

        # Form array with point pairs from test file
        test_pairs = rfn.merge_arrays([test_geo_coords, np.roll(test_geo_coord2, 1)], flatten=True)

        # Calculate Vincenty's Inverse Result using Lat, Long Pairs
        vincinv_result = np.array(list(vincinv(*x) for x in test_pairs[['lat1', 'long1', 'lat2', 'long2']]))

        # Calculate Vincenty's Direct Result using Results from Inverse Function
        vincdir_input = rfn.merge_arrays([test_geo_coords, vincinv_result[:, 1], vincinv_result[:, 0]], flatten=True)
        vincdir_input.dtype.names = ['lat1', 'long1', 'az1to2', 'ell_dist']
        vincdir_result = np.array(list(vincdir(*x) for x in vincdir_input[['lat1', 'long1', 'az1to2', 'ell_dist']]))

        np.testing.assert_almost_equal(test_pairs['lat2'],
                                       vincdir_result[:, 0], decimal=8)
        np.testing.assert_almost_equal(test_pairs['long2'],
                                       vincdir_result[:, 1], decimal=8)
        np.testing.assert_almost_equal(vincinv_result[:, 2],
                                       vincdir_result[:, 2])

    def test_vincinv_edgecases(self):
        lat1 = -32.153892
        lon1 = -15.394827
        lat2 = -31.587369
        lon2 = -13.487739
        gdist, az12, az21 = vincinv(lat1, lon1, lat2, lon2)
        lon1 = lon1 + 14
        lon2 = lon2 + 14
        gdist_2, az12_2, az21_2 = vincinv(lat1, lon1, lat2, lon2)
        self.assertEqual(gdist, gdist_2)
        self.assertEqual(az12, az12_2)
        self.assertEqual(az21, az21_2)


if __name__ == '__main__':
    unittest.main()

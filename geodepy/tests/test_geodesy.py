import unittest
from geodepy.convert import dms2dd, dd2dms
from geodepy.geodesy import vincinv

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


if __name__ == '__main__':
    unittest.main()

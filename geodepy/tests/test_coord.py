import unittest
from geodepy.angles import DMSAngle
from geodepy.constants import utm
from geodepy.coord import CoordCart, CoordGeo, CoordTM

cart_ex1 = CoordCart(-4052052.7379, 4212835.9897, -2545104.5898, -14.269)

geo_ex1 = CoordGeo(DMSAngle(-23, 40, 12.39650), DMSAngle(133, 53, 7.87779),
                   603.2489, 588.9799)

tm_ex1 = CoordTM(53, 386353.2343, 7381852.2986, 603.2489, 588.9799,
                 hemi_north=False, projection=utm)


class TestCoord(unittest.TestCase):
    def test_CoordCart(self):
        self.assertEqual(repr(cart_ex1), 'CoordCart: X: -4052052.7379 '
                                         'Y: 4212835.9897 Z: -2545104.5898 '
                                         'NVal: -14.269')
        # Test converting from CoordCart to CoordGeo
        cart2geo = cart_ex1.geo()
        self.assertAlmostEqual(cart2geo.lat, geo_ex1.lat, 8)
        self.assertAlmostEqual(cart2geo.lon, geo_ex1.lon, 8)
        self.assertAlmostEqual(cart2geo.ell_ht, geo_ex1.ell_ht, 3)
        self.assertAlmostEqual(cart2geo.orth_ht, geo_ex1.orth_ht, 3)

    def test_CoordGeo(self):
        self.assertEqual(repr(geo_ex1), 'CoordGeo: Lat: -23.67011013888889 '
                                        'Lon: 133.88552160833333 '
                                        'Ell_Ht: 603.2489 Orth_Ht: 588.9799')
        # Test converting from CoordGeo to CoordCart
        geo2cart = geo_ex1.cart()
        self.assertAlmostEqual(geo2cart.xaxis, cart_ex1.xaxis, 3)
        self.assertAlmostEqual(geo2cart.yaxis, cart_ex1.yaxis, 3)
        self.assertAlmostEqual(geo2cart.zaxis, cart_ex1.zaxis, 3)
        self.assertAlmostEqual(geo2cart.nval, cart_ex1.nval, 3)

    def test_CoordTM(self):
        self.assertEqual(repr(tm_ex1), 'CoordTM: Zone: 53 East: 386353.2343 '
                                       'North: 7381852.2986 Ell_Ht: 603.2489 '
                                       'Orth_Ht: 588.9799 Hemisphere: South')

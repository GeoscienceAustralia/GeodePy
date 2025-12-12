import unittest
from geodepy.angles import DMSAngle, DECAngle
from geodepy.constants import utm, grs80
from geodepy.coord import CoordCart, CoordGeo, CoordTM

cart_ex1 = CoordCart(-4052052.7379, 4212835.9897, -2545104.5898, 14.269)
cart_ex2 = CoordCart(-4052052.7379, 4212835.9897, -2545104.5898)

geo_ex1 = CoordGeo(
    DMSAngle(-23, 40, 12.39650).deca(),
    DMSAngle(133, 53, 7.87779).deca(),
    603.2489,
    588.9799,
)
geo_ex2 = CoordGeo(geo_ex1.lat, geo_ex1.lon, None, 588.9799)
geo_ex3 = CoordGeo(geo_ex1.lat, geo_ex1.lon, 603.2489, None)
geo_ex4 = CoordGeo(geo_ex1.lat, geo_ex1.lon)

tm_ex1 = CoordTM(
    53, 386353.2343, 7381852.2986, 603.2489, 588.9799, hemi_north=False, projection=utm
)
tm_ex2 = CoordTM(
    tm_ex1.zone,
    tm_ex1.east,
    tm_ex1.north,
    None,
    tm_ex1.orth_ht,
    hemi_north=False,
    projection=utm,
)
tm_ex3 = CoordTM(
    tm_ex1.zone,
    tm_ex1.east,
    tm_ex1.north,
    tm_ex1.ell_ht,
    None,
    hemi_north=False,
    projection=utm,
)
tm_ex4 = CoordTM(
    tm_ex1.zone, tm_ex1.east, tm_ex1.north, None, None, hemi_north=False, projection=utm
)


class TestCoord(unittest.TestCase):
    def test_CoordCart(self):
        # Test Overloads
        self.assertEqual(
            repr(cart_ex1),
            "CoordCart: X: -4052052.7379 "
            "Y: 4212835.9897 Z: -2545104.5898 "
            "NVal: 14.269",
        )

        self.assertEqual(cart_ex1, cart_ex1)
        with self.assertRaises(ValueError):
            cart_ex1 == geo_ex1

        cart_ex1_rounded = round(cart_ex1, 1)
        self.assertEqual(cart_ex1_rounded.xaxis, round(cart_ex1.xaxis, 1))
        self.assertEqual(cart_ex1_rounded.yaxis, round(cart_ex1.yaxis, 1))
        self.assertEqual(cart_ex1_rounded.zaxis, round(cart_ex1.zaxis, 1))
        self.assertEqual(cart_ex1_rounded.nval, round(cart_ex1.nval, 1))
        # No NVal
        cart_ex2_rounded = round(cart_ex2, 1)
        self.assertEqual(cart_ex2_rounded.xaxis, round(cart_ex2.xaxis, 1))
        self.assertEqual(cart_ex2_rounded.yaxis, round(cart_ex2.yaxis, 1))
        self.assertEqual(cart_ex2_rounded.zaxis, round(cart_ex2.zaxis, 1))
        self.assertEqual(cart_ex2_rounded.nval, None)

        # Test converting from CoordCart to CoordGeo
        cart2geo = cart_ex1.geo(ellipsoid=grs80, notation=DECAngle)
        self.assertEqual(round(cart2geo.lat, 8), round(geo_ex1.lat, 8))
        self.assertEqual(round(cart2geo.lon, 8), round(geo_ex1.lon, 8))
        self.assertAlmostEqual(cart2geo.ell_ht, geo_ex1.ell_ht, 3)
        self.assertAlmostEqual(cart2geo.orth_ht, geo_ex1.orth_ht, 3)

        # Test converting from CoordCart to CoordTM
        cart2tm = cart_ex1.tm(ellipsoid=grs80, projection=utm)
        self.assertEqual(cart2tm.zone, tm_ex1.zone)
        self.assertAlmostEqual(cart2tm.east, tm_ex1.east, 3)
        self.assertAlmostEqual(cart2tm.north, tm_ex1.north, 3)
        self.assertAlmostEqual(cart2tm.ell_ht, tm_ex1.ell_ht, 3)
        self.assertAlmostEqual(cart2tm.orth_ht, tm_ex1.orth_ht, 3)
        self.assertEqual(cart2tm.hemi_north, tm_ex1.hemi_north)

    def test_CoordGeo(self):
        # Test Overloads
        self.assertEqual(
            repr(geo_ex1),
            "CoordGeo: Lat: -23.67011013888889 "
            "Lon: 133.88552160833333 "
            "Ell_Ht: 603.2489 Orth_Ht: 588.9799",
        )

        self.assertEqual(geo_ex1, geo_ex1)
        with self.assertRaises(ValueError):
            geo_ex1 == cart_ex1

        geo_ex1_rounded = round(geo_ex1, 1)
        self.assertEqual(geo_ex1_rounded.lat, round(geo_ex1.lat, 1))
        self.assertEqual(geo_ex1_rounded.lon, round(geo_ex1.lon, 1))
        self.assertEqual(geo_ex1_rounded.ell_ht, round(geo_ex1.ell_ht, 1))
        self.assertEqual(geo_ex1_rounded.orth_ht, round(geo_ex1.orth_ht, 1))
        # No Ell Ht
        geo_ex2_rounded = round(geo_ex2, 1)
        self.assertEqual(geo_ex2_rounded.lat, round(geo_ex2.lat, 1))
        self.assertEqual(geo_ex2_rounded.lon, round(geo_ex2.lon, 1))
        self.assertEqual(geo_ex2_rounded.ell_ht, None)
        self.assertEqual(geo_ex2_rounded.orth_ht, round(geo_ex2.orth_ht, 1))
        # No Orth Ht
        geo_ex3_rounded = round(geo_ex3, 1)
        self.assertEqual(geo_ex3_rounded.lat, round(geo_ex3.lat, 1))
        self.assertEqual(geo_ex3_rounded.lon, round(geo_ex3.lon, 1))
        self.assertEqual(geo_ex3_rounded.ell_ht, round(geo_ex3.ell_ht, 1))
        self.assertEqual(geo_ex3_rounded.orth_ht, None)
        # No Ht
        geo_ex4_rounded = round(geo_ex4, 1)
        self.assertEqual(geo_ex4_rounded.lat, round(geo_ex4.lat, 1))
        self.assertEqual(geo_ex4_rounded.lon, round(geo_ex4.lon, 1))
        self.assertEqual(geo_ex4_rounded.ell_ht, None)
        self.assertEqual(geo_ex4_rounded.orth_ht, None)

        # Test converting from CoordGeo to CoordCart
        geo2cart = geo_ex1.cart(ellipsoid=grs80)
        self.assertAlmostEqual(geo2cart.xaxis, cart_ex1.xaxis, 3)
        self.assertAlmostEqual(geo2cart.yaxis, cart_ex1.yaxis, 3)
        self.assertAlmostEqual(geo2cart.zaxis, cart_ex1.zaxis, 3)
        self.assertAlmostEqual(geo2cart.nval, cart_ex1.nval, 3)

        # Test converting from CoordGeo to CoordTM
        geo2tm = geo_ex1.tm(ellipsoid=grs80, projection=utm)
        self.assertEqual(geo2tm.zone, tm_ex1.zone)
        self.assertEqual(geo2tm.east, tm_ex1.east)
        self.assertEqual(geo2tm.north, tm_ex1.north)
        self.assertEqual(geo2tm.ell_ht, tm_ex1.ell_ht)
        self.assertEqual(geo2tm.orth_ht, tm_ex1.orth_ht)
        self.assertEqual(geo2tm.hemi_north, tm_ex1.hemi_north)
        self.assertEqual(geo2tm, tm_ex1)

    def test_CoordTM(self):
        # Test Overloads
        self.assertEqual(
            repr(tm_ex1),
            "CoordTM: Zone: 53 East: 386353.2343 "
            "North: 7381852.2986 Ell_Ht: 603.2489 "
            "Orth_Ht: 588.9799 Hemisphere: South",
        )

        self.assertEqual(tm_ex1, tm_ex1)
        with self.assertRaises(ValueError):
            tm_ex1 == cart_ex1

        tm_ex1_rounded = round(tm_ex1, 1)
        self.assertEqual(tm_ex1_rounded.zone, tm_ex1.zone)
        self.assertEqual(tm_ex1_rounded.east, round(tm_ex1.east, 1))
        self.assertEqual(tm_ex1_rounded.north, round(tm_ex1.north, 1))
        self.assertEqual(tm_ex1_rounded.ell_ht, round(tm_ex1.ell_ht, 1))
        self.assertEqual(tm_ex1_rounded.orth_ht, round(tm_ex1.orth_ht, 1))
        self.assertEqual(tm_ex1_rounded.hemi_north, tm_ex1.hemi_north)
        # No Ell Ht
        tm_ex2_rounded = round(tm_ex2, 1)
        self.assertEqual(tm_ex2_rounded.zone, tm_ex2.zone)
        self.assertEqual(tm_ex2_rounded.east, round(tm_ex2.east, 1))
        self.assertEqual(tm_ex2_rounded.north, round(tm_ex2.north, 1))
        self.assertEqual(tm_ex2_rounded.ell_ht, None)
        self.assertEqual(tm_ex2_rounded.orth_ht, round(tm_ex2.orth_ht, 1))
        self.assertEqual(tm_ex2_rounded.hemi_north, tm_ex2.hemi_north)
        # No Orth Ht
        tm_ex3_rounded = round(tm_ex3, 1)
        self.assertEqual(tm_ex3_rounded.zone, tm_ex3.zone)
        self.assertEqual(tm_ex3_rounded.east, round(tm_ex3.east, 1))
        self.assertEqual(tm_ex3_rounded.north, round(tm_ex3.north, 1))
        self.assertEqual(tm_ex3_rounded.ell_ht, round(tm_ex3.ell_ht, 1))
        self.assertEqual(tm_ex3_rounded.orth_ht, None)
        self.assertEqual(tm_ex3_rounded.hemi_north, tm_ex3.hemi_north)
        # No Ht
        tm_ex4_rounded = round(tm_ex4, 1)
        self.assertEqual(tm_ex4_rounded.zone, tm_ex4.zone)
        self.assertEqual(tm_ex4_rounded.east, round(tm_ex4.east, 1))
        self.assertEqual(tm_ex4_rounded.north, round(tm_ex4.north, 1))
        self.assertEqual(tm_ex4_rounded.ell_ht, None)
        self.assertEqual(tm_ex4_rounded.orth_ht, None)
        self.assertEqual(tm_ex4_rounded.hemi_north, tm_ex4.hemi_north)

        # Test converting from CoordTM to CoordGeo
        tm2geo = tm_ex1.geo(ellipsoid=grs80, notation=DECAngle)
        self.assertEqual(round(tm2geo.lat, 8), round(geo_ex1.lat, 8))
        self.assertEqual(round(tm2geo.lon, 8), round(geo_ex1.lon, 8))
        self.assertAlmostEqual(tm2geo.ell_ht, geo_ex1.ell_ht, 3)
        self.assertAlmostEqual(tm2geo.orth_ht, geo_ex1.orth_ht, 3)

        # Test converting from CoordTM to CoordCart
        tm2cart = tm_ex1.cart(ellipsoid=grs80)
        self.assertAlmostEqual(tm2cart.xaxis, cart_ex1.xaxis, 3)
        self.assertAlmostEqual(tm2cart.yaxis, cart_ex1.yaxis, 3)
        self.assertAlmostEqual(tm2cart.zaxis, cart_ex1.zaxis, 3)
        self.assertAlmostEqual(tm2cart.nval, cart_ex1.nval, 3)

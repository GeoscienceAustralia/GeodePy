import unittest
import os
from math import radians, pi

from geodepy.angles import (
    DECAngle,
    HPAngle,
    GONAngle,
    DMSAngle,
    DDMAngle,
    dec2hp,
    dec2hpa,
    dec2gon,
    dec2gona,
    dec2dms,
    dec2ddm,
    hp2dec,
    hp2deca,
    hp2gon,
    hp2gona,
    hp2dms,
    hp2ddm,
    hp2rad,
    gon2dec,
    gon2deca,
    gon2hp,
    gon2hpa,
    gon2dms,
    gon2ddm,
    gon2rad,
    dd2sec,
    angular_typecheck,
)

rad_exs = [
    radians(123.74875),
    radians(12.575),
    radians(-12.575),
    radians(0.0525),
    radians(0.005),
    radians(-0.005),
]

dec_ex = 123.74875
dec_ex2 = 12.575
dec_ex3 = -12.575
dec_ex4 = 0.0525
dec_ex5 = 0.005
dec_ex6 = -0.005
dec_exs = [dec_ex, dec_ex2, dec_ex3, dec_ex4, dec_ex5, dec_ex6]

deca_ex = DECAngle(123.74875)
deca_ex2 = DECAngle(12.575)
deca_ex3 = DECAngle(-12.575)
deca_ex4 = DECAngle(0.0525)
deca_ex5 = DECAngle(0.005)
deca_ex6 = DECAngle(-0.005)
deca_exs = [deca_ex, deca_ex2, deca_ex3, deca_ex4, deca_ex5, deca_ex6]

hp_ex = 123.44555
hp_ex2 = 12.3430
hp_ex3 = -12.3430
hp_ex4 = 0.0309
hp_ex5 = 0.0018
hp_ex6 = -0.0018
hp_exs = [hp_ex, hp_ex2, hp_ex3, hp_ex4, hp_ex5, hp_ex6]

hpa_ex = HPAngle(123.44555)
hpa_ex2 = HPAngle(12.3430)
hpa_ex3 = HPAngle(-12.3430)
hpa_ex4 = HPAngle(0.0309)
hpa_ex5 = HPAngle(0.0018)
hpa_ex6 = HPAngle(-0.0018)
hpa_exs = [hpa_ex, hpa_ex2, hpa_ex3, hpa_ex4, hpa_ex5, hpa_ex6]

dms_ex = DMSAngle(123, 44, 55.5)
dms_ex2 = DMSAngle(12, 34, 30)
dms_ex3 = DMSAngle(-12, -34, -30)
dms_ex4 = DMSAngle(0, 3, 9)
dms_ex5 = DMSAngle(0, 0, 18)
# dms_ex6 = DMSAngle(-0, 0, -18)
dms_ex6 = DMSAngle(0, 0, 18, positive=False)
dms_exs = [dms_ex, dms_ex2, dms_ex3, dms_ex4, dms_ex5, dms_ex6]

dms_str = "123 44 55.5"
dms_str2 = "12 34 30"
dms_str3 = "-12 34 30"
dms_str4 = "0 3 9"
dms_str5 = "0 0 18"
dms_str6 = "-0 0 18"
dms_strs = [dms_str, dms_str2, dms_str3, dms_str4, dms_str5, dms_str6]

ddm_ex = DDMAngle(123, 44.925)
ddm_ex2 = DDMAngle(12, 34.5)
ddm_ex3 = DDMAngle(-12, -34.5)
ddm_ex4 = DDMAngle(0, 3.15)
ddm_ex5 = DDMAngle(0, 0.3)
ddm_ex6 = DDMAngle(0, 0.3, positive=False)
ddm_exs = [ddm_ex, ddm_ex2, ddm_ex3, ddm_ex4, ddm_ex5, ddm_ex6]

ddm_str = "123 44.925"
ddm_str2 = "12 34.5"
ddm_str3 = "-12 34.5"
ddm_str4 = "0 3.15"
ddm_str5 = "0 0.3"
ddm_str6 = "-0 0.3"
ddm_strs = [ddm_str, ddm_str2, ddm_str3, ddm_str4, ddm_str5, ddm_str6]

gon_ex = 137.4986111111111
gon_ex2 = 13.97222222222222
gon_ex3 = -13.97222222222222
gon_ex4 = 0.05833333333333333
gon_ex5 = 0.00555555555555555
gon_ex6 = -0.00555555555555555
gon_exs = [gon_ex, gon_ex2, gon_ex3, gon_ex4, gon_ex5, gon_ex6]

gona_ex = GONAngle(137.4986111111111)
gona_ex2 = GONAngle(13.97222222222222)
gona_ex3 = GONAngle(-13.97222222222222)
gona_ex4 = GONAngle(0.05833333333333333)
gona_ex5 = GONAngle(0.00555555555555555)
gona_ex6 = GONAngle(-0.00555555555555555)
gona_exs = [gona_ex, gona_ex2, gona_ex3, gona_ex4, gona_ex5, gona_ex6]


class TestConvert(unittest.TestCase):
    def setUp(self):
        self.testData = []
        degreeValues = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256]
        dec_places = 13
        error = 10 ** -(dec_places - 4)
        for deg in degreeValues:
            for min in range(60):
                for sec in range(60):
                    if sec:
                        hp_minus = float(
                            f"{deg:4d}.{min:02d}{sec-1:02d}" + "9" * (dec_places - 4)
                        )
                        dec_minus = deg + (min / 60.0 + (sec - error) / 3600.0)
                        gon_minus = 400.0 / 360.0 * dec_minus
                        rad_minus = pi / 180.0 * dec_minus
                        self.testData.append(
                            [hp_minus, dec_minus, gon_minus, rad_minus]
                        )
                    hp = float(f"{deg:4d}.{min:02d}{sec:02d}")
                    hp_plus = float(
                        f"{deg:4d}.{min:02d}{sec:02d}" + "0" * (dec_places - 5) + "1"
                    )
                    dec = deg + (min / 60.0 + sec / 3600.0)
                    gon = 400.0 / 360.0 * dec
                    rad = pi / 180.0 * dec
                    self.testData.append([hp, dec, gon, rad])
                    dec_plus = deg + (min / 60.0 + (sec + error) / 3600.0)
                    gon_plus = 400.0 / 360.0 * dec_plus
                    rad_plus = pi / 180.0 * dec_plus
                    self.testData.append([hp_plus, dec_plus, gon_plus, rad_plus])

    def test_DECAngle(self):
        # Test DECAngle Methods
        for num, ex in enumerate(deca_exs):
            self.assertEqual(ex.rad(), rad_exs[num])
            self.assertEqual(-ex.rad(), -rad_exs[num])
            self.assertEqual(ex.dec(), dec_exs[num])
            self.assertEqual(-ex.dec(), -dec_exs[num])
            self.assertEqual(round(ex.hpa(), 13), hpa_exs[num])
            self.assertEqual(round(-ex.hpa(), 13), -hpa_exs[num])
            self.assertEqual(ex.hp(), hp_exs[num])
            self.assertEqual(-ex.hp(), -hp_exs[num])
            self.assertAlmostEqual(ex.gon(), gon_exs[num], 13)
            self.assertAlmostEqual(-ex.gon(), -gon_exs[num], 13)
            self.assertEqual(round(ex.gona(), 13), round(gona_exs[num], 13))
            self.assertEqual(round(-ex.gona(), 13), round(-gona_exs[num], 13))
            self.assertEqual(round(ex.dms(), 10), dms_exs[num])
            self.assertEqual(round(-ex.dms(), 10), -dms_exs[num])
            self.assertEqual(round(ex.ddm(), 12), ddm_exs[num])
            self.assertEqual(round(-ex.ddm(), 12), -ddm_exs[num])
            self.assertEqual(str(ex), str(dec_exs[num]))
            self.assertEqual(DECAngle(str(ex)), ex)
            self.assertEqual(int(ex), int(dec_exs[num]))
            self.assertEqual(float(ex), dec_exs[num])
            self.assertEqual(round(ex), DECAngle(round(dec_exs[num])))

        # Test HPAngle Representation
        self.assertEqual(repr(deca_ex), "{DECAngle: +123.74875}")
        self.assertEqual(repr(deca_ex3), "{DECAngle: -12.575}")

        # Test DECAngle Overloads
        self.assertEqual(dec_ex + dec_ex2, (deca_ex + deca_ex2).dec())
        self.assertEqual(dec_ex2 + dec_ex, (deca_ex2 + deca_ex).dec())
        self.assertEqual(dec_ex - dec_ex2, (deca_ex - deca_ex2).dec())
        self.assertEqual(dec_ex2 - dec_ex, (deca_ex2 - deca_ex).dec())
        self.assertEqual(dec_ex * 5, (deca_ex * 5).dec())
        self.assertEqual(5 * dec_ex, (5 * deca_ex).dec())
        self.assertEqual(dec_ex / 3, (deca_ex / 3).dec())
        self.assertEqual(abs(-deca_ex), deca_ex)
        self.assertEqual(-deca_ex2, deca_ex3)
        self.assertEqual(deca_ex2, abs(deca_ex3))
        self.assertTrue(deca_ex == deca_ex)
        self.assertFalse(deca_ex == deca_ex2)
        self.assertTrue(deca_ex != deca_ex2)
        self.assertFalse(deca_ex != deca_ex)
        self.assertTrue(deca_ex > deca_ex2)
        self.assertFalse(deca_ex2 > deca_ex)
        self.assertTrue(deca_ex2 < deca_ex)
        self.assertFalse(deca_ex < deca_ex2)
        with self.assertRaises(TypeError):
            deca_ex * "a"
        with self.assertRaises(TypeError):
            "a" * deca_ex
        with self.assertRaises(TypeError):
            deca_ex / "a"
        with self.assertRaises(TypeError):
            deca_ex + "a"
        with self.assertRaises(TypeError):
            "a" + deca_ex
        with self.assertRaises(TypeError):
            deca_ex - "a"
        with self.assertRaises(TypeError):
            "a" - deca_ex

    def test_HPAngle(self):
        # Test HPAngle Methods
        for num, ex in enumerate(hpa_exs):
            self.assertAlmostEqual(ex.rad(), rad_exs[num], 16)
            self.assertAlmostEqual(-ex.rad(), -rad_exs[num], 16)
            self.assertAlmostEqual(ex.dec(), dec_exs[num], 13)
            self.assertAlmostEqual(-ex.dec(), -dec_exs[num], 13)
            self.assertEqual(round(ex.deca(), 13), deca_exs[num])
            self.assertEqual(round(-ex.deca(), 13), -deca_exs[num])
            self.assertAlmostEqual(ex.hp(), hp_exs[num], 16)
            self.assertAlmostEqual(-ex.hp(), -hp_exs[num], 16)
            self.assertAlmostEqual(ex.gon(), gon_exs[num], 13)
            self.assertAlmostEqual(-ex.gon(), -gon_exs[num], 13)
            self.assertEqual(round(ex.gona(), 13), round(gona_exs[num], 13))
            self.assertEqual(round(-ex.gona(), 13), round(-gona_exs[num], 13))
            self.assertEqual(round(ex.dms(), 10), dms_exs[num])
            self.assertEqual(round(-ex.dms(), 10), -dms_exs[num])
            self.assertEqual(round(ex.ddm(), 12), ddm_exs[num])
            self.assertEqual(round(-ex.ddm(), 12), -ddm_exs[num])
            self.assertEqual(str(ex), str(hp_exs[num]))
            self.assertEqual(HPAngle(str(ex)), ex)
            self.assertEqual(int(ex), int(hp_exs[num]))
            self.assertEqual(float(ex), hp_exs[num])
            self.assertEqual(round(ex), HPAngle(round(hp_exs[num])))

        # Test HPAngle Representation
        self.assertEqual(repr(hpa_ex), "{HPAngle: +123.44555}")
        self.assertEqual(repr(hpa_ex3), "{HPAngle: -12.343}")

        # Test DMSAngle Overloads
        self.assertAlmostEqual(dec_ex + dec_ex2, (hpa_ex + hpa_ex2).dec(), 12)
        self.assertAlmostEqual(dec_ex2 + dec_ex, (hpa_ex2 + hpa_ex).dec(), 12)
        self.assertAlmostEqual(dec_ex - dec_ex2, (hpa_ex - hpa_ex2).dec(), 12)
        self.assertAlmostEqual(dec_ex2 - dec_ex, (hpa_ex2 - hpa_ex).dec(), 12)
        self.assertAlmostEqual(dec_ex * 5, (hpa_ex * 5).dec(), 12)
        self.assertAlmostEqual(5 * dec_ex, (5 * hpa_ex).dec(), 12)
        self.assertAlmostEqual(dec_ex / 3, (hpa_ex / 3).dec(), 12)
        self.assertEqual(abs(-hpa_ex), hpa_ex)
        self.assertEqual(-hpa_ex2, hpa_ex3)
        self.assertEqual(hpa_ex2, abs(hpa_ex3))
        self.assertTrue(hpa_ex == hpa_ex)
        self.assertFalse(hpa_ex == hpa_ex2)
        self.assertTrue(hpa_ex != hpa_ex2)
        self.assertFalse(hpa_ex != hpa_ex)
        self.assertTrue(hpa_ex > hpa_ex2)
        self.assertFalse(hpa_ex2 > hpa_ex)
        self.assertTrue(hpa_ex2 < hpa_ex)
        self.assertFalse(hpa_ex < hpa_ex2)
        with self.assertRaises(TypeError):
            hpa_ex * "a"
        with self.assertRaises(TypeError):
            "a" * hpa_ex
        with self.assertRaises(TypeError):
            hpa_ex / "a"
        with self.assertRaises(TypeError):
            hpa_ex + "a"
        with self.assertRaises(TypeError):
            "a" + hpa_ex
        with self.assertRaises(TypeError):
            hpa_ex - "a"
        with self.assertRaises(TypeError):
            "a" - hpa_ex

        # Test Validity Checking of HP format input
        with self.assertRaises(ValueError):
            HPAngle(123.7)
        with self.assertRaises(ValueError):
            HPAngle(123.007)

    def test_GONAngle(self):
        # Test GONAngle Methods
        for num, ex in enumerate(gona_exs):
            self.assertAlmostEqual(ex.rad(), rad_exs[num], 15)
            self.assertAlmostEqual(-ex.rad(), -rad_exs[num], 15)
            self.assertAlmostEqual(ex.dec(), dec_exs[num], 13)
            self.assertAlmostEqual(-ex.dec(), -dec_exs[num], 13)
            self.assertEqual(round(ex.deca(), 13), deca_exs[num])
            self.assertEqual(round(-ex.deca(), 13), -deca_exs[num])
            self.assertAlmostEqual(ex.hp(), hp_exs[num], 16)
            self.assertAlmostEqual(-ex.hp(), -hp_exs[num], 16)
            self.assertEqual(ex.hpa(), hpa_exs[num])
            self.assertEqual(-ex.hpa(), -hpa_exs[num])
            self.assertAlmostEqual(ex.gon(), gon_exs[num], 13)
            self.assertEqual(round(ex.dms(), 10), dms_exs[num])
            self.assertEqual(round(-ex.dms(), 10), -dms_exs[num])
            self.assertEqual(round(ex.ddm(), 12), ddm_exs[num])
            self.assertEqual(round(-ex.ddm(), 12), -ddm_exs[num])
            self.assertEqual(str(ex), str(gon_exs[num]))
            self.assertEqual(GONAngle(str(ex)), ex)
            self.assertEqual(int(ex), int(gon_exs[num]))
            self.assertEqual(float(ex), gon_exs[num])
            self.assertEqual(round(ex), GONAngle(round(gon_exs[num])))

        # Test GONAngle Representation
        self.assertEqual(repr(gona_ex), "{GONAngle: +137.4986111111111}")
        self.assertEqual(repr(gona_ex3), "{GONAngle: -13.97222222222222}")

        # Test DMSAngle Overloads
        self.assertAlmostEqual(dec_ex + dec_ex2, (gona_ex + gona_ex2).dec(), 12)
        self.assertAlmostEqual(dec_ex2 + dec_ex, (gona_ex2 + gona_ex).dec(), 12)
        self.assertAlmostEqual(dec_ex - dec_ex2, (gona_ex - gona_ex2).dec(), 12)
        self.assertAlmostEqual(dec_ex2 - dec_ex, (gona_ex2 - gona_ex).dec(), 12)
        self.assertAlmostEqual(dec_ex * 5, (gona_ex * 5).dec(), 12)
        self.assertAlmostEqual(5 * dec_ex, (5 * gona_ex).dec(), 12)
        self.assertAlmostEqual(dec_ex / 3, (gona_ex / 3).dec(), 12)
        self.assertEqual(abs(-gona_ex), gona_ex)
        self.assertEqual(-gona_ex2, gona_ex3)
        self.assertEqual(gona_ex2, abs(gona_ex3))
        self.assertTrue(gona_ex == gona_ex)
        self.assertFalse(gona_ex == gona_ex2)
        self.assertTrue(gona_ex != gona_ex2)
        self.assertFalse(gona_ex != gona_ex)
        self.assertTrue(gona_ex > gona_ex2)
        self.assertFalse(gona_ex2 > gona_ex)
        self.assertTrue(gona_ex2 < gona_ex)
        self.assertFalse(gona_ex < gona_ex2)
        with self.assertRaises(TypeError):
            gona_ex * "a"
        with self.assertRaises(TypeError):
            "a" * gona_ex
        with self.assertRaises(TypeError):
            gona_ex / "a"
        with self.assertRaises(TypeError):
            gona_ex + "a"
        with self.assertRaises(TypeError):
            "a" + gona_ex
        with self.assertRaises(TypeError):
            gona_ex - "a"
        with self.assertRaises(TypeError):
            "a" - gona_ex

    def test_DMSAngle(self):
        # Test DMSAngle Methods
        for num, ex in enumerate(dms_exs):
            self.assertAlmostEqual(ex.rad(), rad_exs[num], 16)
            self.assertAlmostEqual(-ex.rad(), -rad_exs[num], 16)
            self.assertAlmostEqual(ex.dec(), dec_exs[num], 16)
            self.assertAlmostEqual(-ex.dec(), -dec_exs[num], 16)
            self.assertEqual(round(ex.deca(), 16), deca_exs[num])
            self.assertEqual(round(-ex.deca(), 16), -deca_exs[num])
            self.assertAlmostEqual(ex.hp(), hp_exs[num], 16)
            self.assertAlmostEqual(-ex.hp(), -hp_exs[num], 16)
            self.assertEqual(ex.hpa(), hpa_exs[num])
            self.assertEqual(-ex.hpa(), -hpa_exs[num])
            self.assertAlmostEqual(ex.gon(), gon_exs[num], 13)
            self.assertAlmostEqual(-ex.gon(), -gon_exs[num], 13)
            self.assertEqual(round(ex.gona(), 13), round(gona_exs[num], 13))
            self.assertEqual(round(-ex.gona(), 13), round(-gona_exs[num], 13))
            self.assertEqual(ex.ddm(), ddm_exs[num])
            self.assertEqual(-ex.ddm(), -ddm_exs[num])
            self.assertEqual((round(ex, 2)).hp(), round(hp_exs[num], 5))
            self.assertEqual((round(-ex, 2)).hp(), round(-hp_exs[num], 5))

        # Test DMSAngle Sign Conventions
        self.assertEqual(
            -dec_ex, DMSAngle(-dms_ex.degree, dms_ex.minute, dms_ex.second).dec()
        )
        self.assertEqual(
            dec_ex, DMSAngle(dms_ex.degree, -dms_ex.minute, -dms_ex.second).dec()
        )
        self.assertAlmostEqual(
            -dec_ex4, DMSAngle(0, -dms_ex4.minute, dms_ex4.second).dec(), 9
        )
        self.assertAlmostEqual(
            dec_ex4, DMSAngle(0, dms_ex4.minute, dms_ex4.second).dec(), 9
        )
        self.assertEqual(-dec_ex5, DMSAngle(0, 0, -dms_ex5.second).dec())
        self.assertEqual(dec_ex5, DMSAngle(0, 0, dms_ex5.second).dec())
        self.assertEqual(-dms_ex3, DMSAngle(12, 34, -30))
        self.assertTrue(dms_ex.positive)
        self.assertFalse((-dms_ex).positive)
        self.assertTrue(dms_ex4.positive)
        self.assertFalse((-dms_ex4).positive)
        self.assertTrue(dms_ex5.positive)
        self.assertFalse((-dms_ex5).positive)
        self.assertFalse(DMSAngle(-1, 2, 3).positive)
        self.assertTrue(DMSAngle(1, -2, 3).positive)
        self.assertTrue(DMSAngle(1, 2, -3).positive)
        self.assertFalse(DMSAngle(0, -1, 2).positive)
        self.assertFalse(DMSAngle(0, 0, -3).positive)
        self.assertFalse(DMSAngle(0, 1, 2, positive=False).positive)
        self.assertFalse(DMSAngle(-0.0, 1, 2).positive)
        self.assertEqual(repr(dms_ex), "{DMSAngle: +123d 44m 55.5s}")
        self.assertEqual(repr(dms_ex3), "{DMSAngle: -12d 34m 30s}")

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
            dms_ex * "a"
        with self.assertRaises(TypeError):
            "a" * dms_ex
        with self.assertRaises(TypeError):
            dms_ex / "a"
        with self.assertRaises(TypeError):
            dms_ex + "a"
        with self.assertRaises(TypeError):
            "a" + dms_ex
        with self.assertRaises(TypeError):
            dms_ex - "a"
        with self.assertRaises(TypeError):
            "a" - dms_ex

        # Test reading in formatted strings
        for num, ex in enumerate(dms_strs):
            self.assertEqual(DMSAngle(ex), dms_exs[num])
            self.assertEqual(ex, str(dms_exs[num]))

    def test_DDMAngle(self):
        # Test DDMAngle Methods
        for num, ex in enumerate(ddm_exs):
            self.assertAlmostEqual(ex.rad(), rad_exs[num], 16)
            self.assertAlmostEqual(-ex.rad(), -rad_exs[num], 16)
            self.assertAlmostEqual(ex.dec(), dec_exs[num], 16)
            self.assertAlmostEqual(-ex.dec(), -dec_exs[num], 16)
            self.assertEqual(round(ex.deca(), 16), deca_exs[num])
            self.assertEqual(round(-ex.deca(), 16), -deca_exs[num])
            self.assertAlmostEqual(ex.hp(), hp_exs[num], 16)
            self.assertAlmostEqual(-ex.hp(), -hp_exs[num], 16)
            self.assertEqual(ex.hpa(), hpa_exs[num])
            self.assertEqual(-ex.hpa(), -hpa_exs[num])
            self.assertAlmostEqual(ex.gon(), gon_exs[num], 13)
            self.assertAlmostEqual(-ex.gon(), -gon_exs[num], 13)
            self.assertEqual(round(ex.gona(), 13), round(gona_exs[num], 13))
            self.assertEqual(round(-ex.gona(), 13), round(-gona_exs[num], 13))
            self.assertEqual(round(ex.dms(), 13), dms_exs[num])
            self.assertEqual(round(-ex.dms(), 13), -dms_exs[num])
            self.assertEqual((round(ex)).hp(), round(hp_exs[num], 2))
            self.assertEqual((round(-ex)).hp(), round(-hp_exs[num], 2))

        # Test DMSAngle Sign Conventions
        self.assertEqual(-dec_ex, DDMAngle(-dms_ex.degree, ddm_ex.minute).dec())
        self.assertEqual(dec_ex, DDMAngle(dms_ex.degree, -ddm_ex.minute).dec())
        self.assertAlmostEqual(-dec_ex4, DDMAngle(0, -ddm_ex4.minute).dec(), 9)
        self.assertAlmostEqual(dec_ex4, DDMAngle(0, ddm_ex4.minute).dec(), 9)
        self.assertEqual(-ddm_ex3, DDMAngle(12, 34.5))
        self.assertTrue(ddm_ex.positive)
        self.assertFalse((-ddm_ex).positive)
        self.assertTrue(ddm_ex4.positive)
        self.assertFalse((-ddm_ex4).positive)
        self.assertTrue(ddm_ex5.positive)
        self.assertFalse((-ddm_ex5).positive)
        self.assertFalse(DDMAngle(-1, 2).positive)
        self.assertTrue(DDMAngle(1, -2).positive)
        self.assertTrue(DDMAngle(1, 2).positive)
        self.assertFalse(DDMAngle(0, -1).positive)
        self.assertFalse(DDMAngle(0, 1, positive=False).positive)
        self.assertFalse(DDMAngle(-0.0, 1).positive)
        self.assertEqual(repr(ddm_ex), "{DDMAngle: +123d 44.925m}")
        self.assertEqual(repr(ddm_ex3), "{DDMAngle: -12d 34.5m}")

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
            ddm_ex * "a"
        with self.assertRaises(TypeError):
            "a" * ddm_ex
        with self.assertRaises(TypeError):
            ddm_ex / "a"
        with self.assertRaises(TypeError):
            ddm_ex + "a"
        with self.assertRaises(TypeError):
            "a" + ddm_ex
        with self.assertRaises(TypeError):
            ddm_ex - "a"
        with self.assertRaises(TypeError):
            "a" - ddm_ex

        # Test reading in formatted strings
        for num, ex in enumerate(ddm_strs):
            self.assertEqual(DDMAngle(ex), ddm_exs[num])
            self.assertEqual(ex, str(ddm_exs[num]))

    def test_angles_interoperability(self):
        self.assertEqual(DMSAngle(1, 2, 3) + DDMAngle(2, 3), DMSAngle(3, 5, 3))
        self.assertEqual(DMSAngle(3, 2, 0) - DDMAngle(2, 2.5), DMSAngle(0, 59, 30))
        self.assertEqual(DDMAngle(2, 3) + DMSAngle(1, 2, 3), DDMAngle(3, 5.05))
        self.assertEqual(DDMAngle(3, 2) - DMSAngle(2, 2, 30), DDMAngle(0, 59.5))

    def test_dec2hp(self):
        for num, ex in enumerate(hp_exs):
            self.assertAlmostEqual(ex, dec2hp(dec_exs[num]), 13)
            self.assertAlmostEqual(-ex, dec2hp(-dec_exs[num]), 13)
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(hp, dec2hp(dec), 13)
            self.assertAlmostEqual(-hp, dec2hp(-dec), 13)

    def test_dec2hpa(self):
        for num, ex in enumerate(dec_exs):
            self.assertEqual(dec2hpa(ex), hpa_exs[num])
            self.assertEqual(dec2hpa(-ex), -hpa_exs[num])

    def test_dec2gon(self):
        for num, ex in enumerate(dec_exs):
            self.assertAlmostEqual(dec2gon(ex), gon_exs[num], 13)
            self.assertAlmostEqual(dec2gon(-ex), -gon_exs[num], 13)
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(gon, dec2gon(dec), 13)
            self.assertAlmostEqual(-gon, dec2gon(-dec), 13)

    def test_dec2gona(self):
        for num, ex in enumerate(dec_exs):
            self.assertEqual(round(dec2gona(ex), 13), round(gona_exs[num], 13))
            self.assertEqual(round(dec2gona(-ex), 13), round(-gona_exs[num], 13))

    def test_dec2dms(self):
        self.assertEqual(dms_ex, dec2dms(dec_ex))
        self.assertEqual(-dms_ex, dec2dms(-dec_ex))

    def test_dec2ddm(self):
        self.assertEqual(ddm_ex, dec2ddm(dec_ex))
        self.assertEqual(-ddm_ex, dec2ddm(-dec_ex))

    def test_hp2dec(self):
        for num, ex in enumerate(dec_exs):
            self.assertAlmostEqual(ex, hp2dec(hp_exs[num]), 13)
            self.assertAlmostEqual(-ex, hp2dec(-hp_exs[num]), 13)
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(dec, hp2dec(hp), 13)
            self.assertAlmostEqual(-dec, hp2dec(-hp), 13)

        self.assertAlmostEqual(0, hp2dec(0), 13)
        self.assertAlmostEqual(258, hp2dec(258), 13)
        self.assertAlmostEqual(
            hp2dec(hp_exs[0]) + hp2dec(hp_exs[1]), dec_exs[0] + dec_exs[1], 13
        )
        # Test that invalid minutes and seconds components raise errors
        # (Minutes and/or Seconds can't be greater than 60, therefore
        #  123.718 - representing 123 degrees, 71 minutes, 80 seconds
        #  should raise an Error)
        with self.assertRaises(ValueError):
            hp2dec(123.7)
        with self.assertRaises(ValueError):
            hp2dec(123.318)
        with self.assertRaises(ValueError):
            hp2dec(123.718)
        with self.assertRaises(ValueError):
            hp2dec(-123.7)
        with self.assertRaises(ValueError):
            hp2dec(-123.318)
        with self.assertRaises(ValueError):
            hp2dec(-123.718)

    def test_hp2deca(self):
        for num, ex in enumerate(hp_exs):
            self.assertEqual(round(hp2deca(ex), 13), deca_exs[num])
            self.assertEqual(round(hp2deca(-ex), 13), -deca_exs[num])

    def test_hp2gon(self):
        for num, ex in enumerate(hp_exs):
            self.assertAlmostEqual(hp2gon(ex), gon_exs[num], 13)
            self.assertAlmostEqual(hp2gon(-ex), -gon_exs[num], 13)
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(gon, hp2gon(hp), 13)
            self.assertAlmostEqual(-gon, hp2gon(-hp), 13)

    def test_hp2gona(self):
        for num, ex in enumerate(hp_exs):
            self.assertEqual(round(hp2gona(ex), 13), round(gona_exs[num], 13))
            self.assertEqual(round(hp2gona(-ex), 13), round(-gona_exs[num], 13))

    def test_hp2rad(self):
        for num, ex in enumerate(hp_exs):
            self.assertAlmostEqual(hp2rad(ex), rad_exs[num], 15)
            self.assertAlmostEqual(hp2rad(-ex), -rad_exs[num], 15)
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(rad, hp2rad(hp), 15)
            self.assertAlmostEqual(-rad, hp2rad(-hp), 15)

    def test_hp2dms(self):
        self.assertEqual(dms_ex.degree, hp2dms(hp_ex).degree)
        self.assertEqual(dms_ex.minute, hp2dms(hp_ex).minute)
        self.assertAlmostEqual(dms_ex.second, hp2dms(hp_ex).second, 10)

        self.assertEqual(-dms_ex, round(hp2dms(-hp_ex), 10))
        self.assertEqual(dms_ex.degree, hp2dms(-hp_ex).degree)
        self.assertEqual(dms_ex.minute, hp2dms(-hp_ex).minute)
        self.assertAlmostEqual(dms_ex.second, hp2dms(-hp_ex).second, 10)

    def test_hp2ddm(self):
        self.assertEqual(ddm_ex, hp2ddm(hp_ex))
        self.assertEqual(-ddm_ex, hp2ddm(-hp_ex))

    def test_gon2dec(self):
        for num, ex in enumerate(gon_exs):
            self.assertAlmostEqual(gon2dec(ex), dec_exs[num], 14)
            self.assertAlmostEqual(gon2dec(-ex), -dec_exs[num], 14)
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(dec, gon2dec(gon), delta=5.8e-14)
            self.assertAlmostEqual(-dec, gon2dec(-gon), delta=5.8e-14)

    def test_gon2deca(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(round(gon2deca(ex), 14), deca_exs[num])
            self.assertEqual(round(gon2deca(-ex), 14), -deca_exs[num])

    def test_gon2hp(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(gon2hp(ex), hp_exs[num])
            self.assertEqual(gon2hp(-ex), -hp_exs[num])
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(hp, gon2hp(gon), 13)
            self.assertAlmostEqual(-hp, gon2hp(-gon), 13)

    def test_gon2hpa(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(gon2hpa(ex), hpa_exs[num])
            self.assertEqual(gon2hpa(-ex), -hpa_exs[num])

    def test_gon2rad(self):
        for num, ex in enumerate(gon_exs):
            self.assertAlmostEqual(gon2rad(ex), rad_exs[num], 15)
            self.assertAlmostEqual(gon2rad(-ex), -rad_exs[num], 15)
        for check in self.testData:
            hp, dec, gon, rad = check
            self.assertAlmostEqual(rad, gon2rad(gon), 13)
            self.assertAlmostEqual(-rad, gon2rad(-gon), 13)

    def test_gon2dms(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(round(gon2dms(ex), 10), dms_exs[num])
            self.assertEqual(round(gon2dms(-ex), 10), -dms_exs[num])

    def test_gon2ddm(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(round(gon2ddm(ex), 12), ddm_exs[num])
            self.assertEqual(round(gon2ddm(-ex), 12), -ddm_exs[num])

    def test_dd2sec(self):
        self.assertEqual(dd2sec(1), 3600)
        self.assertEqual(dd2sec(-1), -3600)
        self.assertAlmostEqual(dd2sec(hp2dec(0.0001)), 1, 12)
        self.assertAlmostEqual(dd2sec(hp2dec(-0.0001)), -1, 12)
        self.assertAlmostEqual(dd2sec(hp2dec(0.00001)), 0.1, 12)
        self.assertEqual(dd2sec(dec_ex4), 189)
        self.assertEqual(dd2sec(-dec_ex4), -189)
        self.assertEqual(dd2sec(dec_ex2), 45270)
        self.assertEqual(dd2sec(-dec_ex2), -45270)

    def test_angular_typecheck(self):
        class_exs = [deca_exs, hpa_exs, gona_exs, dms_exs, ddm_exs]
        for class_ex in class_exs:
            for num, ex in enumerate(class_ex):
                self.assertAlmostEqual(angular_typecheck(ex), dec_exs[num], 13)
                self.assertAlmostEqual(angular_typecheck(-ex), -dec_exs[num], 13)
        for ex in dec_exs:
            self.assertEqual(angular_typecheck(ex), ex)
            self.assertEqual(angular_typecheck(-ex), -ex)


if __name__ == "__main__":
    unittest.main()

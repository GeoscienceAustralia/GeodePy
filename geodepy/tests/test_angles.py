import unittest
from math import radians

from geodepy.angles import (dec2hp, hp2dec,
                            DMSAngle, DDMAngle, HPAngle, DECAngle,
                            dec2dms, dec2ddm, dec2hpa, dec2gon,
                            hp2dms, hp2ddm, hp2rad, hp2deca, hp2gon,
                            gon2dec, gon2deca, gon2rad, gon2hp, gon2hpa,
                            gon2dms, gon2ddm,
                            dd2sec)

rad_exs = [radians(123.74875), radians(12.575), radians(-12.575),
           radians(0.0525), radians(0.005)]

dec_ex = 123.74875
dec_ex2 = 12.575
dec_ex3 = -12.575
dec_ex4 = 0.0525
dec_ex5 = 0.005
dec_exs = [dec_ex, dec_ex2, dec_ex3, dec_ex4, dec_ex5]

deca_ex = DECAngle(123.74875)
deca_ex2 = DECAngle(12.575)
deca_ex3 = DECAngle(-12.575)
deca_ex4 = DECAngle(0.0525)
deca_ex5 = DECAngle(0.005)
deca_exs = [deca_ex, deca_ex2, deca_ex3, deca_ex4, deca_ex5]

hp_ex = 123.44555
hp_ex2 = 12.3430
hp_ex3 = -12.3430
hp_ex4 = 0.0309
hp_ex5 = 0.0018
hp_exs = [hp_ex, hp_ex2, hp_ex3, hp_ex4, hp_ex5]

hpa_ex = HPAngle(123.44555)
hpa_ex2 = HPAngle(12.3430)
hpa_ex3 = HPAngle(-12.3430)
hpa_ex4 = HPAngle(0.0309)
hpa_ex5 = HPAngle(0.0018)
hpa_exs = [hpa_ex, hpa_ex2, hpa_ex3, hpa_ex4, hpa_ex5]

dms_ex = DMSAngle(123, 44, 55.5)
dms_ex2 = DMSAngle(12, 34, 30)
dms_ex3 = DMSAngle(-12, -34, -30)
dms_ex4 = DMSAngle(0, 3, 9)
dms_ex5 = DMSAngle(0, 0, 18)
dms_exs = [dms_ex, dms_ex2, dms_ex3, dms_ex4, dms_ex5]

ddm_ex = DDMAngle(123, 44.925)
ddm_ex2 = DDMAngle(12, 34.5)
ddm_ex3 = DDMAngle(-12, -34.5)
ddm_ex4 = DDMAngle(0, 3.15)
ddm_ex5 = DDMAngle(0, 0.3)
ddm_exs = [ddm_ex, ddm_ex2, ddm_ex3, ddm_ex4, ddm_ex5]

gon_ex = 137.4986111111111
gon_ex2 = 13.97222222222222
gon_ex3 = -13.97222222222222
gon_ex4 = 0.05833333333333333
gon_ex5 = 0.00555555555555555
gon_exs = [gon_ex, gon_ex2, gon_ex3, gon_ex4, gon_ex5]


class TestConvert(unittest.TestCase):
    def test_dec2hp(self):
        for num, ex in enumerate(hp_exs):
            self.assertAlmostEqual(ex, dec2hp(dec_exs[num]), 13)
            self.assertAlmostEqual(-ex, dec2hp(-dec_exs[num]), 13)

    def test_hp2dec(self):
        for num, ex in enumerate(dec_exs):
            self.assertAlmostEqual(ex, hp2dec(hp_exs[num]), 13)
            self.assertAlmostEqual(-ex, hp2dec(-hp_exs[num]), 13)
        self.assertAlmostEqual(hp2dec(hp_exs[0]) + hp2dec(hp_exs[1]),
                               dec_exs[0] + dec_exs[1], 13)
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
            self.assertEqual(round(ex.dms(), 10), dms_exs[num])
            self.assertEqual(round(-ex.dms(), 10), -dms_exs[num])
            self.assertEqual(round(ex.ddm(), 12), ddm_exs[num])
            self.assertEqual(round(-ex.ddm(), 12), -ddm_exs[num])

        # Test HPAngle Representation
        self.assertEqual(repr(deca_ex), '{DECAngle: +123.74875}')
        self.assertEqual(repr(deca_ex3), '{DECAngle: -12.575}')

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
            deca_ex * 'a'
        with self.assertRaises(TypeError):
            'a' * deca_ex
        with self.assertRaises(TypeError):
            deca_ex / 'a'
        with self.assertRaises(TypeError):
            deca_ex + 'a'
        with self.assertRaises(TypeError):
            'a' + deca_ex
        with self.assertRaises(TypeError):
            deca_ex - 'a'
        with self.assertRaises(TypeError):
            'a' - deca_ex

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
            self.assertEqual(ex.ddm(), ddm_exs[num])
            self.assertEqual(-ex.ddm(), -ddm_exs[num])

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
            self.assertEqual(ex.dms(), dms_exs[num])
            self.assertEqual(-ex.dms(), -dms_exs[num])

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
            self.assertEqual(round(ex.dms(), 10), dms_exs[num])
            self.assertEqual(round(-ex.dms(), 10), -dms_exs[num])
            self.assertEqual(round(ex.ddm(), 12), ddm_exs[num])
            self.assertEqual(round(-ex.ddm(), 12), -ddm_exs[num])

        # Test HPAngle Representation
        self.assertEqual(repr(hpa_ex), '{HPAngle: +123.44555}')
        self.assertEqual(repr(hpa_ex3), '{HPAngle: -12.343}')

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
            hpa_ex * 'a'
        with self.assertRaises(TypeError):
            'a' * hpa_ex
        with self.assertRaises(TypeError):
            hpa_ex / 'a'
        with self.assertRaises(TypeError):
            hpa_ex + 'a'
        with self.assertRaises(TypeError):
            'a' + hpa_ex
        with self.assertRaises(TypeError):
            hpa_ex - 'a'
        with self.assertRaises(TypeError):
            'a' - hpa_ex

    def test_angles_interoperability(self):
        self.assertEqual(DMSAngle(1, 2, 3) + DDMAngle(2, 3), DMSAngle(3, 5, 3))
        self.assertEqual(DMSAngle(3, 2, 0) - DDMAngle(2, 2.5),
                         DMSAngle(0, 59, 30))
        self.assertEqual(DDMAngle(2, 3) + DMSAngle(1, 2, 3), DDMAngle(3, 5.05))
        self.assertEqual(DDMAngle(3, 2) - DMSAngle(2, 2, 30),
                         DDMAngle(0, 59.5))

    def test_dec2dms(self):
        self.assertEqual(dms_ex, dec2dms(dec_ex))
        self.assertEqual(-dms_ex, dec2dms(-dec_ex))

    def test_dec2ddm(self):
        self.assertEqual(ddm_ex, dec2ddm(dec_ex))
        self.assertEqual(-ddm_ex, dec2ddm(-dec_ex))

    def test_dec2hpa(self):
        for num, ex in enumerate(dec_exs):
            self.assertEqual(dec2hpa(ex), hpa_exs[num])
            self.assertEqual(dec2hpa(-ex), -hpa_exs[num])

    def test_dec2gon(self):
        for num, ex in enumerate(dec_exs):
            self.assertAlmostEqual(dec2gon(ex), gon_exs[num], 13)
            self.assertAlmostEqual(dec2gon(-ex), -gon_exs[num], 13)

    def test_hp2dms(self):
        self.assertEqual(dms_ex.degree, hp2dms(hp_ex).degree)
        self.assertEqual(dms_ex.minute, hp2dms(hp_ex).minute)
        self.assertAlmostEqual(dms_ex.second, hp2dms(hp_ex).second, 10)

        self.assertEqual(-dms_ex.sign, hp2dms(-hp_ex).sign)
        self.assertEqual(dms_ex.degree, hp2dms(-hp_ex).degree)
        self.assertEqual(dms_ex.minute, hp2dms(-hp_ex).minute)
        self.assertAlmostEqual(dms_ex.second, hp2dms(-hp_ex).second, 10)

    def test_hp2deca(self):
        for num, ex in enumerate(hp_exs):
            self.assertEqual(round(hp2deca(ex), 13), deca_exs[num])
            self.assertEqual(round(hp2deca(-ex), 13), -deca_exs[num])

    def test_hp2gon(self):
        for num, ex in enumerate(hp_exs):
            self.assertAlmostEqual(hp2gon(ex), gon_exs[num], 13)
            self.assertAlmostEqual(hp2gon(-ex), -gon_exs[num], 13)

    def test_hp2rad(self):
        for num, ex in enumerate(hp_exs):
            self.assertEqual(hp2rad(ex), rad_exs[num])
            self.assertEqual(hp2rad(-ex), -rad_exs[num])

    def test_hp2ddm(self):
        self.assertEqual(ddm_ex, hp2ddm(hp_ex))
        self.assertEqual(-ddm_ex, hp2ddm(-hp_ex))

    def test_gon2dec(self):
        for num, ex in enumerate(gon_exs):
            self.assertAlmostEqual(gon2dec(ex), dec_exs[num], 14)
            self.assertAlmostEqual(gon2dec(-ex), -dec_exs[num], 14)

    def test_gon2deca(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(round(gon2deca(ex), 14), deca_exs[num])
            self.assertEqual(round(gon2deca(-ex), 14), -deca_exs[num])

    def test_gon2rad(self):
        for num, ex in enumerate(gon_exs):
            self.assertAlmostEqual(gon2rad(ex), rad_exs[num], 15)
            self.assertAlmostEqual(gon2rad(-ex), -rad_exs[num], 15)

    def test_gon2hp(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(gon2hp(ex), hp_exs[num])
            self.assertEqual(gon2hp(-ex), -hp_exs[num])

    def test_gon2hpa(self):
        for num, ex in enumerate(gon_exs):
            self.assertEqual(gon2hpa(ex), hpa_exs[num])
            self.assertEqual(gon2hpa(-ex), -hpa_exs[num])

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


if __name__ == '__main__':
    unittest.main()

import unittest
from geodepy.survey import first_vel_params, first_vel_corrn, precise_inst_ht, va_conv


class TestSurveyConvert(unittest.TestCase):
    def test_first_vel_params(self):
        temperature = 12
        pressure = 1013.25
        relative_humidity = 60
        edm_wavelength = 0.850
        param_c = 281.781
        param_d = 79.393
        params_new = first_vel_params(edm_wavelength, temperature, pressure, relative_humidity)
        self.assertEqual(round(params_new[0], 3), param_c)
        self.assertEqual(round(params_new[1], 3), param_d)

    def test_first_vel_corrn(self):
        params = first_vel_params(0.85)
        raw_obs_distance = 1117.8517
        obs_temperature = 6.8
        obs_pressure = 960.8
        obs_relative_humidity = 58.6
        correction = first_vel_corrn(raw_obs_distance, params, obs_temperature, obs_pressure, obs_relative_humidity)
        corrected_obs_distance = raw_obs_distance + correction
        self.assertEqual(round(corrected_obs_distance, 4), 1117.8618)

    def test_va_conv(self):
        test1 = va_conv(92.24305555555556, 2116.254)
        self.assertAlmostEqual(test1[0], -2.243055555555556, 14)
        self.assertEqual(test1[1], 2116.254)
        self.assertAlmostEqual(test1[2], 2114.6325, 5)
        self.assertAlmostEqual(test1[3], -82.82744, 5)
        test2 = va_conv(83.18694444444445, 145.145, 1.62, 0.05)
        self.assertEqual(test2[0], 7.427622182644272)
        self.assertAlmostEqual(test2[1], 145.33961, 5)
        self.assertAlmostEqual(test2[2], 144.12006, 5)
        self.assertAlmostEqual(test2[3], 18.78858, 5)
        test3 = va_conv(1.25, 12.23, 2.05)
        self.assertAlmostEqual(test3[0], 88.92943808746813, 14)
        self.assertAlmostEqual(test3[1], 14.27958, 5)
        self.assertAlmostEqual(test3[2], 0.266796, 5)
        self.assertAlmostEqual(test3[3], 14.277089, 5)
        test4 = va_conv(264.70444444444445, 3.25, 1.4)
        self.assertAlmostEqual(test4[0], 27.71315483945665, 14)
        self.assertAlmostEqual(test4[1], 3.65546, 5)
        self.assertAlmostEqual(test4[2], 3.23613, 5)
        self.assertAlmostEqual(test4[3], 1.69995, 5)
        with self.assertRaises(ValueError):
            va_conv(-3.62, 1.5)
            va_conv('brian', 16)

    def test_precise_inst_ht(self):
        va1 = 99.50920833333333
        va2 = 95.67597222222223
        va3 = 91.80145833333333
        va4 = 87.90518055555556
        test1 = precise_inst_ht([va1, va2, va3, va4], 0.4, 0.8)
        self.assertEqual(test1[0], 1.78458)
        self.assertEqual(test1[1], 0.00083)
        test2 = precise_inst_ht([va3, va1, va2], 0.4, 0.8)
        self.assertEqual(test2[0], 1.78441)
        self.assertEqual(test2[1], 0.00109)
        with self.assertRaises(ValueError):
            precise_inst_ht([va1, va2], 0.4, 0.8)


if __name__ == '__main__':
    unittest.main()

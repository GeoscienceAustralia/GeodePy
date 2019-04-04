import unittest
from geodepy.survey import first_vel_params, first_vel_corrn


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


if __name__ == '__main__':
    unittest.main()

import unittest

from geodepy.transform import (
    conform7,
    conform14,
    atrf2014_to_gda2020,
    transform_gda2020_to_atrf2014,
    transform_atrf2014_to_gda2020,
    transform_mga94_to_mga2020,
    transform_mga2020_to_mga94,
)
from geodepy.constants import itrf2014_to_gda2020, gda94_to_gda2020
from datetime import date


class TestTransforms(unittest.TestCase):
    # Tests equality between values produced using grid2geo and geo2grid

    def test_conform7(self):
        # Replication of tests in GDA2020 Tech Manual v1.2 - Sect 3.1.1
        alic_gda1994 = (-4052051.7643, 4212836.2017, -2545106.0245)
        alic_gda2020 = (-4052052.7379, 4212835.9897, -2545104.5898)
        alic_gda2020_comp = conform7(*alic_gda1994, gda94_to_gda2020)
        alic_gda1994_comp = conform7(*alic_gda2020, -gda94_to_gda2020)
        assert abs(alic_gda2020_comp[0] - alic_gda2020[0]) < 5e-5
        assert abs(alic_gda2020_comp[1] - alic_gda2020[1]) < 5e-5
        assert abs(alic_gda2020_comp[2] - alic_gda2020[2]) < 5e-5
        assert abs(alic_gda1994_comp[0] - alic_gda1994[0]) < 5e-5
        assert abs(alic_gda1994_comp[1] - alic_gda1994[1]) < 5e-5
        assert abs(alic_gda1994_comp[2] - alic_gda1994[2]) < 5e-5

    def test_conform14(self):
        # Replication of tests in GDA2020 Tech Manual v1.2 - Sect 3.3.1
        alic_gda2020 = (-4052052.7373, 4212835.9835, -2545104.5867)
        alic_itrf14at2018 = (-4052052.6588, 4212835.9938, -2545104.6946)
        alic_itrf14at2018_comp = conform14(
            *alic_gda2020, date(2018, 1, 1), -itrf2014_to_gda2020
        )
        alic_gda2020_comp = conform14(
            *alic_itrf14at2018, date(2018, 1, 1), itrf2014_to_gda2020
        )
        assert abs(alic_itrf14at2018_comp[0] - alic_itrf14at2018[0]) < 5e-5
        assert abs(alic_itrf14at2018_comp[1] - alic_itrf14at2018[1]) < 5e-5
        assert abs(alic_itrf14at2018_comp[2] - alic_itrf14at2018[2]) < 5e-5
        assert abs(alic_gda2020_comp[0] - alic_gda2020[0]) < 5e-5
        assert abs(alic_gda2020_comp[1] - alic_gda2020[1]) < 5e-5
        assert abs(alic_gda2020_comp[2] - alic_gda2020[2]) < 5e-5

    def test_transform_atrf2014_to_gda2020(self):
        alic_gda2020 = (-4052052.7373, 4212835.9835, -2545104.5867)
        alic_itrf14at2018 = (-4052052.6588, 4212835.9938, -2545104.6946)
        alic_itrf14at2018_comp = transform_gda2020_to_atrf2014(
            *alic_gda2020, date(2018, 1, 1)
        )
        alic_gda2020_comp = transform_atrf2014_to_gda2020(
            *alic_itrf14at2018, date(2018, 1, 1)
        )
        assert abs(alic_itrf14at2018_comp[0] - alic_itrf14at2018[0]) < 5e-5
        assert abs(alic_itrf14at2018_comp[1] - alic_itrf14at2018[1]) < 5e-5
        assert abs(alic_itrf14at2018_comp[2] - alic_itrf14at2018[2]) < 5e-5
        assert abs(alic_gda2020_comp[0] - alic_gda2020[0]) < 5e-5
        assert abs(alic_gda2020_comp[1] - alic_gda2020[1]) < 5e-5
        assert abs(alic_gda2020_comp[2] - alic_gda2020[2]) < 5e-5

    def test_transform_mga94_to_mga2020(self):
        alic_mga94 = (53, 386352.3979, 7381850.7689, 603.3466)
        alic_mga20 = (53, 386353.2343, 7381852.2986, 603.2489)
        # Test with no ellipsoid height supplied
        alic_mga20_noellht_comp = transform_mga94_to_mga2020(
            alic_mga94[0], alic_mga94[1], alic_mga94[2]
        )
        assert (alic_mga20_noellht_comp[0] - alic_mga20[0]) == 0
        assert abs(alic_mga20_noellht_comp[1] - alic_mga20[1]) < 5e-5
        assert abs(alic_mga20_noellht_comp[2] - alic_mga20[2]) < 5e-5
        assert alic_mga20_noellht_comp[3] == 0
        # Test with ellipsoid height supplied
        alic_mga20_ellht_comp = transform_mga94_to_mga2020(
            alic_mga94[0], alic_mga94[1], alic_mga94[2], alic_mga94[3]
        )
        assert (alic_mga20_ellht_comp[0] - alic_mga20[0]) == 0
        assert abs(alic_mga20_ellht_comp[1] - alic_mga20[1]) < 5e-5
        assert abs(alic_mga20_ellht_comp[2] - alic_mga20[2]) < 5e-5
        assert abs(alic_mga20_ellht_comp[3] - alic_mga20[3]) < 5e-5

    def test_transform_mga2020_to_mga94(self):
        alic_mga94 = (53, 386352.3979, 7381850.7689, 603.3466)
        alic_mga20 = (53, 386353.2343, 7381852.2986, 603.2489)
        # Test with no ellipsoid height supplied
        alic_mga94_noellht_comp = transform_mga2020_to_mga94(
            alic_mga20[0], alic_mga20[1], alic_mga20[2]
        )
        assert (alic_mga94_noellht_comp[0] - alic_mga94[0]) == 0
        assert abs(alic_mga94_noellht_comp[1] - alic_mga94[1]) < 5e-5
        assert abs(alic_mga94_noellht_comp[2] - alic_mga94[2]) < 5e-5
        assert alic_mga94_noellht_comp[3] == 0
        # Test with ellipsoid height supplied
        alic_mga94_ellht_comp = transform_mga2020_to_mga94(
            alic_mga20[0], alic_mga20[1], alic_mga20[2], alic_mga20[3]
        )
        assert (alic_mga94_ellht_comp[0] - alic_mga94[0]) == 0
        assert abs(alic_mga94_ellht_comp[1] - alic_mga94[1]) < 5e-5
        assert abs(alic_mga94_ellht_comp[2] - alic_mga94[2]) < 5e-5
        assert abs(alic_mga94_ellht_comp[3] - alic_mga94[3]) < 5e-5


if __name__ == "__main__":
    unittest.main()

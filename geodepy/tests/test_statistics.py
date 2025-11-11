import unittest
from geodepy import statistics
import numpy as np

lat = 19.4792453
lon = 70.69315634
vcv = np.array([[1.44, -1.32, 1.32], [-1.32, 1.22, -1.20], [1.32, -1.20, 1.20]])
var = np.array([[1.44], [1.20], [1.20]])


class TestStatistics(unittest.TestCase):
    def test_rotation_matrix(self):
        expected_result = np.array(
            [
                [-0.94376147, -0.1102527, 0.3117028],
                [0.33062712, -0.31471177, 0.88974278],
                [0.0, 0.94276235, 0.33346538],
            ]
        )

        result = statistics.rotation_matrix(lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3x3(self):
        expected_result = np.array(
            [
                [2.23971834, -1.86955194, 0.3599339],
                [-1.86955194, 1.55096504, -0.29615085],
                [0.3599339, -0.29615085, 0.06931662],
            ]
        )

        result = statistics.vcv_cart2local(vcv, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3x1(self):
        expected_result = np.array([[1.41376457], [1.20291736], [1.22331807]])

        result = statistics.vcv_cart2local(var, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3X2(self):
        v_cart = np.zeros((3, 2))

        with self.assertRaises(ValueError):
            statistics.vcv_cart2local(v_cart, 0.0, 0.0)

    def test_vcv_local2cart_3X3(self):
        expected_result = np.array(
            [
                [0.44517136, -1.15507667, 0.44844663],
                [-1.15507667, 2.95156126, -1.15249271],
                [0.44844663, -1.15249271, 0.46326737],
            ]
        )

        result = statistics.vcv_local2cart(vcv, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_local2cart_3X1(self):
        expected_result = np.array([[1.41376457], [1.20291736], [1.22331807]])

        result = statistics.vcv_cart2local(var, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_local2cart_3X2(self):
        v_cart = np.zeros((3, 2))

        with self.assertRaises(ValueError):
            statistics.vcv_local2cart(v_cart, lat, lon)

    def test_error_ellipse(self):
        expected_result = (1.6292867776015223, 0.07365185899111726, 132.61817915463692)

        result = statistics.error_ellipse(vcv)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_relative_error(self):
        lat1 = -33.371389383333
        lon1 = 145.673034975000
        var1 = np.array(
            [
                [1.7671344090e-04, -8.9986817000e-05, 1.1789951440e-04],
                [-8.9986817000e-05, 1.0963890720e-04, -7.1820721020e-05],
                [1.1789951440e-04, -7.1820721020e-05, 1.4015891560e-04],
            ]
        )
        var2 = np.array(
            [
                [1.7105233790e-04, -8.8389059620e-05, 1.1156043490e-04],
                [-8.8389059620e-05, 1.0754378720e-04, -6.8637225330e-05],
                [1.1156043490e-04, -6.8637225330e-05, 1.3097055280e-04],
            ]
        )
        cov12 = np.array(
            [
                [1.5518691860e-04, -7.8107376300e-05, 1.0312885480e-04],
                [-7.8392614650e-05, 9.5884929150e-05, -6.2159214580e-05],
                [1.0310493330e-04, -6.1905790840e-05, 1.2226743540e-04],
            ]
        )

        expected_results = (
            0.003122110653270988,
            0.0028032035987053208,
            133.65771030508103,
            0.008473125158578584,
        )

        results = statistics.relative_error(lat1, lon1, var1, var2, cov12)
        np.testing.assert_almost_equal(expected_results, results)
        self.assertEqual(type(expected_results), type(results))

    def test_circ_hz_pu(self):
        a = 1
        b = 0

        expected_result = 1.96079

        result = statistics.circ_hz_pu(a, b)
        self.assertEqual(expected_result, result)

    def test_k_val95_typeError(self):
        dof = [[], {}, ""]

        for item in dof:
            with self.assertRaises(TypeError):
                statistics.k_val95(item)

    def test_k_val95_less_1(self):
        dof = -1

        expected_result = statistics.ttable_p95[0]

        result = statistics.k_val95(dof)
        self.assertEqual(expected_result, result)

    def test_k_val95_greater_120(self):
        dof = 121

        expected_result = 1.96

        result = statistics.k_val95(dof)
        self.assertEqual(expected_result, result)

    def test_k_val95_between_1_and_120(self):
        dof = 100

        expected_result = statistics.ttable_p95[dof - 1]

        result = statistics.k_val95(dof)
        self.assertEqual(expected_result, result)


if __name__ == "__main__":
    unittest.main()

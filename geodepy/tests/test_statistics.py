import unittest
from geodepy import statistics
import numpy as np

lat = 19.4792453
lon = 70.69315634
vcv = np.array([
    [1.44, -1.32, 1.32],
    [-1.32, 1.22, -1.20],
    [1.32, -1.20, 1.20]
])
var = np.array([
    [1.44], [1.20], [1.20]
])


class TestStatistics(unittest.TestCase):
    def test_rotation_matrix(self):
        expected_result = np.array([
            [-0.94376147, -0.1102527, 0.3117028],
            [0.33062712, -0.31471177, 0.88974278],
            [0.0,  0.94276235,  0.33346538],
        ])

        result = statistics.rotation_matrix(lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3x3(self):
        expected_result = np.array([
            [2.23753205, -1.8674709, 0.35405044],
            [-1.8674709, 1.54898417, -0.2905506],
            [0.35405044, -0.2905506, 0.05348378]
        ])

        result = statistics.vcv_cart2local(vcv, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3x1(self):
        expected_result = np.array([
            [1.41376457], [1.20291736], [1.22331807]
        ])

        result = vcv_cart2local(var, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3X2(self):
        v_cart = np.zeros((3, 2))

        with self.assertRaises(ValueError):
            statistics.vcv_cart2local(v_cart, 0.0, 0.0)

    def test_vcv_local2cart_3X3(self):
        expected_result = np.array([
            [0.44492825, -1.15577063,  0.45052547],
            [-1.15577063,  2.94958039, -1.14655874],
            [0.45052547, -1.14655874, 0.44549136]
        ])

        result = statistics.vcv_local2cart(vcv, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_local2cart_3X1(self):
        expected_result = np.array([
            [1.41376457], [1.20291736], [1.22331807]
        ])

        result = statistics.vcv_cart2local(var, lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_local2cart_3X2(self):
        v_cart = np.zeros((3, 2))

        with self.assertRaises(ValueError):
            statistics.vcv_local2cart(v_cart, lat, lon)

    def test_error_ellipse(self):
        expected_result = (1.6292867776015223, 0.07365185899111726,
                           132.61817915463692)

        result = statistics.error_ellipse(vcv)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

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
                statistics.k_val95(dof)

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


if __name__ == '__main__':
    unittest.main()

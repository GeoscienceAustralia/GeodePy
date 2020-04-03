import unittest
from geodepy import statistics
from geodepy.statistics import np


class TestStatistics(unittest.TestCase):
    def test_rotation_matrix(self):
        lat = 19.4792
        lon = 70.6931

        expected_result = np.array([
            [-0.94376114, -0.11025276, 0.31170376],
            [0.33062805, -0.31471096, 0.88974272],
            [0.0,  0.94276261,  0.33346463],
        ])

        result = statistics.rotation_matrix(lat, lon)

        np.testing.assert_almost_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3X3(self):
        lat = 19.4792453
        lon = 70.69315634
        v_cart = np.array([
            [1.44, -1.32, 1.32],
            [-1.32, 1.20, -1.20],
            [1.32, -1.20, 1.20]
        ])
        expected_result = np.array([
            [2.23753205, -1.8674709, 0.35405044],
            [-1.8674709, 1.54898417, -0.2905506],
            [0.35405044, -0.2905506, 0.05348378]
        ])

        result = statistics.vcv_cart2local(v_cart, lat, lon)
        print(result)

        self.assertAlmostEqual(expected_result, result, 6)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_and_vcv_local2cart2_3X2(self):
        lat = 0.0
        lon = 0.0
        v_cart = np.array([
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
        ])

        with self.assertRaises(SystemExit):
            statistics.vcv_cart2local(v_cart, lat, lon)

        with self.assertRaises(SystemExit):
            statistics.vcv_local2cart(v_cart, lat, lon)

    def test_vcv_cart2local_and_vcv_local2cart2_2X3(self):
        lat = 0.0
        lon = 0.0
        v_cart = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ])

        with self.assertRaises(SystemExit):
            statistics.vcv_cart2local(v_cart, lat, lon)

        with self.assertRaises(SystemExit):
            statistics.vcv_local2cart(v_cart, lat, lon)

    def test_vcv_cart2local_and_vcv_local2cart2_1X3(self):
        lat = 0.0
        lon = 0.0
        v_cart = np.array([
            [0.0],
            [0.0],
            [0.0],
        ])
        expected_result = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ])

        result_cart2local = statistics.vcv_cart2local(v_cart, lat, lon)
        result_local2cart = statistics.vcv_local2cart(v_cart, lat, lon)

        np.testing.assert_array_equal(expected_result, result_cart2local)
        np.testing.assert_array_equal(expected_result, result_local2cart)
        self.assertEqual(type(expected_result), type(result_cart2local))
        self.assertEqual(type(expected_result), type(result_local2cart))

    def test_error_ellipse(self):
        vcv = np.array([
            [90, 0, 0],
            [0, 90, 0],
            [0, 0, 90],
        ])
        expected_result = (9.486832980505138, 9.486832980505138, 90.0)

        result = statistics.error_ellipse(vcv)

        self.assertEqual(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_circ_hz_pu(self):
        a = 1
        b = 0

        expeted_result = 1.96079

        result = statistics.circ_hz_pu(a, b)
        self.assertEqual(expeted_result, result)

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

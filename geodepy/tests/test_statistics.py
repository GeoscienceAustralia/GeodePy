import sys
sys.path.append("../geodepy")

import unittest
from geodepy import statistics
from geodepy.statistics import np
from math import radians, sin, cos, sqrt, atan2, degrees


class TestStatistics(unittest.TestCase):
    def test_rotation_matrix(self):
        
        lat = 19.4792
        lon =  70.6931

        expected_result = np.array([
            [-0.94376114, -0.11025276, 0.31170376],
            [0.33062805, -0.31471096, 0.88974272],
            [0.0       ,  0.94276261,  0.33346463]
        ])
        
        result =  statistics.rotation_matrix(lat, lon)
        
        np.array_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))


    def test_vcv_cart2local_3X3(self):
        
        lat = 19.4792
        lon = 70.6931
        v_cart = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]
        ])

        expected_result = np.array([
            [ 0.10931491, -0.10405227,  0.2941739 ],
            [-0.10405227,  1.87664567,  0.34874419],
            [ 0.2941739,   0.34874419,  1.01403942]
        ])
        
        result = statistics.vcv_cart2local(v_cart, lat, lon)

        np.array_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))

    def test_vcv_cart2local_3X2(self):
        
        lat = 0.0
        lon = 0.0
        v_cart = np.array([
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0]
        ])

        with self.assertRaises(SystemExit):
            statistics.vcv_cart2local(v_cart, lat, lon)


    def test_vcv_cart2local_2X3(self):
        
        lat = 0.0
        lon = 0.0
        v_cart = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ])

        with self.assertRaises(SystemExit):
            statistics.vcv_cart2local(v_cart, lat, lon)


    def test_vcv_cart2local_1X3(self):
        
        lat = 0.0
        lon = 0.0
        v_cart = np.array([
            [0.0],
            [0.0],
            [0.0]
        ])

        expected_result = np.array([
            [0.0, 0.0, 0.0,],
            [0.0, 0.0, 0.0,],
            [0.0, 0.0, 0.0,]
        ])
        
        result = statistics.vcv_cart2local(v_cart, lat, lon)
        
        np.array_equal(expected_result, result)
        self.assertEqual(type(expected_result), type(result))


    def test_vcv_local2cart2(self):
        pass


    def test_error_ellipse(self):
        pass


    def test_circ_hz_pu(self):
        pass


    def test_k_val95(self):
        pass


if __name__ == '__main__':
    unittest.main()

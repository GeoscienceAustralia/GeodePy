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
        
        np.array_equal(result, expected_result)
        self.assertEqual(type(result), np.ndarray)


    def test_vcv_cart2local(self):
        pass


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
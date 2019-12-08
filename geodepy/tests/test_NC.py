# -*- coding: utf-8 -*-
"""
Created on Fri Dec 06 14:20:43 2019

@author: u84157
"""

 
import unittest
import numpy as np
import geodepy.Normal_Correction

#___________________________________________________________________________#
## Some test Cases to run, check againts ga.gov.au/ausgeoid
# Test Coordinates
Long1=141.478560
Lat1=-31.599542
Height1=368.40
Long2=141.478560
Lat2=-31.599542
Height2=3690.40
# What the output should be
RECOVERED_GRAV=9.79062606
NC=0.80179403

class TestNC(unittest.TestCase):
    def test_NC(self):
        self.assertAlmostEqual(np.asscalar(geodepy.Normal_Correction.normal_correction(Lat1,Long1,Height1,Lat2,Long2,Height2)[0]),NC,7)
    def test_Grav(self):
        self.assertAlmostEqual(np.asscalar(geodepy.Normal_Correction.normal_correction(Lat1,Long1,Height1,Lat2,Long2,Height2)[1]),RECOVERED_GRAV,7)

if __name__ == '__main__':
    unittest.main()
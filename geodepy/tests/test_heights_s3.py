# -*- coding: utf-8 -*-
"""
Created on Fri Nov 08 10:52:07 2019

@author: u84157
"""
import unittest
import numpy as np
import heights_s3

#___________________________________________________________________________#
## Some test Cases to run, check againts ga.gov.au/ausgeoid
# Test Coordinates
Long=148.716
Lat=-25.716
Height=0
# What the output should be
AVWS_H=-40.79437480
AVWS_H_STD=0.05786271
AHD_H=-42.231267
AHD_H_STD=0.10002191
DOVPM=-13.71323159
DOVPV=-2.4971921

class TestHeights(unittest.TestCase):
    def test_AVWS_H(self):
        self.assertAlmostEqual(np.asscalar(heights_s3.GPS_to_AVWS(Lat,Long,Height)[0]),AVWS_H,7)
    def test_AVWS_H_STD(self):
        self.assertAlmostEqual(np.asscalar(heights_s3.GPS_to_AVWS(Lat,Long,Height)[1]),AVWS_H_STD,7)
    def test_DOVPV(self):
        self.assertAlmostEqual(np.asscalar(heights_s3.DOV(Lat,Long)[1]),DOVPV,7)
    def test_DOVPM(self):
        self.assertAlmostEqual(np.asscalar(heights_s3.DOV(Lat,Long)[0]),DOVPM,7)
    def test_AHD_H(self):
        self.assertAlmostEqual(np.asscalar(heights_s3.GPS_to_AHD(Lat,Long,Height)[0]),AHD_H,7)
    def test_AHD_H_STD(self):
        self.assertAlmostEqual(np.asscalar(heights_s3.GPS_to_AHD(Lat,Long,Height)[1]),AHD_H_STD,7)


if __name__ == '__main__':
    unittest.main()
#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Geoid Module

In Development
"""

import numpy as np
from scipy import interpolate
from geodepy.conversions import dms2dd


# Define test grid points
nvals = np.array([(dms2dd(-31.51), dms2dd(133.48), -8.806),
                  (dms2dd(-31.51), dms2dd(133.49), -8.743),
                  (dms2dd(-31.52), dms2dd(133.48), -8.870),
                  (dms2dd(-31.52), dms2dd(133.49), -8.805)])

lat = dms2dd(-31.515996736)
long = dms2dd(133.483540489)



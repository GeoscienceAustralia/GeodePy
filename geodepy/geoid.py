#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Geoid Module

In Development
"""

import numpy as np
from geodepy.convert import hp2dec


# Define test grid points
nvals = np.array([(hp2dec(-31.51), hp2dec(133.48), -8.806),
                  (hp2dec(-31.51), hp2dec(133.49), -8.743),
                  (hp2dec(-31.52), hp2dec(133.48), -8.870),
                  (hp2dec(-31.52), hp2dec(133.49), -8.805)])

lat = hp2dec(-31.515996736)
long = hp2dec(133.483540489)



#!/usr/bin/env python3

from decimal import *

c_vac = 299792.458
k_0 = 0.9996

# Ellipsoid Constants
# Development - separate out UTM parameters from GRS80 constants, add additional ellipsoids
grs80 = [6378137, Decimal('298.25722210088'), 500000,
        10000000, Decimal('0.9996'), 6, -177]


def list():
    """This function lists all the variables contained within this module
    """
    list = 'c_vac\t speed of light in a vacuum (in m/s)\n'
    list += 'k_0\t\t the scale factor on the central meridian (MGA2020)'
    print(list)
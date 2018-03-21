#!/usr/bin/env python3

from decimal import *

c_vac = 299792.458
k_0 = 0.9996

# Ellipsoid Constants
class Ellipsoid(object):
    def __init__(self, semimaj, inversef):
        self.semimaj = semimaj
        self.inversef = inversef

grs80 = Ellipsoid(6378137, Decimal('298.25722210088'))

# Projections
class Projection(object):
    def __init__(self, falseeast, falsenorth, cmscale, zonewidth, initialcm):
        self.falseeast = falseeast
        self.falsenorth = falsenorth
        self.cmscale = cmscale
        self.zonewidth = zonewidth
        self.initialcm = initialcm

utm = Projection(500000, 10000000, Decimal('0.9996'), 6, -177)
# Development - separate out UTM parameters from GRS80 constants, add additional ellipsoids
# grs80 = [6378137, Decimal('298.25722210088'), 500000,
#        10000000, Decimal('0.9996'), 6, -177]

# Helmert 7 Parameter Transformation Parameters
conform_gda94to20 = [0.06155, -0.01087, -0.04019, -0.009994, -0.0394924, -0.0327221, -0.0328979]


def list():
    """This function lists all the variables contained within this module
    """
    list = 'c_vac\t speed of light in a vacuum (in m/s)\n'
    list += 'k_0\t\t the scale factor on the central meridian (MGA2020)'
    print(list)
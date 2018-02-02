from decimal import *
from math import sqrt, sin, cos, degrees, atan, atan2
from constants import grs80
from conversions import dd2dms


# Universal Transverse Mercator Projection Parameters
proj = grs80
# Calculate Projection Constants
f = float(1 / proj[1])
a = proj[0]
b = a * (1 - f)
ecc1sq = f * (2 - f)
ecc2sq = ecc1sq/(1-ecc1sq)

"""
# Test Data
x = -3753473.196
y = 3912741.0310
z = -3347959.6998
"""


def xyz2llh(x, y, z):
    """
    Input: Cartesian XYZ coordinate in metres
    Output: Latitude and Longitude in Decimal
    Degrees and Ellipsoidal Height in Metres
    """
    # Calculate Longitude
    long = atan2(y, x)
    # Calculate Latitude
    p = sqrt(x**2 + y**2)
    latinit = atan((z*(1+ecc2sq))/p)
    lat = latinit
    itercheck = 1
    while abs(itercheck) > 1e-10:
        nu = a/(sqrt(1 - ecc1sq * (sin(lat))**2))
        itercheck = lat - atan((z + nu * ecc1sq * sin(lat))/p)
        lat = atan((z + nu * ecc1sq * sin(lat))/p)
    nu = a/(sqrt(1 - ecc1sq * (sin(lat))**2))
    ellht = p/(cos(lat)) - nu
    # Convert Latitude and Longitude to Degrees
    lat = degrees(lat)
    long = degrees(long)
    return lat, long, ellht

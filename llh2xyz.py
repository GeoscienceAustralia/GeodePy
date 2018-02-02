from decimal import *
from math import sqrt, sin, cos, radians
from constants import grs80
from conversions import dd2dms, dms2dd


"""
# Debug - Test Data
lat = -31.515996736
long = 133.483540489
ellht = 144.7229
lat = dms2dd(lat)
long = dms2dd(long)
"""


# Universal Transverse Mercator Projection Parameters
proj = grs80
# Calculate Projection Constants
f = float(1 / proj[1])
a = proj[0]
b = a * (1 - f)
e2 = f * (2 - f)


def llh2xyz(lat, long, ellht):
    """
    Input: Latitude and Longitude in Decimal Degrees, Ellipsoidal Height in metres
    Output: Cartesian X, Y, Z Coordinates in metres
    """
    # Convert lat & long to radians
    lat = radians(lat)
    long = radians(long)
    # Calculate Ellipsoid Radius of Curvature in the Prime Vertical - nu
    if lat == 0:
        nu = proj[0]
    else:
        nu = a/(sqrt(1 - e2 * (sin(lat)**2)))
    # Calculate x, y, z
    x = Decimal(str((nu + ellht) * cos(lat) * cos(long)))
    y = Decimal(str((nu + ellht) * cos(lat) * sin(long)))
    z = Decimal(str(((b**2 / a**2) * nu + ellht) * sin(lat)))
    return x, y, z


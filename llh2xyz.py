from decimal import *
from math import sqrt, sin, cos, radians


def llh2xyz(lat, long, ellht):
    """
    Input: Latitude and Longitude in Decimal Degrees, Ellipsoidal Height in metres
    Output: Cartesian X, Y, Z Coordinates in metres
    """
    # Convert lat & long to radians
    lat = radians(lat)
    long = radians(long)
    # Universal Transverse Mercator Projection Parameters
    Proj = [6378137, Decimal('298.25722210088'), 500000,
            10000000, Decimal('0.9996'), 6, -177]
    # Calculate Projection Constants
    f = float(1/Proj[1])
    a = Proj[0]
    b = a*(1-f)
    e2 = f*(2 - f)
    # Calculate Ellipsoid Radius of Curvature in the Prime Vertical - nu
    if lat == 0:
        nu = Proj[0]
    else:
        nu = a/(sqrt(1 - e2 * (sin(lat)**2)))
    # Calculate x, y, z
    x = Decimal(str((nu + ellht) * cos(lat) * cos(long)))
    y = Decimal(str((nu + ellht) * cos(lat) * sin(long)))
    z = Decimal(str(((b**2 / a**2) * nu + ellht) * sin(lat)))
    return x, y, z


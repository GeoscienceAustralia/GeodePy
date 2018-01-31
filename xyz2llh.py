from decimal import *
from math import sqrt, sin, cos, degrees, atan, atan2


# Universal Transverse Mercator Projection Parameters
Proj = [6378137, Decimal('298.25722210088'), 500000,
        10000000, Decimal('0.9996'), 6, -177]
# Calculate Projection Constants
f = float(1 / Proj[1])
a = Proj[0]
b = a * (1 - f)
e2 = f * (2 - f)


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
    latcentric = atan2(p, z)
    lat = latcentric
    itercheck = 1
    while abs(itercheck) > 0.0000001:
        nu = a/(sqrt(1 - e2 * (sin(lat))**2))
        ellht = p/cos(lat) - nu
        latcheck = lat
        lat = atan(z / p * sqrt(1 - e2 * (nu / (nu + ellht))))
        itercheck = latcheck - lat
    nu = a/(sqrt(1 - e2 * (sin(lat))**2))
    if abs(degrees(lat)) >= 90:
        l = z + e2 * nu * sin(lat)
        ellht = l/sin(lat) - nu
    else:
        ellht = p/cos(lat) - nu
    # Convert Latitude and Longitude to Degrees
    lat = degrees(lat)
    long = degrees(long)
    return lat, long, ellht, locals()

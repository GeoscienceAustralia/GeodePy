from decimal import *
from math import sqrt, sin, cos, radians, atan, atan2
def xyz2llh(x, y, z):
    """
    Write Desc"""
    # Universal Transverse Mercator Projection Parameters
    Proj = [6378137,Decimal('298.25722210088'),500000,10000000,Decimal('0.9996'),6,-177]
    # Calculate Projection Constants
    f = float(1/Proj[1])
    a = Proj[0]
    b = a * (1 - f)
    e2 = f * (2 - f)
    # Calculate Longitude
    long = degrees(atan2(y,x))
    # Calculate Latitude
    r = sqrt(x**2 + y**2 + z**2)
    p = sqrt(x**2 + y**2)
    latcentric = atan2(p,z)
    lat = latcentric
    itercheck = 1
    while abs(itercheck) > 0.001:
        rn = a/(sqrt(1 - e2 * (sin(lat))**2))
        h = p/cos(lat) - rn
        latcheck = lat
        lat = atan(z/p*sqrt(1-e2*(rn/(rn+h))))
        itercheck = latcheck - lat
    rn = a/(sqrt(1 - e2 * (sin(lat))**2))
    if abs(degrees(lat)) >= 90:
        l = z + e2 * rn * sin(lat)
        h = l/sin(lat) - rn
    else:
        h = p/cos(lat) - rn
    # Convert Latitude and Longitude to Degrees
    lat = degrees(lat)
    long = degrees(long)
    return lat, long, h
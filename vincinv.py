"""
Script to calculate the distance between two points on an ellipsoid
using Vincenty's Inverse Formulae. Reads in a comma separated text
file of the format:
    Latitude1,Longitude1,Latitude2,Longitude2
and outputs a comma separated text file of the format:
    Distance,Azimuth1to2,Azimuth2to1

Distances are in metres and Latitudes, Longitudes and Azimuths are
in DDD.MMSSSSSSSS format.

Ref: http://www.icsm.gov.au/gda/tech.html
Ref: https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
"""

# Author: Josh Batchelor <josh.batchelor@ga.gov.au>

from decimal import *
from math import (pi, degrees, radians, sqrt, sin,
                  cos, tan, asin, atan, atan2)
import os
import csv
from conversions import dd2dms, dms2dd
from constants import grs80

getcontext().prec = 28

"""
# Debug - Test Data
lat1 = -37.57037203
long1 = 144.25295244
lat2 = -37.39101561
long2 = 143.55353839

# Debug - Convert DMS to DD
lat1 = dms2dd(lat1)
long1 = dms2dd(long1)
lat2 = dms2dd(lat2)
long2 = dms2dd(long2)
"""

# Universal Transverse Mercator Projection Parameters
proj = grs80
f = float(1 / proj[1])
a = proj[0]
b = a * (1 - f)


def vincinv(lat1, long1, lat2, long2):
    """
    input: Latitude and Longitude of Points 1 and 2 in Decimal Degreees

    output: Ellipsoid Arc Distance in metres between Points 1 and 2, forward
    and reverse azimuths between Points 1 and 2 in Decimal Degrees
    """
    # Same Point Test
    if lat1 == lat2 and long1 == long2:
        return 0, 0, 0
    # Equation numbering is from GDA2020 Tech Manual v1.0
    # Eq. 71
    tan_u1 = (1 - f) * tan(radians(lat1))
    u1 = atan(tan_u1)
    sin_u1 = sin(u1)
    cos_u1 = cos(u1)
    # Eq. 72
    tan_u2 = (1 - f) * tan(radians(lat2))
    u2 = atan(tan_u2)
    sin_u2 = sin(u2)
    cos_u2 = cos(u2)
    # Initial Approximation
    # Eq. 73
    long_diff = radians(long2 - long1)
    omega = long_diff
    sigma = 1
    sigma_diff = 0.5
    # Sigma Iteration
    while sigma_diff > 1e-10:
        # Eq. 74
        sin_sq_sig = ((cos_u2 * sin(long_diff)) ** 2
                      + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos(long_diff)) ** 2)
        sin_sigma = sqrt(sin_sq_sig)
        # Eq. 75
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos(long_diff)
        # Eq. 76
        tan_sigma = sin_sigma / cos_sigma
        # Iterator
        sigma_diff = abs(sigma - atan(tan_sigma))
        sigma = atan(tan_sigma)
        # Eq. 77
        sin_alpha = (cos_u1 * cos_u2 * sin(long_diff)) / sin_sigma
        alpha = asin(sin_alpha)
        # Eq. 78
        cos_2sigm = cos_sigma - ((2 * sin_u1 * sin_u2) / (cos(alpha) ** 2))
        # Eq. 79
        C = (f / 16) * (cos(alpha) ** 2) * (4 + f * (4 - 3 * (cos(alpha) ** 2)))
        # Eq. 80
        long_diff = (omega + (1 - C)
                     * f * sin_alpha
                     * (sigma + C * sin_sigma
                        * (cos_2sigm + C * sigma
                           * (-1 + 2 * cos_2sigm ** 2))))

    # Eq. 81
    u_sq = (cos(alpha) ** 2 * (a ** 2 - b ** 2)) / b ** 2
    # Eq. 82
    A = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    # Eq. 83
    B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    # Eq. 84
    delta_sigma = (B * sin_sigma
                   * (cos_2sigm + (B / 4)
                      * (cos_sigma * (-1 + 2 * cos_2sigm ** 2)
                         - (B / 6) * cos_2sigm * (-3 + 4 * sin_sigma ** 2)
                         * (-3 + 4 * cos_2sigm ** 2))))
    # Calculate Distance
    if sigma < 0:
        sigma = sigma + pi
    # Eq. 85
    dist = b * A * (sigma - delta_sigma)
    # Calculate Alpha1
    azimuth1to2 = degrees(atan2((cos_u2 * sin(long_diff)),
                                (cos_u1 * sin_u2 - sin_u1
                                 * cos_u2 * cos(long_diff))))
    if azimuth1to2 < 0:
        azimuth1to2 = azimuth1to2 + 360
    # Calculate Alpha2
    azimuth2to1 = degrees(atan2((cos_u1 * sin(long_diff)),
                                (-sin_u1 * cos_u2 + cos_u1
                                 * sin_u2 * cos(long_diff))))
    azimuth2to1 = azimuth2to1 + 180
    # Meridian Critical Case Tests
    if long1 == long2 and lat1 > lat2:
        return dist, 180, 0
    if long1 == long2 and lat1 < lat2:
        return dist, 0, 180
    return dist, azimuth1to2, azimuth2to1


def vincinvio():
    # Enter Filename
    print('Enter co-ordinate file:')
    fn = input()
    # Open Filename
    csvfile = open(fn)
    csvreader = csv.reader(csvfile)
    # Create Output File
    fn_part = (os.path.splitext(fn))
    fn_out = fn_part[0] + '_out' + fn_part[1]
    outfile = open(fn_out, 'w')
    # Write Output
    outfilewriter = csv.writer(outfile)
    outfilewriter.writerow(['Ell_Dist', 'Azimuth1to2', 'Azimuth2to1'])
    for row in csvreader:
        lat1 = dms2dd(float(row[0]))
        long1 = dms2dd(float(row[1]))
        lat2 = dms2dd(float(row[2]))
        long2 = dms2dd(float(row[3]))
        dist, azimuth1to2, azimuth2to1 = vincinv(lat1, long1, lat2, long2)
        azimuth1to2 = dd2dms(azimuth1to2)
        azimuth2to1 = dd2dms(azimuth2to1)
        output = (dist, azimuth1to2, azimuth2to1)
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()

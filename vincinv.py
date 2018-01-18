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

# import os
# import csv
# import pprint
from dd2dms import dd2dms
from decimal import *
from math import (pi, degrees, radians, sqrt, sin,
                  cos, tan, asin, acos, atan, atan2)

getcontext().prec = 28

"""Test Data
lat1 = -37.39101561
long1 = 144.25295244
lat2 = -37.39101561
long2 = 143.55353839
"""

# Universal Transverse Mercator Projection Parameters
Proj = [6378137, Decimal('298.25722210088'), 500000,
        10000000, Decimal('0.9996'), 6, -177]
f = float(1 / Proj[1])
a = Proj[0]
b = a * (1 - f)


def vincinv(lat1, long1, lat2, long2):
    """Vincentys Inverse Tool

    Takes the Latitude and Longitude of two points in Decimal Degrees
    and returns the ellipsoidal distance between the points and the
    forward and reverse azimuths at each point.
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
    omega = radians(long2 - long1)
    long_diff = omega
    sigma = 0.5
    # Sigma Iteration
    while True:
        sin_sigma = sqrt((cos_u2 * sin(long_diff)) ** 2
                         + (cos_u1 * sin_u2 - sin_u1
                            * cos_u2 * cos(long_diff)) ** 2)
        # Eq. 75
        cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(long_diff)
        # Eq. 76
        tan_sigma = sin_sigma / cos_sigma
        sigma_diff = abs(sigma - atan(tan_sigma))
        if sigma_diff < 1e-10:
            break
        sigma = atan(tan_sigma)
        # Eq. 77
        sin_alpha = (cos(u1) * cos(u2) * sin(long_diff)) / sin_sigma
        alpha = asin(sin_alpha)
        # Eq. 78
        cos_2sigm = cos_sigma - ((2 * sin(u1) * sin(u2)) / (cos(alpha)) ** 2)
        # Eq. 79
        C = (f / 16) * cos(alpha) ** 2 * (4 + f * (4 - 3 * cos(alpha) ** 2))
        # Eq. 80
        long_diff = (omega + (1 - C) * f * sin_alpha
                     * (acos(cos_sigma) + C * sin(acos(cos_sigma))
                        * (cos_2sigm + C * cos_sigma * (-1 + 2 * cos_2sigm ** 2))))
        # Eq. 81
        u2 = (cos(alpha) ** 2 * (a ** 2 - b ** 2)) / b ** 2
        # Eq. 82
        A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
        # Eq. 83
        B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
        # Eq. 84
        delta_sigma = (B * sin_sigma
                       * (cos_2sigm + B / 4
                          * (cos_sigma * (-1 + 2 * cos_2sigm ** 2)
                             - B / 6 * cos_2sigm * (-3 + 4 * sin_sigma ** 2)
                             * (-3 + 4 * cos_2sigm ** 2))))
    # Calculate Distance
    if sigma < 0:
        sigma = sigma + pi
    # Eq. 85
    dist = b * A * (sigma - delta_sigma)
    # Calculate Alpha1
    Azimuth1to2 = degrees(atan2((cos_u2 * sin(long_diff)),
                                (cos_u1 * sin_u2 - sin_u1
                                 * cos_u2 * cos(long_diff))))
    if Azimuth1to2 < 0:
        Azimuth1to2 = Azimuth1to2 + 360
    Azimuth1to2 = dd2dms(Azimuth1to2)
    # Calculate Alpha2
    Azimuth2to1 = degrees(atan2((cos_u1 * sin(long_diff)),
                                (-sin_u1 * cos_u2 + cos_u1
                                 * sin_u2 * cos(long_diff))))
    Azimuth2to1 = dd2dms(Azimuth2to1 + 180)
    # Meridian Critical Case Tests
    if long1 == long2 and lat1 > lat2:
        return dist, 180, 0
    if long1 == long2 and lat1 < lat2:
        return dist, 0, 180
    return round(dist, 3), round(Azimuth1to2, 7), round(Azimuth2to1, 7)


"""
# Enter Filename
print('Enter co-ordinate file:')
fn = input()
# Open Filename
csvFile = open(fn)
csvReader = csv.reader(csvFile)
# Create Output File
fn_part = (os.path.splitext(fn))
fn_out = fn_part[0] + '_out' + fn_part[1]
outFile = open(fn_out, 'w')
# Write Output
outFilewriter = csv.writer(outFile)
outFilewriter.writerow(['Ell_Dist', 'Azimuth1to2', 'Azimuth2to1'])
for row in csvReader:
    lat1 = DMS2dd(float(row[0]))
    long1 = DMS2dd(float(row[1]))
    lat2 = DMS2dd(float(row[2]))
    long2 = DMS2dd(float(row[3]))
    output = vincinv(lat1, long1, lat2, long2)
    outFilewriter.writerow(output)
# Close Files
outFile.close()
csvFile.close()
"""

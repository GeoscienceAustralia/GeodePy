"""
This module defines functions for using Vincentys direct
formula. Given the latitude and longitude of a point, the
geodetic azimuth and the ellipsoidal distance to a second
point, Vincenty direct formula can be used to calculate
the latitude and longitude of the second point and the
reverse azimuth.

Functions: vincdir, vincdirio

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

# Test Data
lat1 = -37.57037203
long1 = 144.25295244
azimuth1to2 = 306.5205373
dist = 54972.271


# Universal Transverse Mercator Projection Parameters
proj = grs80
f = float(1 / proj[1])
a = proj[0]
b = a * (1 - f)


def vincdir(lat1, long1, azimuth1to2, dist):
    """
    input: Latitude, Longitude of Point 1 in Decimal Degrees,
    Geodetic Azimuth from Point 1 to 2 in Decimal Degrees,
    Ellipsoidal Distance in metres.
    
    output: Latitude, Longitude of Point 1 in Decimal Degrees,
    Azimuth from Point 2 to 1 in Decimal Degrees.
    """
    azimuth1to2 = radians(azimuth1to2)
    
    # Equation numbering is from GDA2020 Tech Manual v1.0
    tan_u1 = (1 - f) * tan(radians(lat1))  # Eq. 88
    u1 = atan(tan_u1)
    sin_u1 = sin(u1)
    cos_u1 = cos(u1)
    tan_sigma1 = tan_u1 / cos(azimuth1to2)  # Eq. 89
    sigma1 = atan(tan_sigma1)
    sin_alpha = cos_u1 * sin(azimuth1to2)  # Eq. 90
    alpha = asin(sin_alpha)
    u2 = ((cos(alpha)) ** 2
          * (a ** 2 - b ** 2)
          / b ** 2)  # Eq. 91
    A = (1 + (u2 / 16384)
         * (4096 + u2
            * (-768 + u2
               * (320 - 175 * u2))))  # Eq. 92
    B = (u2 / 1024)\
        * (256 + u2
           * (-128 + u2
              * (74 - 47 * u2)))  # Eq. 93
    sigma = dist / (b * A)  # Eq. 94
    # Sigma Iteration
    while True:
        sigm2 = 2 * sigma1 + sigma  # Eq. 95
        sigma_change = (B * sin(sigma) *
                        (cos(sigm2)
                         + (B / 4)
                         * (cos(sigma)
                            * (-1 + 2 * (cos(sigm2)) ** 2)
                            - (B / 6)
                            * cos(sigm2)
                            * (-3 + 4 * sin(sigma) ** 2)
                            * (-3 + 4 * cos(sigm2) ** 2)
                            )))  # Eq. 96
        sigma = (dist / (b * A)) + sigma_change  # Eq. 97
        if abs(sigma_change) < 1e-5:
            break
    sin_sigma = sin(sigma)
    cos_sigma = cos(sigma)

    # Calculate Latitude of Pt. 2
    lat2 = atan2((sin_u1 * cos_sigma + cos_u1 * sin(sigma) * cos(azimuth1to2)),
                 ((1 - f)
                  * sqrt(sin(alpha) ** 2
                         + (sin_u1
                            * sin(sigma)
                            - cos_u1
                            * cos_sigma
                            * cos(azimuth1to2)) ** 2)
                  ))  # Eq. 98
    lat2 = degrees(lat2)

    # Calculate Longitude of Pt. 2
    long = atan2((sin_sigma * sin(azimuth1to2)),
                 (cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos(azimuth1to2)))  # Eq. 99
    C = (f / 16) * cos(alpha) ** 2 * (4 + f * (4 - 3 * cos(alpha) ** 2))  # Eq. 100
    omega = (long
             - (1 - C)
             * f
             * sin_alpha
             * (sigma + C * sin_sigma
                * (cos(sigm2) + C * cos(sigma) * (-1 + 2 * cos(sigm2) ** 2))
                ))  # Eq. 101
    long2 = float(long1) + degrees(omega)  # Eq. 102

    # Calculate Reverse Azimuth
    azimuth2to1 = atan(sin_alpha / (-sin_u1 * sin_sigma + cos_u1 * cos_sigma * cos(azimuth1to2))) + pi
    azimuth2to1 = degrees(azimuth2to1)
    return round(lat2, 11), round(long2, 11), round(azimuth2to1, 9)


def vincdirio():
    """
    No Input:
    Prompts the user for the name of a file in csv format. Data in the file
    must be in the form Latitude, Longitude of Point 1 in Degrees Minutes Seconds,
    Geodetic Azimuth from Point 1 to 2 in Degrees Minutes Seconds and Distance
    in metres with no header line.
    
    No Output:
    Uses the function vincdir to calculate for each row in the csv file the
    geographic coordinate (lat, long) of Point 2 and the Azimuth from Point 2
    to Point 1, all in Degrees Minutes Seconds. This data is written to a new
    file with the name <inputfile>_out.csv
    """
    # Enter Filename
    fn = input('Enter co-ordinate file:\n')
    # Open Filename
    csvfile = open(fn)
    csvreader = csv.reader(csvfile)
    # Create Output File
    fn_part = (os.path.splitext(fn))
    fn_out = fn_part[0] + '_out' + fn_part[1]
    outfile = open(fn_out, 'w')
    # Write Output
    outfilewriter = csv.writer(outfile)
    # outfilewriter.writerow(['Latitude2', 'Longitude2', 'azimuth2to1'])
    for row in csvreader:
        lat1 = dms2dd(float(row[0]))
        long1 = dms2dd(float(row[1]))
        azimuth1to2 = dms2dd(float(row[2]))
        dist = float(row[3])
        lat2, long2, azimuth2to1 = vincdir(lat1, long1, azimuth1to2, dist)
        lat2 = dd2dms(lat2)
        long2 = dd2dms(long2)
        azimuth2to1 = dd2dms(azimuth2to1)
        output = [lat2, long2, azimuth2to1]
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()

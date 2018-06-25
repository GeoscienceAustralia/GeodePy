#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Geodesy Module
"""

import os
import csv
from math import (pi, degrees, radians, sqrt, sin,
                  cos, tan, asin, atan, atan2)
import numpy as np
from geodepy.convert import dd2dms, dms2dd
from geodepy.constants import grs80


def enu2xyz(lat, long, east, north, up):
    """
    function to convert a vector in a local east, north, up reference frame to a vector in a
    cartesian x, y, z reference frame
    :param lat: latitude in decimal degrees
    :param long: longitude in decimal degrees
    :param east: in metres
    :param north: in metres
    :param up: in metres
    :return: x, y, z in metres
    """
    lat = radians(lat)
    long = radians(long)
    # Create ENU Vector
    enu = np.array([[east],
                    [north],
                    [up]])
    # Create Rotation Matrix
    rotate = np.array([[-sin(long), -sin(lat)*cos(long), cos(lat)*cos(long)],
                       [cos(long), -sin(lat)*sin(long), cos(lat)*sin(long)],
                       [0, cos(lat), sin(lat)]])
    xyz = np.dot(rotate, enu)
    # Assign to separate variables
    x = float(xyz[0])
    y = float(xyz[1])
    z = float(xyz[2])
    return x, y, z


def xyz2enu(lat, long, x, y, z):
    """
    function to convert a vector in a cartesian x, y, z reference frame to a vector in a
    local east, north, up reference frame
    :param lat: latitude in decimal degrees
    :param long: longitude in decimal degrees
    :param x: in metres
    :param y: in metres
    :param z: in metres
    :return: east, north, up in metres
    """
    lat = radians(lat)
    long = radians(long)
    # Create XYZ Vector
    xyz = np.array([[x],
                    [y],
                    [z]])
    # Create Rotation Matrix
    rotate = np.array([[-sin(long), cos(long), 0],
                       [-sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat)],
                       [cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)]])
    enu = np.dot(rotate, xyz)
    # Assign to separate variables
    east = float(enu[0])
    north = float(enu[1])
    up = float(enu[2])
    return east, north, up


def vincdir(lat1, long1, azimuth1to2, ell_dist, ellipsoid=grs80):
    """
    Vincentys Direct Formula
    :param lat1: Latitude of Point 1 (Decimal Degrees)
    :param long1: Longitude of Point 1 (Decimal Degrees)
    :param azimuth1to2: Azimuth from Point 1 to 2 (Decimal Degrees)
    :param ell_dist: Ellipsoidal Distance between Points 1 and 2 (m)
    :param ellipsoid: Ellipsoid Object
    :return: Latitude of Point 2 (Decimal Degrees),
             Longitude of Point 2 (Decimal Degrees),
             Azimuth from Point 2 to 1 (Decimal Degrees)
    """
    azimuth1to2 = radians(azimuth1to2)

    # Equation numbering is from GDA2020 Tech Manual v1.0
    # Eq. 88
    tan_u1 = (1 - ellipsoid.f) * tan(radians(lat1))
    u1 = atan(tan_u1)
    sin_u1 = sin(u1)
    cos_u1 = cos(u1)
    # Eq. 89
    tan_sigma1 = tan_u1 / cos(azimuth1to2)
    sigma1 = atan(tan_sigma1)
    # Eq. 90
    sin_alpha = cos_u1 * sin(azimuth1to2)
    alpha = asin(sin_alpha)
    # Eq. 91
    u2 = ((cos(alpha)) ** 2
          * (ellipsoid.semimaj ** 2 - ellipsoid.semimin ** 2)
          / ellipsoid.semimin ** 2)
    # Eq. 92
    A = (1 + (u2 / 16384)
         * (4096 + u2
            * (-768 + u2
               * (320 - 175 * u2))))
    # Eq. 93
    B = ((u2 / 1024)
         * (256 + u2
            * (-128 + u2
               * (74 - 47 * u2))))
    # Eq. 94
    sigma = ell_dist / (ellipsoid.semimin * A)
    # Sigma Iteration
    while True:
        # Eq. 95
        sigm2 = 2 * sigma1 + sigma
        # Eq. 96
        sigma_change = (B * sin(sigma) *
                        (cos(sigm2)
                         + (B / 4)
                         * (cos(sigma)
                            * (-1 + 2 * (cos(sigm2)) ** 2)
                            - (B / 6)
                            * cos(sigm2)
                            * (-3 + 4 * sin(sigma) ** 2)
                            * (-3 + 4 * cos(sigm2) ** 2)
                            )))
        # Eq. 97
        sigma = (ell_dist / (ellipsoid.semimin * A)) + sigma_change
        if abs(sigma_change) < 1e-5:
            break
    sin_sigma = sin(sigma)
    cos_sigma = cos(sigma)

    # Calculate Latitude of Pt. 2
    # Eq. 98
    lat2 = atan2((sin_u1 * cos_sigma + cos_u1 * sin(sigma) * cos(azimuth1to2)),
                 ((1 - ellipsoid.f)
                  * sqrt(sin(alpha) ** 2
                         + (sin_u1
                            * sin(sigma)
                            - cos_u1
                            * cos_sigma
                            * cos(azimuth1to2)) ** 2)
                  ))
    lat2 = degrees(lat2)

    # Calculate Longitude of Pt. 2
    # Eq. 99
    long = atan2((sin_sigma * sin(azimuth1to2)),
                 (cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos(azimuth1to2)))
    # Eq. 100
    C = (ellipsoid.f / 16) * cos(alpha) ** 2 * (4 + ellipsoid.f * (4 - 3 * cos(alpha) ** 2))
    # Eq. 101
    omega = (long
             - (1 - C)
             * ellipsoid.f
             * sin_alpha
             * (sigma + C * sin_sigma
                * (cos(sigm2) + C * cos(sigma) * (-1 + 2 * cos(sigm2) ** 2))
                ))
    # Eq. 102
    long2 = float(long1) + degrees(omega)

    # Calculate Reverse Azimuth
    azimuth2to1 = atan(sin_alpha / (-sin_u1 * sin_sigma + cos_u1 * cos_sigma * cos(azimuth1to2))) + pi
    azimuth2to1 = degrees(azimuth2to1)
    return round(lat2, 11), round(long2, 11), round(azimuth2to1, 9)


def vincinv(lat1, long1, lat2, long2, ellipsoid=grs80):
    """
    Vincenty's Inverse Formula
    :param lat1: Latitude of Point 1 (Decimal Degrees)
    :param long1: Longitude of Point 1 (Decimal Degrees)
    :param lat2: Latitude of Point 2 (Decimal Degrees)
    :param long2: Longitude of Point 2 (Decimal Degrees)
    :param ellipsoid: Ellipsoid Object
    :return: Ellipsoidal Distance between Points 1 and 2 (m),
             Azimuth from Point 1 to 2 (Decimal Degrees),
             Azimuth from Point 2 to 1 (Decimal Degrees)
    """
    if lat1 == lat2 and long1 == long2:
        return 0, 0, 0
    # Equation numbering is from GDA2020 Tech Manual v1.0
    # Eq. 71
    tan_u1 = (1 - ellipsoid.f) * tan(radians(lat1))
    u1 = atan(tan_u1)
    sin_u1 = sin(u1)
    cos_u1 = cos(u1)
    # Eq. 72
    tan_u2 = (1 - ellipsoid.f) * tan(radians(lat2))
    u2 = atan(tan_u2)
    sin_u2 = sin(u2)
    cos_u2 = cos(u2)
    # Initial Approximation
    # Eq. 73
    long_diff = radians(long2 - long1)
    omega = long_diff
    sigma = 1
    sigma_diff = 0.5
    itercount = 0
    # Sigma Iteration
    while sigma_diff > 1e-10 and itercount < 100:
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
        C = (ellipsoid.f / 16) * (cos(alpha) ** 2) * (4 + ellipsoid.f * (4 - 3 * (cos(alpha) ** 2)))
        # Eq. 80
        long_diff = (omega + (1 - C)
                     * ellipsoid.f * sin_alpha
                     * (sigma + C * sin_sigma
                        * (cos_2sigm + C * sigma
                           * (-1 + 2 * cos_2sigm ** 2))))
        itercount += 1

    # Eq. 81
    u_sq = (cos(alpha) ** 2 * (ellipsoid.semimaj ** 2 - ellipsoid.semimin ** 2)) / ellipsoid.semimin ** 2
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
    ell_dist = ellipsoid.semimin * A * (sigma - delta_sigma)
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
        return ell_dist, 180, 0
    if long1 == long2 and lat1 < lat2:
        return ell_dist, 0, 180
    return round(ell_dist, 3), round(azimuth1to2, 9), round(azimuth2to1, 9)


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
        ell_dist = float(row[3])
        lat2, long2, azimuth2to1 = vincdir(lat1, long1, azimuth1to2, ell_dist)
        lat2 = dd2dms(lat2)
        long2 = dd2dms(long2)
        azimuth2to1 = dd2dms(azimuth2to1)
        output = [lat2, long2, azimuth2to1]
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()


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
        ell_dist, azimuth1to2, azimuth2to1 = vincinv(lat1, long1, lat2, long2)
        azimuth1to2 = dd2dms(azimuth1to2)
        azimuth2to1 = dd2dms(azimuth2to1)
        output = (ell_dist, azimuth1to2, azimuth2to1)
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()


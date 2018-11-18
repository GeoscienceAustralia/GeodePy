#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Geodesy Module
"""

import os
import csv
from math import (pi, degrees, radians, sqrt, sin,
                  cos, tan, asin, acos, atan, atan2)
import numpy as np
from geodepy.convert import dec2hp, hp2dec
from geodepy.constants import grs80
from geodepy.transform import grid2geo


def enu2xyz(lat, long, east, north, up):
    """
    function to convert a vector in a local east, north, up reference frame to
    a vector in a cartesian x, y, z reference frame
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
    function to convert a vector in a cartesian x, y, z reference frame to a
    vector in a local east, north, up reference frame
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


def vincdir(lat1, lon1, azimuth1to2, ell_dist, ellipsoid=grs80):
    """
    Vincenty's Direct Formula
    :param lat1: Latitude of Point 1 (Decimal Degrees)
    :param lon1: Longitude of Point 1 (Decimal Degrees)
    :param azimuth1to2: Azimuth from Point 1 to 2 (Decimal Degrees)
    :param ell_dist: Ellipsoidal Distance between Points 1 and 2 (m)
    :param ellipsoid: Ellipsoid Object
    :return: lat2: Latitude of Point 2 (Decimal Degrees),
             lon2: Longitude of Point 2 (Decimal Degrees),
             azimuth2to1: Azimuth from Point 2 to 1 (Decimal Degrees)

    Code review: 14-08-2018 Craig Harrison
    """

    azimuth1to2 = radians(azimuth1to2)

    # Equation numbering is from the GDA2020 Tech Manual v1.0

    # Eq. 88
    u1 = atan((1 - ellipsoid.f) * tan(radians(lat1)))

    # Eq. 89
    sigma1 = atan2(tan(u1), cos(azimuth1to2))

    # Eq. 90
    alpha = asin(cos(u1) * sin(azimuth1to2))

    # Eq. 91
    u_squared = cos(alpha)**2 \
        * (ellipsoid.semimaj**2 - ellipsoid.semimin**2) \
        / ellipsoid.semimin**2

    # Eq. 92
    a = 1 + (u_squared / 16384) \
        * (4096 + u_squared * (-768 + u_squared * (320 - 175 * u_squared)))

    # Eq. 93
    b = (u_squared / 1024) \
        * (256 + u_squared * (-128 + u_squared * (74 - 47 * u_squared)))

    # Eq. 94
    sigma = ell_dist / (ellipsoid.semimin * a)

    # Iterate until the change in sigma, delta_sigma, is insignificant (< 1e-9)
    # or after 1000 iterations have been completed
    two_sigma_m = 0
    for i in range(1000):

        # Eq. 95
        two_sigma_m = 2*sigma1 + sigma

        # Eq. 96
        delta_sigma = b * sin(sigma) * (cos(two_sigma_m) + (b/4)
                                        * (cos(sigma)
                                           * (-1 + 2 * cos(two_sigma_m)**2)
                                           - (b/6) * cos(two_sigma_m)
                                           * (-3 + 4 * sin(sigma)**2)
                                           * (-3 + 4 * cos(two_sigma_m)**2)))
        new_sigma = (ell_dist / (ellipsoid.semimin * a)) + delta_sigma
        sigma_change = new_sigma - sigma
        sigma = new_sigma

        if abs(sigma_change) < 1e-12:
            break

    # Calculate the Latitude of Point 2
    # Eq. 98
    lat2 = atan2(sin(u1)*cos(sigma) + cos(u1)*sin(sigma)*cos(azimuth1to2),
                 (1 - ellipsoid.f)
                 * sqrt(sin(alpha)**2 + (sin(u1)*sin(sigma)
                        - cos(u1)*cos(sigma)*cos(azimuth1to2))**2))
    lat2 = degrees(lat2)

    # Calculate the Longitude of Point 2
    # Eq. 99
    lon = atan2(sin(sigma)*sin(azimuth1to2),
                cos(u1)*cos(sigma) - sin(u1)*sin(sigma)*cos(azimuth1to2))

    # Eq. 100
    c = (ellipsoid.f/16)*cos(alpha)**2 \
        * (4 + ellipsoid.f*(4 - 3*cos(alpha)**2))

    # Eq. 101
    omega = lon - (1-c)*ellipsoid.f*sin(alpha) \
        * (sigma + c*sin(sigma)*(cos(two_sigma_m) + c*cos(sigma)
                                 * (-1 + 2*cos(two_sigma_m)**2)))

    # Eq. 102
    lon2 = float(lon1) + degrees(omega)

    # Calculate the Reverse Azimuth
    azimuth2to1 = degrees(atan2(sin(alpha), -sin(u1)*sin(sigma)
                          + cos(u1)*cos(sigma)*cos(azimuth1to2))) + 180

    return round(lat2, 11), round(lon2, 11), round(azimuth2to1, 9)


def vincinv(lat1, lon1, lat2, lon2, ellipsoid=grs80):
    """
    Vincenty's Inverse Formula
    :param lat1: Latitude of Point 1 (Decimal Degrees)
    :param lon1: Longitude of Point 1 (Decimal Degrees)
    :param lat2: Latitude of Point 2 (Decimal Degrees)
    :param lon2: Longitude of Point 2 (Decimal Degrees)
    :param ellipsoid: Ellipsoid Object
    :return: ell_dist: Ellipsoidal Distance between Points 1 and 2 (m),
             azimuth1to2: Azimuth from Point 1 to 2 (Decimal Degrees),
             azimuth2to1: Azimuth from Point 2 to 1 (Decimal Degrees)

    Code review: 14-08-2018 Craig Harrison
    """

    # Exit if the two input points are the same
    if lat1 == lat2 and lon1 == lon2:
        return 0, 0, 0

    # Equation numbering is from the GDA2020 Tech Manual v1.0

    # Eq. 71
    u1 = atan((1 - ellipsoid.f) * tan(radians(lat1)))

    # Eq. 72
    u2 = atan((1 - ellipsoid.f) * tan(radians(lat2)))

    # Eq. 73; initial approximation
    lon = radians(lon2 - lon1)
    omega = lon

    # Iterate until the change in lambda, lambda_sigma, is insignificant
    # (< 1e-12) or after 1000 iterations have been completed
    alpha = 0
    sigma = 0
    two_sigma_m = 0
    for i in range(1000):

        # Eq. 74
        sin_sigma = sqrt((cos(u2)*sin(lon))**2
                         + (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lon))**2)

        # Eq. 75
        cos_sigma = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(lon)

        # Eq. 76
        sigma = atan2(sin_sigma, cos_sigma)

        # Eq. 77
        alpha = asin((cos(u1)*cos(u2)*sin(lon)) / sin_sigma)

        # Eq. 78
        two_sigma_m = acos(cos(sigma) - 2*sin(u1)*sin(u2) / cos(alpha)**2)

        # Eq. 79
        c = (ellipsoid.f / 16) * cos(alpha)**2 * (4 + ellipsoid.f
                                                  * (4 - 3*cos(alpha)**2))

        # Eq. 80
        new_lon = omega + (1 - c) * ellipsoid.f * sin(alpha) \
            * (sigma + c*sin(sigma) * (cos(two_sigma_m) + c*cos(sigma)
               * (-1 + 2*cos(two_sigma_m)**2)))
        delta_lon = new_lon - lon
        lon = new_lon

        if abs(delta_lon) < 1e-12:
            break

    # Eq. 81
    u_squared = cos(alpha)**2 \
        * (ellipsoid.semimaj**2 - ellipsoid.semimin**2) \
        / ellipsoid.semimin**2

    # Eq. 82
    a = 1 + (u_squared / 16384) \
        * (4096 + u_squared * (-768 + u_squared * (320 - 175 * u_squared)))

    # Eq. 83
    b = (u_squared / 1024) \
        * (256 + u_squared * (-128 + u_squared * (74 - 47 * u_squared)))

    # Eq. 84
    delta_sigma = b*sin(sigma) * (cos(two_sigma_m) + (b / 4)
                                  * (cos(sigma) * (-1 + 2*cos(two_sigma_m)**2)
                                  - (b / 6)*cos(two_sigma_m)
                                  * (-3 + 4*sin(sigma)**2)
                                  * (-3 + 4*cos(two_sigma_m)**2)))
    # Calculate the ellipsoidal distance
    # Eq. 85
    ell_dist = ellipsoid.semimin*a * (sigma - delta_sigma)

    # Calculate the azimuth from point 1 to point 2
    azimuth1to2 = degrees(atan2((cos(u2)*sin(lon)),
                                (cos(u1)*sin(u2)
                                 - sin(u1)*cos(u2)*cos(lon))))

    if azimuth1to2 < 0:
        azimuth1to2 = azimuth1to2 + 360

    # Calculate the azimuth from point 2 to point 1
    azimuth2to1 = degrees(atan2(cos(u1)*sin(lon),
                                (-sin(u1)*cos(u2)
                                 + cos(u1)*sin(u2)*cos(lon)))) + 180

    # Meridian Critical Case Tests
    #if lon1 == lon2 and lat1 > lat2:
    #    azimuth1to2 = 180
    #    azimuth2to1 = 0
    #elif lon1 == lon2 and lat1 < lat2:
    #    azimuth1to2 = 0
    #    azimuth2to1 = 180

    return round(ell_dist, 3), round(azimuth1to2, 9), round(azimuth2to1, 9)


def vincinv_utm(zone1, east1, north1, zone2, east2, north2, hemisphere1='south', hemisphere2='south', ellipsoid=grs80):
    # Convert utm to geographic
    pt1 = grid2geo(zone1, east1, north1, hemisphere1, ellipsoid)
    pt2 = grid2geo(zone2, east2, north2, hemisphere2, ellipsoid)
    # Use vincinv
    return vincinv(pt1[0], pt1[1], pt2[0], pt2[1], ellipsoid)


def vincdirio():
    """
    No Input:
    Prompts the user for the name of a file in csv format. Data in the file
    must be in the form Latitude, Longitude of Point 1 in Degrees Minutes
    Seconds, Geodetic Azimuth from Point 1 to 2 in Degrees Minutes Seconds and
    Distance in metres with no header line.

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
        lat1 = hp2dec(float(row[0]))
        long1 = hp2dec(float(row[1]))
        azimuth1to2 = hp2dec(float(row[2]))
        ell_dist = float(row[3])
        lat2, long2, azimuth2to1 = vincdir(lat1, long1, azimuth1to2, ell_dist)
        lat2 = dec2hp(lat2)
        long2 = dec2hp(long2)
        azimuth2to1 = dec2hp(azimuth2to1)
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
        lat1 = hp2dec(float(row[0]))
        long1 = hp2dec(float(row[1]))
        lat2 = hp2dec(float(row[2]))
        long2 = hp2dec(float(row[3]))
        ell_dist, azimuth1to2, azimuth2to1 = vincinv(lat1, long1, lat2, long2)
        azimuth1to2 = dec2hp(azimuth1to2)
        azimuth2to1 = dec2hp(azimuth2to1)
        output = (ell_dist, azimuth1to2, azimuth2to1)
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()


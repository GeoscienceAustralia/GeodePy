#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Transform Module

Ref1: http://www.icsm.gov.au/gda/tech.html
Ref2: http://www.mygeodesy.id.au/documents/Karney-Krueger%20equations.pdf
"""


import os
import csv
from decimal import *
from math import sqrt, log, degrees, radians, sin, cos, tan, sinh, cosh, atan, atan2
import numpy as np
from constants import grs80, utm
from conversions import dd2dms, dms2dd


getcontext().prec = 28
# Universal Transverse Mercator Projection Parameters
proj = utm

# Ellipsoidal Constants
ellipsoid = grs80

f = 1 / ellipsoid.inversef
semi_min = float(ellipsoid.semimaj * (1 - f))
ecc1sq = float(f * (2 - f))
ecc2sq = float(ecc1sq/(1 - ecc1sq))
ecc1 = sqrt(ecc1sq)
n = f / (2 - f)
n = float(n)
n2 = n ** 2


def rect_radius(inversef):
    """
    Computes the Rectifying Radius of an Ellipsoid with specified Inverse Flattening
    (See Ref 2 Equation 3)
    :param inversef: Ellipsoid Inverse Flattening
    :return: Ellipsoid Rectifying Radius
    """
    nval = (1 / float(inversef)) / (2 - (1 / float(inversef)))
    nval2 = nval ** 2
    return (ellipsoid.semimaj / (1 + nval) * ((nval2 *
                                              (nval2 *
                                               (nval2 *
                                                (25 * nval2 + 64)
                                                + 256)
                                               + 4096)
                                              + 16384)
                                              / 16384.))


def alpha_coeff(inversef):
    """
    Computes the set of Alpha coefficients of an Ellipsoid with specified Inverse Flattening
    (See Ref 2 Equation 5)
    :param inversef: Ellipsoid Inverse Flattening
    :return: Alpha coefficients a2, a4 ... a16
    """
    nval = (1 / float(inversef)) / (2 - (1 / float(inversef)))
    a2 = ((nval *
           (nval *
            (nval *
             (nval *
              (nval *
               (nval *
                ((37884525 - 75900428 * nval)
                 * nval + 42422016)
                - 89611200)
               + 46287360)
              + 63504000)
             - 135475200)
            + 101606400))
          / 203212800.)

    a4 = ((nval ** 2 *
           (nval *
            (nval *
             (nval *
              (nval *
               (nval *
                (148003883 * nval + 83274912)
                - 178508970)
               + 77690880)
              + 67374720)
             - 104509440)
            + 47174400))
          / 174182400.)

    a6 = ((nval ** 3 *
           (nval *
            (nval *
             (nval *
              (nval *
               (318729724 * nval - 738126169)
               + 294981280)
              + 178924680)
             - 234938880)
            + 81164160))
          / 319334400.)

    a8 = ((nval ** 4 *
           (nval *
            (nval *
             ((14967552000 - 40176129013 * nval) * nval + 6971354016)
             - 8165836800)
            + 2355138720))
          / 7664025600.)

    a10 = ((nval ** 5 *
            (nval *
             (nval *
              (10421654396 * nval + 3997835751)
              - 4266773472)
             + 1072709352))
           / 2490808320.)

    a12 = ((nval ** 6 *
            (nval *
             (175214326799 * nval - 171950693600)
             + 38652967262))
           / 58118860800.)

    a14 = ((nval ** 7 *
            (13700311101 - 67039739596 * nval))
           / 12454041600.)

    a16 = (1424729850961 * nval ** 8) / 743921418240.
    return a2, a4, a6, a8, a10, a12, a14, a16


def beta_coeff(inversef):
    """
    Computes the set of Beta coefficients of an Ellipsoid with specified Inverse Flattening
    (See Ref 2 Equation 23)
    :param inversef: Ellipsoid Inverse Flattening
    :return: Alpha coefficients a2, a4 ... a16
    """
    nval = (1 / float(inversef)) / (2 - (1 / float(inversef)))
    b2 = ((nval *
           (nval *
            (nval *
             (nval *
              (nval *
               (nval *
                ((37845269 - 31777436 * nval) - 43097152)
                + 42865200)
               + 752640)
              - 104428800)
             + 180633600)
            - 135475200))
          / 270950400.)

    b4 = ((nval ** 2 *
           (nval *
            (nval *
             (nval *
              (nval *
               ((-24749483 * nval - 14930208) * nval + 100683990)
               - 152616960)
              + 105719040)
             - 23224320)
            - 7257600))
          / 348364800.)

    b6 = ((nval ** 3 *
           (nval *
            (nval *
             (nval *
              (nval *
               (232468668 * nval - 101880889)
               - 39205760)
              + 29795040)
             + 28131840)
            - 22619520))
          / 638668800.)

    b8 = ((nval ** 4 *
           (nval *
            (nval *
             ((-324154477 * nval - 1433121792) * nval + 876745056)
             + 167270400)
            - 208945440))
          / 7664025600.)

    b10 = ((nval ** 5 *
            (nval *
             (nval *
              (312227409 - 457888660 * nval)
              + 67920528)
             - 70779852))
           / 2490808320.)

    b12 = ((nval ** 6 *
            (nval *
             (19841813847 * nval + 3665348512)
             - 3758062126))
           / 116237721600.)

    b14 = ((nval ** 7 *
            (1989295244 * nval - 1979471673))
           / 49816166400.)

    b16 = ((-191773887257 * nval ** 8) / 3719607091200.)
    return b2, b4, b6, b8, b10, b12, b14, b16


A = rect_radius(ellipsoid.inversef)
a = alpha_coeff(ellipsoid.inversef)
b = beta_coeff(ellipsoid.inversef)


def psfandgridconv(xi1, eta1, lat, long, cm, conf_lat):
    lat = radians(lat)
    long_diff = radians(long - cm)

    # Point Scale Factor
    p = 1
    q = 0
    for r in range(1, 9):
        p += 2*r * a[r-1] * cos(2*r * xi1) * cosh(2*r * eta1)
        q += 2*r * a[r-1] * sin(2*r * xi1) * sinh(2*r * eta1)
    q = -q
    psf = (float(proj.cmscale)
           * (A / ellipsoid.semimaj)
           * sqrt(q**2 + p**2)
           * ((sqrt(1 + (tan(lat)**2)) * sqrt(1 - ecc1sq * (sin(lat)**2)))
              / sqrt((tan(conf_lat)**2) + (cos(long_diff)**2))))

    # Grid Convergence
    grid_conv = degrees(atan(abs(q / p))
                        + atan(abs(tan(conf_lat) * tan(long_diff))
                               / sqrt(1 + tan(conf_lat)**2)))
    if cm > long and lat < 0:
        grid_conv = -grid_conv
    elif cm < long and lat > 0:
        grid_conv = -grid_conv

    return psf, grid_conv


def geo2grid(lat, long, zone=0):
    """
    Takes a geographic co-ordinate (latitude, longitude) and returns its corresponding
    Hemisphere, Zone and Projection Easting and Northing, Point Scale Factor and Grid
    Convergence. Default Projection is Universal Transverse Mercator Projection using
    GRS80 Ellipsoid parameters.
    :param lat: Latitude in Decimal Degrees
    :param long: Longitude in Decimal Degrees
    :param zone: Optional Zone Number - Only required if calculating grid co-ordinate
                outside zone boundaries
    :return: hemisphere, zone, east (m), north (m), Point Scale Factor, Grid Convergence (Decimal Degrees)
    """
    # Input Exception Handling - UTM Extents and Values
    try:
        zone = int(zone)
        if zone < 0 or zone > 60:
            raise ValueError
    except ValueError:
        print('ValueError: Invalid Zone - Zones from 1 to 60')
        return
    try:
        if lat < -80 or lat > 84:
            raise ValueError
    except ValueError:
        print('ValueError: Invalid Latitude - Latitudes from -80 to +84')
        return
    try:
        if long < -180 or long > 180:
            raise ValueError
    except ValueError:
        print('ValueError: Invalid Longitude - Longitudes from -180 to +180')
        return

    lat = radians(lat)
    # Calculate Zone
    if zone == 0:
        zone = int((float(long) - (proj.initialcm - (1.5 * proj.zonewidth))) / proj.zonewidth)
    cm = float(zone * proj.zonewidth) + (proj.initialcm - proj.zonewidth)

    # Conformal Latitude
    sigx = (ecc1 * tan(lat)) / sqrt(1 + (tan(lat) ** 2))
    sig = sinh(ecc1 * (0.5 * log((1 + sigx) / (1 - sigx))))
    conf_lat = tan(lat) * sqrt(1 + sig ** 2) - sig * sqrt(1 + (tan(lat) ** 2))
    conf_lat = atan(conf_lat)

    # Longitude Difference
    long_diff = radians(long - cm)

    # Gauss-Schreiber Ratios
    xi1 = atan(tan(conf_lat) / cos(long_diff))
    eta1x = sin(long_diff) / (sqrt(tan(conf_lat) ** 2 + cos(long_diff) ** 2))
    eta1 = log(eta1x + sqrt(1 + eta1x ** 2))

    # Transverse Mercator Ratios
    eta = eta1
    xi = xi1
    for r in range(1, 9):
        eta += a[r-1] * cos(2*r * xi1) * sinh(2*r * eta1)
        xi += a[r-1] * sin(2*r * xi1) * cosh(2*r * eta1)

    # Transverse Mercator Co-ordinates
    x = A * eta
    y = A * xi

    # Hemisphere-dependent UTM Projection Co-ordinates
    east = proj.cmscale * Decimal(str(x)) + proj.falseeast
    if y < 0:
        hemisphere = 'South'
        north = proj.cmscale * Decimal(str(y)) + proj.falsenorth
    else:
        hemisphere = 'North'
        falsenorth = 0
        north = proj.cmscale * Decimal(str(y)) + falsenorth

    # Point Scale Factor and Grid Convergence
    psf, grid_conv = psfandgridconv(xi1, eta1, degrees(lat), long, cm, conf_lat)

    return hemisphere, zone, round(float(east), 3), round(float(north), 3), round(psf, 8), grid_conv


def grid2geo(zone, east, north, hemisphere='south'):
    """
    Takes a Transverse Mercator grid co-ordinate (Zone, Easting, Northing, Hemisphere)
    and returns its corresponding Geographic Latitude and Longitude, Point Scale Factor
    and Grid Convergence. Default Projection is Universal Transverse Mercator Projection
    using GRS80 Ellipsoid parameters.
    :param zone: Zone Number - 1 to 60
    :param east: Easting (m, within 3330km of Central Meridian)
    :param north: Northing (m, 0 to 10,000,000m)
    :param hemisphere: String - 'North' or 'South'(default)
    :return: Latitude and Longitude (Decimal Degrees), Point Scale Factor, Grid Convergence (Decimal Degrees)
    """
    # Input Exception Handling - UTM Extents and Values
    try:
        zone = int(zone)
        if zone < 0 or zone > 60:
            raise ValueError
    except ValueError:
        print('ValueError: Invalid Zone - Zones from 1 to 60')
        return
    try:
        if east < -2830000 or east > 3830000:
            raise ValueError
    except ValueError:
        print('ValueError: Invalid Easting - Must be within 3330km of Central Meridian')
        return
    try:
        if north < 0 or north > 10000000:
            raise ValueError
    except ValueError:
        print('ValueError: Invalid Northing - Must be between 0 and 10,000,000m')
        return
    try:
        h = hemisphere.lower()
        if h != 'north' and h != 'south':
            raise ValueError
    except ValueError:
        print('ValueError: Invalid Hemisphere - String, either North or South')
        return

    # Transverse Mercator Co-ordinates
    x = (east - float(proj.falseeast)) / float(proj.cmscale)
    if hemisphere == 'north':
        y = north / float(proj.cmscale)
    else:
        y = (north - float(proj.falsenorth)) / float(proj.cmscale)

    # Transverse Mercator Ratios
    xi = y / A
    eta = x / A

    # Gauss-Schreiber Ratios
    eta1 = eta
    xi1 = xi
    for r in range(1, 9):
        eta1 += b[r-1] * cos(2*r * xi) * sinh(2*r * eta)
        xi1 += b[r-1] * sin(2*r * xi) * cosh(2*r * eta)

    # Conformal Latitude
    conf_lat = (sin(xi1)) / (sqrt((sinh(eta1)) ** 2 + (cos(xi1)) ** 2))
    t1 = conf_lat
    conf_lat = atan(conf_lat)

    # Finding t using Newtons Method
    def sigma(tn):
        return (sinh(ecc1
                     * 0.5
                     * log((1 + ((ecc1 * tn) / (sqrt(1 + tn ** 2))))
                           / (1 - ((ecc1 * tn) / (sqrt(1 + tn ** 2)))))))

    def ftn(tn):
        return t * sqrt(1 + (sigma(tn)) ** 2) - sigma(tn) * sqrt(1 + tn ** 2) - t1

    def f1tn(tn):
        return ((sqrt(1 + (sigma(tn)) ** 2) * sqrt(1 + tn ** 2) - sigma(tn) * tn)
                * (((1 - float(ecc1sq)) * sqrt(1 + t ** 2))
                   / (1 + (1 - float(ecc1sq)) * t ** 2)))

    diff = 1
    t = t1
    loopcount = 0
    while diff > 1e-50:
        loopcount += 1
        t_before = t
        t = t - (ftn(t) / f1tn(t))
        diff = abs(t - t_before)
    lat = degrees(atan(t))

    # Compute Longitude
    cm = float((zone * proj.zonewidth) + proj.initialcm - proj.zonewidth)
    long_diff = degrees(atan(sinh(eta1) / cos(xi1)))
    long = cm + long_diff

    # Point Scale Factor and Grid Convergence
    psf, grid_conv = psfandgridconv(xi1, eta1, lat, long, cm, conf_lat)

    return round(lat, 11), round(long, 11), round(psf, 8), grid_conv


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
        nu = ellipsoid.semimaj/(sqrt(1 - ecc1sq * (sin(lat))**2))
        itercheck = lat - atan((z + nu * ecc1sq * sin(lat))/p)
        lat = atan((z + nu * ecc1sq * sin(lat))/p)
    nu = ellipsoid.semimaj/(sqrt(1 - ecc1sq * (sin(lat))**2))
    ellht = p/(cos(lat)) - nu
    # Convert Latitude and Longitude to Degrees
    lat = degrees(lat)
    long = degrees(long)
    return lat, long, ellht


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
        nu = grs80.semimaj
    else:
        nu = ellipsoid.semimaj/(sqrt(1 - ecc1sq * (sin(lat)**2)))
    # Calculate x, y, z
    x = Decimal(str((nu + ellht) * cos(lat) * cos(long)))
    y = Decimal(str((nu + ellht) * cos(lat) * sin(long)))
    z = Decimal(str(((semi_min**2 / ellipsoid.semimaj**2) * nu + ellht) * sin(lat)))
    return x, y, z


# # conform7 Debug - Test Values ALIC
# x = -4052051.7643
# y = 4212836.2017
# z = -2545106.0245
# # Debug - Test Helmert Params (GDA94 to GDA2020, scale in ppm and rotations in seconds)
# conform_gda94to20 = [0.06155, -0.01087, -0.04019, -0.009994, -0.0394924, -0.0327221, -0.0328979]


def conform7(x, y, z, trans):
    """
    Performs a Helmert 7 Parameter Transformation using Cartesian point co-ordinates
    and a predefined transformation object.
    :param x: Cartesian X
    :param y: Cartesian Y
    :param z: Cartesian Z
    :param trans: Transformation Object
    :return: Transformed X, Y, Z Cartesian Co-ordinates
    """
    # Create XYZ Vector
    xyz_before = np.array([[x],
                           [y],
                           [z]])
    # Convert Units for Transformation Parameters
    scale = trans.sc / 1000000
    rx = radians(dms2dd(trans.rx / 10000))
    ry = radians(dms2dd(trans.ry / 10000))
    rz = radians(dms2dd(trans.rz / 10000))
    # Create Translation Vector
    translation = np.array([[trans.tx],
                            [trans.ty],
                            [trans.tz]])
    # Create Rotation Matrix
    rotation = np.array([[1., rz, -ry],
                         [-rz, 1., rx],
                         [ry, -rx, 1.]])
    # Conformal Transform Eq
    xyz_after = translation + (1 + scale) * np.dot(rotation, xyz_before)
    # Convert Vector to Separate Variables
    xtrans = float(xyz_after[0])
    ytrans = float(xyz_after[1])
    ztrans = float(xyz_after[2])
    return xtrans, ytrans, ztrans


def grid2geoio():
    """
    No Input:
    Prompts the user for the name of a file in csv format. Data in the file
    must be in the form Point ID, UTM Zone, Easting (m), Northing (m) with
    no header line.

    No Output:
    Uses the function grid2geo to convert each row in the csv file into a
    latitude and longitude in Degrees, Minutes and Seconds. This data is
    written to a new file with the name <inputfile>_out.csv
    """
    # Enter Filename
    print('Enter co-ordinate file (\.csv)\:')
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
    # Optional Header Row
    # outfilewriter.writerow(['Pt', 'Latitude', 'Longitude', 'Point Scale Factor', 'Grid Convergence'])
    for row in csvreader:
        pt_num = row[0]
        zone = float(row[1])
        east = float(row[2])
        north = float(row[3])
        # Calculate Conversion
        lat, long, psf, grid_conv = grid2geo(zone, east, north)
        lat = dd2dms(lat)
        long = dd2dms(long)
        grid_conv = dd2dms(grid_conv)
        output = [pt_num, lat, long, psf, grid_conv]
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()


def geo2gridio():
    """
    No Input:
    Prompts the user for the name of a file in csv format. Data in the file
    must be in the form Point ID, Latitude, Longitude in Decimal Degrees with
    no header line.

    No Output:
    Uses the function geo2grid to convert each row in the csv file into a
    coordinate with UTM Zone, Easting (m), Northing (m). This data is written
    to a new file with the name <inputfile>_out.csv
    """
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
    # Optional Header Row
    # outfilewriter.writerow(['Pt', 'Zone', 'Easting', 'Northing', 'Point Scale Factor', 'Grid Convergence'])
    for row in csvreader:
        pt_num = row[0]
        lat = dms2dd(float(row[1]))
        long = dms2dd(float(row[2]))
        # Calculate Conversion
        hemisphere, zone, east, north, psf, grid_conv = geo2grid(lat, long)
        grid_conv = dms2dd(grid_conv)
        output = [pt_num] + [hemisphere, zone, east, north, psf, grid_conv]
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()

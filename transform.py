"""
Functions: geo2grid, grid2geo, xyz2llh, llh2xyz, geo2gridio, grid2geoio

geo2grid:
    input: Latitude and Longitude in Decimal Degrees.

    output: Zone, Easting and Northing of a point in metres.
    (Default projection is Universal Transverse Mercator.)

grid2geo:
    input: Zone, Easting and Northing of a point in metres.
    (Default projection is Universal Transverse Mercator.)

    output: Latitude and Longitude in Decimal Degrees.

xyz2llh:
    Input: Cartesian XYZ coordinate in metres.

    Output: Latitude and Longitude in Decimal.
    Degrees and Ellipsoidal Height in Metres.


def llh2xyz:
    Input: Latitude and Longitude in Decimal Degrees, Ellipsoidal Height in metres.

    Output: Cartesian X, Y, Z Coordinates in metres.


geo2gridio:
    No Input:
    Prompts the user for the name of a file in csv format. Data in the file
    must be in the form Point ID, Latitude, Longitude in Decimal Degrees with
    no header line.

    No Output:
    Uses the function geo2grid to convert each row in the csv file into a
    coordinate with UTM Zone, Easting (m), Northing (m). This data is written
    to a new file with the name <inputfile>_out.csv

grid2geoio:
    No Input:
    Prompts the user for the name of a file in csv format. Data in the file
    must be in the form Point ID, UTM Zone, Easting (m), Northing (m) with
    no header line.

    No Output:
    Uses the function grid2geo to convert each row in the csv file into a
    latitude and longitude in Degrees, Minutes and Seconds. This data is
    written to a new file with the name <inputfile>_out.csv

Ref: http://www.icsm.gov.au/gda/tech.html
Ref: http://www.mygeodesy.id.au/documents/Karney-Krueger%20equations.pdf
"""

# Author: Josh Batchelor <josh.batchelor@ga.gov.au>

from decimal import *
from math import sqrt, log, degrees, radians, sin, cos, tan, sinh, cosh, atan, atan2
import numpy as np
import os
import csv
from constants import grs80, utm
from conversions import dd2dms, dms2dd


getcontext().prec = 28
# Universal Transverse Mercator Projection Parameters
# proj = grs80
proj = utm
# Ellipsoidal Constants
ellipsoid = grs80
f = 1 / ellipsoid.inversef
semi_maj = ellipsoid.semimaj
semi_min = float(semi_maj * (1 - f))
ecc1sq = float(f * (2 - f))
ecc2sq = float(ecc1sq/(1 - ecc1sq))
ecc1 = sqrt(ecc1sq)
n = f / (2 - f)
n = float(n)
n2 = n ** 2


# Rectifying Radius (Horner Form)
A = ellipsoid.semimaj / (1 + n) * ((n2 *
                          (n2 *
                           (n2 *
                            (25 * n2 + 64)
                            + 256)
                           + 4096)
                          + 16384)
                         / 16384.)

# Alpha Coefficients (Horner Form)
a2 = ((n *
       (n *
        (n *
         (n *
          (n *
           (n *
            ((37884525 - 75900428 * n)
             * n + 42422016)
            - 89611200)
           + 46287360)
          + 63504000)
         - 135475200)
        + 101606400))
      / 203212800.)

a4 = ((n2 *
       (n *
        (n *
         (n *
          (n *
           (n *
            (148003883 * n + 83274912)
            - 178508970)
           + 77690880)
          + 67374720)
         - 104509440)
        + 47174400))
      / 174182400.)

a6 = ((n ** 3 *
       (n *
        (n *
         (n *
          (n *
           (318729724 * n - 738126169)
           + 294981280)
          + 178924680)
         - 234938880)
        + 81164160))
      / 319334400.)

a8 = ((n ** 4 *
       (n *
        (n *
         ((14967552000 - 40176129013 * n) * n + 6971354016)
         - 8165836800)
        + 2355138720))
      / 7664025600.)

a10 = ((n ** 5 *
        (n *
         (n *
          (10421654396 * n + 3997835751)
          - 4266773472)
         + 1072709352))
       / 2490808320.)

a12 = ((n ** 6 *
        (n *
         (175214326799 * n - 171950693600)
         + 38652967262))
       / 58118860800.)

a14 = ((n ** 7 *
        (13700311101 - 67039739596 * n))
       / 12454041600.)

a16 = (1424729850961 * n ** 8) / 743921418240.
a = (a2, a4, a6, a8, a10, a12, a14, a16)

# Beta Coefficients (Horner Form)
b2 = ((n *
       (n *
        (n *
         (n *
          (n *
           (n *
            ((37845269 - 31777436 * n) - 43097152)
            + 42865200)
           + 752640)
          - 104428800)
         + 180633600)
        - 135475200))
      / 270950400.)

b4 = ((n ** 2 *
       (n *
        (n *
         (n *
          (n *
           ((-24749483 * n - 14930208) * n + 100683990)
           - 152616960)
          + 105719040)
         - 23224320)
        - 7257600))
      / 348364800.)

b6 = ((n ** 3 *
       (n *
        (n *
         (n *
          (n *
           (232468668 * n - 101880889)
           - 39205760)
          + 29795040)
         + 28131840)
        - 22619520))
      / 638668800.)

b8 = ((n ** 4 *
       (n *
        (n *
         ((-324154477 * n - 1433121792) * n + 876745056)
         + 167270400)
        - 208945440))
      / 7664025600.)

b10 = ((n ** 5 *
        (n *
         (n *
          (312227409 - 457888660 * n)
          + 67920528)
         - 70779852))
       / 2490808320.)

b12 = ((n ** 6 *
        (n *
         (19841813847 * n + 3665348512)
         - 3758062126))
       / 116237721600.)

b14 = ((n ** 7 *
        (1989295244 * n - 1979471673))
       / 49816166400.)

b16 = ((-191773887257 * n ** 8) / 3719607091200.)
b = (b2, b4, b6, b8, b10, b12, b14, b16)


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
    :return:
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
    long_diff = radians(Decimal(long) - Decimal(str(cm)))

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

    return hemisphere, zone, round(float(east), 3), round(float(north), 3), psf


def grid2geo(zone, easting, northing):
    """
    input: Zone, Easting and Northing of a point in metres.
    (Default projection is Universal Transverse Mercator.)

    output: Latitude and Longitude in Decimal Degrees.
    """
    # Transverse Mercator Co-ordinates
    x = (easting - float(ellipsoid.falseeast)) / float(ellipsoid.cmscale)
    y = (northing - float(ellipsoid.falsenorth)) / float(ellipsoid.cmscale)
    # Transverse Mercator Ratios
    xi = y / A
    eta = x / A
    # Gauss-Schreiber Ratios
    xi2 = b2 * sin(2 * xi) * cosh(2 * eta)
    xi4 = b4 * sin(4 * xi) * cosh(4 * eta)
    xi6 = b6 * sin(6 * xi) * cosh(6 * eta)
    xi8 = b8 * sin(8 * xi) * cosh(8 * eta)
    xi10 = b10 * sin(10 * xi) * cosh(10 * eta)
    xi12 = b12 * sin(12 * xi) * cosh(12 * eta)
    xi14 = b14 * sin(14 * xi) * cosh(14 * eta)
    xi16 = b16 * sin(16 * xi) * cosh(16 * eta)
    eta2 = b2 * cos(2 * xi) * sinh(2 * eta)
    eta4 = b4 * cos(4 * xi) * sinh(4 * eta)
    eta6 = b6 * cos(6 * xi) * sinh(6 * eta)
    eta8 = b8 * cos(8 * xi) * sinh(8 * eta)
    eta10 = b10 * cos(10 * xi) * sinh(10 * eta)
    eta12 = b12 * cos(12 * xi) * sinh(12 * eta)
    eta14 = b14 * cos(14 * xi) * sinh(14 * eta)
    eta16 = b16 * cos(16 * xi) * sinh(16 * eta)
    xi1 = xi + xi2 + xi4 + xi6 + xi8 + xi10 + xi12 + xi14 + xi16
    eta1 = eta + eta2 + eta4 + eta6 + eta8 + eta10 + eta12 + eta14 + eta16
    # Conformal Latitude
    conf_lat = (sin(xi1)) / (sqrt((sinh(eta1)) ** 2 + (cos(xi1)) ** 2))
    t1 = conf_lat
    conf_lat = atan(conf_lat)

    # Finding t using Newtons Method
    def sigma(t):
        sigma = sinh(
            ecc1 * 0.5 * log((1 + ((ecc1 * t) / (sqrt(1 + t ** 2)))) / (1 - ((ecc1 * t) / (sqrt(1 + t ** 2))))))
        return sigma

    def ftn(t):
        ftn = t * sqrt(1 + (sigma(t)) ** 2) - sigma(t) * sqrt(1 + t ** 2) - t1
        return ftn

    def f1tn(t):
        f1tn = (sqrt(1 + (sigma(t)) ** 2) * sqrt(1 + t ** 2) - sigma(t) * t) * (
                ((1 - float(ecc1sq)) * sqrt(1 + t ** 2)) / (1 + (1 - float(ecc1sq)) * t ** 2))
        return f1tn

    t2 = t1 - (ftn(t1)) / (f1tn(t1))
    t3 = t2 - (ftn(t2)) / (f1tn(t2))
    t4 = t3 - (ftn(t3)) / (f1tn(t3))
    # Test No of Iterations Required (this will impact script performance)
    # t5 = t4 - (ftn(t4))/(f1tn(t4))
    # Compute Latitude
    lat = degrees(atan(t4))
    # Compute Longitude
    cm = float((zone * ellipsoid.zonewidth) + ellipsoid.initialcm - ellipsoid.zonewidth)
    long_diff = degrees(atan(sinh(eta1) / cos(xi1)))
    long = cm + long_diff
    return round(lat, 11), round(long, 11)


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
        nu = semi_maj/(sqrt(1 - ecc1sq * (sin(lat))**2))
        itercheck = lat - atan((z + nu * ecc1sq * sin(lat))/p)
        lat = atan((z + nu * ecc1sq * sin(lat))/p)
    nu = semi_maj/(sqrt(1 - ecc1sq * (sin(lat))**2))
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
        nu = semi_maj/(sqrt(1 - ecc1sq * (sin(lat)**2)))
    # Calculate x, y, z
    x = Decimal(str((nu + ellht) * cos(lat) * cos(long)))
    y = Decimal(str((nu + ellht) * cos(lat) * sin(long)))
    z = Decimal(str(((semi_min**2 / semi_maj**2) * nu + ellht) * sin(lat)))
    return x, y, z


"""
# conform7 Debug - Test Values ALIC
x = -4052051.7643
y = 4212836.2017
z = -2545106.0245
# Debug - Test Helmert Params (GDA94 to GDA2020, scale in ppm and rotations in seconds)
conform_gda94to20 = [0.06155, -0.01087, -0.04019, -0.009994, -0.0394924, -0.0327221, -0.0328979]
"""


def conform7(x, y, z, conform7_param):
    """
    input: x, y, z: 3D Cartesian Coordinate X, Y, Z in metres
           conform7_param: list of 7 Helmert Parameters [tx, ty, tz, sc, rx, ry, rz]
            tx, ty, tz: 3 Translations in metres
            sc: Scale factor in parts per million
            rx, ry, rz: 3 Rotations in decimal seconds
    return: xnew, ynew, znew: Transformed 3D Cartesian Coordinate X, Y, Z in metres
    """
    # Create XYZ Vector
    xyz_before = np.array([[x],
                           [y],
                           [z]])
    # Convert Units for Transformation Parameters
    scale = conform7_param[3] / 1000000
    rx = radians(dms2dd(conform7_param[4] / 10000))
    ry = radians(dms2dd(conform7_param[5] / 10000))
    rz = radians(dms2dd(conform7_param[6] / 10000))
    # Create Translation Vector
    translation = np.array([[conform7_param[0]],
                            [conform7_param[1]],
                            [conform7_param[2]]])
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
    # outfilewriter.writerow(['Pt', 'Latitude', 'Longitude'])
    for row in csvreader:
        pt_num = row[0]
        zone = float(row[1])
        east = float(row[2])
        north = float(row[3])
        # Calculate Conversion
        lat, long = grid2geo(zone, east, north)
        lat = dd2dms(lat)
        long = dd2dms(long)
        output = [pt_num, lat, long]
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
    # outfilewriter.writerow(['Pt', 'Zone', 'Easting', 'Northing'])
    for row in csvreader:
        pt_num = row[0]
        lat = dms2dd(float(row[1]))
        long = dms2dd(float(row[2]))
        # Calculate Conversion
        output = geo2grid(lat, long)
        output = [pt_num] + list(output)
        outfilewriter.writerow(output)
    # Close Files
    outfile.close()
    csvfile.close()

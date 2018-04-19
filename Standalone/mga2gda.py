# Author: Josh Batchelor <josh.batchelor@ga.gov.au>


from decimal import *
from math import sqrt, log, degrees, sin, cos, sinh, cosh, atan
import os
import argparse
import csv

# Universal Transverse Mercator Projection Parameters
proj = [6378137, Decimal('298.25722210088'), 500000,
        10000000, Decimal('0.9996'), 6, -177]

# Ellipsoidal Constants
f = 1 / proj[1]
semi_maj = proj[0]
semi_min = float(semi_maj * (1 - f))
ecc1sq = float(f * (2 - f))
ecc2sq = float(ecc1sq/(1 - ecc1sq))
ecc1 = sqrt(ecc1sq)
n = f / (2 - f)
n = float(n)
n2 = n ** 2

# Rectifying Radius (Horner Form)
A = proj[0] / (1 + n) * ((n2 *
                          (n2 *
                           (n2 *
                            (25 * n2 + 64)
                            + 256)
                           + 4096)
                          + 16384)
                         / 16384.)

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


def dd2dms(dd):
    minutes, seconds = divmod(abs(dd) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    dms = degrees + (minutes / 100) + (seconds / 10000)
    return dms if dd >= 0 else -dms


def grid2geo(zone, easting, northing):
    """
    input: Zone, Easting and Northing of a point in metres.
    (Default projection is Universal Transverse Mercator.)

    output: Latitude and Longitude in Decimal Degrees.
    """
    # Transverse Mercator Co-ordinates
    x = (easting - float(proj[2])) / float(proj[4])
    y = (northing - float(proj[3])) / float(proj[4])
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
    cm = float((zone * proj[5]) + proj[6] - proj[5])
    long_diff = degrees(atan(sinh(eta1) / cos(xi1)))
    long = cm + long_diff
    return round(lat, 11), round(long, 11)


def grid2geoio(fn):
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
    # Open Filename
    csvfile = open(fn)
    csvreader = csv.reader(csvfile)
    # Create Output File
    fn_part = (os.path.splitext(fn))
    fn_out = fn_part[0] + '_out' + fn_part[1]
    outfile = open(fn_out, 'w', newline='')
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
    return 'Output saved as ' + str(fn_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Batch Converter of Map Grid of Australia grid co-ordinates to '
                                                 'Geodetic Datum of Australia geographic co-ordinates. Files must '
                                                 'be .csv and of the format Pt ID, Zone, Easting, Northing with no '
                                                 'header line.')
    parser.add_argument('f', metavar='--file', type=str, help='Input Filename')
    args = parser.parse_args()
    fn = args.f

    grid2geoio(fn)

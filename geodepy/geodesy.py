#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Geodesy Module
"""

from math import degrees, radians, sqrt, sin, cos, tan, asin, acos, atan, atan2
import numpy as np
from geodepy.constants import grs80, utm
from geodepy.convert import geo2grid, grid2geo, angular_typecheck
from geodepy.statistics import rotation_matrix
from geodepy.survey import radiations


def enu2xyz(lat, lon, east, north, up):
    """Convert a column vector in the local reference frame to a column vector
    in the Cartesian reference frame.
    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :param east: in metres
    :param north: in metres
    :param up: in metres
    :return: x, y, z in metres
    """
    rot_matrix = rotation_matrix(angular_typecheck(lat), angular_typecheck(lon))
    enu = np.array([[east], [north], [up]])
    xyz = rot_matrix @ enu
    x = xyz[0, 0]
    y = xyz[1, 0]
    z = xyz[2, 0]
    return x, y, z


def xyz2enu(lat, lon, x, y, z):
    """Convert a column vector in the Cartesian reference frame to a column
    vector in the local reference frame.
    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :param x: in metres
    :param y: in metres
    :param z: in metres
    :return: east, north, up in metres
    """
    rot_matrix = rotation_matrix(angular_typecheck(lat), angular_typecheck(lon))
    xyz = np.array([[x], [y], [z]])
    enu = rot_matrix.transpose() @ xyz
    east = enu[0, 0]
    north = enu[1, 0]
    up = enu[2, 0]
    return east, north, up


def vincdir(lat1, lon1, azimuth1to2, ell_dist, ellipsoid=grs80):
    """
    Vincenty's Direct Formula
    :param lat1: Latitude of Point 1 (decimal degrees)
    :type lat1: float (decimal degrees), DMSAngle or DDMAngle
    :param lon1: Longitude of Point 1 (decimal degrees)
    :type lon1: float (decimal degrees), DMSAngle or DDMAngle
    :param azimuth1to2: Azimuth from Point 1 to 2 (decimal degrees)
    :type azimuth1to2: float (decimal degrees), DMSAngle or DDMAngle
    :param ell_dist: Ellipsoidal Distance between Points 1 and 2 (metres)
    :param ellipsoid: Ellipsoid Object
    :return: lat2: Latitude of Point 2 (Decimal Degrees),
             lon2: Longitude of Point 2 (Decimal Degrees),
             azimuth2to1: Azimuth from Point 2 to 1 (Decimal Degrees)

    Code review: 14-08-2018 Craig Harrison
    """

    # Convert Angles to Decimal Degrees (if required)
    lat1 = angular_typecheck(lat1)
    lon1 = angular_typecheck(lon1)
    azimuth1to2 = radians(angular_typecheck(azimuth1to2))

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
    :param lat1: Latitude of Point 1 (decimal degrees)
    :type lat1: float (decimal degrees), DMSAngle or DDMAngle
    :param lon1: Longitude of Point 1 (decimal degrees)
    :type lon1: float (decimal degrees), DMSAngle or DDMAngle
    :param lat2: Latitude of Point 2 (decimal degrees)
    :type lat2: float (decimal degrees), DMSAngle or DDMAngle
    :param lon2: Longitude of Point 2 (decimal degrees)
    :type lon2: float (decimal degrees), DMSAngle or DDMAngle
    :param ellipsoid: Ellipsoid Object
    :return: ell_dist: Ellipsoidal Distance between Points 1 and 2 (m),
             azimuth1to2: Azimuth from Point 1 to 2 (Decimal Degrees),
             azimuth2to1: Azimuth from Point 2 to 1 (Decimal Degrees)
    Code review: 14-08-2018 Craig Harrison
    """

    # Convert Angles to Decimal Degrees (if required)
    lat1 = angular_typecheck(lat1)
    lon1 = angular_typecheck(lon1)
    lat2 = angular_typecheck(lat2)
    lon2 = angular_typecheck(lon2)

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
    cos_two_sigma_m = 0
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
        cos_two_sigma_m = cos(sigma) - (2*sin(u1)*sin(u2) / cos(alpha)**2)

        # Eq. 79
        c = (ellipsoid.f / 16) * cos(alpha)**2 * (4 + ellipsoid.f
                                                  * (4 - 3*cos(alpha)**2))

        # Eq. 80
        new_lon = omega + (1 - c) * ellipsoid.f * sin(alpha) * (sigma + c*sin(sigma)
                                                                * (cos_two_sigma_m + c * cos(sigma)
                                                                   * (-1 + 2*(cos_two_sigma_m**2))))

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
    delta_sigma = b*sin(sigma) * (cos_two_sigma_m + (b / 4)
                                  * (cos(sigma) * (-1 + 2*cos_two_sigma_m**2)
                                  - (b / 6)*cos_two_sigma_m
                                  * (-3 + 4*sin(sigma)**2)
                                  * (-3 + 4*cos_two_sigma_m**2)))

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


def vincdir_utm(zone1, east1, north1, grid1to2, grid_dist,
                hemisphere='south', ellipsoid=grs80):
    """
    Perform Vincenty's Direct Computation using UTM Grid Coordinates, a
    grid bearing and grid distance.
    Note: Point 2 UTM Coordinates use the Zone specified for Point 1, even if
    Point 2 would typically be computed in a different zone. This keeps the grid
    bearings and line scale factor all relative to the same UTM Zone.
    :param zone1: Point 1 Zone Number - 1 to 60
    :param east1: Point 1 Easting (m, within 3330km of Central Meridian)
    :param north1: Point 1 Northing (m, 0 to 10,000,000m)
    :param grid1to2: Grid Bearing from Point 1 to 2 (decimal degrees),
    :param grid_dist: UTM Grid Distance between Points 1 and 2 (m)
    :param hemisphere: String - 'North' or 'South'(default)
    :param ellipsoid: Ellipsoid Object (default: GRS80)
    :return: zone2: Point 2 Zone Number - 1 to 60
             east2: Point 2 Easting (m, within 3330km of Central Meridian)
             north2: Point 2 Northing (m, 0 to 10,000,000m)
             grid2to1: Grid Bearing from Point 2 to 1 (decimal degrees)
             lsf: Line Scale Factor (for Point 1 Zone)
    """
    # Convert angle to decimal degrees (if required)
    grid1to2 = angular_typecheck(grid1to2)

    # Convert UTM Coords to Geographic
    lat1, lon1, psf1, gridconv1 = grid2geo(zone1, east1, north1,
                                           hemisphere, ellipsoid)

    # Convert Grid Bearing to Geodetic Azimuth
    az1to2 = grid1to2 - gridconv1

    # Estimate Line Scale Factor (LSF)
    zone2, east2, north2 = (zone1, *radiations(east1, north1,
                                               grid1to2, grid_dist))
    lsf = line_sf(zone1, east1, north1, zone2, east2, north2)

    # Iteratively estimate Pt 2 Coordinates, refining LSF each time
    lsf_diff = 1
    max_iter = 10
    while lsf_diff > 1e-9:
        lsf_previous = lsf
        lat2, lon2, az2to1 = vincdir(lat1, lon1, az1to2,
                                     grid_dist / lsf, ellipsoid)
        (hemisphere2, zone2, east2,
         north2, psf2, gridconv2) = geo2grid(lat2, lon2,
                                             zone1, ellipsoid)
        lsf = line_sf(zone1, east1, north1,
                      zone2, east2, north2,
                      hemisphere, ellipsoid)
        lsf_diff = abs(lsf_previous - lsf)

    # Compute Grid Bearing for Station 2
    lat2, lon2, psf2, gridconv2 = grid2geo(zone2, east2, north2,
                                           hemisphere, ellipsoid)
    grid2to1 = az2to1 + gridconv2

    return zone2, east2, north2, grid2to1, lsf


def vincinv_utm(zone1, east1, north1, zone2, east2, north2,
                hemisphere='south', ellipsoid=grs80):
    """
    Perform Vincentys Inverse Computation using UTM Grid Coordinates
    Note: Where coordinates from different zones are used, UTM Grid Distance
    is relative to Point 1's Zone. Grid Bearings of Points 1 and 2 are
    relative to each of their respective Zones.
    :param zone1: Point 1 Zone Number - 1 to 60
    :param east1: Point 1 Easting (m, within 3330km of Central Meridian)
    :param north1: Point 1 Northing (m, 0 to 10,000,000m)
    :param zone2: Point 2 Zone Number - 1 to 60
    :param east2: Point 2 Easting (m, within 3330km of Central Meridian)
    :param north2: Point 2 Northing (m, 0 to 10,000,000m)
    :param hemisphere: String - 'North' or 'South'(default)
    :param ellipsoid: Ellipsoid Object (default: GRS80)
    :return: grid_dist: UTM Grid Distance between Points 1 and 2 (m),
             grid1to2: Grid Bearing from Point 1 to 2 (decimal degrees),
             grid2to1: Grid Bearing from Point 2 to 1 (decimal degrees)
             lsf: Line Scale Factor (for Point 1 Zone)
    """
    # Convert utm to geographic
    pt1 = grid2geo(zone1, east1, north1, hemisphere, ellipsoid)
    pt2 = grid2geo(zone2, east2, north2, hemisphere, ellipsoid)
    ell_dist, az1to2, az2to1 = vincinv(pt1[0], pt1[1],
                                       pt2[0], pt2[1], ellipsoid)

    # Compute Grid Distance using Line Scale Factor
    lsf = line_sf(zone1, east1, north1,
                  zone2, east2, north2, hemisphere, ellipsoid)
    grid_dist = ell_dist * lsf

    # Compute Grid Bearings using Grid Convergence
    grid1to2 = az1to2 + pt1[3]
    grid2to1 = az2to1 + pt2[3]

    return grid_dist, grid1to2, grid2to1, lsf


def line_sf(zone1, east1, north1, zone2, east2, north2,
            hemisphere='south', ellipsoid=grs80, projection=utm):
    """
    Computes Line Scale Factor for a pair of Transverse Mercator Coordinates
    Ref: Deakin 2010, Traverse Computations on the Ellipsoid and on the
    Universal Transverse Mercator Projection, pp 35
    http://www.mygeodesy.id.au/documents/Trav_Comp_V2.1.pdf
    :param zone1: Station 1 Zone Number - 1 to 60
    :param east1: Station 1 Easting (m, within 3330km of Central Meridian)
    :param north1: Station 1 Northing (m, 0 to 10,000,000m)
    :param zone2: Station 2 Zone Number - 1 to 60
    :param east2: Station 2 Easting (m, within 3330km of Central Meridian)
    :param north2: Station 2 Northing (m, 0 to 10,000,000m)
    :param hemisphere: String - 'North' or 'South'(default)
    :param ellipsoid: Ellipsoid Object (default GRS80)
    :param projection: Projection Object (default Universal Transverse Mercator)
    :return: Line Scale Factor (relative to Zone 1 if different zones
    are entered)
    """
    # Re-project cross-zone coordinate to same UTM zone
    if zone1 != zone2:
        # Re-project Station 2 Coordinates to same zone as Station 1
        stn2_geo = grid2geo(zone2, east2, north2, hemisphere, ellipsoid)
        stn2_zone1 = geo2grid(stn2_geo[0], stn2_geo[1], zone1, ellipsoid)
        zone2 = stn2_zone1[1]
        east2 = stn2_zone1[2]
        north2 = stn2_zone1[3]

    # Comute easting distances from Central Meridian
    eastofcm1 = east1 - projection.falseeast
    eastofcm2 = east2 - projection.falseeast

    # Compute Mean Latitude
    lat1 = grid2geo(zone1, east1, north1, hemisphere, ellipsoid)[0]
    lat2 = grid2geo(zone2, east2, north2, hemisphere, ellipsoid)[0]
    lat_mean = (lat1 + lat2) / 2

    # Compute Line Scale Factor (Deakin 2010 Eq. 13)
    r_sq_m = (rho(lat_mean, ellipsoid) *
              nu(lat_mean, ellipsoid) *
              projection.cmscale ** 2)
    k1 = ((eastofcm1 ** 2 + eastofcm1 * eastofcm2 + eastofcm2 ** 2) /
          (6 * r_sq_m))
    k2 = ((eastofcm1 ** 2 + eastofcm1 * eastofcm2 + eastofcm2 ** 2) /
          (36 * r_sq_m))
    return projection.cmscale * (1 + k1 * (1 + k2))


def rho(lat, ellipsoid=grs80):
    """
    Return the radius of curvature of the ellipsoid in the meridian plane
    (rho) at a given latitude
    :param lat: latitude in decimal degrees
    :param ellipsoid: Ellipsoid Object
    :return: rho at specified latitude
    """
    return ((ellipsoid.semimaj * (1 - ellipsoid.ecc1sq)) /
            (1 - ellipsoid.ecc1sq * (sin(radians(lat)) ** 2)) ** 1.5)


def nu(lat, ellipsoid=grs80):
    """
    Return the radius of curvature of the ellipsoid in the prime vertical plane
    (nu) at a given latitude
    :param lat: latitude in decimal degrees
    :param ellipsoid: Ellipsoid Object
    :return: nu at specified latitude
    """
    return (ellipsoid.semimaj /
            sqrt(1 - ellipsoid.ecc1sq * (sin(radians(lat)) ** 2)))

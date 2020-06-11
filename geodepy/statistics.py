#!/usr/bin/env python3

from math import radians, sin, cos, sqrt, atan2, degrees
import numpy as np


def rotation_matrix(lat, lon):
    """Returns the rotation matrix between the local and Cartesian reference
    systems at a given latitude and longitude

    See Section 4.2.3 of the DynAdjust User's Guide v1.0

    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :return: rot_matrix (3x3)
    """
    (rlat, rlon) = (radians(lat), radians(lon))
    rot_matrix = np.array(
        [[-sin(rlon), -sin(rlat) * cos(rlon), cos(rlat) * cos(rlon)],
         [cos(rlon), -sin(rlat) * sin(rlon), cos(rlat) * sin(rlon)],
         [0.0, cos(rlat), sin(rlat)]]
    )
    return rot_matrix


def vcv_cart2local(vcv_cart, lat, lon):
    """Transforms a VCV from the Cartesian to the local reference frame. If
    only a column vector of variances is supplied then they are assumed to be
    uncorrelated and only the column vector of variances is returned.

    See Section 4.4.1 of the DynAdjust User's Guide v1.0

    :param vcv_cart: a VCV in the Cartesian reference frame (3x1 or 3x3)
    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :return: vcv_local (3x1 or 3x3)
    """
    column_vector = False
    if vcv_cart.shape[0] == 3:
        if vcv_cart.shape[1] == 1:
            column_vector = True
            vcv_cart = np.array([[vcv_cart[0, 0], 0.0, 0.0],
                                 [0.0, vcv_cart[1, 0], 0.0],
                                 [0.0, 0.0, vcv_cart[2, 0]]])
        elif vcv_cart.shape[1] == 3:
            pass
        else:
            raise ValueError('Matrix must be either 3x1 or 3x3')
    else:
        raise ValueError('Matrix must be either 3x1 or 3x3')
    rot_matrix = rotation_matrix(lat, lon)
    vcv_local = rot_matrix.transpose() @ vcv_cart @ rot_matrix
    if column_vector:
        vcv_local = np.array([[vcv_local[0, 0]], [vcv_local[1, 1]],
                              [vcv_local[2, 2]]])
    return vcv_local


def vcv_local2cart(vcv_local, lat, lon):
    """Transform a VCV from the local to the Cartesian reference frame. If
    only a column vector of variances is supplied then they are assumed to be
    uncorrelated and only the column vector of variances is returned.

    See Section 4.4.1 of the DynAdjust User's Guide v1.0

    :param vcv_local: a VCV in the local reference frame (3x1 or 3x3)
    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :return: vcv_cart (3x1 or 3x3)
    """
    column_vector = False
    if vcv_local.shape[0] == 3:
        if vcv_local.shape[1] == 1:
            column_vector = True
            vcv_local = np.array([[vcv_local[0, 0], 0.0, 0.0],
                                 [0.0, vcv_local[1, 0], 0.0],
                                 [0.0, 0.0, vcv_local[2, 0]]])
        elif vcv_local.shape[1] == 3:
            pass
        else:
            raise ValueError('Matrix must be either 3x1 or 3x3')
    else:
        raise ValueError('Matrix must be either 3x1 or 3x3')
    rot_matrix = rotation_matrix(lat, lon)
    vcv_cart = rot_matrix @ vcv_local @ rot_matrix.transpose()
    if column_vector:
        vcv_cart = np.array([[vcv_cart[0, 0]], [vcv_cart[1, 1]],
                             [vcv_cart[2, 2]]])
    return vcv_cart


def error_ellipse(vcv):
    """Calculate the semi-major axis, semi-minor axis, and the orientation of
    the error ellipse defined by a VCV

    See Section 7.3.3.1 of the DynaNet User's Guide v3.3

    :param vcv: a VCV (3x3)
    :return: a, semi-major axis
    :return: b, semi-minor axis
    :return: orientation, the orientation of the error ellipse
    """
    z = sqrt((vcv[0, 0] - vcv[1, 1])**2 + 4 * vcv[0, 1]**2)
    a = sqrt(0.5 * (vcv[0, 0] + vcv[1, 1] + z))
    b = sqrt(0.5 * (vcv[0, 0] + vcv[1, 1] - z))
    orientation = 90 - degrees(0.5 * atan2((2 * vcv[0, 1]),
                                           (vcv[0, 0] - vcv[1, 1])))

    return a, b, orientation


def relative_error(lat, lon, var1, var2, cov12):
    """
    Function to compute relative error between two 3D stations:
        - 2D relative error ellipse [semi-major axis, semi-minor axis, bearing]
        - 1D relative 'up' error
        
    Adapted from Harvey B.R. (1998) Practical least squares and statistics for surveyors,
    Monograph 13, Section 4.4, p.135

    :param lat: latitude at Stn1 in decimal degrees
    :param lon: longitude at Stn1 in decimal degrees
    :param var1: 3x3 Cartesian XYZ variance matrix of Stn1
    :param var2: 3x3 Cartesian XYZ variance matrix of Stn2
    :param cov12: 3x3 Cartesian XYZ covariance block between Stn1 and Stn2
    :return: Relative error ellipse components [smaj, smin, brg] and relative 'up' error
    """

    # rotate matrices from cartesian to local
    var1_enu = vcv_cart2local(var1, lat, lon)
    var2_enu = vcv_cart2local(var2, lat, lon)
    cov12_enu = vcv_cart2local(cov12, lat, lon)

    # form relative variance matrix
    rel_var = np.zeros((3, 3))

    # diagonal terms
    rel_var[0, 0] = var1_enu[0, 0] + var2_enu[0, 0] - 2 * cov12_enu[0, 0]
    rel_var[1, 1] = var1_enu[1, 1] + var2_enu[1, 1] - 2 * cov12_enu[1, 1] 
    rel_var[2, 2] = var1_enu[2, 2] + var2_enu[2, 2] - 2 * cov12_enu[2, 2]

    # east-north covariance
    rel_var[0, 1] = var1_enu[0, 1] + var2_enu[0, 1] - cov12_enu[0, 1] - cov12_enu[1, 0]
    rel_var[1, 0] = rel_var[0, 1]

    # east-up covariance
    rel_var[0, 2] = var1_enu[0, 2] + var2_enu[0, 2] - cov12_enu[0, 2] - cov12_enu[2, 0]
    rel_var[2, 0] = rel_var[0, 2]

    # north-up covariance
    rel_var[1, 2] = var1_enu[1, 2] + var2_enu[1, 2] - cov12_enu[1, 2] - cov12_enu[2, 1]
    rel_var[2, 1] = rel_var[1, 2]

    # relative error ellipse [smaj, smin, brg]
    ree = error_ellipse(rel_var)

    # relative up error
    rue = rel_var[2, 2] ** 0.5

    return ree[0], ree[1], ree[2], rue


def circ_hz_pu(a, b):
    """Calculate the circularised horizontal PU from the semi-major and
    semi-minor axes of an error ellipse

    :param a: semi-major axis
    :param b: semi-minor axis
    :return: r the radius of the circularised error
    """
    q0 = 1.960790
    q1 = 0.004071
    q2 = 0.114276
    q3 = 0.371625
    c = b / a
    k = q0 + q1 * c + q2 * c**2 + q3 * c**3
    r = a * k

    return r


def k_val95(dof):
    """
    Returns the Coverage Factor k for a given 1 sigma (68.27%) Standard
    Deviation to allow conversion to a 95% Standard Deviation. This
    uses a simplified table of k values rounded to 5 decimal places
    for Degrees of Freedom (DOF) in the range 1 to 120. For DOF above
    120, returns k value of 1.96 and for DOF below 1, returns k value
    for DOF = 1. Coverage Factor produced using following scipy stats
    code: stats.t.ppf(1-0.025,dof)

    :param dof: Degrees of Freedom (Number of Measurements Minus 1)
    :return: Coverage Factor k
    """
    if not isinstance(dof, int):
        raise TypeError('Degrees of Freedom must be Int')
    if dof < 1:
        return ttable_p95[0]
    elif dof > 120:
        return 1.96
    else:
        return ttable_p95[dof - 1]


ttable_p95 = ([12.7062, 4.30265, 3.18245, 2.77645, 2.57058, 2.44691,
              2.36462, 2.30600, 2.26216, 2.22814, 2.20099, 2.17881,
              2.16037, 2.14479, 2.13145, 2.11991, 2.10982, 2.10092,
              2.09302, 2.08596, 2.07961, 2.07387, 2.06866, 2.06390,
              2.05954, 2.05553, 2.05183, 2.04841, 2.04523, 2.04227,
              2.03951, 2.03693, 2.03452, 2.03224, 2.03011, 2.02809,
              2.02619, 2.02439, 2.02269, 2.02108, 2.01954, 2.01808,
              2.01669, 2.01537, 2.01410, 2.01290, 2.01174, 2.01063,
              2.00958, 2.00856, 2.00758, 2.00665, 2.00575, 2.00488,
              2.00404, 2.00324, 2.00247, 2.00172, 2.00100, 2.00030,
              1.99962, 1.99897, 1.99834, 1.99773, 1.99714, 1.99656,
              1.99601, 1.99547, 1.99495, 1.99444, 1.99394, 1.99346,
              1.99300, 1.99254, 1.99210, 1.99167, 1.99125, 1.99085,
              1.99045, 1.99006, 1.98969, 1.98932, 1.98896, 1.98861,
              1.98827, 1.98793, 1.98761, 1.98729, 1.98698, 1.98667,
              1.98638, 1.98609, 1.98580, 1.98552, 1.98525, 1.98498,
              1.98472, 1.98447, 1.98422, 1.98397, 1.98373, 1.98350,
              1.98326, 1.98304, 1.98282, 1.98260, 1.98238, 1.98217,
              1.98197, 1.98177, 1.98157, 1.98137, 1.98118, 1.98099,
              1.98081, 1.98063, 1.98045, 1.98027, 1.98010, 1.97993])

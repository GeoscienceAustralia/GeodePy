#!/usr/bin/env python3

from math import radians, sin, cos, sqrt, atan2, degrees
import numpy as np


def rotation_matrix(lat, lon):
    """Returns the 3x3 rotation matrix for a given latitude and longitude
    (given in decimal degrees)
    See Section 4.2.3 of the DynaNet User's Guide v3.3
    """
    (rlat, rlon) = (radians(lat), radians(lon))
    rot_matrix = np.array(
        [[-sin(rlon), -sin(rlat) * cos(rlon), cos(rlat) * cos(rlon)],
         [cos(rlon), -sin(rlat) * sin(rlon), cos(rlat) * sin(rlon)],
         [0.0, cos(rlat), sin(rlat)]]
    )
    return rot_matrix


def vcv_cart2local(vcv_cart, lat, lon):
    """Transform a 3x3 VCV from the Cartesian to the local reference frame. If
    only a column vector of variances is supplied (3x1) then the full VCV is
    padded out with zeros for the transformation. In these cases, only the
    column vector of variances (3x1) is returned.

    See Section 4.4.1 of the DynaNet User's Guide v3.3
    """
    if vcv_cart.shape[0] == 3:
        if vcv_cart.shape[1] == 1:
            vcv_cart = np.array([[vcv_cart[0,0], 0.0, 0.0],
                                 [0.0, vcv_cart[1,0], 0.0],
                                 [0.0, 0.0, vcv_cart[2,0]]])
        elif vcv_cart.shape[1] == 3:
            pass
        else:
            sys.exit('Matrix must be either 3x1 or 3x3')
    else:
         sys.exit('Matrix must be either 3x1 or 3x3')
    rot_matrix = rotation_matrix(lat, lon)
    vcv_local = rot_matrix.transpose() @ vcv_cart @ rot_matrix
    return vcv_local


def vcv_local2cart(vcv_local, lat, lon):
    """Transform a 3x3 VCV from the local to the Cartesian reference frame. If
    only a column vector of variances is supplied (3x1) then the full VCV is
    padded out with zeros for the transformation. In these cases, only the
    column vector of variances (3x1) is returned.

    See Section 4.4.1 of the DynaNet User's Guide v3.3
    """
    if vcv_local.shape[0] == 3:
        if vcv_local.shape[1] == 1:
            vcv_local = np.array([[vcv_local[0, 0], 0.0, 0.0],
                                 [0.0, vcv_local[1, 0], 0.0],
                                 [0.0, 0.0, vcv_local[2, 0]]])
        elif vcv_local.shape[1] == 3:
            pass
        else:
            sys.exit('Matrix must be either 3x1 or 3x3')
    else:
        sys.exit('Matrix must be either 3x1 or 3x3')
    rot_matrix = rotation_matrix(lat, lon)
    vcv_cart = rot_matrix @ vcv_local @ rot_matrix.transpose()
    return vcv_cart


def error_ellipse(vcv):
    """Calculate the semi-major axis, semi-minor axis, and the orientation of
    the error ellipse calculated from a 3x3 VCV
    See Section 7.3.3.1 of the DynaNet User's Guide v3.3
    """
    z = sqrt((vcv[0, 0] - vcv[1, 1])**2 + 4 * vcv[0, 1]**2)
    a = sqrt(0.5 * (vcv[0, 0] + vcv[1, 1] + z))
    b = sqrt(0.5 * (vcv[0, 0] + vcv[1, 1] - z))
    orientation = 90 - degrees(0.5 * atan2((2 * vcv[0, 1]),
                                           (vcv[0, 0] - vcv[1, 1])))

    return a, b, orientation


def circ_hz_pu(a, b):
    """Calculate the circularised horizontal PU form the semi-major and
    semi-minor axes
    """
    q0 = 1.960790
    q1 = 0.004071
    q2 = 0.114276
    q3 = 0.371625
    c = b / a
    k = q0 + q1 * c + q2 * c**2 + q3 * c**3
    r = a * k

    return r

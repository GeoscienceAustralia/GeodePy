import numpy as np
from math import radians
from conversions import dms2dd


# copied from transform.py


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


def conform14(x, y, z, conform14_array, datefrom='refepoch', dateto='refepoch'):
    """

    :param x:
    :param y:
    :param z:
    :param conform14_array:
    :param datefrom:
    :param dateto:
    :return: xtrans, ytrans, ztrans [, conform7_param]
    """
    # Test for changes in epoch
    if datefrom != 'refepoch':
        return True
        # do stuff
    if dateto != 'refepoch':
        return True
        # do stuff
    # apply datefrom and dateto to conform14_array rates
    # generate final set of transformation parameters
    # apply these to x, y, z using conform7() method
    return x, y, z, conform14_array, datefrom, dateto

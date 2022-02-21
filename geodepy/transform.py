#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Transform Module

Ref1
http://www.icsm.gov.au/sites/default/files/GDA2020TechnicalManualV1.1.1.pdf

Ref2
http://www.mygeodesy.id.au/documents/Karney-Krueger%20equations.pdf
"""

import datetime
from math import radians
import numpy as np
from geodepy.constants import (Transformation, TransformationSD,
                               atrf2014_to_gda2020, gda94_to_gda2020)
from geodepy.statistics import vcv_local2cart, vcv_cart2local
from geodepy.convert import (hp2dec, geo2grid,
                             grid2geo, xyz2llh, llh2xyz)
from geodepy.ntv2reader import NTv2Grid, interpolate_ntv2


def conform7(x, y, z, trans, vcv=None):
    """
    Performs a Helmert 7 Parameter Conformal Transformation using Cartesian point co-ordinates
    and a predefined transformation object.
    :param x: Cartesian X (m)
    :param y: Cartesian Y (m)
    :param z: Cartesian Z (m)
    :param trans: Transformation Object (note: this function ignores all time-dependent variables)
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :return: Transformed X, Y, Z Cartesian Co-ordinates, vcv matrix
    """
    if type(trans) != Transformation:
        raise ValueError('trans must be a Transformation Object')
    # Create XYZ Vector
    xyz_before = np.array([[x],
                           [y],
                           [z]])
    # Convert Units for Transformation Parameters
    scale = 1 + trans.sc / 1000000
    rx = radians(hp2dec(trans.rx / 10000))
    ry = radians(hp2dec(trans.ry / 10000))
    rz = radians(hp2dec(trans.rz / 10000))
    # Create Translation Vector
    translation = np.array([[trans.tx],
                            [trans.ty],
                            [trans.tz]])
    # Create Rotation Matrix
    rotation = np.array([[1., rz, -ry],
                         [-rz, 1., rx],
                         [ry, -rx, 1.]])

    rot_xyz = rotation @ xyz_before

    # Conformal Transform Eq
    xyz_after = translation + scale * rot_xyz
    # Convert Vector to Separate Variables
    xtrans = float(xyz_after[0])
    ytrans = float(xyz_after[1])
    ztrans = float(xyz_after[2])

    # Transformation uncertainty propagation
    # Adapted from Harvey B.R. (1998) Practical least squares and statistics for surveyors,
    # Monograph 13 Section 8.7.2, p.274
    if (type(trans.tf_sd) == TransformationSD) and (vcv is not None):
        # Q matrix:
        q_mat = np.zeros((10, 10))
        # xyz_before vcv
        for i in range(3):
            for j in range(3):
                q_mat[i, j] = vcv[i, j]

        # transformation variances
        q_mat[3, 3] = (trans.tf_sd.sd_sc / 1000000)**2
        q_mat[4, 4] = radians(trans.tf_sd.sd_rx/3600)**2
        q_mat[5, 5] = radians(trans.tf_sd.sd_ry/3600)**2
        q_mat[6, 6] = radians(trans.tf_sd.sd_rz/3600)**2
        q_mat[7, 7] = trans.tf_sd.sd_tx**2
        q_mat[8, 8] = trans.tf_sd.sd_ty**2
        q_mat[9, 9] = trans.tf_sd.sd_tz**2

        # Jacobian matrix:
        j_mat = np.zeros((3, 10))

        # scaled rotations
        j_mat[0, 0] = scale
        j_mat[0, 1] = scale * rz
        j_mat[0, 2] = -scale * ry
        j_mat[1, 0] = -j_mat[0, 1]
        j_mat[1, 1] = scale
        j_mat[1, 2] = scale * rx
        j_mat[2, 0] = -j_mat[0, 2]
        j_mat[2, 1] = -j_mat[1, 2]
        j_mat[2, 2] = scale

        # XYZ rotated
        j_mat[0, 3] = rot_xyz[0]
        j_mat[1, 3] = rot_xyz[1]
        j_mat[2, 3] = rot_xyz[2]

        # scaled XYZ
        j_mat[0, 5] = -scale * xyz_before[2]
        j_mat[0, 6] = scale * xyz_before[1]
        j_mat[1, 4] = scale * xyz_before[2]
        j_mat[1, 6] = -scale * xyz_before[0]
        j_mat[2, 4] = -scale * xyz_before[1]
        j_mat[2, 5] = scale * xyz_before[0]

        # Identity
        j_mat[0, 7] = 1
        j_mat[1, 8] = 1
        j_mat[2, 9] = 1

        # multiply J Q J_trans
        vcv_after = j_mat @ q_mat @ j_mat.transpose()

        return xtrans, ytrans, ztrans, vcv_after

    else:
        return xtrans, ytrans, ztrans, None


def conform14(x, y, z, to_epoch, trans, vcv=None):
    """
    Performs a Helmert 14 Parameter Conformal Transformation using Cartesian point co-ordinates
    and a predefined transformation object. The transformation parameters are projected from
    the transformation objects reference epoch to a specified epoch.
    :param x: Cartesian X (m)
    :param y: Cartesian Y (m)
    :param z: Cartesian Z (m)
    :param to_epoch: Epoch co-ordinate transformation is performed at (datetime.date Object)
    :param trans: Transformation Object
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :return: Cartesian X, Y, Z co-ordinates and vcv matrix transformed using Transformation parameters at desired epoch
    """
    if type(trans) != Transformation:
        raise ValueError('trans must be a Transformation Object')
    if type(to_epoch) != datetime.date:
        raise ValueError('to_epoch must be a datetime.date Object')
    # Calculate 7 Parameters from 14 Parameter Transformation Object
    timetrans = trans + to_epoch

    # Perform Transformation
    xtrans, ytrans, ztrans, trans_vcv = conform7(x, y, z, timetrans, vcv=vcv)
    return xtrans, ytrans, ztrans, trans_vcv


def transform_mga94_to_mga2020(zone, east, north, ell_ht=False, vcv=None):
    """
    Performs conformal transformation of Map Grid of Australia 1994 to Map Grid of Australia 2020 Coordinates
    using the GDA2020 Tech Manual v1.2 7 parameter similarity transformation parameters
    :param zone: Zone Number - 1 to 60
    :param east: Easting (m, within 3330km of Central Meridian)
    :param north: Northing (m, 0 to 10,000,000m)
    :param ell_ht: Ellipsoid Height (m) (optional)
    :param vcv: Optional 3*3 numpy array in local enu units to propagate tf uncertainty
    :return: MGA2020 Zone, Easting, Northing, Ellipsoid Height (if none provided, returns 0), and vcv matrix
    """
    lat, lon, psf, gridconv = grid2geo(zone, east, north)
    if ell_ht is False:
        ell_ht_in = 0
    else:
        ell_ht_in = ell_ht
    if vcv is not None:
        vcv = vcv_local2cart(vcv, lat, lon)
    x94, y94, z94 = llh2xyz(lat, lon, ell_ht_in)
    x20, y20, z20, vcv20 = conform7(x94, y94, z94, gda94_to_gda2020, vcv)
    lat, lon, ell_ht_out = xyz2llh(x20, y20, z20)
    if vcv20 is not None:
        vcv20 = vcv_cart2local(vcv20, lat, lon)
    if ell_ht is False:
        ell_ht_out = 0
    hemisphere, zone20, east20, north20, psf, gridconv = geo2grid(lat, lon)
    return zone20, east20, north20, round(ell_ht_out, 4), vcv20


def transform_mga2020_to_mga94(zone, east, north, ell_ht=False, vcv=None):
    """
    Performs conformal transformation of Map Grid of Australia 2020 to Map Grid of Australia 1994 Coordinates
    using the reverse form of the GDA2020 Tech Manual v1.2 7 parameter similarity transformation parameters
    :param zone: Zone Number - 1 to 60
    :param east: Easting (m, within 3330km of Central Meridian)
    :param north: Northing (m, 0 to 10,000,000m)
    :param ell_ht: Ellipsoid Height (m) (optional)
    :param vcv: Optional 3*3 numpy array in local enu units to propagate tf uncertainty
    :return: MGA1994 Zone, Easting, Northing, Ellipsoid Height (if none provided, returns 0), and vcv matrix
    """
    lat, lon, psf, gridconv = grid2geo(zone, east, north)
    if ell_ht is False:
        ell_ht_in = 0
    else:
        ell_ht_in = ell_ht
    if vcv is not None:
        vcv = vcv_local2cart(vcv, lat, lon)
    x94, y94, z94 = llh2xyz(lat, lon, ell_ht_in)
    x20, y20, z20, vcv94 = conform7(x94, y94, z94, -gda94_to_gda2020, vcv=vcv)
    lat, lon, ell_ht_out = xyz2llh(x20, y20, z20)
    if vcv94 is not None:
        vcv94 = vcv_cart2local(vcv94, lat, lon)
    if ell_ht is False:
        ell_ht_out = 0
    hemisphere, zone20, east20, north20, psf, gridconv = geo2grid(lat, lon)
    return zone20, east20, north20, round(ell_ht_out, 4), vcv94


def transform_atrf2014_to_gda2020(x, y, z, epoch_from, vcv=None):
    """
    Transforms Cartesian (x, y, z) Coordinates in terms of the Australian Terrestrial Reference Frame (ATRF) at
    a specified epoch to coordinates in terms of Geocentric Datum of Australia 2020 (GDA2020 - reference epoch 2020.0)
    :param x: ATRF Cartesian X Coordinate (m)
    :param y: ATRF Cartesian Y Coordinate (m)
    :param z: ATRF Cartesian Z Coordinate (m)
    :param epoch_from: ATRF Coordinate Epoch (datetime.date Object)
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :return: Cartesian X, Y, Z Coordinates and vcv matrix in terms of GDA2020
    """
    return conform14(x, y, z, epoch_from, atrf2014_to_gda2020, vcv=vcv)


def transform_gda2020_to_atrf2014(x, y, z, epoch_to, vcv=None):
    """
    Transforms Cartesian (x, y, z) Coordinates in terms of Geocentric Datum of Australia 2020
    (GDA2020 - reference epoch 2020.0) to coordinates in terms of the Australian Terrestrial Reference Frame (ATRF) at
    a specified epoch
    :param x: GDA2020 Cartesian X Coordinate (m)
    :param y: GDA2020 Cartesian Y Coordinate (m)
    :param z: GDA2020 Cartesian Z Coordinate (m)
    :param epoch_to: ATRF Coordinate Epoch (datetime.date Object)
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :return: Cartesian X, Y, Z Coordinates and vcv matrix in terms of ATRF at the specified Epoch
    """
    return conform14(x, y, z, epoch_to, -atrf2014_to_gda2020, vcv=vcv)


def ntv2_2d(ntv2_grid, lat, lon, forward_tf=True, method='bicubic'):
    """
    Performs a 2D transformation based on ntv2 grid shifts.
    :param ntv2_grid: Ntv2Grid object (create with read_ntv2_file() function in geodepy.ntv2reader module)
    :param lat: latitude in decimal degrees
    :param lon: longitude in decimal degrees
    :param forward_tf: True/False:
                       - True applies the shifts in the direction given in the grid.
                       - False applies the shifts in the opposite direction of the grid
    :param method: Interpolation strategy - either 'bicubic' or 'bilinear'
    :return: Transformed latitude and longitude
    """

    # validate input data
    if not isinstance(ntv2_grid, NTv2Grid):
        raise TypeError('ntv2_grid must be Ntv2Grid object')
    if method != 'bicubic' and method != 'bilinear':
        raise ValueError(f'interpolation strategy "{method}" not supported')

    # interrogate grid
    shifts = interpolate_ntv2(ntv2_grid, lat, lon, method=method)

    # null results are outside of grid extents.
    if shifts[0] is None:
        raise ValueError('Coordinate outside of grid extents')

    if forward_tf:
        tf_lat = lat + shifts[0] / 3600
        tf_lon = lon - shifts[1] / 3600
    else:
        tf_lat = lat - shifts[0] / 3600
        tf_lon = lon + shifts[1] / 3600

    return tf_lat, tf_lon


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
from geodepy.constants import Transformation, atrf_gda2020,\
    gda94_to_gda2020
from geodepy.convert import hp2dec, geo2grid, \
    grid2geo, xyz2llh, llh2xyz


def conform7(x, y, z, trans):
    """
    Performs a Helmert 7 Parameter Conformal Transformation using Cartesian point co-ordinates
    and a predefined transformation object.
    :param x: Cartesian X (m)
    :param y: Cartesian Y (m)
    :param z: Cartesian Z (m)
    :param trans: Transformation Object (note: this function ignores all time-dependent variables)
    :return: Transformed X, Y, Z Cartesian Co-ordinates
    """
    if type(trans) != Transformation:
        raise ValueError('trans must be a Transformation Object')
    # Create XYZ Vector
    xyz_before = np.array([[x],
                           [y],
                           [z]])
    # Convert Units for Transformation Parameters
    scale = trans.sc / 1000000
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
    # Conformal Transform Eq
    xyz_after = translation + (1 + scale) * np.dot(rotation, xyz_before)
    # Convert Vector to Separate Variables
    xtrans = float(xyz_after[0])
    ytrans = float(xyz_after[1])
    ztrans = float(xyz_after[2])
    return xtrans, ytrans, ztrans


def conform14(x, y, z, to_epoch, trans):
    """
    Performs a Helmert 14 Parameter Conformal Transformation using Cartesian point co-ordinates
    and a predefined transformation object. The transformation parameters are projected from
    the transformation objects reference epoch to a specified epoch.
    :param x: Cartesian X (m)
    :param y: Cartesian Y (m)
    :param z: Cartesian Z (m)
    :param to_epoch: Epoch co-ordinate transformation is performed at (datetime.date Object)
    :param trans: Transformation Object
    :return: Cartesian X, Y, Z co-ordinates transformed using Transformation parameters at desired epoch
    """
    if type(trans) != Transformation:
        raise ValueError('trans must be a Transformation Object')
    if type(to_epoch) != datetime.date:
        raise ValueError('to_epoch must be a datetime.date Object')
    # Calculate 7 Parameters from 14 Parameter Transformation Object
    timetrans = trans + to_epoch
    # Perform Transformation
    xtrans, ytrans, ztrans = conform7(x, y, z, timetrans)
    return xtrans, ytrans, ztrans


def mga94_to_mga2020(zone, east, north, ell_ht=False):
    """
    Performs conformal transformation of Map Grid of Australia 1994 to Map Grid of Australia 2020 Coordinates
    using the GDA2020 Tech Manual v1.2 7 parameter similarity transformation parameters
    :param zone: Zone Number - 1 to 60
    :param east: Easting (m, within 3330km of Central Meridian)
    :param north: Northing (m, 0 to 10,000,000m)
    :param ell_ht: Ellipsoid Height (m) (optional)
    :return: MGA2020 Zone, Easting, Northing and Ellipsoid Height (if none provided, returns 0)
    """
    lat, lon, psf, gridconv = grid2geo(zone, east, north)
    if ell_ht is False:
        ell_ht_in = 0
    else:
        ell_ht_in = ell_ht
    x94, y94, z94 = llh2xyz(lat, lon, ell_ht_in)
    x20, y20, z20 = conform7(x94, y94, z94, gda94_to_gda2020)
    lat, lon, ell_ht_out = xyz2llh(x20, y20, z20)
    if ell_ht is False:
        ell_ht_out = 0
    hemisphere, zone20, east20, north20, psf, gridconv = geo2grid(lat, lon)
    return zone20, east20, north20, round(ell_ht_out, 4)


def mga2020_to_mga94(zone, east, north, ell_ht=False):
    """
    Performs conformal transformation of Map Grid of Australia 2020 to Map Grid of Australia 1994 Coordinates
    using the reverse form of the GDA2020 Tech Manual v1.2 7 parameter similarity transformation parameters
    :param zone: Zone Number - 1 to 60
    :param east: Easting (m, within 3330km of Central Meridian)
    :param north: Northing (m, 0 to 10,000,000m)
    :param ell_ht: Ellipsoid Height (m) (optional)
    :return: MGA1994 Zone, Easting, Northing and Ellipsoid Height (if none provided, returns 0)
    """
    lat, lon, psf, gridconv = grid2geo(zone, east, north)
    if ell_ht is False:
        ell_ht_in = 0
    else:
        ell_ht_in = ell_ht
    x94, y94, z94 = llh2xyz(lat, lon, ell_ht_in)
    x20, y20, z20 = conform7(x94, y94, z94, -gda94_to_gda2020)
    lat, lon, ell_ht_out = xyz2llh(x20, y20, z20)
    if ell_ht is False:
        ell_ht_out = 0
    hemisphere, zone20, east20, north20, psf, gridconv = geo2grid(lat, lon)
    return zone20, east20, north20, round(ell_ht_out, 4)


def atrf2014_to_gda2020(x, y, z, epoch_from):
    """
    Transforms Cartesian (x, y, z) Coordinates in terms of the Australian Terrestrial Reference Frame (ATRF) at
    a specified epoch to coordinates in terms of Geocentric Datum of Australia 2020 (GDA2020 - reference epoch 2020.0)
    :param x: ATRF Cartesian X Coordinate (m)
    :param y: ATRF Cartesian Y Coordinate (m)
    :param z: ATRF Cartesian Z Coordinate (m)
    :param epoch_from: ATRF Coordinate Epoch (datetime.date Object)
    :return: Cartesian X, Y, Z Coordinates in terms of GDA2020
    """
    return conform14(x, y, z, epoch_from, atrf_gda2020)


def gda2020_to_atrf2014(x, y, z, epoch_to):
    """
    Transforms Cartesian (x, y, z) Coordinates in terms of Geocentric Datum of Australia 2020
    (GDA2020 - reference epoch 2020.0) to coordinates in terms of the Australian Terrestrial Reference Frame (ATRF) at
    a specified epoch
    :param x: GDA2020 Cartesian X Coordinate (m)
    :param y: GDA2020 Cartesian Y Coordinate (m)
    :param z: GDA2020 Cartesian Z Coordinate (m)
    :param epoch_to: ATRF Coordinate Epoch (datetime.date Object)
    :return: Cartesian X, Y, Z Coordinate in terms of ATRF at the specified Epoch
    """
    return conform14(x, y, z, epoch_to, -atrf_gda2020)

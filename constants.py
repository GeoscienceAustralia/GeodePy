#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Constants Module
"""

from decimal import *

c_vac = 299792.458
k_0 = 0.9996


# Ellipsoid Constants
class Ellipsoid(object):
    def __init__(self, semimaj, inversef):
        self.semimaj = semimaj
        self.inversef = inversef


grs80 = Ellipsoid(6378137, Decimal('298.25722210088'))


# Projections
class Projection(object):
    def __init__(self, falseeast, falsenorth, cmscale, zonewidth, initialcm):
        self.falseeast = falseeast
        self.falsenorth = falsenorth
        self.cmscale = cmscale
        self.zonewidth = zonewidth
        self.initialcm = initialcm


utm = Projection(500000, 10000000, Decimal('0.9996'), 6, -177)


# Helmert 7 Parameter Transformation Parameters
class Transformation(object):
    def __init__(self, from_datum, to_datum, ref_epoch,
                 tx, ty, tz, sc, rx, ry, rz,
                 d_tx=0, d_ty=0, d_tz=0, d_sc=0, d_rx=0, d_ry=0, d_rz=0):
        self.from_datum = from_datum    # Datum Transforming From
        self.to_datum = to_datum        # Datum Transforming To
        self.ref_epoch = ref_epoch      # Reference Epoch (YYYY.DOY)
        self.tx = tx                    # Translation in x (m)
        self.ty = ty                    # Translation in y (m)
        self.tz = tz                    # Translation in z (m)
        self.sc = sc                    # Scale Change (parts per million)
        self.rx = rx                    # Rotation about x (arcseconds)
        self.ry = ry                    # Rotation about y (arcseconds)
        self.rz = rz                    # Rotation about z (arcseconds)
        self.d_tx = d_tx                # Rate of change in Translation in x (m per year)
        self.d_ty = d_ty                # Rate of change in Translation in y (m per year)
        self.d_tz = d_tz                # Rate of change in Translation in z (m per year)
        self.d_sc = d_sc                # Rate of change in Scale Change (parts per million per year)
        self.d_rx = d_rx                # Rate of change in Rotation about x (arcseconds per year)
        self.d_ry = d_ry                # Rate of change in Rotation about y (arcseconds per year)
        self.d_rz = d_rz                # Rate of change in Rotation about z (arcseconds per year)


gda94to20 = Transformation('GDA94', 'GDA2020', 0,
                           0.06155, -0.01087, -0.04019, -0.009994, -0.0394924, -0.0327221, -0.0328979)

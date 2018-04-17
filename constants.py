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


# Helmert 14 Parameter Transformation Parameters
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

    def __repr__(self):
        return ('Transformation: '
                + 'From ' + repr(self.from_datum) + ' to ' + repr(self.to_datum) + '\n'
                + 'Reference Epoch: ' + repr(self.ref_epoch) + '\n'
                + '  tx: ' + repr(self.tx) + 'm + ' + repr(self.d_tx) + 'm/yr' + '\n'
                + '  ty: ' + repr(self.ty) + 'm + ' + repr(self.d_ty) + 'm/yr' + '\n'
                + '  tz: ' + repr(self.tz) + 'm + ' + repr(self.d_tz) + 'm/yr' + '\n'
                + '  sc: ' + repr(self.sc) + 'ppm + ' + repr(self.d_sc) + 'ppm/yr' + '\n'
                + '  rx: ' + repr(self.rx) + '\" + ' + repr(self.d_rx) + '\"/yr' + '\n'
                + '  ry: ' + repr(self.ry) + '\" + ' + repr(self.d_ry) + '\"/yr' + '\n'
                + '  rz: ' + repr(self.rz) + '\" + ' + repr(self.d_rz) + '\"/yr' + '\n')

    def __neg__(self):
        return Transformation(self.to_datum,
                              self.from_datum,
                              self.ref_epoch,
                              -self.tx, -self.ty, -self.tz,
                              -self.sc,
                              -self.rx, -self.ry, -self.rz,
                              -self.d_tx, -self.d_ty, -self.d_tz,
                              -self.d_sc,
                              -self.d_rx, -self.d_ry, -self.d_rz)


def iers2trans(itrf_from, itrf_to, ref_epoch, tx, ty, tz, sc, rx, ry, rz, d_tx, d_ty, d_tz, d_sc, d_rx, d_ry, d_rz):
    """
    Used to convert IERS transformation parameters into Transformation Class parameters.
    Note: All rotation and delta rotation terms have sign change applied
    :param itrf_from: ITRF Realization Transforming From
    :param itrf_to: ITRF Realization Transforming To
    :param ref_epoch: Reference Epoch (YYYY.DOY)
    :param tx: Translation in x (mm)
    :param ty: Translation in y (mm)
    :param tz: Translation in z (mm)
    :param sc: Scale Change (parts per billion)
    :param rx: Rotation about x (milliarcseconds)
    :param ry: Rotation about y (milliarcseconds)
    :param rz: Rotation about z (milliarcseconds)
    :param d_tx: Rate of change in Translation in x (mm per year)
    :param d_ty: Rate of change in Translation in y (mm per year)
    :param d_tz: Rate of change in Translation in z (mm per year)
    :param d_sc: Rate of change in Scale Change (parts per billion per year)
    :param d_rx: Rate of change in Rotation about x (milliarcseconds per year)
    :param d_ry: Rate of change in Rotation about y (milliarcseconds per year)
    :param d_rz: Rate of change in Rotation about z (milliarcseconds per year)
    :return: Transformation Object following Australian Convention
    """
    return Transformation(itrf_from, itrf_to, ref_epoch,
                          round(tx / 1000, 8), round(ty / 1000, 8), round(tz / 1000, 8),
                          round(sc / 1000, 8),
                          round(-rx / 1000, 8), round(-ry / 1000, 8), round(-rz / 1000, 8),
                          round(d_tx / 1000, 8), round(d_ty / 1000, 8), round(d_tz / 1000, 8),
                          round(d_sc / 1000, 8),
                          round(-d_rx / 1000, 8), round(-d_ry / 1000, 8), round(-d_rz / 1000, 8))


# GDA1994 to GDA2020 Transformation Parameters from GDA2020 Tech Manual v1.1.1

gda94to20 = Transformation('GDA94', 'GDA2020', 0,
                           0.06155, -0.01087, -0.04019, -0.009994, -0.0394924, -0.0327221, -0.0328979)

# ITRF2014 Parameters
# link: http://itrf.ign.fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt

itrf14to08 = iers2trans('ITRF2014', 'ITRF2008', 2010.0,
                        1.6, 1.9, 2.4, -0.02, 0, 0, 0,
                        0.0, 0.0, -0.1, 0.03, 0, 0, 0)

itrf14to05 = iers2trans('ITRF2014', 'ITRF2005', 2010.0,
                        2.6, 1.0, -2.3, 0.92, 0, 0, 0,
                        0.3, 0.0, -0.1, 0.03, 0, 0, 0)

itrf14to00 = iers2trans('ITRF2014', 'ITRF2000', 2010.0,
                        0.7, 1.2, -26.1, 2.12, 0, 0, 0,
                        0.1, 0.1, -1.9, 0.11, 0, 0, 0)

itrf14to97 = iers2trans('ITRF2014', 'ITRF1997', 2010.0,
                        7.4, -0.5, -62.8, 3.80, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)

itrf14to96 = iers2trans('ITRF2014', 'ITRF1996', 2010.0,
                        7.4, -0.5, -62.8, 3.80, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)

itrf14to94 = iers2trans('ITRF2014', 'ITRF1994', 2010.0,
                        7.4, -0.5, -62.8, 3.80, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)

itrf14to93 = iers2trans('ITRF2014', 'ITRF1993', 2010.0,
                        -50.4, 3.3, -60.2, 4.29, -2.81, -3.38, 0.40,
                        -2.8, -0.1, -2.5, 0.12, -0.11, -0.19, 0.07)

itrf14to92 = iers2trans('ITRF2014', 'ITRF1992', 2010.0,
                        15.4, 1.5, -70.8, 3.09, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)

itrf14to91 = iers2trans('ITRF2014', 'ITRF1991', 2010.0,
                        27.4, 15.5, -76.8, 4.49, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)

itrf14to90 = iers2trans('ITRF2014', 'ITRF1990', 2010.0,
                        25.4, 11.5, -92.8, 4.79, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)

itrf14to89 = iers2trans('ITRF2014', 'ITRF1989', 2010.0,
                        30.4, 35.5, -130.8, 8.19, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)

itrf14to88 = iers2trans('ITRF2014', 'ITRF1988', 2010.0,
                        25.4, -0.5, -154.8, 11.29, 0.1, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.2)
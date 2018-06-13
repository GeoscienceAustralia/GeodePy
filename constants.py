#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Constants Module
"""

from math import sqrt

c_vac = 299792.458
k_0 = 0.9996


# Ellipsoid Constants
class Ellipsoid(object):
    def __init__(self, semimaj, inversef):
        self.semimaj = semimaj
        self.inversef = inversef
        self.f = 1 / self.inversef
        self.semimin = float(self.semimaj * (1 - self.f))
        self.ecc1sq = float(self.f * (2 - self.f))
        self.ecc2sq = float(self.ecc1sq / (1 - self.ecc1sq))
        self.ecc1 = sqrt(self.ecc1sq)
        self.n = float(self.f / (2 - self.f))
        self.n2 = self.n ** 2


# Geodetic Reference System 1980
grs80 = Ellipsoid(6378137, 298.25722210088)

# Australian National Spheroid (See GDA94 Tech Manual v2.4 pp 7)
ans = Ellipsoid(6378160, 298.25)


# Projections
class Projection(object):
    def __init__(self, falseeast, falsenorth, cmscale, zonewidth, initialcm):
        self.falseeast = falseeast
        self.falsenorth = falsenorth
        self.cmscale = cmscale
        self.zonewidth = zonewidth
        self.initialcm = initialcm


utm = Projection(500000, 10000000, 0.9996, 6, -177)


# Helmert 14 Parameter Transformation Parameters
class Transformation(object):
    def __init__(self, from_datum, to_datum, ref_epoch,
                 tx, ty, tz, sc, rx, ry, rz,
                 d_tx=0.0, d_ty=0.0, d_tz=0.0, d_sc=0.0, d_rx=0.0, d_ry=0.0, d_rz=0.0):
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
        """
        Reverses Direction of Transformation Object
        i.e. ITRF2014 to ITRF2000 transformation becomes ITRF2000 to ITRF2014 transformation
        :return: Reversed Direction Transformation Object
        """
        return Transformation(self.to_datum,
                              self.from_datum,
                              self.ref_epoch,
                              -self.tx, -self.ty, -self.tz,
                              -self.sc,
                              -self.rx, -self.ry, -self.rz,
                              -self.d_tx, -self.d_ty, -self.d_tz,
                              -self.d_sc,
                              -self.d_rx, -self.d_ry, -self.d_rz)

    def __add__(self, other):
        """
        Change Reference Epoch (add decimal year).
        Advances all transformation parameters by their respective rates of change.
        :param other: Decimal Year
        :return: Transformation object with parameters and ref epoch advanced by input year
        """
        if type(other) == int or type(other) == float:
            return Transformation(self.to_datum,
                                  self.from_datum,
                                  self.ref_epoch + other,
                                  round(self.tx + (self.d_tx * other), 8),
                                  round(self.ty + (self.d_ty * other), 8),
                                  round(self.tz + (self.d_tz * other), 8),
                                  round(self.sc + (self.d_sc * other), 8),
                                  round(self.rx + (self.d_rx * other), 8),
                                  round(self.ry + (self.d_ry * other), 8),
                                  round(self.rz + (self.d_rz * other), 8),
                                  self.d_tx,
                                  self.d_ty,
                                  self.d_tz,
                                  self.d_sc,
                                  self.d_rx,
                                  self.d_ry,
                                  self.d_rz)
        """ Experimental - Not Yet Functional
        elif type(other) == Transformation:
            epoch_diff = other.ref_epoch - self.ref_epoch
            return Transformation(other.to_datum,
                                  self.from_datum,
                                  self.ref_epoch,
                                  round(self.tx + (other.tx * epoch_diff), 8),
                                  round(self.ty + (other.ty * epoch_diff), 8),
                                  round(self.tz + (other.tz * epoch_diff), 8),
                                  round(self.sc + (other.sc * epoch_diff), 8),
                                  round(self.rx + (other.rx * epoch_diff), 8),
                                  round(self.ry + (other.ry * epoch_diff), 8),
                                  round(self.rz + (other.rz * epoch_diff), 8),
                                  round(self.d_tx + (other.d_tx * epoch_diff), 8),
                                  round(self.d_ty + (other.d_ty * epoch_diff), 8),
                                  round(self.d_tz + (other.d_tz * epoch_diff), 8),
                                  round(self.d_sc + (other.d_sc * epoch_diff), 8),
                                  round(self.d_rx + (other.d_rx * epoch_diff), 8),
                                  round(self.d_ry + (other.d_ry * epoch_diff), 8),
                                  round(self.d_rz + (other.d_rz * epoch_diff), 8))
            """

    def __sub__(self, other):
        """
        Change Reference Epoch (subtract decimal year).
        Retracts all transformation parameters by their respective rates of change.
        :param other: Decimal Year
        :return: Transformation object with parameters and ref epoch advanced by input year
        """
        return Transformation(self.to_datum,
                              self.from_datum,
                              self.ref_epoch - other,
                              round(self.tx - (self.d_tx * other), 8),
                              round(self.ty - (self.d_ty * other), 8),
                              round(self.tz - (self.d_tz * other), 8),
                              round(self.sc - (self.d_sc * other), 8),
                              round(self.rx - (self.d_rx * other), 8),
                              round(self.ry - (self.d_ry * other), 8),
                              round(self.rz - (self.d_rz * other), 8),
                              self.d_tx,
                              self.d_ty,
                              self.d_tz,
                              self.d_sc,
                              self.d_rx,
                              self.d_ry,
                              self.d_rz)


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

gda94to20 = Transformation('GDA1994', 'GDA2020', 0,
                           0.06155, -0.01087, -0.04019, -0.009994, -0.0394924, -0.0327221, -0.0328979)

itrf14togda20 = Transformation('ITRF2014', 'GDA2020', 2020,
                               0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0.00150379, 0.00118346, 0.00120716)

# GDA1994 to ITRF Transformation Parameters from Dawson and Woods (2010)
# AGD66 and AGD84 to GDA94 Transformation Parameters from GDA94 Tech Manual v2.4
# link: http://www.icsm.gov.au/datum/gda2020-and-gda94-technical-manuals

itrf08togda94 = Transformation('ITRF2008', 'GDA1994', 1994.0,
                               -0.08468, -0.01942, 0.03201, 0.00971, -0.0004254, 0.0022578, 0.0024015,
                               0.00142, 0.00134, 0.00090, 0.000109, 0.0015461, 0.001820, 0.0011551)

itrf05togda94 = Transformation('ITRF2005', 'GDA1994', 1994.0,
                               -0.07973, -0.00686, 0.03803, 0.006636, -0.0000351, 0.0021211, 0.0021411,
                               0.00225, -0.00062, -0.00056, 0.000294, 0.0014707, 0.0011443, 0.0011701)

itrf00togda94 = Transformation('ITRF2000', 'GDA1994', 1994.0,
                               -0.04591, -0.02985, -0.02037, 0.00707, -0.0016705, 0.0004594, 0.0019356,
                               -0.00466, 0.00355, 0.01124, 0.000249, 0.0017454, 0.0014868, 0.001224)

itrf97togda94 = Transformation('ITRF1997', 'GDA1994', 1994.0,
                               -0.01463, -0.02762, -0.02532, 0.006695, -0.0017893, -0.0006047, 0.0009962,
                               -0.00860, 0.00036, 0.01125, 0.000007, 0.0016394, 0.0015198, 0.0013801)

itrf96togda94 = Transformation('ITRF1996', 'GDA1994', 1994.0,
                               0.02454, -0.03643, -0.06812, 0.006901, -0.0027359, -0.0020431, 0.0003731,
                               -0.02180, 0.00471, 0.02627, 0.000388, 0.0020203, 0.0021735, 0.0016290)

agd84togda94 = Transformation('AGD84', 'GDA94', 0,
                              -117.763, -51.510, 139.061, -0.191, -0.292, -0.443, -0.277)

agd66togda94 = Transformation('AGD1966', 'GDA1994', 0,
                              -117.808, -51.536, 137.784, -0.290, -0.303, -0.446, -0.234)

agd66togda94_act = Transformation('AGD66', 'GDA94', 0,
                                  -129.193, -41.212, 130.730, -2.955, -0.246, -0.374, -0.329)

agd66togda94_tas = Transformation('AGD66', 'GDA94', 0,
                                  -120.271, -64.543, 161.632, 2.499, -0.217, 0.067, 0.129)

agd66togda94_vicnsw = Transformation('AGD66', 'GDA94', 0,
                                     -119.353, -48.301, 139.484, -0.613, -0.415, -0.260, -0.437)

agd66togda94_nt = Transformation('AGD66', 'GDA94', 0,
                                 -124.133, -42.003, 137.400, -1.854, 0.008, -0.557, -0.178)


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
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

itrf14to96 = iers2trans('ITRF2014', 'ITRF1996', 2010.0,
                        7.4, -0.5, -62.8, 3.80, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

itrf14to94 = iers2trans('ITRF2014', 'ITRF1994', 2010.0,
                        7.4, -0.5, -62.8, 3.80, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

itrf14to93 = iers2trans('ITRF2014', 'ITRF1993', 2010.0,
                        -50.4, 3.3, -60.2, 4.29, -2.81, -3.38, 0.40,
                        -2.8, -0.1, -2.5, 0.12, -0.11, -0.19, 0.07)

itrf14to92 = iers2trans('ITRF2014', 'ITRF1992', 2010.0,
                        15.4, 1.5, -70.8, 3.09, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

itrf14to91 = iers2trans('ITRF2014', 'ITRF1991', 2010.0,
                        27.4, 15.5, -76.8, 4.49, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

itrf14to90 = iers2trans('ITRF2014', 'ITRF1990', 2010.0,
                        25.4, 11.5, -92.8, 4.79, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

itrf14to89 = iers2trans('ITRF2014', 'ITRF1989', 2010.0,
                        30.4, 35.5, -130.8, 8.19, 0, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

itrf14to88 = iers2trans('ITRF2014', 'ITRF1988', 2010.0,
                        25.4, -0.5, -154.8, 11.29, 0.1, 0, 0.26,
                        0.1, -0.5, -3.3, 0.12, 0, 0, 0.02)

# ITRF2008 Parameters
# link: http://itrf.ign.fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt

itrf08to05 = iers2trans('ITRF2008', 'ITRF2005', 2000.0,
                        -2.0, -0.9, -4.7, 0.94, 0, 0, 0,
                        0.3, 0.0, 0.0, 0.0, 0, 0, 0)

itrf08to00 = iers2trans('ITRF2008', 'ITRF2000', 2000.0,
                        -1.9, -1.7, -10.5, 1.34, 0, 0, 0,
                        0.1, 0.1, -1.8, 0.08, 0, 0, 0)

itrf08to97 = iers2trans('ITRF2008', 'ITRF1997', 2000.0,
                        4.8, 2.6, -33.2, 2.92, 0, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

itrf08to96 = iers2trans('ITRF2008', 'ITRF1996', 2000.0,
                        4.8, 2.6, -33.2, 2.92, 0, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

itrf08to94 = iers2trans('ITRF2008', 'ITRF1994', 2000.0,
                        4.8, 2.6, -33.2, 2.92, 0, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

itrf08to93 = iers2trans('ITRF2008', 'ITRF1993', 2000.0,
                        -24.0, 2.4, -38.6, 3.41, -1.71, -1.48, -0.30,
                        -2.8, -0.1, -2.4, 0.09, -0.11, -0.19, 0.07)

itrf08to92 = iers2trans('ITRF2008', 'ITRF1992', 2000.0,
                        12.8, 4.6, -41.2, 2.21, 0, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

itrf08to91 = iers2trans('ITRF2008', 'ITRF1991', 2000.0,
                        24.8, 18.6, -47.2, 3.61, 0, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

itrf08to90 = iers2trans('ITRF2008', 'ITRF1990', 2000.0,
                        22.8, 14.6, -63.2, 3.91, 0, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

itrf08to89 = iers2trans('ITRF2008', 'ITRF1989', 2000.0,
                        27.8, 38.6, -101.2, 7.31, 0, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

itrf08to88 = iers2trans('ITRF2008', 'ITRF1988', 2000.0,
                        22.8, 2.6, -125.2, 10.41, 0.10, 0, 0.06,
                        0.1, -0.5, -3.2, 0.09, 0, 0, 0.02)

# ITRF2005 Parameters
# link: http://itrf.ensg.ign.fr/ITRF_solutions/2005/tp_05-00.php

itrf05to00 = iers2trans('ITRF2005', 'ITRF2000', 2000.0,
                        0.1, -0.8, -5.8, 0.40, 0, 0, 0,
                        -0.2, 0.1, -1.8, 0.08, 0, 0, 0)

# ITRF2000 Parameters
# link: ftp://itrf.ensg.ign.fr/pub/itrf/ITRF.TP
# NOTE: This ref lists translations in centimetres. All other ITRF transformations are shown in millimetres.
# NOTE: All translations and rates of translation shown below have been converted to millimetres.

itrf00to97 = iers2trans('ITRF2000', 'ITRF1997', 1997.0,
                        6.7, 6.1, -18.5, 1.55, 0, 0, 0,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

itrf00to96 = iers2trans('ITRF2000', 'ITRF1996', 1997.0,
                        6.7, 6.1, -18.5, 1.55, 0, 0, 0,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

itrf00to94 = iers2trans('ITRF2000', 'ITRF1994', 1997.0,
                        6.7, 6.1, -18.5, 1.55, 0, 0, 0,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

itrf00to93 = iers2trans('ITRF2000', 'ITRF1993', 1988.0,
                        12.7, 6.5, -20.9, 1.95, -0.39, 0.80, -1.14,
                        -2.9, -0.2, -0.6, 0.01, -0.11, -0.19, 0.07)

itrf00to92 = iers2trans('ITRF2000', 'ITRF1992', 1988.0,
                        14.7, 13.5, -13.9, 0.75, 0, 0, -0.18,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

itrf00to91 = iers2trans('ITRF2000', 'ITRF1991', 1988.0,
                        26.7, 27.5, -19.9, 2.15, 0, 0, -0.18,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

itrf00to90 = iers2trans('ITRF2000', 'ITRF1990', 1988.0,
                        14.7, 13.5, -13.9, 0.75, 0, 0, -0.18,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

itrf00to89 = iers2trans('ITRF2000', 'ITRF1989', 1988.0,
                        29.7, 47.5, -73.9, 5.85, 0, 0, -0.18,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

itrf00to88 = iers2trans('ITRF2000', 'ITRF1988', 1988.0,
                        24.7, 11.5, -97.9, 8.95, 0, 0, -0.18,
                        0.0, -0.6, -1.4, 0.01, 0, 0, 0.02)

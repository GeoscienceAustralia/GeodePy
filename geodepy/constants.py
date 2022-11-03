#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Constants Module
"""

from math import sqrt
from datetime import date

c_vac = 299792.458
k_0 = 0.9996


# Ellipsoid constants
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
        self.meanradius = (2 * self.semimaj + self.semimin)/3


# Geodetic Reference System 1980
# www.epsg-registry.org/export.htm?gml=urn:ogc:def:ellipsoid:EPSG::7019
grs80 = Ellipsoid(6378137, 298.257222101)


# World Geodetic System 1984
# www.epsg-registry.org/export.htm?gml=urn:ogc:def:ellipsoid:EPSG::7030
wgs84 = Ellipsoid(6378137, 298.257223563)


# Australian National Spheroid
# www.epsg-registry.org/export.htm?gml=urn:ogc:def:ellipsoid:EPSG::7003
ans = Ellipsoid(6378160, 298.25)


# International (Hayford) 1924
# www.epsg-registry.org/export.htm?gml=urn:ogc:def:ellipsoid:EPSG::7022
intl24 = Ellipsoid(6378388, 297)


# Projections
class Projection(object):
    def __init__(self, falseeast, falsenorth, cmscale, zonewidth, initialcm):
        """
        Transverse Mercator Projection Parameters
        :param falseeast: Easting (m) assigned to Central Meridian
        :param falsenorth: Northing (m) assigned to Equator
        :param cmscale: Central Meridian Scale Factor (unitless, 1 is no scale)
        :param zonewidth: Width (decimal degrees) of each TM Zone
        :param initialcm: Longitude (decimal degrees) of TM Zone 1
        """
        self.falseeast = falseeast
        self.falsenorth = falsenorth
        self.cmscale = cmscale
        self.zonewidth = zonewidth
        self.initialcm = initialcm


utm = Projection(500000, 10000000, 0.9996, 6, -177)

# Integrated Survey Grid - used in NSW as the projection for AGD66
# Spatial Services projections page - https://www.spatial.nsw.gov.au/surveying/geodesy/projections
# ISG Technical Manual - https://www.spatial.nsw.gov.au/__data/assets/pdf_file/0017/25730/ISG.pdf
isg = Projection(300000, 5000000, 0.99994, 2, -177)


# Helmert 14 parameter transformation
class Transformation(object):
    def __init__(self, from_datum, to_datum, ref_epoch,
                 tx, ty, tz, sc, rx, ry, rz,
                 d_tx=0.0, d_ty=0.0, d_tz=0.0, d_sc=0.0,
                 d_rx=0.0, d_ry=0.0, d_rz=0.0, tf_sd=None):
        self.from_datum = from_datum   # Datum transforming from
        self.to_datum = to_datum       # Datum transforming to
        self.ref_epoch = ref_epoch     # Ref. epoch (datetime.date object)
        self.tx = tx                   # Translation in X (m)
        self.ty = ty                   # Translation in Y (m)
        self.tz = tz                   # Translation in Z (m)
        self.sc = sc                   # Scale factor (ppm)
        self.rx = rx                   # Rotation about X (arcsec)
        self.ry = ry                   # Rotation about Y (arcsec)
        self.rz = rz                   # Rotation about Z (arcsec)
        self.d_tx = d_tx               # Translation rate of change in X (m/yr)
        self.d_ty = d_ty               # Translation rate of change in Y (m/yr)
        self.d_tz = d_tz               # Translation rate of change in Z (m/yr)
        self.d_sc = d_sc               # Scale factor rate of change (ppm/yr)
        self.d_rx = d_rx               # Rot rate of change about X (arcsec/yr)
        self.d_ry = d_ry               # Rot rate of change about Y (arcsec/yr)
        self.d_rz = d_rz               # Rot rate of change about Z (arcsec/yr)
        self.tf_sd = tf_sd             # TransformationSD object

    def __repr__(self):
        return ('Transformation: '
                + 'From ' + repr(self.from_datum) +
                ' to ' + repr(self.to_datum) + '\n'
                + 'Reference Epoch: ' + repr(self.ref_epoch) + '\n'
                + '  tx: ' + repr(self.tx) + 'm + ' + repr(self.d_tx) +
                'm/yr' + '\n'
                + '  ty: ' + repr(self.ty) + 'm + ' + repr(self.d_ty) +
                'm/yr' + '\n'
                + '  tz: ' + repr(self.tz) + 'm + ' + repr(self.d_tz) +
                'm/yr' + '\n'
                + '  sc: ' + repr(self.sc) + 'ppm + ' + repr(self.d_sc) +
                'ppm/yr' + '\n'
                + '  rx: ' + repr(self.rx) + '\" + ' + repr(self.d_rx) +
                '\"/yr' + '\n'
                + '  ry: ' + repr(self.ry) + '\" + ' + repr(self.d_ry) +
                '\"/yr' + '\n'
                + '  rz: ' + repr(self.rz) + '\" + ' + repr(self.d_rz) +
                '\"/yr' + '\n')

    def __neg__(self):
        """
        Reverses the direction of the transformation object, i.e.,
        the ITRF2014 to ITRF2000 transformation becomes the ITRF2000 to
        ITRF2014 transformation.

        :return: reversed direction transformation object
        """
        return Transformation(self.to_datum,
                              self.from_datum,
                              self.ref_epoch,
                              -self.tx, -self.ty, -self.tz,
                              -self.sc,
                              -self.rx, -self.ry, -self.rz,
                              -self.d_tx, -self.d_ty, -self.d_tz,
                              -self.d_sc,
                              -self.d_rx, -self.d_ry, -self.d_rz,
                              tf_sd=self.tf_sd
                              )

    def __add__(self, other):
        """
        Change the transformation epoch to a specified date. Advances all
        transformation parameters by their respective rates of change.

        :param other: datetime.date object
        :return: transformation object with parameters and ref epoch moved to
        the specified date
        """
        if type(other) == date:
            timediff = (other - self.ref_epoch).days/365.25

            if type(self.tf_sd) == TransformationSD:
                self.tf_sd.sd_tx = (self.tf_sd.sd_tx**2 + (self.tf_sd.sd_d_tx * timediff)**2) ** 0.5
                self.tf_sd.sd_ty = (self.tf_sd.sd_ty**2 + (self.tf_sd.sd_d_ty * timediff)**2) ** 0.5
                self.tf_sd.sd_tz = (self.tf_sd.sd_tz**2 + (self.tf_sd.sd_d_tz * timediff)**2) ** 0.5
                self.tf_sd.sd_sc = (self.tf_sd.sd_sc**2 + (self.tf_sd.sd_d_sc * timediff)**2) ** 0.5
                self.tf_sd.sd_rx = (self.tf_sd.sd_rx**2 + (self.tf_sd.sd_d_rx * timediff)**2) ** 0.5
                self.tf_sd.sd_ry = (self.tf_sd.sd_ry**2 + (self.tf_sd.sd_d_ry * timediff)**2) ** 0.5
                self.tf_sd.sd_rz = (self.tf_sd.sd_rz**2 + (self.tf_sd.sd_d_rz * timediff)**2) ** 0.5

            return Transformation(self.to_datum,
                                  self.from_datum,
                                  other,
                                  round(self.tx + (self.d_tx * timediff), 8),
                                  round(self.ty + (self.d_ty * timediff), 8),
                                  round(self.tz + (self.d_tz * timediff), 8),
                                  round(self.sc + (self.d_sc * timediff), 8),
                                  round(self.rx + (self.d_rx * timediff), 8),
                                  round(self.ry + (self.d_ry * timediff), 8),
                                  round(self.rz + (self.d_rz * timediff), 8),
                                  self.d_tx,
                                  self.d_ty,
                                  self.d_tz,
                                  self.d_sc,
                                  self.d_rx,
                                  self.d_ry,
                                  self.d_rz,
                                  self.tf_sd
                                  )
        else:
            ValueError('supports adding datetime.date objects only')


class TransformationSD(object):
    def __init__(self, sd_tx=None, sd_ty=None, sd_tz=None, sd_sc=None,
                 sd_rx=None, sd_ry=None, sd_rz=None,
                 sd_d_tx=None, sd_d_ty=None, sd_d_tz=None, sd_d_sc=None,
                 sd_d_rx=None, sd_d_ry=None, sd_d_rz=None):
        self.sd_tx = sd_tx         # one-sigma uncertainty of tx (m)
        self.sd_ty = sd_ty         # one-sigma uncertainty of ty (m)
        self.sd_tz = sd_tz         # one-sigma uncertainty of tz (m)
        self.sd_sc = sd_sc         # one-sigma uncertainty of sc (ppm)
        self.sd_rx = sd_rx         # one-sigma uncertainty of rx (arcsec/yr)
        self.sd_ry = sd_ry         # one-sigma uncertainty of ry (arcsec/yr)
        self.sd_rz = sd_rz         # one-sigma uncertainty of rz (arcsec/yr)
        self.sd_d_tx = sd_d_tx     # one-sigma uncertainty of d_tx (m/yr)
        self.sd_d_ty = sd_d_ty     # one-sigma uncertainty of d_ty (m/yr)
        self.sd_d_tz = sd_d_tz     # one-sigma uncertainty of d_tz (m/yr)
        self.sd_d_sc = sd_d_sc     # one-sigma uncertainty of d_sc (ppm/yr)
        self.sd_d_rx = sd_d_rx     # one-sigma uncertainty of d_rx (arcsec/yr)
        self.sd_d_ry = sd_d_ry     # one-sigma uncertainty of d_ry (arcsec/yr)
        self.sd_d_rz = sd_d_rz     # one-sigma uncertainty of d_rz (arcsec/yr)


def iers2trans(itrf_from, itrf_to, ref_epoch, tx, ty, tz, sc, rx, ry, rz,
               d_tx, d_ty, d_tz, d_sc, d_rx, d_ry, d_rz):
    """
    Used to convert IERS transformation parameters into GeodePy Transformation
    class parameters.
    Note: All rotation and delta rotation terms have the sign change applied.
    :param itrf_from: ITRF realization transforming from
    :param itrf_to: ITRF realization transforming to
    :param ref_epoch: Reference epoch (YYYY.DOY)
    :param tx: Translation in X (mm)
    :param ty: Translation in Y (mm)
    :param tz: Translation in Z (mm)
    :param sc: Scale factor (ppb)
    :param rx: Rotation about X (milliarcsec)
    :param ry: Rotation about Y (milliarcsec)
    :param rz: Rotation about Z (milliarcsec)
    :param d_tx: Translation rate of change in X (mm/yr)
    :param d_ty: Translation rate of change in Y (mm/yr)
    :param d_tz: Translation rate of change in Z (mm/yr)
    :param d_sc: Scale factor rate of change (ppb/yr)
    :param d_rx: Rate of change in Rotation about X (milliarcsec/yr)
    :param d_ry: Rate of change in Rotation about X (milliarcsec/yr)
    :param d_rz: Rate of change in Rotation about X (milliarcsec/yr)
    :return: Transformation object following the Australian convention
    """
    return Transformation(itrf_from, itrf_to, ref_epoch,
                          round(tx / 1000, 8), round(ty / 1000, 8),
                          round(tz / 1000, 8), round(sc / 1000, 8),
                          round(-rx / 1000, 8), round(-ry / 1000, 8),
                          round(-rz / 1000, 8), round(d_tx / 1000, 8),
                          round(d_ty / 1000, 8), round(d_tz / 1000, 8),
                          round(d_sc / 1000, 8), round(-d_rx / 1000, 8),
                          round(-d_ry / 1000, 8), round(-d_rz / 1000, 8))


"""
Transformation Parameter Definitions

Please note transformation parameter definitions should align to the following
GeodePy conventions:

Transformation Variable Naming:
datumYYYY_to_datumYYYY[_suffix]
where:
datum - common abbreviation of datum name
YYYY - year (last two digits of year if before 2000, 4 digits if 2000 or later)
Suffix - optional extension used to differentiate between regional
transformation parameters

Transformation Standard Deviation Variable Naming:
datumYYYY_to_datumYYYY[_suffix]_sd

All parameters in a transformation variable definition should be explicitly
named rather than positionally defined only, as this eases manual checking
of parameters against sources. For example, the from datum should be defined
as from_datum='datumYYYY' rather than just 'datumYYYY'.

Reverse transformations should be explicitly defined using the negative
overload, i.e:
datum2YYYY_to_datum1YYYY = -datum1YYYY_to_datum2YYYY

References to the source of the transformation(s) parameters should be given,
including a web URL
"""
# GDA94 to GDA2020 transformation parameters [GDA2020 Tech Manual v1.2]
# Ref: https://www.icsm.gov.au/gda2020-and-gda94-technical-manuals

gda94_to_gda2020_sd = TransformationSD(
    sd_tx=0.0007, sd_ty=0.0006, sd_tz=0.0007,
    sd_sc=0.00010,
    sd_rx=0.000011, sd_ry=0.000010, sd_rz=0.000011,)

gda94_to_gda2020 = Transformation(
    from_datum='GDA94', to_datum='GDA2020', ref_epoch=0,
    tx=0.06155, ty=-0.01087, tz=-0.04019,
    sc=-0.009994,
    rx=-0.0394924, ry=-0.0327221, rz=-0.0328979,
    tf_sd=gda94_to_gda2020_sd)

gda2020_to_gda94 = -gda94_to_gda2020


# ITRF2014 to GDA2020 transformation parameters [GDA2020 Tech Manual v1.2].
# This transformation is also called the Australian Plate Motion Model and it
# was derived using the 109 ARGN and AuScope GNSS CORS that were used to define
# the National Measurement (Recognized-Value Standard of Measurement of
# Position)Determination 2017.
# Ref: https://www.legislation.gov.au/Details/F2017L01352

itrf2014_to_gda2020_sd = TransformationSD(
    sd_tx=0.00, sd_ty=0.00, sd_tz=0.00,
    sd_sc=0.00,
    sd_rx=0.00, sd_ry=0.00, sd_rz=0.00,
    sd_d_tx=0.00, sd_d_ty=0.00, sd_d_tz=0.00,
    sd_d_sc=0.00,
    sd_d_rx=0.00000417, sd_d_ry=0.00000401, sd_d_rz=0.00000370)

itrf2014_to_gda2020 = Transformation(
    from_datum='ITRF2014', to_datum='GDA2020', ref_epoch=date(2020, 1, 1),
    tx=0, ty=0, tz=0, sc=0, rx=0, ry=0, rz=0,
    d_tx=0, d_ty=0, d_tz=0, d_sc=0,
    d_rx=0.00150379, d_ry=0.00118346, d_rz=0.00120716,
    tf_sd=itrf2014_to_gda2020_sd)

gda2020_to_itrf2014 = -itrf2014_to_gda2020


# ATRF2014 to GDA2020 transformation parameters. Note this transformation is
# identical to the above Australian Plate Motion Model
# Ref: https://www.icsm.gov.au/australian-terrestrial-reference-frame

atrf2014_to_gda2020_sd = TransformationSD(
    sd_tx=0.00, sd_ty=0.00, sd_tz=0.00,
    sd_sc=0.00,
    sd_rx=0.00, sd_ry=0.00, sd_rz=0.00,
    sd_d_tx=0.00, sd_d_ty=0.00, sd_d_tz=0.00,
    sd_d_sc=0.00,
    sd_d_rx=0.00000417, sd_d_ry=0.00000401, sd_d_rz=0.00000370)

atrf2014_to_gda2020 = Transformation(
    from_datum='ATRF2014', to_datum='GDA2020', ref_epoch=date(2020, 1, 1),
    tx=0, ty=0, tz=0, sc=0, rx=0, ry=0, rz=0,
    d_tx=0, d_ty=0, d_tz=0, d_sc=0,
    d_rx=0.00150379, d_ry=0.00118346, d_rz=0.00120716,
    tf_sd=atrf2014_to_gda2020_sd)


# GDA94 to ITRF transformation parameters [Dawson and Woods (2010)]
# AGD66 and AGD84 to GDA94 transformation parameters [GDA94 Tech Manual v2.4]
# http://www.icsm.gov.au/datum/gda2020-and-gda94-technical-manuals

# Transformation Standard Deviations
itrf2008_to_gda94_sd = TransformationSD(
    sd_tx=0.00091, sd_ty=0.00078, sd_tz=0.00106,
    sd_sc=0.000126,
    sd_rx=0.0000221, sd_ry=0.0000236, sd_rz=0.0000194,
    sd_d_tx=0.00008, sd_d_ty=0.00007, sd_d_tz=0.00011,
    sd_d_sc=0.000013,
    sd_d_rx=0.0000028, sd_d_ry=0.0000030, sd_d_rz=0.0000023)

itrf2005_to_gda94_sd = TransformationSD(
    sd_tx=0.00256, sd_ty=0.00187, sd_tz=0.00337,
    sd_sc=0.000227,
    sd_rx=0.0000883, sd_ry=0.0000972, sd_rz=0.0000600,
    sd_d_tx=0.00028, sd_d_ty=0.00020, sd_d_tz=0.00036,
    sd_d_sc=0.000022,
    sd_d_rx=0.0000096, sd_d_ry=0.0000106, sd_d_rz=0.0000070)

itrf2000_to_gda94_sd = TransformationSD(
    sd_tx=0.00874, sd_ty=0.00412, sd_tz=0.01105,
    sd_sc=0.000423,
    sd_rx=0.0002852, sd_ry=0.0003602, sd_rz=0.0001567,
    sd_d_tx=0.00166, sd_d_ty=0.00083, sd_d_tz=0.00209,
    sd_d_sc=0.000076,
    sd_d_rx=0.0000537, sd_d_ry=0.0000677, sd_d_rz=0.0000317)

itrf97_to_gda94_sd = TransformationSD(
    sd_tx=0.01107, sd_ty=0.00574, sd_tz=0.01427,
    sd_sc=0.000425,
    sd_rx=0.0003757, sd_ry=0.0004642, sd_rz=0.0001955,
    sd_d_tx=0.00276, sd_d_ty=0.00136, sd_d_tz=0.00371,
    sd_d_sc=0.000101,
    sd_d_rx=0.0000967, sd_d_ry=0.0001196, sd_d_rz=0.0000453)

itrf96_to_gda94_sd = TransformationSD(
    sd_tx=0.02033, sd_ty=0.01059, sd_tz=0.02866,
    sd_sc=0.000653,
    sd_rx=0.0007627, sd_ry=0.0009049, sd_rz=0.0003225,
    sd_d_tx=0.00753, sd_d_ty=0.00391, sd_d_tz=0.01074,
    sd_d_sc=0.000228,
    sd_d_rx=0.0002850, sd_d_ry=0.0003382, sd_d_rz=0.0001169)

# Transformations
itrf2008_to_gda94 = Transformation(
    from_datum='ITRF2008', to_datum='GDA94', ref_epoch=date(1994, 1, 1),
    tx=-0.08468, ty=-0.01942, tz=0.03201,
    sc=0.00971,
    rx=-0.0004254, ry=0.0022578, rz=0.0024015,
    d_tx=0.00142, d_ty=0.00134, d_tz=0.00090,
    d_sc=0.000109,
    d_rx=0.0015461, d_ry=0.0011820, d_rz=0.0011551,
    tf_sd=itrf2008_to_gda94_sd)

itrf2005_to_gda94 = Transformation(
    from_datum='ITRF2005', to_datum='GDA94', ref_epoch=date(1994, 1, 1),
    tx=-0.07973, ty=-0.00686, tz=0.03803,
    sc=0.006636,
    rx=-0.0000351, ry=0.0021211, rz=0.0021411,
    d_tx=0.00225, d_ty=-0.00062, d_tz=-0.00056,
    d_sc=0.000294,
    d_rx=0.0014707, d_ry=0.0011443, d_rz=0.0011701,
    tf_sd=itrf2005_to_gda94_sd)

itrf2000_to_gda94 = Transformation(
    from_datum='ITRF2000', to_datum='GDA94', ref_epoch=date(1994, 1, 1),
    tx=-0.04591, ty=-0.02985, tz=-0.02037,
    sc=0.00707,
    rx=-0.0016705, ry=0.0004594, rz=0.0019356,
    d_tx=-0.00466, d_ty=0.00355, d_tz=0.01124,
    d_sc=0.000249,
    d_rx=0.0017454, d_ry=0.0014868, d_rz=0.001224,
    tf_sd=itrf2000_to_gda94_sd)

itrf97_to_gda94 = Transformation(
    from_datum='ITRF97', to_datum='GDA94', ref_epoch=date(1994, 1, 1),
    tx=-0.01463, ty=-0.02762, tz=-0.02532,
    sc=0.006695,
    rx=-0.0017893, ry=-0.0006047, rz=0.0009962,
    d_tx=-0.00860, d_ty=0.00036, d_tz=0.01125,
    d_sc=0.000007,
    d_rx=0.0016394, d_ry=0.0015198, d_rz=0.0013801,
    tf_sd=itrf97_to_gda94_sd)

itrf96_to_gda94 = Transformation(
    from_datum='ITRF96', to_datum='GDA94', ref_epoch=date(1994, 1, 1),
    tx=0.02454, ty=-0.03643, tz=-0.06812,
    sc=0.006901,
    rx=-0.0027359, ry=-0.0020431, rz=0.0003731,
    d_tx=-0.02180, d_ty=0.00471, d_tz=0.02627,
    d_sc=0.000388,
    d_rx=0.0020203, d_ry=0.0021735, d_rz=0.0016290,
    tf_sd=itrf96_to_gda94_sd)

agd84_to_gda94 = Transformation(
    from_datum='AGD84', to_datum='GDA94', ref_epoch=0,
    tx=-117.763, ty=-51.510, tz=139.061,
    sc=-0.191,
    rx=-0.292, ry=-0.443, rz=-0.277)

agd66_to_gda94 = Transformation(
    from_datum='AGD66', to_datum='GDA94', ref_epoch=0,
    tx=-117.808, ty=-51.536, tz=137.784,
    sc=-0.290,
    rx=-0.303, ry=-0.446, rz=-0.234)

agd66_to_gda94_act = Transformation(
    from_datum='AGD66', to_datum='GDA94', ref_epoch=0,
    tx=-129.193, ty=-41.212, tz=130.730,
    sc=-2.955,
    rx=-0.246, ry=-0.374, rz=-0.329)

agd66_to_gda94_tas = Transformation(
    from_datum='AGD66', to_datum='GDA94', ref_epoch=0,
    tx=-120.271, ty=-64.543, tz=161.632,
    sc=2.499,
    rx=-0.217, ry=0.067, rz=0.129)

agd66_to_gda94_vicnsw = Transformation(
    from_datum='AGD66', to_datum='GDA94', ref_epoch=0,
    tx=-119.353, ty=-48.301, tz=139.484,
    sc=-0.613,
    rx=-0.415, ry=-0.260, rz=-0.437)

agd66_to_gda94_nt = Transformation(
    from_datum='AGD66', to_datum='GDA94', ref_epoch=0,
    tx=-124.133, ty=-42.003, tz=137.400,
    sc=-1.854,
    rx=0.008, ry=-0.557, rz=-0.178)

gda94_to_itrf2008 = -itrf2008_to_gda94
gda94_to_itrf2005 = -itrf2005_to_gda94
gda94_to_itrf2000 = -itrf2000_to_gda94
gda94_to_itrf97 = -itrf97_to_gda94
gda94_to_itrf96 = -itrf96_to_gda94
gda94_to_agd84 = -agd84_to_gda94
gda94_to_agd66 = -agd66_to_gda94
gda94_to_agd66_act = -agd66_to_gda94_act
gda94_to_agd66_tas = -agd66_to_gda94_tas
gda94_to_agd66_vicnsw = -agd66_to_gda94_vicnsw
gda94_to_agd66_nt = -agd66_to_gda94_nt

# ITRF2020 parameters
# Note: These parameter definitions use the units, rotation and delta-rotation
# sign convention used by the IERS. Units are converted and rotations flipped in
# the transformation object to allow compatibility with GeodePy functions.
# Ref: https://itrf.ign.fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt

itrf2020_to_itrf2014 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF2014', ref_epoch=date(2015, 1, 1),
    tx=-1.4, ty=-0.9, tz=1.4,
    sc=-0.42,
    rx=0.0, ry=0.0, rz=0.0,
    d_tx=0.0, d_ty=-0.1, d_tz=0.2,
    d_sc=0.0,
    d_rx=0.0, d_ry=0.0, d_rz=0.0)

itrf2020_to_itrf2008 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF2008', ref_epoch=date(2015, 1, 1),
    tx=0.2, ty=1.0, tz=3.3,
    sc=-0.29,
    rx=0.0, ry=0.0, rz=0.0,
    d_tx=0.0, d_ty=-0.1, d_tz=0.1,
    d_sc=0.03,
    d_rx=0.0, d_ry=0.0, d_rz=0.0)

itrf2020_to_itrf2005 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF2005', ref_epoch=date(2015, 1, 1),
    tx=2.7, ty=0.1, tz=-1.4,
    sc=0.65,
    rx=0.0, ry=0.0, rz=0.0,
    d_tx=0.3, d_ty=-0.1, d_tz=0.1,
    d_sc=0.03,
    d_rx=0.0, d_ry=0.0, d_rz=0.0)

itrf2020_to_itrf2000 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF2000', ref_epoch=date(2015, 1, 1),
    tx=-0.2, ty=0.8, tz=-34.2,
    sc=2.25,
    rx=0.0, ry=0.0, rz=0.0,
    d_tx=0.1, d_ty=0.0, d_tz=-1.7,
    d_sc=0.11,
    d_rx=0.0, d_ry=0.0, d_rz=0.0)

itrf2020_to_itrf97 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF97', ref_epoch=date(2015, 1, 1),
    tx=6.5, ty=-3.9, tz=-77.9,
    sc=3.98,
    rx=0.0, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)

itrf2020_to_itrf96 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF96', ref_epoch=date(2015, 1, 1),
    tx=6.5, ty=-3.9, tz=-77.9,
    sc=3.98,
    rx=0.0, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)

itrf2020_to_itrf94 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF94', ref_epoch=date(2015, 1, 1),
    tx=6.5, ty=-3.9, tz=-77.9,
    sc=3.98,
    rx=0.0, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)

itrf2020_to_itrf93 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF93', ref_epoch=date(2015, 1, 1),
    tx=-65.8, ty=1.9, tz=-71.3,
    sc=4.47,
    rx=-3.36, ry=-4.33, rz=0.75,
    d_tx=-2.8, d_ty=-0.2, d_tz=-2.3,
    d_sc=0.12,
    d_rx=-0.11, d_ry=-0.19, d_rz=0.07)

itrf2020_to_itrf92 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF92', ref_epoch=date(2015, 1, 1),
    tx=14.5, ty=-1.9, tz=-85.9,
    sc=3.27,
    rx=0.0, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)

itrf2020_to_itrf91 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF91', ref_epoch=date(2015, 1, 1),
    tx=26.5, ty=12.1, tz=-91.9,
    sc=4.67,
    rx=0.0, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)

itrf2020_to_itrf90 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF90', ref_epoch=date(2015, 1, 1),
    tx=24.5, ty=8.1, tz=-107.9,
    sc=4.97,
    rx=0.0, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)

itrf2020_to_itrf89 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF89', ref_epoch=date(2015, 1, 1),
    tx=29.5, ty=32.1, tz=-145.9,
    sc=8.37,
    rx=0.0, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)

itrf2020_to_itrf88 = iers2trans(
    itrf_from='ITRF2020', itrf_to='ITRF88', ref_epoch=date(2015, 1, 1),
    tx=24.5, ty=-3.9, tz=-169.9,
    sc=11.47,
    rx=0.1, ry=0.0, rz=0.36,
    d_tx=0.1, d_ty=-0.6, d_tz=-3.1,
    d_sc=0.12,
    d_rx=0.0, d_ry=0.0, d_rz=0.02)
    
itrf2014_to_itrf2020 = -itrf2020_to_itrf2014
itrf2008_to_itrf2020 = -itrf2020_to_itrf2008
itrf2005_to_itrf2020 = -itrf2020_to_itrf2005
itrf2000_to_itrf2020 = -itrf2020_to_itrf2000
itrf97_to_itrf2020 = -itrf2020_to_itrf97
itrf96_to_itrf2020 = -itrf2020_to_itrf96
itrf94_to_itrf2020 = -itrf2020_to_itrf94
itrf93_to_itrf2020 = -itrf2020_to_itrf93
itrf92_to_itrf2020 = -itrf2020_to_itrf92
itrf91_to_itrf2020 = -itrf2020_to_itrf91
itrf90_to_itrf2020 = -itrf2020_to_itrf90
itrf89_to_itrf2020 = -itrf2020_to_itrf89
itrf88_to_itrf2020 = -itrf2020_to_itrf88
    

# ITRF2014 parameters
# Note: These parameter definitions use the units, rotation and delta-rotation
# sign convention used by the IERS. Units are converted and rotations flipped in
# the transformation object to allow compatibility with GeodePy functions.
# For more information, see the GDA2020 tech manual section 2.2.1 available
# here: https://www.icsm.gov.au/gda2020-and-gda94-technical-manuals
# Ref: http://itrf.ign.fr/doc_ITRF/Transfo-ITRF2014_ITRFs.txt

itrf2014_to_itrf2008 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF2008', ref_epoch=date(2010, 1, 1),
    tx=1.6, ty=1.9, tz=2.4,
    sc=-0.02,
    rx=0, ry=0, rz=0,
    d_tx=0.0, d_ty=0.0, d_tz=-0.1,
    d_sc=0.03,
    d_rx=0, d_ry=0, d_rz=0)

itrf2014_to_itrf2005 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF2005', ref_epoch=date(2010, 1, 1),
    tx=2.6, ty=1.0, tz=-2.3,
    sc=0.92,
    rx=0, ry=0, rz=0,
    d_tx=0.3, d_ty=0.0, d_tz=-0.1,
    d_sc=0.03,
    d_rx=0, d_ry=0, d_rz=0)

itrf2014_to_itrf2000 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF2000', ref_epoch=date(2010, 1, 1),
    tx=0.7, ty=1.2, tz=-26.1,
    sc=2.12,
    rx=0, ry=0, rz=0,
    d_tx=0.1, d_ty=0.1, d_tz=-1.9,
    d_sc=0.11,
    d_rx=0, d_ry=0, d_rz=0)

itrf2014_to_itrf97 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF97', ref_epoch=date(2010, 1, 1),
    tx=7.4, ty=-0.5, tz=-62.8,
    sc=3.80,
    rx=0, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2014_to_itrf96 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF96', ref_epoch=date(2010, 1, 1),
    tx=7.4, ty=-0.5, tz=-62.8,
    sc=3.80,
    rx=0, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2014_to_itrf94 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF94', ref_epoch=date(2010, 1, 1),
    tx=7.4, ty=-0.5, tz=-62.8,
    sc=3.80,
    rx=0, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2014_to_itrf93 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF93', ref_epoch=date(2010, 1, 1),
    tx=-50.4, ty=3.3, tz=-60.2,
    sc=4.29,
    rx=-2.81, ry=-3.38, rz=0.40,
    d_tx=-2.8, d_ty=-0.1, d_tz=-2.5,
    d_sc=0.12,
    d_rx=-0.11, d_ry=-0.19, d_rz=0.07)

itrf2014_to_itrf92 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF92', ref_epoch=date(2010, 1, 1),
    tx=15.4, ty=1.5, tz=-70.8,
    sc=3.09,
    rx=0, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2014_to_itrf91 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF91', ref_epoch=date(2010, 1, 1),
    tx=27.4, ty=15.5, tz=-76.8,
    sc=4.49,
    rx=0, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2014_to_itrf90 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF90', ref_epoch=date(2010, 1, 1),
    tx=25.4, ty=11.5, tz=-92.8,
    sc=4.79,
    rx=0, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2014_to_itrf89 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF89', ref_epoch=date(2010, 1, 1),
    tx=30.4, ty=35.5, tz=-130.8,
    sc=8.19,
    rx=0, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2014_to_itrf88 = iers2trans(
    itrf_from='ITRF2014', itrf_to='ITRF88', ref_epoch=date(2010, 1, 1),
    tx=25.4, ty=-0.5, tz=-154.8,
    sc=11.29,
    rx=0.1, ry=0, rz=0.26,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.3,
    d_sc=0.12,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf2014 = -itrf2014_to_itrf2008
itrf2005_to_itrf2014 = -itrf2014_to_itrf2005
itrf2000_to_itrf2014 = -itrf2014_to_itrf2000
itrf97_to_itrf2014 = -itrf2014_to_itrf97
itrf96_to_itrf2014 = -itrf2014_to_itrf96
itrf94_to_itrf2014 = -itrf2014_to_itrf94
itrf93_to_itrf2014 = -itrf2014_to_itrf93
itrf92_to_itrf2014 = -itrf2014_to_itrf92
itrf91_to_itrf2014 = -itrf2014_to_itrf91
itrf90_to_itrf2014 = -itrf2014_to_itrf90
itrf89_to_itrf2014 = -itrf2014_to_itrf89
itrf88_to_itrf2014 = -itrf2014_to_itrf88


# ITRF2008 parameters
# Note: These parameter definitions use the units, rotation and delta-rotation
# sign convention used by the IERS. Units are converted and rotations flipped in
# the transformation object to allow compatibility with GeodePy functions.
# For more information, see the GDA2020 tech manual section 2.2.1 available
# here: https://www.icsm.gov.au/gda2020-and-gda94-technical-manuals
# Ref: http://itrf.ign.fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt

itrf2008_to_itrf2005 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF2005', ref_epoch=date(2000, 1, 1),
    tx=-2.0, ty=-0.9, tz=-4.7,
    sc=0.94,
    rx=0, ry=0, rz=0,
    d_tx=0.3, d_ty=0.0, d_tz=0.0,
    d_sc=0.0,
    d_rx=0, d_ry=0, d_rz=0)

itrf2008_to_itrf2000 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF2000', ref_epoch=date(2000, 1, 1),
    tx=-1.9, ty=-1.7, tz=-10.5,
    sc=1.34,
    rx=0, ry=0, rz=0,
    d_tx=0.1, d_ty=0.1, d_tz=-1.8,
    d_sc=0.08,
    d_rx=0, d_ry=0, d_rz=0)

itrf2008_to_itrf97 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF97', ref_epoch=date(2000, 1, 1),
    tx=4.8, ty=2.6, tz=-33.2,
    sc=2.92,
    rx=0, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf96 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF96', ref_epoch=date(2000, 1, 1),
    tx=4.8, ty=2.6, tz=-33.2,
    sc=2.92,
    rx=0, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf94 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF94', ref_epoch=date(2000, 1, 1),
    tx=4.8, ty=2.6, tz=-33.2,
    sc=2.92,
    rx=0, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf93 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF93', ref_epoch=date(2000, 1, 1),
    tx=4-24.0, ty=2.4, tz=-38.6,
    sc=3.41,
    rx=-1.71, ry=-1.48, rz=-0.30,
    d_tx=-2.8, d_ty=-0.1, d_tz=-2.4,
    d_sc=0.09,
    d_rx=-0.11, d_ry=-0.19, d_rz=0.07)

itrf2008_to_itrf92 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF92', ref_epoch=date(2000, 1, 1),
    tx=12.8, ty=4.6, tz=-41.2,
    sc=2.21,
    rx=0, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf91 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF91', ref_epoch=date(2000, 1, 1),
    tx=24.8, ty=18.6, tz=-47.2,
    sc=3.61,
    rx=0, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf90 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF90', ref_epoch=date(2000, 1, 1),
    tx=22.8, ty=14.6, tz=-63.2,
    sc=3.91,
    rx=0, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf89 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF89', ref_epoch=date(2000, 1, 1),
    tx=27.8, ty=38.6, tz=-101.2,
    sc=7.31,
    rx=0, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2008_to_itrf88 = iers2trans(
    itrf_from='ITRF2008', itrf_to='ITRF88', ref_epoch=date(2000, 1, 1),
    tx=22.8, ty=2.6, tz=-125.2,
    sc=10.41,
    rx=0.10, ry=0, rz=0.06,
    d_tx=0.1, d_ty=-0.5, d_tz=-3.2,
    d_sc=0.09,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2005_to_itrf2008 = -itrf2008_to_itrf2005
itrf2000_to_itrf2008 = -itrf2008_to_itrf2000
itrf97_to_itrf2008 = -itrf2008_to_itrf97
itrf96_to_itrf2008 = -itrf2008_to_itrf96
itrf94_to_itrf2008 = -itrf2008_to_itrf94
itrf93_to_itrf2008 = -itrf2008_to_itrf93
itrf92_to_itrf2008 = -itrf2008_to_itrf92
itrf91_to_itrf2008 = -itrf2008_to_itrf91
itrf90_to_itrf2008 = -itrf2008_to_itrf90
itrf89_to_itrf2008 = -itrf2008_to_itrf89
itrf88_to_itrf2008 = -itrf2008_to_itrf88


# ITRF2005 parameters
# Note: These parameter definitions use the units, rotation and delta-rotation
# sign convention used by the IERS. Units are converted and rotations flipped in
# the transformation object to allow compatibility with GeodePy functions.
# For more information, see the GDA2020 tech manual section 2.2.1 available
# here: https://www.icsm.gov.au/gda2020-and-gda94-technical-manuals
# Ref: http://itrf.ensg.ign.fr/ITRF_solutions/2005/tp_05-00.php

itrf2005_to_itrf2000 = iers2trans(
    itrf_from='ITRF2005', itrf_to='ITRF2000', ref_epoch=date(2000, 1, 1),
    tx=0.1, ty=-0.8, tz=-5.8,
    sc=0.40,
    rx=0, ry=0, rz=0,
    d_tx=-0.2, d_ty=0.1, d_tz=-1.8,
    d_sc=0.08,
    d_rx=0, d_ry=0, d_rz=0)

itrf2000_to_itrf2005 = -itrf2005_to_itrf2000

# ITRF2000 parameters
# ftp://ftp.iers.org/products/reference-systems/terrestrial/itrf/ITRF.TP
# NOTE: This ref lists translations in centimetres. All other ITRF
# transformations are shown in millimetres.
# NOTE: All translations and rates of translation shown below have been
# converted to millimetres.
itrf2000_to_itrf97 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF97', ref_epoch=date(1997, 1, 1),
    tx=6.7, ty=6.1, tz=-18.5,
    sc=1.55,
    rx=0, ry=0, rz=0,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2000_to_itrf96 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF96', ref_epoch=date(1997, 1, 1),
    tx=6.7, ty=6.1, tz=-18.5,
    sc=1.55,
    rx=0, ry=0, rz=0,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2000_to_itrf94 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF94', ref_epoch=date(1997, 1, 1),
    tx=6.7, ty=6.1, tz=-18.5,
    sc=1.55,
    rx=0, ry=0, rz=0,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2000_to_itrf93 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF93', ref_epoch=date(1988, 1, 1),
    tx=12.7, ty=6.5, tz=-20.9,
    sc=1.95,
    rx=-0.39, ry=0.80, rz=-1.14,
    d_tx=-2.9, d_ty=-0.2, d_tz=-0.6,
    d_sc=0.01,
    d_rx=-0.11, d_ry=-0.19, d_rz=0.07)

itrf2000_to_itrf92 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF92', ref_epoch=date(1988, 1, 1),
    tx=14.7, ty=13.5, tz=-13.9,
    sc=0.75,
    rx=0, ry=0, rz=-0.18,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2000_to_itrf91 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF91', ref_epoch=date(1988, 1, 1),
    tx=26.7, ty=27.5, tz=-19.9,
    sc=2.15,
    rx=0, ry=0, rz=-0.18,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2000_to_itrf90 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF90', ref_epoch=date(1988, 1, 1),
    tx=24.7, ty=23.5, tz=-35.9,
    sc=2.45,
    rx=0, ry=0, rz=-0.18,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2000_to_itrf89 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF89', ref_epoch=date(1988, 1, 1),
    tx=29.7, ty=47.5, tz=-73.9,
    sc=5.85,
    rx=0, ry=0, rz=-0.18,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf2000_to_itrf88 = iers2trans(
    itrf_from='ITRF2000', itrf_to='ITRF88', ref_epoch=date(1988, 1, 1),
    tx=24.7, ty=11.5, tz=-97.9,
    sc=8.95,
    rx=0.10, ry=0, rz=-0.18,
    d_tx=0.0, d_ty=-0.6, d_tz=-1.4,
    d_sc=0.01,
    d_rx=0, d_ry=0, d_rz=0.02)

itrf97_to_itrf2000 = -itrf2000_to_itrf97
itrf96_to_itrf2000 = -itrf2000_to_itrf96
itrf94_to_itrf2000 = -itrf2000_to_itrf94
itrf93_to_itrf2000 = -itrf2000_to_itrf93
itrf92_to_itrf2000 = -itrf2000_to_itrf92
itrf91_to_itrf2000 = -itrf2000_to_itrf91
itrf90_to_itrf2000 = -itrf2000_to_itrf90
itrf89_to_itrf2000 = -itrf2000_to_itrf89
itrf88_to_itrf2000 = -itrf2000_to_itrf88


# The locations of files used in the height module
aws_server = '/vsicurl/https://geoid.s3-ap-southeast-2.amazonaws.com/'
file_DOV_PV = aws_server + 'AGQG/DOV_PV.tif'
file_DOV_PM = aws_server + 'AGQG/DOV_PM.tif'
file_AG2020 = aws_server + 'AUSGeoid/AUSGeoid2020_RELEASEV20170908.tif'
file_AG2020_STD = aws_server + 'AUSGeoid/AUSGeoid2020_RELEASEV20170908_err.tif'
file_AVWS = aws_server + 'AGQG/AGQG_20201120.tif'
file_AVWS_STD = aws_server + 'AGQG/AGQG_uncertainty_20201120.tif'
file_GRAV_BA = aws_server + 'AGQG/Bouguer_Grav_RELEASE20191107.tif'
file_AG98=aws_server+'AUSGeoid/AUSGeoid98.tif'
file_AG09=aws_server+'AUSGeoid/AUSGeoid09_V1.01.tif'
file_AG98_DOV_PV=aws_server+'AUSGeoid/AUSGeoid98_DOV_PV.tif'
file_AG98_DOV_PM=aws_server+'AUSGeoid/AUSGeoid98_DOV_PM.tif'
file_AG09_DOV_PV=aws_server+'AUSGeoid/AUSGeoid09_DOV_PV_V1.01.tif'
file_AG09_DOV_PM=aws_server+'AUSGeoid/AUSGeoid09_DOV_PM_V1.01.tif'


# GRS80 normal gravity flattening [Moritz, 2000 Section 4]
grs80_ngf = 0.005302440112

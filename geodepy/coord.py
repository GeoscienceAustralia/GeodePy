#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Coordinate Module
"""

from geodepy.constants import Projection, utm, grs80
from geodepy.angles import DECAngle, angular_typecheck
from geodepy.convert import xyz2llh, llh2xyz


class CoordCart(object):
    """
    Cartesian Coordinate Class
    Used for working with coordinates representing points in a 3 dimensional
    earth-centred earth-fixed system. N Value (optional) is the distance (m)
    at the coordinate that a reference surface (typically a geoid) is above
    or below an ellipsoid
    """
    def __init__(self, xaxis, yaxis, zaxis, nval=None):
        """
        :param xaxis: Distance (m) along X Axis (positive when passing through
        intersection of equator and prime meridian)
        :type xaxis: float
        :param yaxis: Distance (m) along Y Axis (positive when passing through
        intersection of equator and 90 degrees east longitude)
        :type yaxis: float
        :param zaxis: Distance (m) along Z Axis (positive when passing through
        north pole)
        :type zaxis: float
        :param nval: Distance (m) between a reference surface (typically a
        geoid) and an ellipsoid (positive when surface is above ellipsoid)
        :type nval: float
        """
        self.xaxis = float(xaxis)
        self.yaxis = float(yaxis)
        self.zaxis = float(zaxis)
        if not nval:
            self.nval = None
        else:
            self.nval = float(nval)

    def __repr__(self):
        return (f'CoordCart: X: {self.xaxis} Y: {self.yaxis} Z: {self.zaxis} '
                f'NVal: {self.nval}')

    def __eq__(self, other):
        if isinstance(other, CoordCart):
            return vars(self) == vars(other)
        else:
            raise ValueError(f"Can't compare {self} to {other}. If other object"
                             f" is coord, convert to CoordCart first.")

    def __round__(self, n=None):
        if self.nval is None:
            return CoordCart(round(self.xaxis, n), round(self.yaxis, n),
                             round(self.zaxis, n))
        else:
            return CoordCart(round(self.xaxis, n), round(self.yaxis, n),
                             round(self.zaxis, n), round(self.nval, n))

    def geo(self, ellipsoid=grs80):
        """
        Convert coordinates to Geographic
        Note: If no N Value, no Orthometric Height set.
        :param ellipsoid: geodepy.constants.Ellipsoid Object (default: grs80)
        :return: Geographic Coordinate
        :rtype: CoordGeo
        """
        lat, lon, ell_ht = xyz2llh(self.xaxis, self.yaxis,
                                   self.zaxis, ellipsoid)
        if self.nval is None:
            return CoordGeo(DECAngle(lat), DECAngle(lon), ell_ht)
        else:
            return CoordGeo(DECAngle(lat), DECAngle(lon), ell_ht,
                            ell_ht + self.nval)


class CoordGeo(object):
    """
    Geographic Coordinate Class
    Used for working with coordinates representing points in an ellipsoidal
    system with ellipsoid heights (m) relative to the surface of the ellipsoid.
    Orthometric heights (m) are relative to a difference reference surface
    (typically a geoid).
    """
    def __init__(self, lat, lon, ell_ht=None, orth_ht=None):
        """
        :param lat: Geographic Latitude (angle between +180 and -180 degrees,
        positive values are east of central meridian)
        :type lat: float (decimal degrees) or any Angle object in geodepy.angle
        :param lon: Geographic Longitude (angle between +90 and -90 degrees,
        positive values are north of equator)
        :type lon: float (decimal degrees) or any Angle object in geodepy.angle
        :param ell_ht: Ellipsoid Height (m, positive is outside ellipsoid)
        :type ell_ht: float
        :param orth_ht: Orthometric Height (m)
        :type orth_ht: float
        """
        self.lat = angular_typecheck(lat)
        self.lon = angular_typecheck(lon)
        if ell_ht is None:
            self.ell_ht = None
        else:
            self.ell_ht = float(ell_ht)
        if orth_ht is None:
            self.orth_ht = None
        else:
            self.orth_ht = float(orth_ht)

    def __repr__(self):
        return (f'CoordGeo: Lat: {self.lat} Lon: {self.lon} '
                f'Ell_Ht: {self.ell_ht} Orth_Ht: {self.orth_ht}')

    def __eq__(self, other):
        if isinstance(other, CoordGeo):
            return vars(self) == vars(other)
        else:
            raise ValueError(f"Can't compare {self} to {other}. If other object"
                             f" is coord, convert to CoordGeo first.")

    def __round__(self, n=None):
        if self.ell_ht is None and self.orth_ht is None:
            return CoordGeo(round(self.lat, n), round(self.lon, n))
        elif self.ell_ht is None:
            return CoordGeo(round(self.lat, n), round(self.lon, n),
                            None, round(self.orth_ht, n))
        elif self.orth_ht is None:
            return CoordGeo(round(self.lat, n), round(self.lon, n),
                            round(self.ell_ht, n), None)
        else:
            return CoordGeo(round(self.lat, n), round(self.lon, n),
                            round(self.ell_ht, n), round(self.orth_ht, n))

    def cart(self, ellipsoid=grs80):
        """
        Convert coordinates to Cartesian
        Note: If no ellipsoid height set, uses 0m. No N Value output
        :param ellipsoid: geodepy.constants.Ellipsoid Object (default: grs80)
        :return: Cartesian Coordinate
        :rtype: CoordCart
        """
        if self.ell_ht:
            x, y, z = llh2xyz(self.lat, self.lon, self.ell_ht, ellipsoid)
            if self.orth_ht:  # Only N Value if both Ellipsoid and Ortho Heights
                return CoordCart(x, y, z, self.orth_ht - self.ell_ht)
            else:  # No Ortho Height -> No N Value
                return CoordCart(x, y, z)
        else:  # No Ellipsoid Height - set to 0m for conversion
            x, y, z = llh2xyz(self.lat, self.lon, 0, ellipsoid)
            return CoordCart(x, y, z)  # No Ellipsoid Height -> No N Value


class CoordTM(object):
    """
    Transverse Mercator Coordinate Class
    Used for working with coordinates representing points in a Transverse
    Mercator projected planar system. Coordinates are related to an ellipsoid
    via a projection (default is Universal Transverse Mercator, see
    geodepy.constants.Projection for more details). Ellipsoid heights (m)
    are relative to the surface of the ellipsoid. Orthometric heights (m) are
    relative to a difference reference surface (typically a geoid).
    """
    def __init__(self, zone, east, north, ell_ht=None, orth_ht=None,
                 hemi_north=False, projection=utm):
        """
        :param zone: Transverse Mercator Zone Number - 1 to 60
        :type zone: int
        :param east: Easting (m, within 3330km of Central Meridian)
        :type east: float
        :param north: Northing (m, 0 to 10,000,000m)
        :type north: float
        :param ell_ht: Ellipsoid Height (m, positive is outside ellipsoid)
        :type ell_ht: float
        :param orth_ht: Orthometric Height (m)
        :type orth_ht: float
        :param hemi_north: True if coordinate in Northern Hemisphere, False if
        coordinate in Southern Hemisphere (default)
        :type hemi_north: bool
        :param projection: Information defining relationship between projected
        coordinate and the ellipsoid
        :type projection: geodepy.constants.Projection object
        (default: geodepy.constants.utm)
        """
        self.zone = int(zone)
        self.east = float(east)
        self.north = float(north)
        self.ell_ht = float(ell_ht)
        self.orth_ht = float(orth_ht)
        if isinstance(hemi_north, bool):
            self.hemi_north = hemi_north
        else:
            raise TypeError(f'CoordTM: invalid hemi_north, must be bool (True '
                            f'for Northern Hemisphere, False for Southern). '
                            f'hemi_north entered: {hemi_north}')
        if isinstance(projection, Projection):
            self.projection = projection
        else:
            raise TypeError(f'CoordTM: invalid projection, must be geodepy.'
                            f'constants.Projection object (default is geodepy.'
                            f'constants.utm. projection entered: {projection}')

    def __repr__(self):
        if self.hemi_north:
            return (f'CoordTM: Zone: {self.zone} East: {self.east} '
                    f'North: {self.north} Ell_Ht: {self.ell_ht} '
                    f'Orth_Ht: {self.orth_ht} Hemisphere: North')
        else:  # not self.hemi_north
            return (f'CoordTM: Zone: {self.zone} East: {self.east} '
                    f'North: {self.north} Ell_Ht: {self.ell_ht} '
                    f'Orth_Ht: {self.orth_ht} Hemisphere: South')

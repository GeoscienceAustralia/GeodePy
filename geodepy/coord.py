#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Coordinate Module
"""

from geodepy.constants import Projection, utm, grs80
from geodepy.angles import (DECAngle, HPAngle, GONAngle, DMSAngle, DDMAngle,
                            dec2hpa, dec2gona, dec2dms, dec2ddm,
                            angular_typecheck)
from geodepy.convert import xyz2llh, llh2xyz, grid2geo, geo2grid


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

    def geo(self, ellipsoid=grs80, notation=DECAngle):
        """
        Convert coordinates to Geographic
        Note: If no N Value, no Orthometric Height set.
        :param ellipsoid: geodepy.constants.Ellipsoid Object (default: grs80)
        :param notation: Latitude and Longitude Angle Notation format
        :type notation: geodepy.angle class or float
        :return: Geographic Coordinate
        :rtype: CoordGeo
        """
        lat, lon, ell_ht = xyz2llh(self.xaxis, self.yaxis,
                                   self.zaxis, ellipsoid)
        if notation is DECAngle:
            lat = DECAngle(lat)
            lon = DECAngle(lon)
        elif notation is HPAngle:
            lat = DECAngle(lat).hpa()
            lon = DECAngle(lon).hpa()
        elif notation is GONAngle:
            lat = DECAngle(lat).gona()
            lon = DECAngle(lon).gona()
        elif notation is DMSAngle:
            lat = DECAngle(lat).dms()
            lon = DECAngle(lon).dms()
        elif notation is DDMAngle:
            lat = DECAngle(lat).ddm()
            lon = DECAngle(lon).ddm()
        elif notation is float:
            pass  # geodepy.convert.grid2geo returns float dec degrees
        else:
            raise ValueError(f'CoordCart.geo() notation requires class float or'
                             f' class from geodepy.angles module. '
                             f'Supplied: {notation}')
        if self.nval is None:
            return CoordGeo(lat, lon, ell_ht)
        else:
            return CoordGeo(lat, lon, ell_ht,
                            ell_ht - self.nval)

    # TODO: Add functionality to utilise different TM projections

    def tm(self, ellipsoid=grs80, projection=utm):
        """
        Convert coordinates to Transverse Mercator
        :param ellipsoid: geodepy.constants.Ellipsoid Object (default: grs80)
        :param projection: geodepy.constants.Projection Object (default: utm)
        :return: Transverse Mercator Projection Coordinate
        :rtype: CoordTM
        """
        return self.geo(ellipsoid).tm(ellipsoid, projection)


class CoordGeo(object):
    """
    Geographic Coordinate Class
    Used for working with coordinates representing points in an ellipsoidal
    system with ellipsoid heights (m) relative to the surface of the ellipsoid.
    Orthometric heights (m) are relative to a different reference surface
    (typically a geoid).
    """
    def __init__(self, lat, lon, ell_ht=None, orth_ht=None):
        """
        :param lat: Geographic Latitude (angle between +90 and -90 degrees,
        positive values are north of equator)
        :type lat: float (decimal degrees) or any Angle object in geodepy.angles
        :param lon: Geographic Longitude (angle between +180 and -180 degrees,
        positive values are east of central meridian)
        :type lon: float (decimal degrees) or any Angle object in geodepy.angles
        :param ell_ht: Ellipsoid Height (m, positive is outside ellipsoid)
        :type ell_ht: float
        :param orth_ht: Orthometric Height (m)
        :type orth_ht: float
        """
        # Check latitude and longitude are supported types, both same type
        type_list = [float, DECAngle, HPAngle, GONAngle, DMSAngle, DDMAngle]
        if not all(x in type_list for x in [type(lat), type(lon)]):
            raise TypeError(f'CoordGeo Latitude and Longitude must be type: '
                            f'float (decimal degrees) or geodepy.angles class '
                            f'Lat: {type(lat)}, Lon: {type(lon)}')
        if type(lat) != type(lon):
            raise TypeError(f'CoordGeo Latitude and Longitude must have the '
                            f'same notation type: float (decimal degrees) or '
                            f'geodepy.angles class. Lat: {type(lat)}, '
                            f'Lon: {type(lon)}')
        self.lat = lat
        self.lon = lon
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

    def notation(self, notation):
        if type(self.lat) == float:  # Decimal Degrees (float)
            # Use functions to convert from Decimal Degrees (float)
            if notation == float:
                pass
            elif notation == DECAngle:
                new_lat = DECAngle(self.lat)
                new_lon = DECAngle(self.lon)
            elif notation == HPAngle:
                new_lat = dec2hpa(self.lat)
                new_lon = dec2hpa(self.lon)
            elif notation == GONAngle:
                new_lat = dec2gona(self.lat)
                new_lon = dec2gona(self.lon)
            elif notation == DMSAngle:
                new_lat = dec2dms(self.lat)
                new_lon = dec2dms(self.lon)
            elif notation == DDMAngle:
                new_lat = dec2ddm(self.lat)
                new_lon = dec2ddm(self.lon)
            else:
                raise ValueError(
                    f'CoordGeo.notation() notation requires class float or '
                    f'class from geodepy.angles module. '
                    f'Supplied: {notation}')
        elif type(self.lat) in [DECAngle, HPAngle, GONAngle,
                                DMSAngle, DDMAngle]:
            # Use methods to convert from geodepy.angles classes
            if notation == float:
                new_lat = self.lat.dec()
                new_lon = self.lon.dec()
            elif notation == DECAngle:
                new_lat = self.lat.deca()
                new_lon = self.lon.deca()
            elif notation == HPAngle:
                new_lat = self.lat.hpa()
                new_lon = self.lon.hpa()
            elif notation == GONAngle:
                new_lat = self.lat.gona()
                new_lon = self.lon.gona()
            elif notation == DMSAngle:
                new_lat = self.lat.dms()
                new_lon = self.lon.dms()
            elif notation == DDMAngle:
                new_lat = self.lat.ddm()
                new_lon = self.lon.ddm()
            else:
                raise ValueError(
                    f'CoordGeo.notation() notation requires class float or '
                    f'class from geodepy.angles module. '
                    f'Supplied: {notation}')
        return CoordGeo(new_lat, new_lon, self.ell_ht, self.orth_ht)

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
                return CoordCart(x, y, z, self.ell_ht - self.orth_ht)
            else:  # No Ortho Height -> No N Value
                return CoordCart(x, y, z)
        else:  # No Ellipsoid Height - set to 0m for conversion
            x, y, z = llh2xyz(self.lat, self.lon, 0, ellipsoid)
            return CoordCart(x, y, z)  # No Ellipsoid Height -> No N Value

    # TODO: Add functionality to utilise different TM projections

    def tm(self, ellipsoid=grs80, projection=utm):
        """
        Convert coordinates to Universal Transverse Mercator Projection
        Note: Heights are not used in calculation
        :param ellipsoid: geodepy.constants.Ellipsoid Object (default: grs80)
        :param projection: geodepy.constants.Projection Object (default: utm)
        :return: Transverse Mercator Projection Coordinate
        :rtype: CoordTM
        """
        hemi, zone, east, north, psf, gc = geo2grid(self.lat, self.lon,
                                                    0, ellipsoid)
        if hemi == 'North':
            hemi_north = True
        else:  # hemi == 'South'
            hemi_north = False

        return CoordTM(zone, east, north,
                       self.ell_ht, self.orth_ht,
                       hemi_north, projection)


class CoordTM(object):
    """
    Transverse Mercator Coordinate Class
    Used for working with coordinates representing points in a Transverse
    Mercator projected planar system. Coordinates are related to an ellipsoid
    via a projection (default is Universal Transverse Mercator, see
    geodepy.constants.Projection for more details). Ellipsoid heights (m)
    are relative to the surface of the ellipsoid. Orthometric heights (m) are
    relative to a different reference surface (typically a geoid).
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
        if ell_ht is None:
            self.ell_ht = None
        else:
            self.ell_ht = float(ell_ht)
        if orth_ht is None:
            self.orth_ht = None
        else:
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

    def __eq__(self, other):
        if isinstance(other, CoordTM):
            return vars(self) == vars(other)
        else:
            raise ValueError(f"Can't compare {self} to {other}. If other object"
                             f" is coord, convert to CoordTM first.")

    def __round__(self, n=None):
        if self.ell_ht is None and self.orth_ht is None:
            return CoordTM(self.zone, round(self.east, n), round(self.north, n),
                           None, None,
                           self.hemi_north, self.projection)
        elif self.ell_ht is None:
            return CoordTM(self.zone, round(self.east, n), round(self.north, n),
                           None, round(self.orth_ht, n),
                           self.hemi_north, self.projection)
        elif self.orth_ht is None:
            return CoordTM(self.zone, round(self.east, n), round(self.north, n),
                           round(self.ell_ht, n), None,
                           self.hemi_north, self.projection)
        else:
            return CoordTM(self.zone, round(self.east, n), round(self.north, n),
                           round(self.ell_ht, n), round(self.orth_ht, n),
                           self.hemi_north, self.projection)

    def geo(self, ellipsoid=grs80, notation=DECAngle):
        """

        :param ellipsoid: geodepy.constants.Ellipsoid Object (default: grs80)
        :param notation: Latitude and Longitude Angle Notation format
        :type notation: geodepy.angle class or float
        :return:
        """
        if self.hemi_north:
            hemi_str = 'north'
        else:
            hemi_str = 'south'
        lat, lon, psf, grid_conv = grid2geo(self.zone, self.east, self.north,
                                            hemi_str, ellipsoid)
        if notation is DECAngle:
            lat = DECAngle(lat)
            lon = DECAngle(lon)
        elif notation is HPAngle:
            lat = DECAngle(lat).hpa()
            lon = DECAngle(lon).hpa()
        elif notation is GONAngle:
            lat = DECAngle(lat).gona()
            lon = DECAngle(lon).gona()
        elif notation is DMSAngle:
            lat = DECAngle(lat).dms()
            lon = DECAngle(lon).dms()
        elif notation is DDMAngle:
            lat = DECAngle(lat).ddm()
            lon = DECAngle(lon).ddm()
        elif notation is float:
            pass  # geodepy.convert.grid2geo returns float dec degrees
        else:
            raise ValueError(f'CoordTM.geo() notation requires class float or '
                             f'class from geodepy.angles module. '
                             f'Supplied: {notation}')
        return CoordGeo(lat, lon, self.ell_ht, self.orth_ht)

    # TODO: Add functionality to utilise different TM projections

    def cart(self, ellipsoid=grs80):
        """
        Convert coordinates to Cartesian
        Note: If no ellipsoid height set, uses 0m. No N Value output
        :param ellipsoid: geodepy.constants.Ellipsoid Object (default: grs80)
        :return: Cartesian Coordinate
        :rtype: CoordCart
        """
        return self.geo(ellipsoid).cart(ellipsoid)

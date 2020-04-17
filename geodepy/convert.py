#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Convert Module
"""

from math import modf, sin, cos, atan2, radians, degrees, sqrt


class DMSAngle(object):
    def __init__(self, degree=0, minute=0, second=0.0):
        # Set sign of object based on sign of any variable
        if degree == 0:
            if str(degree)[0] == '-':
                self.sign = -1
            elif minute < 0:
                self.sign = -1
            elif second < 0:
                self.sign = -1
            else:
                self.sign = 1
        elif degree > 0:
            self.sign = 1
        else:  # degree < 0
            self.sign = -1
        self.degree = abs(int(degree))
        self.minute = abs(int(minute))
        self.second = abs(second)

    def __repr__(self):
        if self.sign == -1:
            signsymbol = '-'
        elif self.sign == 1:
            signsymbol = '+'
        else:
            signsymbol = 'error'
        return '{DMSAngle: ' + signsymbol + str(self.degree) + 'd ' + \
               str(self.minute) + 'm ' + str(self.second) + 's}'

    def __add__(self, other):
        return dec2dms(self.dec() + other.dec())

    def __radd__(self, other):
        return dec2dms(other.dec() + self.dec())

    def __sub__(self, other):
        return dec2dms(self.dec() - other.dec())

    def __rsub__(self, other):
        return dec2dms(other.dec() - self.dec())

    def __mul__(self, other):
        try:
            return dec2dms(self.dec() * other)
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __rmul__(self, other):
        try:
            return dec2dms(other * self.dec())
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __truediv__(self, other):
        try:
            return dec2dms(self.dec() / other)
        except TypeError:
            raise TypeError('Division only defined between DMSAngle Object '
                            'and Int or Float')

    def __abs__(self):
        return DMSAngle(self.degree, self.minute, self.second)

    def __neg__(self):
        if self.sign == 1:
            return DMSAngle(-self.degree, -self.minute, -self.second)
        else:  # sign == -1
            return DMSAngle(self.degree, self.minute, self.second)

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def __lt__(self, other):
        return self.dec() < other.dec()

    def __gt__(self, other):
        return self.dec() > other.dec()

    def dec(self):
        return self.sign * (self.degree + (self.minute / 60) +
                            (self.second / 3600))

    def hp(self):
        return self.sign * (self.degree + (self.minute / 100) +
                            (self.second / 10000))

    def ddm(self):
        if self.sign == 1:
            return DDMAngle(self.degree, self.minute + (self.second/60))
        else:  # sign == -1
            return -DDMAngle(self.degree, self.minute + (self.second/60))


class DDMAngle(object):
    def __init__(self, degree=0, minute=0.0):
        # Set sign of object based on sign of any variable
        if degree == 0:
            if str(degree)[0] == '-':
                self.sign = -1
            elif minute < 0:
                self.sign = -1
            else:
                self.sign = 1
        elif degree > 0:
            self.sign = 1
        else:  # degree < 0
            self.sign = -1
        self.degree = abs(int(degree))
        self.minute = abs(minute)

    def __repr__(self):
        if self.sign == -1:
            signsymbol = '-'
        elif self.sign == 1:
            signsymbol = '+'
        else:
            signsymbol = 'error'
        return '{DDMAngle: ' + signsymbol + str(self.degree) + 'd ' + \
               str(self.minute) + 'm}'

    def __add__(self, other):
        return dec2ddm(self.dec() + other.dec())

    def __sub__(self, other):
        return dec2ddm(self.dec() - other.dec())

    def __mul__(self, other):
        try:
            return dec2ddm(self.dec() * other)
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __rmul__(self, other):
        try:
            return dec2ddm(other * self.dec())
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __truediv__(self, other):
        try:
            return dec2ddm(self.dec() / other)
        except TypeError:
            raise TypeError('Division only defined between DMSAngle Object '
                            'and Int or Float')

    def __abs__(self):
        return DDMAngle(self.degree, self.minute)

    def __neg__(self):
        if self.sign == 1:
            return DDMAngle(-self.degree, -self.minute)
        else:  # sign == -1
            return DDMAngle(self.degree, self.minute)

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def __lt__(self, other):
        return self.dec() < other.dec()

    def __gt__(self, other):
        return self.dec() > other.dec()

    def dec(self):
        return self.sign * (self.degree + (self.minute / 60))

    def hp(self):
        minute_int, second = divmod(self.minute, 1)
        return self.sign * (self.degree + (minute_int / 100) +
                            (second * 0.006))

    def dms(self):
        minute_int, second = divmod(self.minute, 1)
        return self.sign * DMSAngle(self.degree, minute_int, second * 60)


class DNACoord(object):
    def __init__(self, pointid, const, easting, northing, zone, lat,
                 long, ortho_ht, ell_ht, x, y, z, x_sd, y_sd, z_sd, desc):
        self.pointid = pointid
        self.const = const
        self.easting = easting
        self.northing = northing
        self.zone = zone
        self.lat = lat
        self.long = long
        self.ortho_ht = ortho_ht
        self.ell_ht = ell_ht
        self.x = x
        self.y = y
        self.z = z
        self.x_sd = x_sd
        self.y_sd = y_sd
        self.z_sd = z_sd
        self.desc = desc

    def converthptodd(self):
        self.lat = hp2dec(self.lat)
        self.long = hp2dec(self.long)


def dec2hp(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    hp = degree + (minute / 100) + (second / 10000)
    return hp if dec >= 0 else -hp


def hp2dec(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    dec = degree + (minute / 60) + (second / 360)
    return dec if hp >= 0 else -dec


def dec2dms(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    return DMSAngle(degree, minute, second) if dec >= 0 \
        else DMSAngle(-degree, minute, second)


def dec2ddm(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    minute = minute + (second / 60)
    return DDMAngle(degree, minute) if dec >= 0 else DDMAngle(-degree, minute)


def hp2dms(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    return DMSAngle(degree, minute, second * 10) if hp >= 0 \
        else DMSAngle(-degree, minute, second * 10)


def hp2ddm(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    minute = minute + (second / 6)
    return DDMAngle(degree, minute) if hp >= 0 else DDMAngle(-degree, minute)


def polar2rect(r, theta):
    x = r * sin(radians(theta))
    y = r * cos(radians(theta))
    return x, y


def rect2polar(x, y):
    r = sqrt(x ** 2 + y ** 2)
    theta = atan2(x, y)
    if theta < 0:
        theta = degrees(theta) + 360
    else:
        theta = degrees(theta)
    return r, theta


def dd2sec(dd):
    minutes, seconds = divmod(abs(dd) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    sec = (degrees * 3600) + (minutes * 60) + seconds
    return sec if dd >= 0 else -sec


def dec2sex(lon, lat):
    """Convert decimal degrees to sexagesimal format

    Longitudes go from -180 to 180
    Latitudes go from -90 to 90
    """
    def fmt_sex(coord):
        """Function to create a sexagesimal-formatted coordinate as a string
        """
        coord = float(coord)
        if coord < 0:
            flag = -1
            coord = abs(coord)
        else:
            flag = 1
        min_sec, deg = modf(coord)
        deg = int(deg)
        deg *= flag  # deal with negatives
        sec, min = modf(min_sec * 60)
        min = round(min)
        sec *= 60
        sex_coord = '{} {:02d} {:05.2f}'.format(deg, min, sec)

        return sex_coord

    # Convert the coordinates
    sex_lon = fmt_sex(lon)
    sex_lat = fmt_sex(lat)

    return sex_lon, sex_lat


def sex2dec(lon, lat):
    """Convert a sexagesimal coordinate to decimal degrees

    Longitudes go from -180 to 180
    Latitudes go from -90 to 90
    """
    def fmt_dec(coord):
        """Function to piece together a decimal coordinate
        """
        coord = str(coord)
        if coord[:1] == '-':
            flag = -1
            coord = coord[1:]
        else:
            flag = 1
        deg, min, sec = coord.split()
        dec_coord = int(deg) + float(min) / 60 + float(sec) / 3600
        dec_coord = float('{:.6f}'.format(dec_coord))
        dec_coord *= flag  # deal with negatives

        return dec_coord

    # Convert the coordinates
    dec_lon = fmt_dec(lon)
    dec_lat = fmt_dec(lat)

    return dec_lon, dec_lat


def sex2hp(lon, lat):
    """Convert a sexagesimal coordinate to HP notation

    Longitudes go from -180 to 180
    Latitudes go from -90 to 90
    """
    def fmt_hp(coord):
        """Function to piece together a coordinate in HP notation
        """
        coord = str(coord)
        if coord[:1] == '-':
            flag = -1
            coord = coord[1:]
        else:
            flag = 1
        deg, min, sec = coord.split()
        sec = sec.replace('.', '')
        hp_coord = '{}.{:02d}{:04d}'.format(deg, int(min), int(sec))
        if flag == -1:
            hp_coord = '-' + hp_coord  # deal with negatives

        return hp_coord

    # Convert the coordinates
    hp_lon = fmt_hp(lon)
    hp_lat = fmt_hp(lat)

    return hp_lon, hp_lat


def hp2sex(lon, lat):
    """Convert HP notation to a sexagesimal coordinate

    Longitudes go from -180 to 180
    Latitudes go from -90 to 90
    """
    def fmt_sex(coord):
        """Function to piece together a coordinate in HP notation
        """
        coord = str(coord)
        if coord[:1] == '-':
            flag = -1
            coord = coord[1:]
        else:
            flag = 1
        deg, min, sec = coord.split()
        sec = sec.replace('.', '')
        hp_coord = '{}.{:02d}{:04d}'.format(deg, int(min), int(sec))
        if flag == -1:
            hp_coord = '-' + hp_coord  # deal with negatives

        return hp_coord

    # Convert the coordinates
    hp_lon = fmt_sex(lon)
    hp_lat = fmt_sex(lat)

    return hp_lon, hp_lat


def dd2dms_v(dd):
    minutes, seconds = divmod(abs(dd) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    dms = degrees + (minutes / 100) + (seconds / 10000)
    dms[dd <= 0] = -dms[dd <= 0]
    return dms


def dms2dd_v(dms):
    degmin, seconds = divmod(abs(dms) * 1000, 10)
    degrees, minutes = divmod(degmin, 100)
    dd = degrees + (minutes / 60) + (seconds / 360)
    dd[dms <= 0] = -dd[dms <= 0]
    return dd

#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Conversions Module
"""

from math import modf


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


class DMSAngle(object):
    def __init__(self, degree, minute, second):
        self.degree = int(degree)
        self.minute = abs(int(minute))
        self.second = abs(second)

    def __repr__(self):
        return '{DMSAngle: ' + str(self.degree) + 'd ' + str(self.minute) + 'm ' + str(self.second) + 's}'

    def __add__(self, other):
        return dec2dms(self.dec() + other.dec())

    def __sub__(self, other):
        return dec2dms(self.dec() - other.dec())

    def __mul__(self, other):
        try:
            return dec2dms(self.dec() * other)
        except ValueError:
            raise ValueError('Multiply only defined between DMSAngle Object and Int or Float')

    def __rmul__(self, other):
        try:
            return dec2dms(other * self.dec())
        except ValueError:
            raise ValueError('Multiply only defined between DMSAngle Object and Int or Float')

    def __truediv__(self, other):
        return dec2dms(self.dec() / other)

    def __abs__(self):
        return DMSAngle(abs(self.degree), self.minute, self.second)

    def __neg__(self):
        return DMSAngle(-self.degree, self.minute, self.second)

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def dec(self):
        if self.degree >= 0:
            return self.degree + (self.minute / 60) + (self.second / 3600)
        else:
            return -(abs(self.degree) + (self.minute / 60) + (self.second / 3600))

    def hp(self):
        if self.degree >= 0:
            return self.degree + (self.minute / 100) + (self.second / 10000)
        else:
            return -(abs(self.degree) + (self.minute / 100) + (self.second / 10000))

    def ddm(self):
        return DDMAngle(self.degree, self.minute + (self.second/60))


class DDMAngle(object):
    def __init__(self, degree, minute):
        self.degree = int(degree)
        self.minute = abs(minute)

    def __repr__(self):
        return '{DDMAngle: ' + str(self.degree) + 'd ' + str(self.minute) + 'm}'

    def __add__(self, other):
        return dec2ddm(self.dec() + other.dec())

    def __sub__(self, other):
        return dec2ddm(self.dec() - other.dec())

    def __mul__(self, other):
        try:
            return dec2ddm(self.dec() * other)
        except ValueError:
            raise ValueError('Multiply only defined between DMSAngle Object and Int or Float')

    def __rmul__(self, other):
        try:
            return dec2ddm(other * self.dec())
        except ValueError:
            raise ValueError('Multiply only defined between DMSAngle Object and Int or Float')

    def __truediv__(self, other):
        return dec2ddm(self.dec() / other)

    def __abs__(self):
        return DDMAngle(abs(self.degree), self.minute)

    def __neg__(self):
        return DDMAngle(-self.degree, self.minute)

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def dec(self):
        if self.degree >= 0:
            return self.degree + (self.minute / 60)
        else:
            return -(abs(self.degree) + (self.minute / 60))

    def hp(self):
        minute_int, second = divmod(self.minute, 1)
        if self.degree >= 0:
            return self.degree + (minute_int / 100) + (second * 0.006)
        else:
            return -(abs(self.degree) + (minute_int / 100) + (second * 0.006))

    def dms(self):
        minute_int, second = divmod(self.minute, 1)
        return DMSAngle(self.degree, minute_int, second * 60)


def dec2dms(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    return DMSAngle(degree, minute, second) if dec >= 0 else DMSAngle(-degree, minute, second)


def dec2ddm(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    minute = minute + (second / 60)
    return DDMAngle(degree, minute) if dec >= 0 else DDMAngle(-degree, minute)


def hp2dms(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    return DMSAngle(degree, minute, second * 10) if hp >= 0 else DMSAngle(-degree, minute, second * 10)


def hp2ddm(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    minute = minute + (second / 6)
    return DDMAngle(degree, minute) if hp >= 0 else DDMAngle(-degree, minute)


# ----------------
# Legacy Functions
# dd2dms refactored as dec2hp
def dd2dms(dd):
    minutes, seconds = divmod(abs(dd) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    dms = degrees + (minutes / 100) + (seconds / 10000)
    return dms if dd >= 0 else -dms


# dms2dd refactored as hp2dec
def dms2dd(dms):
    degmin, seconds = divmod(abs(dms) * 1000, 10)
    degrees, minutes = divmod(degmin, 100)
    dd = degrees + (minutes / 60) + (seconds / 360)
    return dd if dms >= 0 else -dd

  
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


def read_dnacoord(fn):
    coord_list = []
    with open(fn, 'r') as file:
        dnadata = file.readlines()
        for line in dnadata:
            pointid = line[0:20]
            const = line[21:25]
            easting = float(line[28:40])
            northing = float(line[41:58])
            zone = int(line[60:63])
            lat = float(line[63:78])
            long = float(line[78:92])
            ortho_ht = float(line[93:103])
            ell_ht = float(line[103:114])
            x = float(line[115:129])
            y = float(line[130:144])
            z = float(line[145:159])
            x_sd = float(line[160:171])
            y_sd = float(line[172:181])
            z_sd = float(line[182:191])
            desc = line[192:-1]
            record = DNACoord(pointid.strip(), const.strip(), easting, northing, zone, lat,
                              long, ortho_ht, ell_ht, x, y, z, x_sd, y_sd, z_sd, desc.strip())
            coord_list.append(record)
    return coord_list
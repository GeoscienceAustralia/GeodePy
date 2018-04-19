#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Conversions Module
"""

from math import modf


def dec2hp(lon, lat):
    """Convert decimal degrees to HP notation

    Longitudes go from -180 to 180
    Latitudes go from -90 to 90
    """
    def fmt_hp(coord):
        """Function to piece together a coordinate in HP notation
        """
        coord = float(coord)
        if coord < 0:
            flag = -1
            coord = abs(coord)
        else:
            flag = 1
        deg = int(coord)
        min = int(60 * (coord - deg))
        sec = 60 * (60 * (coord - deg) - min)
        # Move decimal over two places and round to an integer
        sec = round(sec * 100)
        deg *= flag  # deal with negatives
        hp_coord = '{}.{:02d}{:04d}'.format(deg, min, sec)

        return hp_coord

# Convert the coordinates
    hp_lon = fmt_hp(lon)
    hp_lat = fmt_hp(lat)

    return hp_lon, hp_lat


def hp2dec(lon, lat):
    """Convert HP notation to decimal degrees

    Longitudes go from -180 to 180
    Latitudes go from -90 to 90
    """
    def fmt_dec(coord):
        """Function to piece together a decimal coordinate
        """
        coord = float(coord)
        if coord < 0:
            flag = -1
            coord = abs(coord)
        else:
            flag = 1
        min_sec, deg = modf(coord)
        sec, min = modf(min_sec * 100)
        sec *= 100
        dec_coord = deg + min / 60 + sec / 3600
        dec_coord = float('{:.6f}'.format(dec_coord))
        dec_coord *= flag  # deal with negatives

        return dec_coord

# Convert the coordinates
    dec_lon = fmt_dec(lon)
    dec_lat = fmt_dec(lat)

    return dec_lon, dec_lat


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

def dd2dms(dd):
    minutes, seconds = divmod(abs(dd) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    dms = degrees + (minutes / 100) + (seconds / 10000)
    return dms if dd >= 0 else -dms


def dms2dd_v(dms):
    degmin, seconds = divmod(abs(dms) * 1000, 10)
    degrees, minutes = divmod(degmin, 100)
    dd = degrees + (minutes / 60) + (seconds / 360)
    dd[dms <= 0] = -dd[dms <= 0]
    return dd


def dms2dd(dms):
    degmin, seconds = divmod(abs(dms) * 1000, 10)
    degrees, minutes = divmod(degmin, 100)
    dd = degrees + (minutes / 60) + (seconds / 360)
    return dd if dms >= 0 else -dd


class Angle(object):
    def __init__(self, decdeg):
        try:
            self.decdeg = float(decdeg)
        except ValueError:
            print('ValueError: Angle must be numeric')

    def __repr__(self):
        return '{Decimal Angle: ' + str(self.decdeg) + '}'

    def __add__(self, other):
        return self.decdeg + float(other)

    def __sub__(self, other):
        return self.decdeg - float(other)

    def __abs__(self):
        return -self.decdeg if self.decdeg < 0 else self.decdeg

    def __float__(self):
        return self.decdeg

    def hp(self):
        minutes, seconds = divmod(abs(self) * 3600, 60)
        degrees, minutes = divmod(minutes, 60)
        hp = degrees + (minutes / 100) + (seconds / 10000)
        return hp if float(self) >= 0 else -hp

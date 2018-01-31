#!/usr/bin/env python3

import math


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
        min_sec, deg = math.modf(coord)
        sec, min = math.modf(min_sec * 100)
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
        min_sec, deg = math.modf(coord)
        deg = int(deg)
        deg *= flag  # deal with negatives
        sec, min = math.modf(min_sec * 60)
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

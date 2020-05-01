#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
FileIO Module
"""

from geodepy.convert import hp2dec


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
            record = DNACoord(pointid.strip(), const.strip(), easting,
                              northing, zone, lat, long, ortho_ht, ell_ht, x,
                              y, z, x_sd, y_sd, z_sd, desc.strip())
            coord_list.append(record)
    return coord_list

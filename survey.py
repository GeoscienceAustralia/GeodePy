#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Module
"""

import os
from math import sqrt, degrees, radians, sin, cos, asin
from conversions import dd2dms, dms2dd


# Defines a bunch of classes required to convert GSI (or any other format) to DynaNet v3 Format


class AngleObs(object):
    def __init__(self, degree=0, minute=0, second=0):
        self.degree = int(degree)
        self.minute = int(minute)
        self.second = float(second)

    def __repr__(self):
        return repr(self.degree) + 'd ' + repr(self.minute) + '\' ' + repr(self.second) + '\"'

    def decimal(self):
        dd = abs(self.degree + self.minute / 60 + self.second / 3600)
        return dd


class Coordinate(object):
    def __init__(self, pt_id, system, hz_datum, vert_datum, epoch, x=0.0, y=0.0, z=0.0):
        self.pt_id = pt_id
        self.system = system
        self.hz_datum = hz_datum
        self.vert_datum = vert_datum
        self.epoch = epoch
        self.x = x
        self.y = y
        self.z = z


class InstSetup(object):
    def __init__(self, pt_id, coordinate, observation=None):
        self.pt_id = pt_id
        self.coordinate = coordinate
        if observation is None:
            self.observation = []
        else:
            self.observation = []
            self.observation.append(observation)

    def __repr__(self):
        return ('{InstSetup: ' + repr(self.pt_id)
                + ' ' + repr(self.coordinate)
                + '}\n Observations:\n'
                + repr(self.observation))

    def addobs(self, observation):
        return self.observation.append(observation)


class Observation(object):
    def __init__(self, from_id, to_id, inst_height=0.0, target_height=0.0,
                 hz_obs=AngleObs(0, 0, 0), va_obs=AngleObs(0, 0, 0), sd_obs=0.0, hz_dist=0.0, vert_dist=0.0):
        self.from_id = from_id
        self.to_id = to_id
        self.inst_height = inst_height
        self.target_height = target_height
        self.hz_obs = hz_obs
        self.va_obs = va_obs
        self.sd_obs = sd_obs
        self.hz_dist = hz_dist
        self.vert_dist = vert_dist

    def __repr__(self):
        return ('{to: ' + repr(self.to_id)
                + '; target_height ' + repr(self.target_height)
                + '; hz_obs ' + repr(self.hz_obs)
                + '; va_obs ' + repr(self.va_obs)
                + '; sd_obs ' + repr(self.sd_obs)
                + '; target_height ' + repr(self.target_height)
                + '}')


# Functions to read in data to classes from Leica GSI format file (GA_Survey2.frt)


def readgsi(filepath):
    """
    Takes in a gsi file (GA_Survey2.frt) and returns a list
    of stations with their associated observations.
    :param filepath:
    :return:
    """
    # Check file extension, throw except if not .gsi
    ext = os.path.splitext(filepath)[-1].lower()
    try:
        if ext != '.gsi':
            raise ValueError
    except ValueError:
        print('ValueError: file must have .gsi extension')
        return
    # Read data from gsi file
    with open(filepath, 'r') as file:
        gsidata = file.readlines()
        stn_index = [0]
        for i in gsidata:
            # Create list of line numbers of station records (Only stations have '84..' string)
            if '84..' in i:
                lnid = i[3:7]
                lnid = int(lnid.lstrip('0'))
                stn_index.append(lnid)
        gsi_listbystation = []
        # Create lists of gsi data with station records as first element
        for i in range(0, (len(stn_index) - 1)):
            gsi_listbystation.append(gsidata[(stn_index[i]) - 1:(stn_index[i + 1])])
        del gsi_listbystation[0]
    return gsi_listbystation


def gsi2class(gsi_list):
    """
    Takes a list where first entry is station record and
    all remaining records are observations and creates
    a InstSetup Object with Observation Objects included.
    :param gsi_list:
    :return:
    """
    def readgsiword16(linestring, word_id):
        wordstart = str.find(linestring, word_id)
        word_val = linestring[(wordstart + 7):(wordstart + 23)]
        word_val = int(word_val.lstrip('0'))
        return word_val

    def parse_ptid(gsi_line):
        ptid = gsi_line[8:24]
        ptid = ptid.lstrip('0')
        return ptid

    def parse_easting(gsi_line):
        return (readgsiword16(gsi_line, '84..')) / 10000

    def parse_northing(gsi_line):
        return (readgsiword16(gsi_line, '85..')) / 10000

    def parse_elev(gsi_line):
        return (readgsiword16(gsi_line, '86..')) / 10000

    def parse_hz(gsi_line):
        dms = readgsiword16(gsi_line, '21.324') / 100000
        degmin, second = divmod(abs(dms) * 1000, 10)
        degree, minute = divmod(degmin, 100)
        return AngleObs(degree, minute, second * 10)

    def parse_slope(gsi_line):
        return (readgsiword16(gsi_line, '31..')) / 10000

    def parse_dist(gsi_line):
        return (readgsiword16(gsi_line, '32..')) / 10000

    def parse_tgtht(gsi_line):
        return (readgsiword16(gsi_line, '87..')) / 10000

    def parse_instht(gsi_line):
        return (readgsiword16(gsi_line, '88..')) / 10000

    def parse_vert(gsi_line):
        dms = readgsiword16(gsi_line, '22.324') / 100000
        degmin, seconds = divmod(abs(dms) * 1000, 10)
        degrees, minutes = divmod(degmin, 100)
        return AngleObs(degrees, minutes, seconds * 10)

    for record in gsi_list:
        from_stn = parse_ptid(record[0])
        obs_list = []
        for line in record:
            if '31..' in line:
                to_stn = parse_ptid(line)
                hz = parse_hz(line)
                vert = parse_vert(line)
                slope = parse_slope(line)
                tgtht = parse_tgtht(line)
                instht = parse_instht(line)
                obs = Observation(from_stn, to_stn, instht, tgtht, hz, vert, slope)
                obs_list.append(obs)
        if '84..' in record[0]:
            pt_id = parse_ptid(record[0])
            easting = parse_easting(record[0])
            northing = parse_northing(record[0])
            elev = parse_elev(record[0])
            coord = Coordinate(pt_id, 'utm', 'gda', 'gda', '2018', easting, northing, elev)
            setup = InstSetup(pt_id, coord, obs_list)
    return setup


"""
    # Create Coordinate and Instrument Setup Objects
    coord = Coordinate(pt_id, 'utm', 'gda', 'gda', '2018', easting, northing, elev)
    setup = InstSetup(pt_id, coord)
    # Add Instrument Setup to Project
    project.update({'InstSetup_' + str(stncount): setup})
        # Add Observation Records to InstSetup Objects

return project
"""


def readgsiword16(linestring, word_id):
    wordstart = str.find(linestring, word_id)
    word_val = linestring[(wordstart + 7):(wordstart + 23)]
    word_val = 0.0 if word_val.lstrip('0') == '' else int(word_val.lstrip('0'))
    return word_val


# Functions to write out station and obs data to DNA format


def dnaout_sd(observation):
    return ('S '
            + observation.from_name.ljust(20)
            + observation.to_name.ljust(20)
            + ''.ljust(20)
            + str(observation.sd_obs).ljust(14)  # 76
            + ''.ljust(14)
            + '0.0010'.ljust(9)         # add standard deviation
            + str(1.7960).ljust(7)      # add intrument height
            + str(observation.target_height).ljust(7))


def dnaout_va(observation):
    return ('V '
            + observation.from_name.ljust(20)
            + observation.to_name.ljust(20)
            + ''.ljust(34)
            + str(observation.va_obs.degrees).rjust(3)
            + ' '
            + str('%02d' % observation.va_obs.minutes)
            + ' '
            + str(observation.va_obs.seconds).ljust(8)
            + '1.0000'.ljust(9)         # add standard deviation
            + str(1.7960).ljust(7)      # add intrument height
            + str(observation.target_height).ljust(7))


def va_conv(verta_hp, slope_dist, height_inst=0, height_tgt=0):
    """
    Function to convert vertical angles (zenith distances) and slope distances into horizontal
    distances and changes in height. Instrument and Target heights can be entered to allow
    computation of zenith and slope distances between ground points.

    :param verta_hp:        Vertical Angle from Instrument to Target, expressed in HP Format (DDD.MMSSSSSS)
    :param slope_dist:      Slope Distance from Instrument to Target in metres
    :param height_inst:     Height of Instrument. Optional - Default Value of 0m
    :param height_tgt:      Height of Target. Optional - Default Value of 0m

    :return: verta_pt_hp:   Vertical Angle between Ground Points, expressed in HP Format (DDD.MMSSSSSS)
    :return: slope_dist_pt: Slope Distance between Ground Points in metres
    :return: hz_dist:       Horizontal Distance
    :return: delta_ht:      Change in height between Ground Points in metres
    """
    # Convert Zenith Angle to Vertical Angle
    try:
        if verta_hp == 0 or verta_hp == 180:
            raise ValueError
        elif 0 < verta_hp < 180:
            verta = radians(90 - dms2dd(verta_hp))
        elif 180 < verta_hp < 360:
            verta = radians(270 - dms2dd(verta_hp))
        else:
            raise ValueError
    except ValueError:
        print('ValueError: Vertical Angle Invalid')
        return
    # Calculate Horizontal Dist and Delta Height
    hz_dist = slope_dist * cos(verta)
    delta_ht = slope_dist * sin(verta)
    # Account for Target and Instrument Heights
    if height_inst == 0 and height_tgt == 0:
        verta_pt_hp = verta_hp
        slope_dist_pt = slope_dist
    else:
        delta_ht = height_inst + delta_ht - height_tgt
        slope_dist_pt = sqrt(delta_ht ** 2 + hz_dist ** 2)
        verta_pt = asin(delta_ht / slope_dist)
        verta_pt_hp = dd2dms(degrees(verta_pt) + 90)
    return verta_pt_hp, slope_dist_pt, hz_dist, delta_ht


"""
to_stn = ['GA03', 'GA02', 'SYM2', 'ATWR', 'MAHON', 'MAHON', 'ATWR', 'SYM2', 'GA02', 'GA03']
flfr = [359.59597, 14.48289, 25.35515, 215.57043, 300.59046, 120.59060, 35.57045, 205.35530, 194.48282, 179.59581]
flfr_vert = [90.30228, 90.09547, 89.50396, 88.03600, 87.58014, 272.01587, 271.55593, 270.09227, 269.50091, 269.29414]
"""


def hz_round(brg_list):
    """
    Input: an even palindromic list of horizontal angle observations in two faces
    (e.g. [fl1, fl2, fl3, fr3, fr2, fr1])
    Output: an averaged face-left sense list of horizontal angles.
    """
    brg_avg = []
    obs = int((len(brg_list))/2)
    for i in range(0, obs):
        hz_avg = (dms2dd(brg_list[i]) + (dms2dd(brg_list[-(i+1)])-180))/2
        brg_avg.append(round(dd2dms(hz_avg), 7))
    return brg_avg


def va_round(va_list):
    """
    Input: an even palindromic list of vertical angle observations in two faces
    (e.g. [fl1, fl2, fl3, fr3, fr2, fr1])
    Output: an averaged face-left sense list of vertical angles.
    """
    va_avg = []
    obs = int((len(va_list))/2)
    for i in range(0, obs):
        fl_ang = dms2dd(va_list[i]) - 90
        fr_ang = 270 - dms2dd(va_list[-(i+1)])
        ang_avg = (fl_ang + fr_ang)/2 + 90
        va_avg.append(round(dd2dms(ang_avg), 7))
    return va_avg

"""
# Test for round of obs
if to_stn == to_stn[::-1] and len(to_stn) % 2 == 0:
    brg_avg = hz_round(flfr)
    va_avg = va_round(flfr_vert)
    to_stn_avg = to_stn[0:int(len(to_stn)/2)]
else:
    brg_avg = ['nope']
"""
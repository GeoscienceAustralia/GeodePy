#!/usr/bin/env python3

import os
from math import sqrt, degrees, radians, sin, cos, asin
from conversions import dd2dms, dms2dd


# Defines a bunch of classes required to convert GSI (or any other format) to DynaNet v3 Format


class Coordinate(object):
    def __init__(self, pt_id, system, hz_datum, vert_datum, epoch, x=0, y=0, z=0):
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
    def __init__(self, from_id, to_id, inst_height=0, target_height=0, hz_obs=0, va_obs=0, sd_obs=0, hz_dist=0, vert_dist=0):
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
        return ('{to: ' + repr(self.to_name)
                + '; target_height ' + repr(self.target_height)
                + '; hz_obs ' + repr(self.hz_obs)
                + '; va_obs ' + repr(self.va_obs)
                + '; sd_obs ' + repr(self.sd_obs)
                + '; target_height ' + repr(self.target_height)
                + '}')


class AngleObs(object):
    def __init__(self, degrees=0, minutes=0, seconds=0):
        self.degrees = int(degrees)
        self.minutes = int(minutes)
        self.seconds = float(seconds)

    def __repr__(self):
        return repr(self.degrees) + 'd ' + repr(self.minutes) + '\' ' + repr(self.seconds) + '\"'

    def decimal(self):
        dd = abs(self.degrees + self.minutes / 60 + self.seconds / 3600)
        return dd


# Functions to read in data to classes from Leica GSI format file (GA_Survey2.frt)

def readgsiword16(linestring, word_id):
    wordstart = str.find(linestring, word_id)
    word_val = linestring[(wordstart + 7):(wordstart + 23)]
    word_val = int(word_val.lstrip('0'))
    return word_val


def readgsi(filepath):
    # check file extension, throw except if not .gsi
    ext = os.path.splitext(filepath)[-1].lower()
    try:
        if ext != '.gsi':
            raise ValueError
    except ValueError:
        print('ValueError: file must have .gsi extension')
        return
    # Open file and read data line-by-line
    with open(filepath, 'r') as file:
        project = dict()
        stncount = 0
        gsilines = file.readlines()
        for line in gsilines:
            ln_id = int(line[3:7])
            # Create Station Record
            if '84..' in line:
                stncount += 1
                # Parse Pt ID
                pt_id = line[8:24]
                pt_id = pt_id.lstrip('0')
                # Parse Easting
                easting = readgsiword16(line, '84..')
                easting = easting / 10000
                # Parse Northing
                northing = readgsiword16(line, '85..')
                northing = northing / 10000
                # Parse Elev
                elev = readgsiword16(line, '86..')
                elev = elev / 10000
                # Create Coordinate Object
                coord = Coordinate(pt_id, 'utm', 'gda', 'gda', '2018', easting, northing, elev)
                # Create and Add Instrument Setup to Project
                project.update({'InstSetup_' + str(stncount): InstSetup(pt_id, coord)})
    # Wrap this in a while loop for all InstSetups
        # Create InstSetup Object
        # Read all obs into Observation Object until next InstSetup

    return project

# Functions to write out station and obs data to DNA format


def dnaout_sd(observation):
    return ('S '
            + observation.from_name.ljust(20)
            + observation.to_name.ljust(20)
            + ''.ljust(20)
            + str(observation.sd_obs).ljust(14) #76
            + ''.ljust(14)
            + '0.0010'.ljust(9)         #add standard deviation
            + str(1.7960).ljust(7)      #add intrument height
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
            + '1.0000'.ljust(9)         #add standard deviation
            + str(1.7960).ljust(7)      #add intrument height
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


to_stn = ['GA03', 'GA02', 'SYM2', 'ATWR', 'MAHON', 'MAHON', 'ATWR', 'SYM2', 'GA02', 'GA03']
flfr = [359.59597, 14.48289, 25.35515, 215.57043, 300.59046, 120.59060, 35.57045, 205.35530, 194.48282, 179.59581]
flfr_vert = [90.30228, 90.09547, 89.50396, 88.03600, 87.58014, 272.01587, 271.55593, 270.09227, 269.50091, 269.29414]


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


# Test for round of obs
if to_stn == to_stn[::-1] and len(to_stn) % 2 == 0:
    brg_avg = hz_round(flfr)
    va_avg = va_round(flfr_vert)
    to_stn_avg = to_stn[0:int(len(to_stn)/2)]
else:
    brg_avg = ['nope']

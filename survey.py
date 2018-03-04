#!/usr/bin/env python3

from math import sqrt, degrees, radians, sin, cos, asin
from conversions import dd2dms, dms2dd

# Defines a bunch of classes required to convert GSI (or any other format) to DynaNet v3 Format


class Coordinate(object):
    def __init__(self, name, system, hz_datum, vert_datum, epoch, x=0, y=0, z=0):
        self.name = name
        self.system = system
        self.hz_datum = hz_datum
        self.vert_datum = vert_datum
        self.epoch = epoch
        self.x = x
        self.y = y
        self.z = z


class InstSetup(object):
    def __init__(self, name, coordinate, inst_height, observation=None):
        self.name = name
        self.coordinate = coordinate
        self.inst_height = inst_height
        if observation is None:
            self.observation = []
        else:
            self.observation.append(observation)

    def __repr__(self):
        return ('{InstSetup: ' + repr(self.name)
                + ' ' + repr(self.coordinate)
                + '; inst_height ' + repr(self.inst_height)
                + '}\n Observations:\n'
                + repr(self.observation))


class Observation(object):
    def __init__(self, to_name, coordinate, target_height=0, hz_obs=0, va_obs=0, sd_obs=0, hz_dist=0, vert_dist=0):
        self.to_name = to_name
        self.coordinate = coordinate
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
        return

# Functions to write out station and obs data to DNA format






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
    except:
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


# Probably not going to use this stuff

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

#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Module
"""

from math import sqrt, sin, cos, asin, radians, degrees
from geodepy.convert import hp2dec, dec2hp


def va_conv(verta_hp, slope_dist, height_inst=0, height_tgt=0):
    """
    Function to convert vertical angles (zenith distances) and slope distances
    into horizontal distances and changes in height. Instrument and Target
    heights can be entered to allow computation of zenith and slope distances
    between ground points.

    :param verta_hp:        Vertical Angle from Instrument to Target, expressed
                            in HP Format (DDD.MMSSSSSS)
    :param slope_dist:      Slope Distance from Instrument to Target in metres
    :param height_inst:     Height of Instrument. Optional - Default Value of 0m
    :param height_tgt:      Height of Target. Optional - Default Value of 0m

    :return: verta_pt_hp:   Vertical Angle between Ground Points, expressed
                            in HP Format (DDD.MMSSSSSS)
    :return: slope_dist_pt: Slope Distance between Ground Points in metres
    :return: hz_dist:       Horizontal Distance
    :return: delta_ht:      Change in height between Ground Points in metres
    """
    # Convert Zenith Angle to Vertical Angle
    try:
        if verta_hp == 0 or verta_hp == 180:
            raise ValueError
        elif 0 < verta_hp < 180:
            verta = radians(90 - hp2dec(verta_hp))
        elif 180 < verta_hp < 360:
            verta = radians(270 - hp2dec(verta_hp))
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
        verta_pt_hp = dec2hp(degrees(verta_pt) + 90)
    return verta_pt_hp, slope_dist_pt, hz_dist, delta_ht

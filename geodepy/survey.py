#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Module
"""

from math import sqrt, sin, cos, asin, radians, degrees, exp
from geodepy.convert import hp2dec, dec2hp


def first_vel_params(wavelength, temp=12, pressure=1013.25, rel_humidity=60):
    """
    Calculates the constant First Velocity Correction Parameters C and D
    for a given set of standard instrument settings
    :param wavelength: Instrument Carrier Wavelength (micrometres)
    :param temp: Standard Temperature (degrees Celsius)
    :param pressure: Standard Pressure (hectopascals or millibars)
    :param rel_humidity: Standard Relative Humidity (percentage)
    :return: First Velocity Correction Parameters C and D
    """
    group_refractivity = 287.6155 + (4.8866 / wavelength ** 2) + (0.068 / wavelength ** 4)
    param_d = (273.15 / 1013.25) * group_refractivity
    saturation_pressure = ((1.0007 + ((3.46 * pressure) * (10 ** -6)))
                           * 6.1121 * exp((17.502 * temp) / (240.94 + temp)))
    param_c = ((((group_refractivity * pressure) / (273.15 + temp)) * (273.15 / 1013.25))
               - (11.27 * ((saturation_pressure * (rel_humidity / 100)) / (273.15 + temp))))
    return param_c, param_d


def first_vel_corrn(dist, first_vel_param, temp, pressure, rel_humidity):
    """
    Carries out First Velocity Correction of Electronic Distance Measurement,
    given correction parameters and atmospheric observations
    :param dist: Uncorrected Observed Slope Distance
    :param first_vel_param: Tuple of First Velocity Parameters C and D (see function first_vel_params)
    :param temp: Observed Temperature (degrees Celsius)
    :param pressure: Observed Pressure (hectopascals or millibars)
    :param rel_humidity: Observed Relative Humidity (percentage)
    :return: Slope Distance with First Velocity Correction applied
    """
    param_c = first_vel_param[0]
    param_d = first_vel_param[1]
    part_h2o_press = exp((17.269 * temp) / (237.3 + temp))
    first_vel_corrn_ppm = param_c - (
                (param_d * pressure) / (273.15 + temp) + ((11.27 * part_h2o_press) * (0.061078 * rel_humidity)) / (
                    273.15 + temp))
    first_vel_corrn_metres = dist * first_vel_corrn_ppm * (10 ** -6)
    return first_vel_corrn_metres


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

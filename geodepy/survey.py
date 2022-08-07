#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Module
"""

from math import sqrt, sin, cos, atan, radians, degrees, exp
from statistics import mean, stdev
from geodepy.convert import rect2polar, polar2rect


def first_vel_params(wavelength, frequency, n_REF=None, unit_length=None):
    """
    Calculates the constant First Velocity Correction Parameters C and D
    for a given set of standard instrument settings
    :param wavelength: Instrument Carrier Wavelength (Nanometers) - Mandatory
    :param frequency: Instrument modulation frequency (Hz) - Mandatory
    :param n_REF: manufacturers reference refractive index - Recommended
    :param unit_length: unit length of instrument - Optional
    :return: First Velocity Correction Parameters C and D
    Reference Rueger, J.M., 2012, Electronic Distance Measurement – An Introduction, 4rd edition, Springer, Berlin
    """
    
    if not n_REF and not unit_length:
        raise ValueError('Error - n_REF and unit_length cannot both be None')
        
    if not n_REF:
        # Rueger eq 6.3
        n_REF = 299792458 / (2 * unit_length * frequency)
    
    # Rueger eq 6.12
    param_c = (n_REF-1)*10**6
    
    # Rueger eq 5.12
    # https://office.iag-aig.org/doc/5d7b8fda0c032.pdf Resolution 3
    nG_1_10_6 = 287.6155 + (4.8866 / wavelength ** 2) + (0.068 / wavelength ** 4)
    # Rueger eq 6.13
    param_d = (273.15 / 1013.25) * nG_1_10_6
    
    return param_c, param_d


def part_h2o_vap_press(dry_temp, pressure, rel_humidity=None, wet_temp=None):
    """
    :param dry_temp: Observed Dry Temperature (degrees Celsius)
    :param pressure: Observed Pressure (hectopascals or millibars)
    :param rel_humidity: Observed Relative Humidity (percentage) - Optional if wet_temp supplied
    :param wet_temp: Observed Wet Temperature (degrees Celsius) - Optional if rel_humidity supplied
    Reference Rueger, J.M., 2012, Electronic Distance Measurement – An Introduction, 4rd edition, Springer, Berlin
    """
    
    if not rel_humidity and not wet_temp:
        raise ValueError('Error - rel_humidity and wet_temp cannot both be None')
                         
    if rel_humidity:
        wet_temp = dry_temp
    # Rueger eq 5.27
    E_w = ((1.0007 + ((3.46 * pressure) * (10 ** -6)))
         * 6.1121 * exp((17.502 * wet_temp) / (240.94 + wet_temp)))
    
    if rel_humidity:
        # Rueger eq 5.29
        e = (E_w*rel_humidity)/100
    else:
        e = E_w - 0.000662 * pressure * (dry_temp - wet_temp)
    
    return e
    

def first_vel_corrn(dist, first_vel_param, temp, pressure,
                    rel_humidity=None, wet_temp=None):
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
    e = part_h2o_vap_press(temp, pressure, rel_humidity, wet_temp)
    
    # Rueger eq 6.11
    t_273_15 = temp+273.15
    first_vel_corrn_ppm = (param_c-((param_d*pressure)/t_273_15)+((11.27*e)/t_273_15))
    first_vel_corrn_metres = dist * first_vel_corrn_ppm * (10 ** -6)
    
    return first_vel_corrn_metres


def mets_partial_differentials(group_ref_Index=1.00028, 
                               temp=15, 
                               pressure=1013.25,
                               rel_humidity=60):
    """ Returns the sensitivity coefficients for temp, pressure and humidity in ppm
    :param group_ref_Index: manufacturers group refractive index of light
    :param temp: Observed Dry Temperature (degrees Celsius)
    :param pressure: Observed Pressure (hectopascals or millibars)
    :param rel_humidity: Observed Relative Humidity (percentage) - Optional if wet_temp supplied
    """
    ng_1 = group_ref_Index - 1
    t_273_15 = temp + 273.15
    e = part_h2o_vap_press(temp, pressure, rel_humidity)
    
    param_k = (((((ng_1 * 273.15 * pressure)
                  /1013.25)
                 -(11.27 * e * 10**-6))
                /t_273_15**2)
               *(10**6))
    param_l = (0.26957809 * (ng_1 / t_273_15)
               *(10**6))
    param_m = (((11.27 * 10**-6)
               /t_273_15)
               *(10**6))
    
    return param_k, param_l, param_m


def precise_inst_ht(vert_list, spacing, offset):
    """
    Uses a set of Vertical Angle Observations taken to a
    levelling staff at regular intervals to determine the
    height of the instrument above a reference mark
    :param vert_list: List of Vertical (Zenith) Angle Observations (minimum of 3) in Decimal Degrees format
    :param spacing: Distance in metres between each vertical angle observation
    :param offset: Lowest observed height above reference mark
    :return: Instrument Height above reference mark and its standard deviation
    """
    if len(vert_list) < 3:
        raise ValueError('ValueError: 3 or more vertical angles required')
    vert_list.sort(reverse=True)
    vert_pairs = [(va1, va2) for va1, va2 in zip(vert_list, vert_list[1:])]
    base_ht = []
    height_comp = []
    for num, pair in enumerate(vert_pairs):
        base_ht_pair = offset + num * spacing
        base_ht.append(base_ht_pair)
        dist_a = sin(radians(pair[1])) * (spacing / (sin(radians(pair[0] - pair[1]))))
        delta_ht = dist_a * (sin(radians(pair[0] - 90)))
        height_comp.append(delta_ht + base_ht[num])
    return round(mean(height_comp), 5), round(stdev(height_comp), 5)


def joins(east1, north1, east2, north2):
    return rect2polar(east2 - east1, north2 - north1)


def radiations(east1, north1, brg1to2, dist, rotation=0, psf=1):
    delta_east, delta_north = polar2rect(dist * psf, brg1to2 + rotation)
    return east1 + delta_east, north1 + delta_north


def va_conv(zenith_angle, slope_dist, height_inst=0, height_tgt=0):
    """
    Function to convert vertical angles (zenith distances) and slope distances
    into horizontal distances and changes in height. Instrument and Target
    heights can be entered to allow computation of zenith and slope distances
    between ground points.

    :param zenith_angle:    Zenith Angle from Instrument to Target, expressed
                            in decimal degrees
    :param slope_dist:      Slope Distance from Instrument to Target in metres
    :param height_inst:     Height of Instrument. Optional - Default Value of 0m
    :param height_tgt:      Height of Target. Optional - Default Value of 0m

    :return: vert_angle_pt: Vertical Angle between Ground Points, expressed
                            in decimal degrees
    :return: slope_dist_pt: Slope Distance between Ground Points in metres
    :return: hz_dist:       Horizontal Distance
    :return: delta_ht:      Change in height between Ground Points in metres
    """
    # Convert Zenith Angle to Vertical Angle
    try:
        if zenith_angle == 0 or zenith_angle == 180:
            raise ValueError
        elif 0 < zenith_angle < 180:
            zenith_angle = radians(90 - zenith_angle)
        elif 180 < zenith_angle < 360:
            zenith_angle = radians(270 - zenith_angle)
        else:
            raise ValueError
    except ValueError:
        raise ValueError('ValueError: Vertical Angle Invalid')
    # Calculate Horizontal Dist and Delta Height
    hz_dist = slope_dist * cos(zenith_angle)
    delta_ht = slope_dist * sin(zenith_angle)
    # Account for Target and Instrument Heights
    delta_ht = height_inst + delta_ht - height_tgt
    slope_dist_pt = sqrt(delta_ht ** 2 + hz_dist ** 2)
    vert_angle_pt = atan(delta_ht / hz_dist)
    vert_angle_pt = degrees(vert_angle_pt)
    return vert_angle_pt, slope_dist_pt, hz_dist, delta_ht

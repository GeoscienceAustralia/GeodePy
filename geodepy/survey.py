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

    Reference Rueger, J.M., 2012, Electronic Distance Measurement: An Introduction, 4rd edition, Springer, Berlin

    :param wavelength: Instrument Carrier Wavelength (Micrometers) - Mandatory
    :param frequency: Instrument modulation frequency (Hz) - Optional
    :param n_REF: manufacturers reference refractive index - Recommended
    :param unit_length: unit length of instrument - Optional
    :return: First Velocity Correction Parameters C and D

    """

    if not n_REF:
        if not unit_length:
            raise ValueError("Error - n_REF and unit_length cannot both be None")
        if not frequency:
            raise ValueError("Error - n_REF and frequency cannot both be None")
        # Rueger eq 6.3
        n_REF = 299792458 / (2 * unit_length * frequency)

    # Rueger eq 6.12
    param_c = (n_REF - 1) * 10**6

    # Rueger eq 5.12
    # https://office.iag-aig.org/doc/5d7b8fda0c032.pdf Resolution 3
    nG_1_10_6 = 287.6155 + (4.8866 / wavelength**2) + (0.068 / wavelength**4)
    # Rueger eq 6.13
    param_d = (273.15 / 1013.25) * nG_1_10_6

    return param_c, param_d


def part_h2o_vap_press(dry_temp, pressure, rel_humidity=None, wet_temp=None):
    """
    Reference Rueger, J.M., 2012, Electronic Distance Measurement: An Introduction, 4rd edition, Springer, Berlin

    :param dry_temp: Observed Dry Temperature (degrees Celsius)
    :param pressure: Observed Pressure (hectopascals or millibars)
    :param rel_humidity: Observed Relative Humidity (percentage) - Optional if wet_temp supplied
    :param wet_temp: Observed Wet Temperature (degrees Celsius) - Optional if rel_humidity supplied

    """

    if not rel_humidity and not wet_temp:
        raise ValueError("Error - rel_humidity and wet_temp cannot both be None")

    if rel_humidity:
        wet_temp = dry_temp
    # Rueger eq 5.27
    E_w = (
        (1.0007 + ((3.46 * pressure) * (10**-6)))
        * 6.1121
        * exp((17.502 * wet_temp) / (240.94 + wet_temp))
    )

    if rel_humidity:
        # Rueger eq 5.29
        e = (E_w * rel_humidity) / 100
    else:
        e = E_w - 0.000662 * pressure * (dry_temp - wet_temp)

    return e


def first_vel_corrn(
    dist,
    first_vel_param,
    temp,
    pressure,
    rel_humidity=None,
    wet_temp=None,
    CO2_ppm=None,
    wavelength=None,
):
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
    if not CO2_ppm:
        param_d = first_vel_param[1]
        e = part_h2o_vap_press(temp, pressure, rel_humidity, wet_temp)

        # Rueger eq 6.11
        t_273_15 = temp + 273.15
        first_vel_corrn_ppm = (
            param_c - ((param_d * pressure) / t_273_15) + ((11.27 * e) / t_273_15)
        )
        first_vel_corrn_metres = dist * first_vel_corrn_ppm * (10**-6)

    else:
        if all([temp, pressure, rel_humidity, wavelength]):
            e = humidity2part_water_vapour_press(rel_humidity, temp)
            NPROPG_1 = group_refractivity(wavelength, temp, pressure, e, CO2_ppm)
        else:
            raise ValueError(
                "Error - temp, pressure, rel_humidity and wavelength "
                + "are all mandatory when CO2_ppm is not None"
            )

        n_ref = 1 + (param_c / 1.0e6)
        n_g = 1 + (NPROPG_1 / 1.0e8)
        first_vel_corrn_metres = ((n_ref / n_g) - 1) * dist

    return first_vel_corrn_metres


def mets_partial_differentials(
    group_ref_Index=1.00028, temp=15, pressure=1013.25, rel_humidity=60
):
    """
    Calculates the sensitivity coefficients for temp, pressure and humidity in ppm

    :param group_ref_Index: manufacturers group refractive index of light
    :param temp: Observed Dry Temperature (degrees Celsius)
    :param pressure: Observed Pressure (hectopascals or millibars)
    :param rel_humidity: Observed Relative Humidity (percentage) - Optional if wet_temp supplied
    :return: Sensitivity Coefficients K, L and M in ppm
    """
    ng_1 = group_ref_Index - 1
    t_273_15 = temp + 273.15
    e = part_h2o_vap_press(temp, pressure, rel_humidity)

    param_k = (
        (((ng_1 * 273.15 * pressure) / 1013.25) - (11.27 * e * 10**-6)) / t_273_15**2
    ) * (10**6)
    param_l = 0.26957809 * (ng_1 / t_273_15) * (10**6)
    param_m = ((11.27 * 10**-6) / t_273_15) * (10**6)

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
        raise ValueError("ValueError: 3 or more vertical angles required")
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
    """
    Calculates the bearing and distance from point 1 to point 2

    :param east1: Easting of Point 1
    :param north1: Northing of Point 1
    :param east2: Easting of Point 2
    :param north2: Northing of Point 2
    :return: Distance and Bearing from Point 1 to Point 2
    """
    return rect2polar(east2 - east1, north2 - north1)


def radiations(east1, north1, brg1to2, dist, rotation=0, psf=1):
    """
    Calculates the coordinates of Point 2 given Point 1, bearing and distance

    :param east1: Easting of Point 1
    :param north1: Northing of Point 1
    :param brg1to2: Bearing from Point 1 to Point 2 (decimal degrees)
    :param dist: Distance from Point 1 to Point 2 (metres)
    :param rotation: Rotation to be applied to bearing (decimal degrees) - Optional
    :param psf: Prism Scale Factor to be applied to distance - Optional
    :return: Easting and Northing of Point 2
    """
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

    :return: 
        - vert_angle_pt - Vertical Angle between Ground Points, expressed in decimal degrees
        - slope_dist_pt - Slope Distance between Ground Points in metres
        - hz_dist - Horizontal Distance
        - delta_ht - Change in height between Ground Points in metres
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
        raise ValueError("ValueError: Vertical Angle Invalid")
    # Calculate Horizontal Dist and Delta Height
    hz_dist = slope_dist * cos(zenith_angle)
    delta_ht = slope_dist * sin(zenith_angle)
    # Account for Target and Instrument Heights
    delta_ht = height_inst + delta_ht - height_tgt
    slope_dist_pt = sqrt(delta_ht**2 + hz_dist**2)
    vert_angle_pt = atan(delta_ht / hz_dist)
    vert_angle_pt = degrees(vert_angle_pt)
    return vert_angle_pt, slope_dist_pt, hz_dist, delta_ht


"""
The following has been translated from UNSW FORTRAN code
"""
"""
    PROGRAM RINDEX

    -----------------------------------------------------------------
    Date :  6 July 1995

    Version 4.0 with CF = 1.022

    This program calculates the phase and group refractive index
    using equations for the visible and near infrared developed by
    P.E. Ciddor of the National Measurement Laboratory, CSIRO,
    Division of the Applied Physics, Australia. (see reference below)

    Program originally written by S.K. Johnson and J.M. R?eger of the
    School of Geomatic Engineering, University of New South
    Wales, Australia.
    Program converted from FORTRAN to python by Kent Wheeler of Landgate.

    Main References:
    P.E. Ciddor, 1996, Refractive Index of Air: New Equations for the 
          Visible and Near Infrared, Applied Optics (Lasers, Phototonics,
          and Environmental optics), 35(9): 1566-73.
    R.S. Davies, 1992, Equations for Determination of the Density of Moist 
          Air (1981/91), Metrologia, 29: 67-70.
    P. Giacomo, 1982, Equation for Determination of the Density of Moist 
          Air, Metrologia, 18: 33-40.
    E.P. Peck and K. Reeder, 1972, Dispersion of Air, Journal of the
          Optical Society of America, 62(8): 958-962.
    J.C. Owens, 1967, Optical Refractive Index of Air: Dependence on
          Pressure, Temperature and Composition, Applied Optics, 6(1): 51-59
    K.P Birch and M.J. Downs, 1994, Correction to the updated Edlen Equation
          for the Refractive Index of Air, Metrologia, 31: 315-316.
    K.P. Birch and M.J. Downs, 1993, An Updated Edlen Equation for the
          Refractive Index of Air, Metrologia, 30: 155-162.

    Variables used:

    A,B,C,D = constants used to calculate SVP
    A0,A1,A2,B0,B1,C0,C1,DLC,E = constants used to calculate Z
    ALPHA, BETA, GAMMA = constants to calculate enhancement factor of
                          water vapour in air
    CF = correction factor
    CTERM = reference refractivity
    DIST_UC = uncorrected distance (m)
    DIST_CO = corrected distance (m)
    F = enhanvement factor of water vapour in air
    H = humidity (%)
    K0,K1,K2,K3 = constants used to calculate NASX
    LAMDA = wavelength (micrometre)
    MA = molar mass of dry air containing a fraction XC
          of Carbon Dioxide (kg/mol)
    MV = molar mass of water vapour (kg/mol)
    NAS_1 = refractivity of standard air with 450 ppm = (NAS-1)10E8
    NASX_1 = refractivity of standard air with XC ppm = (NASX-1)10E8
    NPROPG_1 = group refractivity = (NPROPgroup-1)10E8
    NPROPP_1 = phase refractivity = (NPROPphase-1)10E8
    NWS_1 = refractivity of water vapour at standard
            conditions = (NWS-1)10E8
    P = pressure (hPa)
    PV = partial water vapour pressure (Pa)
    RHOA = density of dry air (kg per cubic metre)
    RHOASX = density of standard air (kg per cubic metre)
    RHOW = density of pure water vapour (kg per cubic metre)
    RHOWS = density of standard water vapour (kg per cubic metre)
    SIGMA = wavenumber (1/micrometre)
    SVP = saturation vapour pressure of water vapour (Pa)
          in air at temperature (TK)
    TC = temperature (degrees Celcius)
    TK = temperature (Kelvin)
    W0,W1,W2,W3 = constants used to calculate NWS
    XC = carbon dioxide content (ppm)
    Z = compressibility of air
    ZA = compressibility of dry air
    ZWV = compressibility of water vapour
    -----------------------------------------------------------------
"""


def refractivity_constants():
    """
    :return: Refractivity constants used in the refractivity calculations.
    """
    # PECK & REEDER (1972) AS AMMENDED BY CIDDOR
    # (DRY AIR REFRACTIVITY)
    (K0, K1, K2, K3) = (238.0185, 5792105.0, 57.362, 167917.0)

    # OWENS (1967) WATER VAPOUR REFRACTIVITY
    (W0, W1, W2, W3) = (295.235, 2.6422, -0.032380, 0.004028)

    # ENHANCEMENT FACTOR (DAVIES, 1992)
    (ALPHA, BETA, GAMMA) = (1.00062, 3.14e-8, 5.6e-7)

    # COMPRESSIBILITY (DAVIES 1992)
    (A0, A1, A2, B0, B1, C0, C1, DLC, E) = (
        1.58123e-6,
        -2.9331e-8,
        1.1043e-10,
        5.707e-6,
        -2.051e-8,
        1.9898e-4,
        -2.376e-6,
        1.83e-11,
        -0.765e-8,
    )

    # CIDDOR 1995 AND DAVIES 1992.  NEW CF ON 28/06/95
    (CF, R, MV) = (1.022, 8.314510, 0.018015)

    # STANDARD CONDITION DRY AIR (PECK & REEDER 1972)
    (TCSTDAIR, TKSTDAIR, PSTDAIR) = (15.0, 288.15, 101325.0)
    (TCSTDWV, TKSTDWV, PSTDWV) = (20.0, 293.15, 1333.0)

    return (
        (K0, K1, K2, K3),
        (W0, W1, W2, W3),
        (ALPHA, BETA, GAMMA),
        (A0, A1, A2, B0, B1, C0, C1, DLC, E),
        (CF, R, MV),
        (TCSTDAIR, TKSTDAIR, PSTDAIR),
        (TCSTDWV, TKSTDWV, PSTDWV),
    )


def phase_refractivity(LAMDA, TC, P, PV, XC=420):
    """
    Calculates the phase refractivity of moist air using Ciddor's equations

    :param LAMDA: wavelength (micrometre)
    :param TC: temperature (degrees Celcius)
    :param P: pressure (hPa)
    :param PV: partial water vapour pressure (Pa)
    :param XC: carbon dioxide content (ppm)
    
    :return: NPROPP_1 - phase refractivity = (NPROPphase-1)10E8
    """

    (
        (K0, K1, K2, K3),
        (W0, W1, W2, W3),
        (ALPHA, BETA, GAMMA),
        (A0, A1, A2, B0, B1, C0, C1, DLC, E),
        (CF, R, MV),
        (TCSTDAIR, TKSTDAIR, PSTDAIR),
        (TCSTDWV, TKSTDWV, PSTDWV),
    ) = refractivity_constants()

    # CONVERT LAMDA TO SIGMA, THE PRESSURE UNITS TO PASCAL AND CREATE TC AND PP
    SIGMA = 1.0 / LAMDA

    # CONVERSION FROM hPa TO Pa AND C TO K
    P = P * 100.0
    PV = PV * 100.0
    TK = TC + 273.15

    # CALCULATE REFRACTIVITY OF STANDARD AIR (TC=15, P=101325)
    TEMP1 = SIGMA * SIGMA

    # CIDDOR EQ.(1)
    NAS_1 = (K1 / (K0 - TEMP1)) + (K3 / (K2 - TEMP1))

    # CIDDOR EQ.(2) FOR AMBIENT C02 CONTENT (XC)
    NASX_1 = NAS_1 * (1 + 0.534e-6 * (XC - 450.0))
    STD_AIR = {"(Nas-1)10E8": NAS_1, "(Nasx-1)10E8": NASX_1}

    # USING EQ. 2, CALCULATE REFRACTIVITY OF WATER VAPOUR AT STANDARD CONDITIONS
    # (TC=20, P=1333)
    TEMP2 = TEMP1 * TEMP1
    TEMP3 = TEMP1 * TEMP2

    # FOLLOWING OWENS (1967) AS AMENDED TO FIT BIRCH & DOWNS (1988)
    NWS_1 = CF * (W0 + W1 * TEMP1 + W2 * TEMP2 + W3 * TEMP3)
    STD_WATER_VAPOUR = {"(Nws-1)10E8": NWS_1}

    # CALCULATE THE DENSITY OF STANDARD AIR (XV=0, TC=15, P=101325)
    # EQ.(2) DAVIS (1992) XC = MOLECULAR FRACTION OF C02
    # XC IN CIDDOR'S FORMULA IS IN PPM, CONVERTED TO A FRACTION
    MA = (28.9635 + 12.011e-6 * (XC - 400.0)) * 1.0e-3

    # USING CIDDOR EQ.(A1) WITH XV = 0 : DRY AIR COMPRESSIBILITY
    ZA = (
        1.0
        - (PSTDAIR / TKSTDAIR) * (A0 + (A1 * TCSTDAIR) + (A2 * TCSTDAIR * TCSTDAIR))
        + ((PSTDAIR / TKSTDAIR) * (PSTDAIR / TKSTDAIR) * DLC)
    )

    # USING EQ. 3 FOR DRY AIR
    RHOASX = (PSTDAIR * MA) / (ZA * R * TKSTDAIR)
    DENSITY_OF_STD_AIR = {"Ma": MA, "Z_(DRY AIR)": ZA, "RHOasx": RHOASX}

    # CALCULATE THE DENSITY OF WATER VAPOUR AT STANDARD CONDITIONS
    # (XV=1, TC=20, P=1333)

    # USING CIDDOR EQUATION A1 WITH XV = 1 AND PARTIAL
    # WATER VAPOUR PRESSURE
    ZWV = (
        1.0
        - (PSTDWV / TKSTDWV)
        * (
            A0
            + (A1 * TCSTDWV)
            + (A2 * TCSTDWV * TCSTDWV)
            + (B0 + B1 * TCSTDWV)
            + (C0 + C1 * TCSTDWV)
        )
        + ((PSTDWV / TKSTDWV) * (PSTDWV / TKSTDWV) * (DLC + E))
    )

    # USING EQ. 3
    RHOWS = (PSTDWV * MV) / (ZWV * R * TKSTDWV)

    DENSITY_OF_STD_WATER_VAPOUR = {"Z_(WATER VAPOUR)": ZWV, "RHOws": RHOWS}

    # CALCULATE THE DENSITY OF DRY AMBIENT AIR (TC= USER, P=USER)
    F = ALPHA + (BETA * P) + (GAMMA * TC * TC)
    XV = F * PV / P

    # USING CIDDOR EQUATION A1
    Z = (
        1.0
        - (P / TK)
        * (
            A0
            + (A1 * TC)
            + (A2 * TC * TC)
            + ((B0 + B1 * TC) * XV)
            + ((C0 + C1 * TC) * XV * XV)
        )
        + ((P / TK) * (P / TK) * (DLC + (E * XV * XV)))
    )

    RHOA = (P * MA * (1.0 - XV)) / (Z * R * TK)

    DENSITY_OF_DRY_AMBIENT_AIR = {"f": F, "Xv": XV, "Z_(AMB. AIR)": Z, "RHOa": RHOA}

    # CALCULATE THE DENSITY OF AMBIENT WATER VAPOUR
    # USING EQ. 3
    RHOW = (P * MV * XV) / (Z * R * TK)

    # CALCULATE THE PHASE REFRACTIVITY OF MOIST AIR FROM CIDDOR EQ.(4)
    TEMP1 = RHOA / RHOASX
    TEMP2 = RHOW / RHOWS

    DENSITY_OF_AMBIENT_WATER_VAPOUR = {"RHOw": RHOW, "Da": TEMP1, "Dv": TEMP2}

    TEMP1 = TEMP1 * NASX_1
    TEMP2 = TEMP2 * NWS_1

    DENSITY_OF_AMBIENT_WATER_VAPOUR["DaNasx"] = TEMP1
    DENSITY_OF_AMBIENT_WATER_VAPOUR["DvNvs"] = TEMP2

    NPROPP_1 = TEMP1 + TEMP2

    REFRACTIVITY_OF_MOIST_AIR = {"(NpropP-1)10E8": NPROPP_1}

    return NPROPP_1


def group_refractivity(LAMDA, TC, P, PV, XC=420):
    """
    Calculates the group refractivity of moist air using Ciddor's equations

    :param LAMDA: wavelength (micrometre)
    :param TC: temperature (degrees Celcius)
    :param P: pressure (hPa)
    :param PV: partial water vapour pressure (Pa)
    :param XC: carbon dioxide content (ppm)
    
    :return: NPROPG_1 - group refractivity = (NPROPgroup-1)10E8
    """

    (
        (K0, K1, K2, K3),
        (W0, W1, W2, W3),
        (ALPHA, BETA, GAMMA),
        (A0, A1, A2, B0, B1, C0, C1, DLC, E),
        (CF, R, MV),
        (TCSTDAIR, TKSTDAIR, PSTDAIR),
        (TCSTDWV, TKSTDWV, PSTDWV),
    ) = refractivity_constants()

    # CONVERT LAMDA TO SIGMA, THE PRESSURE UNITS TO PASCAL AND CREATE TC AND PP
    SIGMA = 1.0 / LAMDA

    # CONVERSION FROM hPa TO Pa AND C TO K
    P = P * 100.0
    PV = PV * 100.0
    TK = TC + 273.15

    # CALCULATE REFRACTIVITY OF STANDARD AIR (TC=15, P=101325)
    TEMP1 = SIGMA * SIGMA

    # CIDDOR EQ.(9)
    NGAS_1 = K1 * ((K0 + TEMP1) / ((K0 - TEMP1) * (K0 - TEMP1))) + K3 * (
        (K2 + TEMP1) / ((K2 - TEMP1) * (K2 - TEMP1))
    )

    # CIDDOR EQ.(2) FOR AMBIENT C02 CONTENT (XC)
    NGASX_1 = NGAS_1 * (1 + 0.534e-6 * (XC - 450.0))
    STD_AIR = {"(Ngas-1)10E8": NGAS_1, "(Ngasx-1)10E8": NGASX_1}

    # USING EQ. 2, CALCULATE REFRACTIVITY OF WATER VAPOUR AT STANDARD CONDITIONS
    # (TC=20, P=1333)
    TEMP2 = TEMP1 * TEMP1
    TEMP3 = TEMP1 * TEMP2

    # FOLLOWING OWENS (1967) AS AMENDED TO FIT BIRCH & DOWNS (1988)
    NGWS_1 = CF * (W0 + 3.0 * W1 * TEMP1 + 5.0 * W2 * TEMP2 + 7.0 * W3 * TEMP3)
    STD_WATER_VAPOUR = {"(Ngws-1)10E8": NGWS_1}

    # CALCULATE THE DENSITY OF STANDARD AIR (XV=0, TC=15, P=101325)
    # EQ.(2) DAVIS (1992) XC = MOLECULAR FRACTION OF C02
    # XC IN CIDDOR'S FORMULA IS IN PPM, CONVERTED TO A FRACTION
    MA = (28.9635 + 12.011e-6 * (XC - 400.0)) * 1.0e-3

    # USING CIDDOR EQ.(A1) WITH XV = 0 : DRY AIR COMPRESSIBILITY
    ZA = (
        1.0
        - (PSTDAIR / TKSTDAIR) * (A0 + (A1 * TCSTDAIR) + (A2 * TCSTDAIR * TCSTDAIR))
        + ((PSTDAIR / TKSTDAIR) * (PSTDAIR / TKSTDAIR) * DLC)
    )

    # USING EQ. 3 FOR DRY AIR
    RHOASX = (PSTDAIR * MA) / (ZA * R * TKSTDAIR)
    DENSITY_OF_STD_AIR = {"Ma": MA, "Z_(DRY AIR)": ZA, "RHOasx": RHOASX}

    # CALCULATE THE DENSITY OF WATER VAPOUR AT STANDARD CONDITIONS
    # (XV=1, TC=20, P=1333)

    # USING CIDDOR EQUATION A1 WITH XV = 1 AND PARTIAL
    # WATER VAPOUR PRESSURE
    ZWV = (
        1.0
        - (PSTDWV / TKSTDWV)
        * (
            A0
            + (A1 * TCSTDWV)
            + (A2 * TCSTDWV * TCSTDWV)
            + (B0 + B1 * TCSTDWV)
            + (C0 + C1 * TCSTDWV)
        )
        + ((PSTDWV / TKSTDWV) * (PSTDWV / TKSTDWV) * (DLC + E))
    )

    # USING EQ. 3
    RHOWS = (PSTDWV * MV) / (ZWV * R * TKSTDWV)

    DENSITY_OF_STD_WATER_VAPOUR = {"Z_(WATER VAPOUR)": ZWV, "RHOws": RHOWS}

    # CALCULATE THE DENSITY OF DRY AMBIENT AIR (TC= USER, P=USER)
    F = ALPHA + (BETA * P) + (GAMMA * TC * TC)
    XV = F * PV / P

    # USING CIDDOR EQUATION A1
    Z = (
        1.0
        - (P / TK)
        * (
            A0
            + (A1 * TC)
            + (A2 * TC * TC)
            + ((B0 + B1 * TC) * XV)
            + ((C0 + C1 * TC) * XV * XV)
        )
        + ((P / TK) * (P / TK) * (DLC + (E * XV * XV)))
    )

    # USING EQ. 3
    RHOA = (P * MA * (1.0 - XV)) / (Z * R * TK)

    DENSITY_OF_DRY_AMBIENT_AIR = {"f": F, "Xv": XV, "Z_(AMB. AIR)": Z, "RHOa": RHOA}

    # CALCULATE THE DENSITY OF AMBIENT WATER VAPOUR
    # USING EQ. 3
    RHOW = (P * MV * XV) / (Z * R * TK)

    # CALCULATE THE GROUP REFRACTIVITY OF MOIST AIR FROM CIDDOR EQ.(4)
    TEMP1 = RHOA / RHOASX
    TEMP2 = RHOW / RHOWS

    DENSITY_OF_AMBIENT_WATER_VAPOUR = {"RHOw": RHOW, "Da": TEMP1, "Dv": TEMP2}

    TEMP1 = TEMP1 * NGASX_1
    TEMP2 = TEMP2 * NGWS_1

    DENSITY_OF_AMBIENT_WATER_VAPOUR["DaNgasx"] = TEMP1
    DENSITY_OF_AMBIENT_WATER_VAPOUR["DvNgvs"] = TEMP2

    NPROPG_1 = TEMP1 + TEMP2

    REFRACTIVITY_OF_MOIST_AIR = {"(NpropG-1)10E8": NPROPG_1}

    return NPROPG_1


def humidity2part_water_vapour_press(H, TC):
    """
    Calculates the partial water vapour pressure from relative humidity and temperature

    :param H: humidity (%)
    :param TC: temperature (degrees Celcius)
    
    :return: PV - partial water vapour pressure (Pa)
    """

    (A, B, C, D) = (1.2378847e-5, -1.9121316e-2, 33.93711047, -6.3431645e3)
    TK = TC + 273.15

    # NOW CONVERT THE RELATIVE HUMIDITY TO PATIAL WATER VAPOUR PRESSURE (hPa)
    # SEE EQ.(22) GIACOMO 1982
    SVP = exp(A * TK * TK + B * TK + C + D / TK)
    PV = (H / 100.0) * SVP
    PV = PV / 100.0

    return PV

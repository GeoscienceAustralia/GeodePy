from math import sqrt, degrees, radians, sin, cos, asin
from conversions import dd2dms, dms2dd

def va_conv(verta_hp, slope_dist, height_inst=0, height_tgt=0):
    """
    Function to convert vertical angles (zenith distances) and slope distances into horizontal
    distances and changes in height. Instrument and Target heights can be entered to allow
    computation of zenith and slope distances between ground points.

    :param verta_hp: Vertical Angle from Instrument to Target, expressed in HP Format (DDD.MMSSSSSS)
    :param slope_dist: Slope Distance in metres
    :param height_inst: Height of Instrument. Optional - Default Value of 0m
    :param height_tgt: Height of Target. Optional - Default Value of 0m

    :return: verta_pt_hp: Vertical Angle between Ground Points, expressed in HP Format (DDD.MMSSSSSS)
    :return: slope_dist_pt: Slope Distance between Ground Points in metres
    :return: hz_dist: Horizontal Distance
    :return: delta_ht: Change in height between Ground Points in metres
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

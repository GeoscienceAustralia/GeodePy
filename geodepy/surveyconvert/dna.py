#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Data Converter - DynAdjust DNA v3 Output Tools
"""

from math import atan, degrees, sqrt
from geodepy.convert import dd2sec
from geodepy.statistics import k_val95


def dnaout_dirset(obslist, same_stdev=True, stdev_angular=1, stdev_pointing=0.001):
    # Test for Single Observation
    if len(obslist) < 2:
        return []
    fromlist = []
    for ob in obslist:
        fromlist.append(ob.from_id)
    fromlist = list(set(fromlist))
    if len(fromlist) != 1:
        raise ValueError('Direction Set requires obs with same from_id')
    else:
        pass
    dnaobs = []
    if same_stdev:
        # create first line using obslist[0]
        line1 = ('D '
                 + obslist[0].from_id.ljust(20)
                 + obslist[0].to_id.ljust(20)
                 + str(len(obslist)-1).ljust(20)
                 + str(obslist[0].hz_obs.degree).rjust(17)
                 + ' '
                 + str('%02d' % obslist[0].hz_obs.minute)
                 + ' '
                 + str('{:.3f}'.format(obslist[0].hz_obs.second).rjust(6, '0').ljust(8))
                 + ('{:.4f}'.format(stdev_angular)))  # add standard deviation
        dnaobs.append(line1)
        # create other lines using range(1:)
        for num in range(1, len(obslist)):
            line = ('D '
                    + ''.ljust(40)
                    + obslist[num].to_id.ljust(20)
                    + str(obslist[num].hz_obs.degree).rjust(17)
                    + ' '
                    + str('%02d' % obslist[num].hz_obs.minute)
                    + ' '
                    + str('{:.3f}'.format(obslist[num].hz_obs.second).rjust(6, '0').ljust(8))
                    + ('{:.4f}'.format(stdev_angular)))  # add standard deviation
            dnaobs.append(line)
    else:
        # calc stdev for first ob
        if obslist[0].sd_obs == 0:
            stdev_pt = str(stdev_angular)
        else:
            # TODO create GUM Std Dev Model for Dirset
            pointing_err_pt = dd2sec(degrees(atan(stdev_pointing / obslist[0].sd_obs)))
            stdev_pt = str((sqrt((pointing_err_pt ** 2) + (stdev_angular ** 2)))/sqrt(obslist[0].rounds))
        # create first line using obslist[0]
        line1 = ('D '
                 + obslist[0].from_id.ljust(20)
                 + obslist[0].to_id.ljust(20)
                 + str(len(obslist)-1).ljust(20)
                 + str(obslist[0].hz_obs.degree).rjust(17)
                 + ' '
                 + str('%02d' % obslist[0].hz_obs.minute)
                 + ' '
                 + str('{:.3f}'.format(obslist[0].hz_obs.second).rjust(6, '0').ljust(8))
                 + stdev_pt[0:5])  # add standard deviation
        dnaobs.append(line1)
        # create other lines using range(1:)
        for num in range(1, len(obslist)):
            if obslist[num].sd_obs == 0:  # keep
                stdev_pt = str(stdev_angular)
            else:
                pointing_err_pt = dd2sec(degrees(atan(stdev_pointing / obslist[num].sd_obs)))
                stdev_pt = str((sqrt((pointing_err_pt ** 2) + (stdev_angular ** 2)))/sqrt(obslist[num].rounds))
            line = ('D '
                    + ''.ljust(40)
                    + obslist[num].to_id.ljust(20)
                    + str(obslist[num].hz_obs.degree).rjust(17)
                    + ' '
                    + str('%02d' % obslist[num].hz_obs.minute)
                    + ' '
                    + str('{:.3f}'.format(obslist[num].hz_obs.second).rjust(6, '0').ljust(8))
                    + stdev_pt[0:5])  # add standard deviation
            dnaobs.append(line)
    return dnaobs


def dnaout_sd(obslist, stdev_distconst=0.001, stdev_distppm=1):
    dnaobs = []
    for observation in obslist:
        # GUM Uncertainty Estimation Model
        # Characterised Components (Uncertainty, k Value, Population (no of measurements - 1)
        type_a_esdm = observation.sd_stdev / sqrt(observation.rounds)
        if observation.rounds > 2:
            type_a_uncertainty = (type_a_esdm, 1, observation.rounds - 1)
        else:
            type_a_uncertainty = (type_a_esdm, 1, 1)
        inst_centring = (0.001, 2, 30)
        tgt_centring = (0.001, 2, 30)
        inst_uncertainty_const = (stdev_distconst, 2, 100)
        inst_uncertainty_ppm = (stdev_distppm, 2, 100)
        # Standard Error 1 sigma
        sterr_typea = type_a_uncertainty[0] / type_a_uncertainty[1]
        sterr_inst_centring = inst_centring[0] / inst_centring[1]
        sterr_tgt_centring = tgt_centring[0] / tgt_centring[1]
        sterr_inst_const = inst_uncertainty_const[0] / inst_uncertainty_const[1]
        sterr_inst_ppm = observation.sd_obs * ((inst_uncertainty_ppm[0] / inst_uncertainty_ppm[1]) / 1000000)
        combined_stdev1sigma = sqrt(sterr_typea ** 2 + sterr_inst_centring ** 2 + sterr_tgt_centring ** 2 +
                           sterr_inst_const ** 2 + sterr_inst_ppm ** 2)
        # Effective Degrees of Freedom (Welch-Satterthwaite) and
        vi_eff = (combined_stdev1sigma ** 4) / ((sterr_typea ** 4) / type_a_uncertainty[2] +
                                                (sterr_inst_centring ** 4) / inst_centring[2] +
                                                (sterr_tgt_centring ** 4) / tgt_centring[2] +
                                                (sterr_inst_const ** 4) / inst_uncertainty_const[2] +
                                                (sterr_inst_ppm ** 4) / inst_uncertainty_ppm[2])
        stdev = combined_stdev1sigma * k_val95(int(vi_eff))
        # Exclude slope distances of 0m
        if observation.sd_obs == 0:
            pass
        else:
            # Write line output in DNA v3 Format
            line = ('S '
                    + observation.from_id.ljust(20)
                    + observation.to_id.ljust(20)
                    + ''.ljust(19)
                    + ('{:.4f}'.format(observation.sd_obs)).rjust(9)
                    + ''.ljust(21)
                    + ('{:.4f}'.format(stdev)).ljust(8)
                    + ('{:.4f}'.format(observation.inst_height).ljust(7))
                    + ('{:.4f}'.format(observation.target_height)))
            dnaobs.append(line)
    return dnaobs


def dnaout_va(obslist, same_stdev=True, stdev=1, stdev_pointing=0.001):
    dnaobs = []
    if same_stdev:
        for observation in obslist:
            # Format Line
            line = ('V '
                    + observation.from_id.ljust(20)
                    + observation.to_id.ljust(20)
                    + ''.ljust(34)
                    + str(observation.va_obs.degree).rjust(3) + ' '
                    + str(observation.va_obs.minute).rjust(2, '0') + ' '
                    + ('{:.3f}'.format(observation.va_obs.second).rjust(6, '0')).ljust(8)
                    + ('{:.4f}'.format(stdev)).ljust(8)
                    + ('{:.4f}'.format(observation.inst_height)).ljust(7)
                    + ('{:.4f}'.format(observation.target_height)))
            dnaobs.append(line)
    else:
        for observation in obslist:
            # Calc Standard Deviation (in seconds) based on pointing error (in metres)
            # and angular standard deviation (in seconds)
            if observation.sd_obs == 0:
                stdev_pt = str(stdev)
            else:
                # TODO create GUM Std Dev Model for VA
                pointing_err_pt = dd2sec(degrees(atan(stdev_pointing / observation.sd_obs)))
                stdev_pt = str((sqrt((pointing_err_pt ** 2) + (stdev ** 2)))/sqrt(observation.rounds))
            # Format Line
            line = ('V '
                    + observation.from_id.ljust(20)
                    + observation.to_id.ljust(20)
                    + ''.ljust(34)
                    + str(observation.va_obs.degree).rjust(3)
                    + ' '
                    + str(observation.va_obs.minute).rjust(2, '0')
                    + ' '
                    + str(('{:.3f}'.format(observation.va_obs.second).rjust(6, '0')).ljust(8))
                    + stdev_pt[0:5].ljust(8)  # add standard deviation
                    + ('{:.4f}'.format(observation.inst_height)).ljust(7)
                    + ('{:.4f}'.format(observation.target_height)))
            dnaobs.append(line)
    return dnaobs

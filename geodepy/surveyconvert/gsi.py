#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Data Converter - Leica GSI Format Import Tools
Note: These tools are designed to work with the Leica GSI Format File xxxxxxxxxxxx.frt
"""

import os
from datetime import datetime
from geodepy.convert import DMSAngle
from geodepy.survey import first_vel_params
from geodepy.surveyconvert.config import readconfig, renameobs, removeobs, first_vel_cfg, stdev_cfg
from geodepy.surveyconvert.classtools import Coordinate, InstSetup, Observation, reducesetup, first_vel_observations
from geodepy.surveyconvert.dna import dnaout_dirset, dnaout_va, dnaout_sd


def gsi2msr(path, cfg_path=None):
    """
    Converts .gsi format survey observations to DNA v3 .msr for use with DynAdjust
    :param path: .gsi file path
    :param cfg_path: .gpy conversion configuration file
    :return: DNA v3 .msr file (same directory as source .fbk file)
    """
    gsi_project = gsi2class(readgsi(path))
    # Read config file
    if cfg_path is not None:
        cfg = readconfig(cfg_path)
        # Rename obs as per config file
        gsi_project = renameobs(cfg, gsi_project)
        # Remove obs as per config file
        gsi_project = removeobs(cfg, gsi_project)
        # Get First Velocity Correction Observations
        first_vel_obs = first_vel_cfg(cfg)
        # Get Standard Deviation Parameters
        stdev_params = stdev_cfg(cfg)
    else:
        first_vel_obs = None
        stdev_params = None
    # Reduce observations in setups
    for setup in gsi_project:
        reduced_obs = reducesetup(setup.observation, strict=False, zerodist=True)
        setup.observation = reduced_obs
    # Perform First Velocity Correction
    if first_vel_obs is not None:
        # Use wavelength and standard atmospheric parameters to get Parameters C and D
        params = first_vel_params(first_vel_obs[0])
        for setup in gsi_project:
            corrected_obs = first_vel_observations(setup.observation,
                                                   params,
                                                   first_vel_obs[1],  # Observed Temperature
                                                   first_vel_obs[2],  # Observed Pressure
                                                   first_vel_obs[3])  # Observed Relative Humidity
            setup.observation = corrected_obs
    # Produce Measurement format data from setups
    msr_raw = []
    if stdev_params is None:
        stdev_params = (1, 0.001, 0.001, 1)  # Default standard deviation parameters
    for setup in gsi_project:
        dna_dirset = dnaout_dirset(setup.observation, False, stdev_params[0], stdev_params[1])
        dna_va = dnaout_va(setup.observation, False, stdev_params[0], stdev_params[1])
        dna_sd = dnaout_sd(setup.observation, stdev_params[2], stdev_params[3])
        msr_raw.append(dna_dirset + dna_va + dna_sd)
    # Build msr header
    dircount = 0
    vacount = 0
    sdcount = 0
    for group in msr_raw:
        for line in group:
            if line.startswith('D'):
                dircount = 1
            elif line.startswith('V'):
                vacount += 1
            elif line.startswith('S'):
                sdcount += 1
    obscount = dircount + vacount + sdcount
    now = datetime(2020, 1, 1)
    date = (str(now.day).rjust(2, '0') + '.'
            + str(now.month).rjust(2, '0') + '.'
            + str(now.year))
    header = ('!#=DNA 3.01 MSR'.ljust(19)
              + date.ljust(19)
              + 'GDA94'.ljust(9)
              + date.ljust(18)
              + str(obscount))
    msr = [line for sublist in msr_raw for line in sublist]
    msr = [header] + msr
    # Output MSR File
    if cfg_path is not None:
        fn, ext = os.path.splitext(cfg_path)
        msr_fn = fn + '.msr'
    else:
        fn, ext = os.path.splitext(path)
        msr_fn = fn + '.msr'
    with open(msr_fn, 'w+') as msr_file:
        for line in msr:
            msr_file.write(line + '\n')
    # output will be dna msr file
    return msr


def gsi2stn(path, utmzone, cfg_path=None):
    """
    Converts coordinate list file (.txt) associated with .fbk file into DNA v3 stn file
    :param path: .txt co-ordinate list associated with .fbk file
    :param utmzone: UTM Coordinate Zone as string i.e. 'S56' for Southern Hemisphere Zone 56
    :param cfg_path: .gpy DNA conversion config file
    :return: DNA v3 .stn file (same directory as source .fbk file)
    """
    # Read Data from file
    with open(path) as raw:
        ptlist = raw.readlines()
    # Split comma separated values, replace blanks with zeroes
    for num, line in enumerate(ptlist):
        ptlist[num] = ptlist[num].strip()
        ptlist[num] = ptlist[num].split(',')
        for place, item in enumerate(ptlist[num]):
            if item == '0':
                ptlist[num][place] = '0'
    # Read Config file contents
    if cfg_path is not None:
        cfg_list = readconfig(cfg_path)
        constrain_list = []
        rename_list = []
        remove_list = []
        for group in cfg_list:
            group_header = group[0].lower()
            if group_header.startswith('constrain'):
                constrain_list = group[1:]
            elif group_header.startswith('rename'):
                rename_list = group[1:]
            elif group_header.startswith('remove'):
                remove_list = group[1:]

            # Perform Block Shift of Coordinates as Specified in Config
            elif group_header.startswith('blockshift'):
                shift_list = group[1:]
                delta_east = shift_list[0]
                delta_north = shift_list[1]
                delta_up = shift_list[2]
                for pt in ptlist:
                    pt[1] = str('{:.4f}'.format(float(pt[1]) + float(delta_east)))
                    pt[2] = str('{:.4f}'.format(float(pt[2]) + float(delta_north)))
                    pt[3] = str('{:.4f}'.format(float(pt[3]) + float(delta_up)))

        # Rename Points as per Config file
        for num, i in enumerate(rename_list):
            rename_list[num] = i.split(',')
        for pt_rename in rename_list:
            for num, pt in enumerate(ptlist):
                if pt[0] == pt_rename[0]:
                    ptlist[num][0] = pt_rename[1]

        # Remove Points as per Config file
        for pt_rem in remove_list:
            for num, pt in enumerate(ptlist):
                if pt[0] == pt_rem:
                    del ptlist[num]

        # Set Points in Config Constrains list to 'CCC'
        for num, pt in enumerate(ptlist):
            ptlist[num] = ['FFF'] + pt
        for pt_constrain in constrain_list:
            for num, pt in enumerate(ptlist):
                if pt[1] == pt_constrain:
                    ptlist[num][0] = 'CCC'

    # Remove Duplicate Station Points
    cleanlist = []
    deduplist = []
    dedupnum = []
    for num, pt in enumerate(ptlist):
        if pt[1] not in deduplist:
            deduplist.append(pt[1])
            dedupnum.append(num)
    for num in dedupnum:
        cleanlist.append(ptlist[num])

    # Write header line
    stn = []
    now = datetime(2020, 1, 1)
    header = ('!#=DNA 3.01 STN    '
              + str(now.day).rjust(2, '0') + '.'
              + str(now.month).rjust(2, '0') + '.'
              + str(now.year)
              + 'GDA94'.rjust(14)
              + (str(len(cleanlist))).rjust(25))
    stn.append(header)

    # Write line strings in stn format
    for pt in cleanlist:
        line = (pt[1].ljust(20)  # Pt ID
                + pt[0]  # Constraint
                + ' UTM '  # Projection
                + pt[2].rjust(13)  # Easting
                + pt[3].rjust(18)  # Northing
                + pt[4].rjust(16)  # Elevation
                + utmzone.rjust(15))  # Hemisphere/Zone input
        stn.append(line)
    # Write line strings to file
    if cfg_path is not None:
        fn, ext = os.path.splitext(cfg_path)
        stn_fn = fn + '.stn'
    else:
        fn, ext = os.path.splitext(path)
        stn_fn = fn + '.stn'
    with open(stn_fn, 'w+') as stn_file:
        for line in stn:
            stn_file.write(line + '\n')
    return stn


def readgsi(filepath):
    """
    Takes in a gsi file (GA_Survey2.frt) and returns a list
    of stations with their associated observations.
    :param filepath: full directory of .gsi file
    :return: gsi data in list form
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
        stn_index = []
        for num, line in enumerate(gsidata):
            # Create list of line numbers of station records
            # (Only stations have '84..' string)
            if '84..' in line:
                stn_index.append(num + 1)
        stn_index.append(len(gsidata))
        gsi_listbystation = []
        # Create lists of gsi data with station records as first element
        for j in range(0, len(stn_index)):
            gsi_listbystation.append(gsidata[(stn_index[j - 1] - 1):(stn_index[j] - 1)])
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
        try:
            wordstart = str.find(linestring, word_id)
            if wordstart == -1:
                raise ValueError
        except ValueError:
            print('ValueError: GSI record type ' + word_id + ' not found\n'
                  'Line Data: ' + linestring)
            return None
        word_val = linestring[(wordstart + 7):(wordstart + 23)]
        word_val = word_val.lstrip('0')
        if word_val == '':
            return 0
        else:
            return int(word_val)

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
        return DMSAngle(degree, minute, second * 10)

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
        if 0 < degrees <= 180:
            face = 'FL'
        elif 180 < degrees <= 360:
            face = 'FR'
        else:
            face = 'FL'
        return face, DMSAngle(degrees, minutes, seconds * 10)

    gsi_project = []
    for record in gsi_list:
        from_stn = parse_ptid(record[0])
        obs_list = []
        for line in record:
            if '31..' in line:
                to_stn = parse_ptid(line)
                hz = parse_hz(line)
                face, vert = parse_vert(line)
                slope = parse_slope(line)
                tgtht = parse_tgtht(line)
                instht = parse_instht(line)
                obs = Observation(from_stn, to_stn, instht, tgtht,
                                  face, hz, vert, slope)
                obs_list.append(obs)
        if '84..' in record[0]:
            pt_id = parse_ptid(record[0])
            easting = parse_easting(record[0])
            northing = parse_northing(record[0])
            elev = parse_elev(record[0])
            coord = Coordinate(pt_id, 'utm', 'gda94', 'gda94',
                               '2018.1', easting, northing, elev)
            setup = InstSetup(pt_id, coord)
            for i in range(0, len(obs_list)):
                setup.addobs(obs_list[i])
        gsi_project.append(setup)
    return gsi_project


def readgsiword16(linestring, word_id):
    wordstart = str.find(linestring, word_id)
    word_val = linestring[(wordstart + 7):(wordstart + 23)]
    word_val = 0.0 if word_val.lstrip('0') == '' else int(word_val.lstrip('0'))
    return word_val

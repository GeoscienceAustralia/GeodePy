#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Data Converter - Leica GSI Format Import Tools
Note: These tools are designed to work with the Leica GSI Format File xxxxxxxxxxxx.frt
"""

import os
from datetime import datetime
from geodepy.convert import DMSAngle
from geodepy.surveyconvert.config import readconfig, renameobs, removeobs
from geodepy.surveyconvert.classtools import Coordinate, InstSetup, Observation, reducesetup
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
    # Reduce observations in setups
    for setup in gsi_project:
        reduced_obs = reducesetup(setup.observation, strict=False, zerodist=False)
        setup.observation = reduced_obs
    # Produce Measurement format data from setups
    msr_raw = []
    for setup in gsi_project:
        dna_dirset = dnaout_dirset(setup.observation, same_stdev=False)
        dna_va = dnaout_va(setup.observation, same_stdev=False)
        dna_sd = dnaout_sd(setup.observation)
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
    fn, ext = os.path.splitext(path)
    msr_fn = fn + '.msr'
    with open(msr_fn, 'w+') as msr_file:
        for line in msr:
            msr_file.write(line + '\n')
    # output will be dna msr file
    return msr


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

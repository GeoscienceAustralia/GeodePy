#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Data Converter - Geomax Zoom90 Theodelite FBK Format Import Tools
"""

import os
import numpy as np
from geodepy.convert import DMSAngle
from geodepy.survey import first_vel_params
from geodepy.surveyconvert.config import readconfig, renameobs, removeobs, first_vel_cfg
from geodepy.surveyconvert.classtools import Coordinate, InstSetup, Observation, reducesetup, first_vel_observations
from geodepy.surveyconvert.dna import dnaout_dirset, dnaout_va, dnaout_sd


def fbk2msr(path, cfg_path, strict=False, zerodist=False, same_stdev=False):
    """
    Converts .fbk format survey observations to DNA v3 .msr for use with DynAdjust
    :param path: .fbk file path
    :param cfg_path: .gpy DNA conversion config file
    :param strict: If True, all single-face Obs are ignored. If False, all
    single-face Obs are included and converted to Face Left
    :param zerodist: If True, Obs with Slope Distance of Zero are included.
    If False, these are ignored
    :param same_stdev:
    :return: DNA v3 .msr file (same directory as source .fbk file)
    """
    fbk_project = fbk2class(readfbk(path))
    # Read config file
    cfg = readconfig(cfg_path)
    # Rename obs as per config file
    fbk_project = renameobs(cfg, fbk_project)
    # Remove obs as per config file
    fbk_project = removeobs(cfg, fbk_project)
    # Get First Velocity Correction Observations
    first_vel_obs = first_vel_cfg(cfg)
    # Reduce observations in setups
    for setup in fbk_project:
        reduced_obs = reducesetup(setup.observation, strict, zerodist)
        setup.observation = reduced_obs
    # Perform First Velocity Correction
    if first_vel_obs is not None:
        # Use wavelength and standard atmospheric parameters to get Parameters C and D
        params = first_vel_params(first_vel_obs[0])
        for setup in fbk_project:
            corrected_obs = first_vel_observations(setup.observation,
                                                   params,
                                                   first_vel_obs[1],  # Observed Temperature
                                                   first_vel_obs[2],  # Observed Pressure
                                                   first_vel_obs[3])  # Observed Relative Humidity
            setup.observation = corrected_obs
    # Produce Measurement format data from setups
    msr_raw = []
    for setup in fbk_project:
        dna_dirset = dnaout_dirset(setup.observation, same_stdev)
        dna_va = dnaout_va(setup.observation, same_stdev)
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
    day, month, year = fbkdate(path)
    date = day + '.' + month + '.' + year
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


def writestn(file, utmzone):
    """
    Converts coordinate list file (.txt) associated with .fbk file into DNA v3 stn file
    :param file: .txt co-ordinate list associated with .fbk file
    :return: DNA v3 .stn file (same directory as source .fbk file)
    """
    # Read Data from file
    with open(file) as raw:
        ptlist = raw.readlines()
    # Split comma separated values
    for num, line in enumerate(ptlist):
        ptlist[num] = ptlist[num].strip()
        ptlist[num] = ptlist[num].split(',')
    # Get Date from last point
    for i in ptlist[-1]:
        if i.startswith('DATE'):
            month = i[5:7]
            day = i[8:10]
            year = i[11:15]

    # Read Config file contents
    fn, ext = os.path.splitext(file)
    cfg_list = readconfig(fn + '.gpy')
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
    header = ('!#=DNA 3.01 STN    '
              + day + '.'
              + month + '.'
              + year
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
                + utmzone.rjust(15)  # Hemisphere/Zone input
                + ' '
                + pt[5])  # Pt Description
        stn.append(line)
    # Write line strings to file
    fn, ext = os.path.splitext(file)
    stn_fn = fn + '.stn'
    with open(stn_fn, 'w+') as stn_file:
        for line in stn:
            stn_file.write(line + '\n')
    return stn


def stripfile(filedata, listofterms):
    """
    Creates a list with lines starting with strings from a list
    :param filedata: File Data in list form
    :param listofterms: list of strings to use as search terms
    :return: list of file lines starting with strings from list of terms
    """
    datawithterms = []
    for line in filedata:
        if type(line) == str:
            for i in listofterms:
                if line.startswith(i):
                    datawithterms.append(line)
        elif type(line) == list:
            for i in listofterms:
                if line[0] == i:
                    datawithterms.append(line)
    return datawithterms


def addprismht(fbklist):
    prismlist = []
    rowto = []
    numlines = 0
    # Build array with prism height, startline and endline
    for num, line in enumerate(fbklist, 1):
        numlines += 1
        if line.startswith('PRISM'):
            prismlist.append([line[6:11], num])
    prismarray = np.asarray(prismlist)
    prismarrayrows = prismarray.shape[0]
    for num, row in enumerate(prismarray, 1):
        if num > prismarrayrows:
            break
        elif num == prismarrayrows:
            rowto.append(numlines)
        else:
            rowto.append(prismarray[num, 1])
    rowtoarray = np.asarray(rowto)
    prismrange = np.append(prismarray,
                           rowtoarray.reshape(rowtoarray.size, 1),
                           axis=1)
    # Add Prism Heights to each observation in file
    filecontents = [x.strip() for x in fbklist]
    for prism in prismrange:
        for i in range(int(prism[1]), int(prism[2])):
            if not (not filecontents[i].startswith('F1')
                    and not filecontents[i].startswith('F2')):
                filecontents[i] = filecontents[i] + ' ' + prism[0]
    filecontents = [x + '\n' for x in filecontents]
    return filecontents


def readfbk(filepath):
    with open(filepath) as f:
        # Remove non-obs data
        stage1 = stripfile(f, ['PRISM', 'NEZ', 'STN', 'F1', 'F2'])
        # Add prism heights to each observation
        stage2 = addprismht(stage1)
        # Remove Spaces from Obs Descriptions (inside "")

        def wscommstrip(string):
            string_list = list(string)
            commrange = [pos for pos, char in enumerate(string) if char == '\"']
            if len(commrange) != 2:
                return string
            else:
                for char in range(commrange[0] + 1, commrange[1]):
                    if string_list[char] == ' ':
                        string_list[char] = '_'
                string_us = ''.join(string_list)
                return string_us

        stage3 = []
        for line in stage2:
            stage3.append(wscommstrip(line))
        # Split obs
        stage4 = []
        for num, i in enumerate(stage3):
            stage4.append(stage3[num].split())
        # Add coordinates in file to stations
        coordlist = []
        for i in stage4:
            if i[0] == 'NEZ':
                coordlist.append(i)
        stage5 = stripfile(stage4, ['STN', 'F1', 'F2'])
        # Add Coord to setup without coord info (broken)
        # for coord in coordlist:
        #     for num, line in enumerate(stage4):
        #         if line[1] == coord[1]:
        #             stage5[num] = line + coord[2:]
        # Group by Setup
        stn_index = [0]
        for num, i in enumerate(stage5, 1):
            # Create list of line numbers of station records
            # (Only stations have 'STN' string)
            if 'STN' in i:
                lnid = num
                stn_index.append(lnid)
        stn_index.append(len(stage5))
        fbk_listbystation = []
        # Create lists of fbk data with station records as first element
        for i in range(0, len(stn_index) - 1):
            fbk_listbystation.append(stage5[stn_index[i] - 1:stn_index[i + 1] - 1])
        del fbk_listbystation[0]
        # for i in range(0, (len(stn_index) - 1)):
        #     fbk_listbystation.append(stage5[(stn_index[i]) - 1:(stn_index[i + 1])])
        # del fbk_listbystation[0]
    return fbk_listbystation


def fbkdate(filepath):
    with open(filepath) as f:
        day = None
        month = None
        year = None
        for line in f:
            if line.startswith('! DT'):
                month = line[4:6]
                day = line[7:9]
                year = line[10:14]
                break
    return day, month, year


def fbk2class(fbk_list):

    def parse_angle(obs_string):
        dms = float(obs_string)
        degmin, second = divmod(abs(dms) * 1000, 10)
        degree, minute = divmod(degmin, 100)
        return DMSAngle(degree, minute, second * 10)

    fbk_project = []
    for setup_list in fbk_list:
        obs_list = []
        if setup_list[0][0] == 'STN' and len(setup_list[0]) <= 3:
            # This is the station information part
            from_id = setup_list[0][1]
            inst_height = float(setup_list[0][2])
            coord = Coordinate(from_id, 'utm', 'gda94', 'gda94',
                               '2018.1', 0, 0, 0)
            setup = InstSetup(from_id, coord)
        elif setup_list[0][0] == 'STN' and len(setup_list[0]) > 3:
            from_id = setup_list[0][1]
            inst_height = float(setup_list[0][2])
            east = float(setup_list[0][3])
            north = float(setup_list[0][4])
            elev = float(setup_list[0][5])
            coord = Coordinate(from_id, 'utm', 'gda94', 'gda94',
                               '2018.1', east, north, elev)
            setup = InstSetup(from_id, coord)
        for record in setup_list:
            if record[0] == 'F1' or record[0] == 'F2':
                """
                This is the obs information part
                from_id, to_id, inst_height, target_height
                face
                hz_obs
                va_obs
                sd_obs
                hz_dist
                vert_dist
                """
                # Read Face
                if int(float(record[5])) in range(0, 180):
                    face = 'FL'
                elif int(float(record[5])) in range(180, 360):
                    face = 'FR'
                else:
                    ValueError('Invalid Vertical Angle in ' + record)
                to_id = record[2]  # Read To ID
                hz_obs = parse_angle(record[3])  # 3 is Hz Ob (HP)
                sd_obs = float(record[4])
                va_obs = parse_angle(record[5])  # 5 is Vert Ob (HP)
                target_height = float(record[7])
                obs = Observation(from_id, to_id,
                                  inst_height, target_height,
                                  face, hz_obs, va_obs, sd_obs)
                obs_list.append(obs)
            # else:
            #     raise ValueError('Unexpected format found')
        for i in obs_list:
                setup.addobs(i)
        fbk_project.append(setup)
    return fbk_project

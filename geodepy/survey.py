#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Module
"""

import os
import numpy as np
import itertools
from math import sqrt, degrees, radians, sin, cos, asin, atan
from datetime import datetime
from geodepy.convert import dec2hp, hp2dec, dd2sec


# Defines a bunch of classes required to convert GSI (or any other format) to DynaNet v3 Format


class AngleObs(object):
    def __init__(self, degree=0, minute=0, second=0):
        self.degree = abs(int(degree))
        self.minute = abs(int(minute))
        self.second = abs(round(float(second), 3))

    def __repr__(self):
        return repr(self.degree) + 'd '\
               + repr(self.minute) + '\' '\
               + repr(self.second) + '\"'

    def decimal(self):
        dd = abs(self.degree + self.minute / 60 + self.second / 3600)
        return dd if self.degree >= 0 else -dd

    def hp(self):
        hp = abs(self.degree + self.minute / 100 + self.second / 10000)
        return hp if self.degree >= 0 else -hp

    def __add__(self, other):
        degreeadd = self.degree + other.degree
        minuteadd = abs(self.minute) + abs(other.minute)
        secondadd = abs(self.second) + abs(other.second)
        while secondadd >= 60:
            minuteadd += 1
            secondadd -= 60
        while minuteadd >= 60:
            degreeadd += 1
            minuteadd -= 60
        return AngleObs(degreeadd, minuteadd, secondadd)

    def __sub__(self, other):
        degreesub = self.degree - other.degree
        minutesub = abs(self.minute) - abs(other.minute)
        secondsub = abs(self.second) - abs(other.second)
        while secondsub < 0:
            secondsub += 60
            minutesub -= 1
        while minutesub < 0:
            minutesub += 60
            degreesub -= 1
        return AngleObs(degreesub, minutesub, secondsub)

    def __truediv__(self, other):
        degreediv, degreerem = divmod(self.degree, other)
        minutediv, minuterem = divmod(self.minute, other)
        seconddiv, secondrem = divmod(self.second, other)
        minutediv = minutediv + ((degreerem / other) * 60)
        seconddiv = seconddiv + ((minuterem / other) * 60) + (secondrem / other)
        return AngleObs(degreediv, minutediv, seconddiv)


class Coordinate(object):
    def __init__(self, pt_id, system, hz_datum,
                 vert_datum, epoch, x=0.0, y=0.0, z=0.0):
        self.pt_id = pt_id
        self.system = system
        self.hz_datum = hz_datum
        self.vert_datum = vert_datum
        self.epoch = epoch
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return ('{' + str(self.hz_datum) + ', '
                + str(self.x) + ', '
                + str(self.y) + ', '
                + str(self.z) + '}')


class InstSetup(object):
    def __init__(self, pt_id, coordinate, observation=None):
        self.pt_id = pt_id
        self.coordinate = coordinate
        if observation is None:
            self.observation = []
        else:
            self.observation = []
            self.observation.append(observation)

    def __repr__(self):
        return ('{InstSetup: ' + repr(self.pt_id)
                + ' ' + repr(self.coordinate)
                + '}\n Observations:\n'
                + repr(self.observation)) + '\n\n'

    def addobs(self, observation):
        return self.observation.append(observation)

    def __iter__(self):
        return self


class Observation(object):
    def __init__(self, from_id, to_id, inst_height=0.0, target_height=0.0,
                 face='FL', hz_obs=AngleObs(0, 0, 0), va_obs=AngleObs(0, 0, 0),
                 sd_obs=0.0, hz_dist=0.0, vert_dist=0.0):
        self.from_id = from_id
        self.to_id = to_id
        self.inst_height = inst_height
        self.target_height = target_height
        self.face = face
        self.hz_obs = hz_obs
        self.va_obs = va_obs
        self.sd_obs = sd_obs
        self.hz_dist = hz_dist
        self.vert_dist = vert_dist

    def __repr__(self):
        return ('{from: ' + repr(self.from_id)
                + 'to: ' + repr(self.to_id)
                + '; inst_height ' + repr(self.inst_height)
                + '; target_height ' + repr(self.target_height)
                + '; face ' + repr(self.face)
                + '; hz_obs ' + repr(self.hz_obs)
                + '; va_obs ' + repr(self.va_obs)
                + '; sd_obs ' + repr(self.sd_obs)
                + '}')

    def changeface(self):
        # Change Horizontal Angle
        if 0 <= self.hz_obs.degree < 180:
            hz_switch = self.hz_obs + AngleObs(180)
        elif 180 <= self.hz_obs.degree < 360:
            hz_switch = self.hz_obs - AngleObs(180)
        else:
            raise ValueError('Horizontal Angle out of range (0 to 360 degrees)')
        # Change Vertical Angle
        if 0 <= self.va_obs.degree < 360:
            va_switch = AngleObs(360) - self.va_obs
        else:
            raise ValueError('Vertical Angle out of range (0 to 360 degrees)')
        # Change Face Label
        newface = None
        if self.face == 'FL':
            newface = 'FR'
        elif self.face == 'FR':
            newface = 'FL'
        return Observation(self.from_id,
                           self.to_id,
                           self.inst_height,
                           self.target_height,
                           newface,
                           hz_switch,
                           va_switch,
                           self.sd_obs,
                           self.hz_obs,
                           self.vert_dist)


# Functions to configure data in a project using a DNA conversion config file
def readconfig(path):
    """
    Read data in from a DNA conversion config file to a list
    :param path: .gpy file path
    :return: config data in a list
    """
    with open(path) as f:
        fstring = f.read()
        cfg_list = fstring.split('\n\n')
        for num, linegroup in enumerate(cfg_list):
            cfg_list[num] = cfg_list[num].lstrip('\n')
            cfg_list[num] = cfg_list[num].rstrip('\n')
            cfg_list[num] = cfg_list[num].splitlines()
    return cfg_list


def renameobs(cfg_list, project):
    # Find entry in cfg_list startswith rename
    rename_list = []
    for group in cfg_list:
        group_header = group[0].lower()
        if group_header.startswith('rename'):
            rename_list = group[1:]
    for num, i in enumerate(rename_list):
        rename_list[num] = i.split(',')
    # if old in obs.to_id or obs.from_id, replace with new
    for setup in project:
        for obs in setup.observation:
            for rename_pair in rename_list:
                if rename_pair[0] == obs.from_id:
                    obs.from_id = rename_pair[1]
                elif rename_pair[0] == obs.to_id:
                    obs.to_id = rename_pair[1]
    # if old in setup info, replace with new
    for setup in project:
        for rename_pair in rename_list:
            if setup.pt_id == rename_pair[0]:
                setup.pt_id = rename_pair[1]
    return project


def removeobs(cfg_list, project):
    # Find entry in cfg_list startswith remove
    remove_list = []
    for group in cfg_list:
        group_header = group[0].lower()
        if group_header.startswith('remove'):
            remove_list = group[1:]
    for remove_id in remove_list:
        for setup in project:
            if remove_id == setup.pt_id:
                del setup
        for setup in project:
            for num, obs in enumerate(setup.observation):
                if remove_id == obs.from_id:
                    del setup.observation[num]
                elif remove_id == obs.to_id:
                    del setup.observation[num]
    return project


# def dist_sd():
#     pass
#
#
# def pointing_sd():
#     pass


# Functions to read in data from fbk format (Geomax Zoom90 Theodolite used
# in Surat Survey

def fbk2msr(path, cfg_path):
    """
    Converts .fbk format survey observations to DNA v3 .msr for use with DynAdjust
    :param path: .fbk file path
    :return: DNA v3 .msr file (same directory as source .fbk file)
    """
    fbk_project = fbk2class(readfbk(path))
    # Read config file
    cfg = readconfig(cfg_path)
    # Rename obs as per config file
    fbk_project = renameobs(cfg, fbk_project)
    # Remove obs as per config file
    fbk_project = removeobs(cfg, fbk_project)
    # Reduce observations in setups
    for setup in fbk_project:
        reduced_obs = reducesetup(setup.observation)
        setup.observation = reduced_obs
    # Produce Measurement format data from setups
    msr_raw = []
    for setup in fbk_project:
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

# Debug - Example fbk file
# testfbk = '\\geodepy\\tests\\resources\\Site01-152.fbk'


def writestn(file):
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

    # Write header line
    stn = []
    header = ('!#=DNA 3.01 STN    '
              + day + '.'
              + month + '.'
              + year
              + 'GDA94'.rjust(14)
              + (str(len(ptlist))).rjust(25))
    stn.append(header)

    # Write line strings in stn format
    for pt in ptlist:
        line = (pt[1].ljust(20)  # Pt ID
                + pt[0]  # Constraint
                + ' UTM'  # Projection
                + pt[2].rjust(13)  # Easting
                + pt[3].rjust(18)  # Northing
                + pt[4].rjust(16)  # Elevation
                + 'S56'.rjust(16)  # Hemisphere/Zone input
                + ' '
                + pt[4])  # Pt Description
        stn.append(line)
    # Write line strings to file
    fn, ext = os.path.splitext(file)
    stn_fn = fn + '.stn'
    with open(stn_fn, 'w+') as stn_file:
        for line in stn:
            stn_file.write(line + '\n')
    return stn

# Debug - Example txt file
# testcoord = '\\geodepy\\tests\\resources\\Site01-152.txt'


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
        for i in range(0,len(stn_index) - 1):
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
        return AngleObs(degree, minute, second * 10)

    fbk_project = []
    for setup_list in fbk_list:
        obs_list = []
        if setup_list[0][0] == 'STN' and len(setup_list[0]) <= 3:
            # This is the station information part
            from_id = setup_list[0][1]
            inst_height = setup_list[0][2]
            coord = Coordinate(from_id, 'utm', 'gda94', 'gda94',
                               '2018.1', 0, 0, 0)
            setup = InstSetup(from_id, coord)
        elif setup_list[0][0] == 'STN' and len(setup_list[0]) > 3:
            from_id = setup_list[0][1]
            inst_height = setup_list[0][2]
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


# Functions to read data to classes from Leica GSI format file (GA_Survey2.frt)


def gsi2msr(path):
    gsi_class = gsi2class(readgsi(path))
    # Reduce observations in setups
    for setup in gsi_class:
        reduced_obs = reducesetup(setup.observation)
        setup.observation = reduced_obs
    # Produce Measurement format data from setups
    msr_raw = []
    for setup in gsi_class:
        dna_dirset = dnaout_dirset(setup.observation, same_stdev=True)
        dna_va = dnaout_va(setup.observation, same_stdev=True)
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
    now = datetime.now()
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
        return AngleObs(degree, minute, second * 10)

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
        return face, AngleObs(degrees, minutes, seconds * 10)

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


# Functions to write out station and obs data to DNA format


def meanfaces(ob1, ob2):
    """
    Take two Observations and return their mean Face Left Sense Observation
    If one Observation and one None, return Face Left Sense Observation
    :param ob1: Observation Object (or None)
    :param ob2: Observation Object (or None)
    :return: Meaned Observation Object of ob1 and ob2 (Face Left Sense)
    """
    if type(ob1) != Observation and type(ob1) != type(None):
        raise TypeError('Invalid Input Type (ob1)')
    elif type(ob2) != Observation and type(ob2) != type(None):
        raise TypeError('Invalid Input Type (ob2)')
    elif type(ob1) == type(None) and type(ob2) == type(None):
        raise TypeError('Ob1 and Ob2 cannot both be None')
    elif type(ob1) == type(None) and type(ob2) == Observation:
        if ob2.face == 'FL':
            return ob2
        else:
            return ob2.changeface()
    elif type(ob1) == Observation and type(ob2) == type(None):
        if ob1.face == 'FL':
            return ob1
        else:
            return ob1.changeface()
    elif (ob1.from_id != ob2.from_id
            or ob1.to_id != ob2.to_id
            or ob1.inst_height != ob2.inst_height
            or ob1.target_height != ob2.target_height):
        raise ValueError('Inputs must be Observation Objects. '
                         'From, To, Instrument and Target Height must be the same to mean two observations')
    else:
        # Check Face Inputs, convert all to Face Left Sense
        if ob1.face == 'FL' and ob2.face == 'FR':
            ob2 = ob2.changeface()
        elif ob1.face == 'FR' and ob2.face == 'FL':
            ob1 = ob1.changeface()
        elif ob1.face == 'FR' and ob2.face == 'FR':
            ob1 = ob1.changeface()
            ob2 = ob2.changeface()
        elif ob1.face == 'FL' and ob2.face == 'FL':
            pass
        else:
            raise ValueError('Incompatible Face Values')
        # Calculate means, return new observation
        meaned_hz = (ob1.hz_obs + ob2.hz_obs) / 2
        meaned_va = (ob1.va_obs + ob2.va_obs) / 2
        meaned_sd = round(((ob1.sd_obs + ob2.sd_obs) / 2), 4)
        return Observation(ob1.from_id,
                           ob1.to_id,
                           ob1.inst_height,
                           ob1.target_height,
                           ob1.face,
                           meaned_hz,
                           meaned_va,
                           meaned_sd)


def reducesetup(obslist, strict=False):
    """
    Takes a list of Observations from one setup and
    means together FL, FR pairs of Observations.
    :param obslist: List of Observations (i.e. from one InstSetup)
    :param strict: If True, all single-face Obs are ignored. If False, all
    single-face Obs are included and converted to Face Left
    :return: a reduced list of Observations
    """
    # Group obs numbers by to_id
    uniqueto = []
    for ob in obslist:
        uniqueto.append(ob.to_id)
    uniqueto = list(set(uniqueto))
    # Sort Obs by to_id and face
    meanedobs = []
    for unique_id in uniqueto:
        fl_list = []
        fr_list = []
        for ob in obslist:
            if ob.to_id == unique_id and ob.face == 'FL':
                fl_list.append(ob)
            elif ob.to_id == unique_id and ob.face == 'FR':
                fr_list.append(ob)
            elif ob.to_id != unique_id and (ob.face == 'FL' or ob.face == 'FR'):
                pass
            else:
                raise ValueError('Invalid Face')
        obsdict = {unique_id: {'FL': fl_list, 'FR': fr_list}}
        # Group Obs into FL, FR pairs and mean (Remove all non-paired obs)
        if strict:
            for key in obsdict:
                pairedlist = list(zip(obsdict[key]['FL'], obsdict[key]['FR']))
                for pair in pairedlist:
                    meanob = meanfaces(pair[0], pair[1])
                    meanedobs.append(meanob)
        # Group Obs into FL, FR pairs and mean (Keep all non-paired obs)
        elif not strict:
            for key in obsdict:
                pairedlist = list(itertools.zip_longest(obsdict[key]['FL'], obsdict[key]['FR'], fillvalue=None))
                for pair in pairedlist:
                    meanob = meanfaces(pair[0], pair[1])
                    meanedobs.append(meanob)
    return meanedobs


def dnaout_dirset(obslist, same_stdev=True):
    # Test for Single Observation
    if len(obslist) < 2:
        return []
    fromlist = []
    pointing_err = 0.001  # 0.001m
    stdev = '1.0000'  # 1sec
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
                 + stdev.ljust(9))  # add standard deviation
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
                    + stdev.ljust(9))  # add standard deviation
            dnaobs.append(line)
    else:
        stdev = str(dd2sec(degrees(atan(pointing_err / obslist[0].sd_obs))))
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
                 + stdev[0:6].ljust(9))  # add standard deviation
        dnaobs.append(line1)
        # create other lines using range(1:)
        for num in range(1, len(obslist)):
            stdev = str(dd2sec(degrees(atan(pointing_err / obslist[num].sd_obs))))
            line = ('D '
                    + ''.ljust(40)
                    + obslist[num].to_id.ljust(20)
                    + str(obslist[num].hz_obs.degree).rjust(17)
                    + ' '
                    + str('%02d' % obslist[num].hz_obs.minute)
                    + ' '
                    + str('{:.3f}'.format(obslist[num].hz_obs.second).rjust(6, '0').ljust(8))
                    + stdev[0:6].ljust(9))  # add standard deviation
            dnaobs.append(line)
    return dnaobs


def dnaout_sd(obslist):
    dnaobs = []
    for observation in obslist:
        line = ('S '
                + observation.from_id.ljust(20)
                + observation.to_id.ljust(20)
                + ''.ljust(19)
                + ('{:.4f}'.format(observation.sd_obs)).rjust(9)
                + ''.ljust(21)
                + '0.0010'.ljust(8)         # add standard deviation
                + str(observation.inst_height).ljust(7)      # add intrument height
                + str(observation.target_height).ljust(7))
        dnaobs.append(line)
    return dnaobs


def dnaout_va(obslist, same_stdev=True):
    dnaobs = []
    pointing_err = 0.001
    stdev = '1.0000'
    if same_stdev:
        for observation in obslist:
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
                    + stdev[0:6].ljust(8)         # add standard deviation
                    + str(observation.inst_height).ljust(7)      # add intrument height
                    + str(observation.target_height).ljust(7))
            dnaobs.append(line)
    else:
        for observation in obslist:
            # Calc Standard Deviation (in seconds) based on pointing error (in metres)
            stdev = str(dd2sec(degrees(atan(pointing_err / observation.sd_obs))))
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
                    + stdev[0:6].ljust(8)  # add standard deviation
                    + str(observation.inst_height).ljust(7)  # add intrument height
                    + str(observation.target_height).ljust(7))
            dnaobs.append(line)
    return dnaobs


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


"""
to_stn = ['GA03', 'GA02', 'SYM2', 'ATWR', 'MAHON', 'MAHON', 'ATWR', 'SYM2', 'GA02', 'GA03']
flfr = [359.59597, 14.48289, 25.35515, 215.57043, 300.59046, 120.59060, 35.57045, 205.35530, 194.48282, 179.59581]
flfr_vert = [90.30228, 90.09547, 89.50396, 88.03600, 87.58014, 272.01587, 271.55593, 270.09227, 269.50091, 269.29414]
"""


def hz_round(brg_list):
    """
    Input: an even palindromic list of horizontal angle observations in two faces
    (e.g. [fl1, fl2, fl3, fr3, fr2, fr1])
    Output: an averaged face-left sense list of horizontal angles.
    """
    brg_avg = []
    obs = int((len(brg_list))/2)
    for i in range(0, obs):
        hz_avg = (hp2dec(brg_list[i]) + (hp2dec(brg_list[-(i + 1)]) - 180)) / 2
        brg_avg.append(round(dec2hp(hz_avg), 7))
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
        fl_ang = hp2dec(va_list[i]) - 90
        fr_ang = 270 - hp2dec(va_list[-(i + 1)])
        ang_avg = (fl_ang + fr_ang)/2 + 90
        va_avg.append(round(dec2hp(ang_avg), 7))
    return va_avg


"""
# Test for round of obs
if to_stn == to_stn[::-1] and len(to_stn) % 2 == 0:
    brg_avg = hz_round(flfr)
    va_avg = va_round(flfr_vert)
    to_stn_avg = to_stn[0:int(len(to_stn)/2)]
else:
    brg_avg = ['nope']
"""

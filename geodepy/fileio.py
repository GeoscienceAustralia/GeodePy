#!/usr/bin/env python3

"""
Geoscience Australia - GeodePy Package
File input/output module

In Development
"""

from numpy import zeros
import sys

def rd_snx_est(file):
    """This function reads in the SOLUTION/ESTIMATE block of a SINEX file. It
    returns estimate, a list of tuples:

    estimate = [(code, soln, refEpoch, staX, staY, staZ, staX_sd, staY_sd,
                 staZ_sd [, velX, velY, velZ, velX_sd, velY_sd, velZ_sd])...]

    where:
        * code is the stations's 4-character ID
        * soln is the segment of the stations's time series
        * refEpoch is the epoch of the solution in the form YY:DOY:SSSSS (YY
            is the two digit year, DOY is day of year, and SSSSS is the time of
            day in seconds
        * sta[XYZ] is the station coordinates in the Cartesian reference frame
        * sta[XYZ]_sd
        * vel[XYZ] is the station velocity in the Cartesian reference frame
        * vel[XYZ]_sd

    Velocities are not included in all SINEX files and so are only returned if
    present.

    :param file: the input SINEX file
    :return: estimate
    """

    # Create data structures and set variables
    lines = []
    estimate = []
    check_types = True
    velocities = False
    go = False
    code = ''
    soln = ''
    epoch = ''
    stax = ''
    stay = ''
    staz = ''
    stax_sd = ''
    stay_sd = ''
    staz_sd = ''
    velx = ''
    vely = ''
    velz = ''
    velx_sd = ''
    vely_sd = ''
    velz_sd = ''

    # Read the SOLUTION/ESTIMATE block into a list and determine if there is
    # any velocity information
    with open(file) as f:
        for line in f:
            if line[:18] == '-SOLUTION/ESTIMATE':
                break
            if go and line[:11] == '*INDEX TYPE':
                pass
            elif go:
                if check_types:
                    if line[7:10] == 'VEL':
                        velocities = True
                        check_types = False
                lines.append(line)
            if line[:18] == '+SOLUTION/ESTIMATE':
                go = True

    if velocities:
        for line in lines:
            typ = line[7:11]
            if typ == 'STAX':
                code = line[14:18]
                soln = line[23:26].lstrip()
                epoch = line[27:39]
                stax = float(line[47:68])
                stax_sd = float(line[69:80])
            elif typ == 'STAY':
                stay = float(line[47:68])
                stay_sd = float(line[69:80])
            elif typ == 'STAZ':
                staz = float(line[47:68])
                staz_sd = float(line[69:80])
            elif typ == 'VELX':
                velx = float(line[47:68])
                velx_sd = float(line[69:80])
            elif typ == 'VELY':
                vely = float(line[47:68])
                vely_sd = float(line[69:80])
            elif typ == 'VELZ':
                velz = float(line[47:68])
                velz_sd = float(line[69:80])
                info = (code, soln, epoch, stax, stay, staz, stax_sd,
                        stay_sd, staz_sd, velx, vely, velz, velx_sd,
                        vely_sd, velz_sd)
                estimate.append(info)
    else:
        for line in lines:
            typ = line[7:11]
            if typ == 'STAX':
                code = line[14:18]
                soln = line[23:26].lstrip()
                epoch = line[27:39]
                stax = float(line[47:68])
                stax_sd = float(line[69:80])
            elif typ == 'STAY':
                stay = float(line[47:68])
                stay_sd = float(line[69:80])
            elif typ == 'STAZ':
                staz = float(line[47:68])
                staz_sd = float(line[69:80])
                info = (code, soln, epoch, stax, stay, staz, stax_sd,
                        stay_sd, staz_sd)
                estimate.append(info)

    return estimate

def rd_snx_mat(file):

    # Read in the codes (station names) and solutions, and check for velocities
    data = rd_snx_est(file)
    code_solns = []
    for stat in data:
        code_solns.append((stat[0], stat[1]))
    if len(data[0]) == 15:
        velocities = True

    # Read the SOLUTION/MATRIX_ESTIMATE block into a list and determine if the
    # matrix is upper or lower triangular
    lines = []
    go = False
    with open(file) as f:
        for line in f:
            if line[:25] == '-SOLUTION/MATRIX_ESTIMATE':
                break
            if go and line[:12] == '*PARA1 PARA2':
                pass
            elif go:
                lines.append(line)
            if line[:25] == '+SOLUTION/MATRIX_ESTIMATE':
                low_up = line[26]
                go = True

    # Create the VCV matrix
    if velocities:
        n = 6 * int(len(code_solns))
    else:
        n = 3 * int(len(code_solns))
    vcv = zeros((n, n))
    for line in lines:
        print(line)

    return code_solns

def rd_dyna_adj_coords(file):

    stat = ''
    const = ''
    easting = ''
    northing = ''
    zone = ''
    lat = ''
    lon = ''
    h_ortho = ''
    h_ellipse = ''
    x = ''
    y = ''
    z = ''
    sd_e = ''
    sd_n = ''
    sd_u = ''

    lines = []
    go = False
    with open(file) as f:
        for line in f:
            line = line.rstrip()
            if line[:20] == 'Adjusted Coordinates':
                go = True
            if go and line != '':
                lines.append(line)
    for line in lines[5:]:
        print(line)

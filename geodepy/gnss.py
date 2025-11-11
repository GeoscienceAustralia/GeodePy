#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
GNSS Module

In Development

Specific utility functions for SINEX:

    - remove_stns_sinex()
    - remove_velocity_sinex()
    - remove_matrixzeros_sinex()

General functions for reading from SINEX:

    - list_sinex_blocks()
    - read_sinex_comments()
    - read_sinex_header_line()
    - read_sinex_custom()
    - read_sinex_estimate()
    - read_sinex_matrix()
    - read_sinex_sites()
    - read_disconts()
    - read_solution_epochs()
    - read_sinex_header_block()
    - read_sinex_file_reference_block()
    - read_sinex_input_acknowledgments_block()
    - read_sinex_solution_statistics_block()
    - read_sinex_site_receiver_block()
    - read_sinex_site_antenna_block()
    - read_sinex_site_gps_phase_center_block()
    - read_sinex_site_eccentricity_block()
    - read_sinex_site_id_block()
    - read_sinex_solution_epochs_block()
    - read_sinex_solution_estimate_block()
    - read_sinex_solution_apriori_block()
    - read_sinex_solution_matrix_estimate_block()
    - read_sinex_solution_matrix_apriori_block()
    - matrix2dataframe_solution_matrix_estimate()
    - sinex2dataframe_solution_estimate()
    - sinex2dataframe_solution_apriori()
    - sinex2dataframe_solution_matrix_estimate()

General functions for writing to SINEX:

    - print_sinex_comments()
    - set_creation_time()
    - dataframe2sinex_solution_estimate()
    - dataframe2sinex_solution_apriori()
    - dataframe2matrix_snx_vcv()
    - dataframe2matrix_solution_matrix_estimate()
    - dataframe2sinex_solution_matrix_estimate()
    - writeSINEX()

"""

from datetime import datetime
import numpy as np
import pandas as pd
from numpy import zeros, delete
from geodepy.angles import DMSAngle
import re


def list_sinex_blocks(file):
    """
    This script lists the blocks in a SINEX file

    :param str file: the input SINEX file
    """
    blocks = []
    with open(file) as f:
        for line in f.readlines():
            if line.startswith("+"):
                col = line.split(" ")
                block = col[0].replace("+", "")
                blocks.append(block.strip())
    for block in blocks:
        print(block)


def print_sinex_comments(file):
    """This script prints comments in a SINEX file.

    :param str file: the input SINEX file
    """
    go = False
    with open(file) as f:
        for line in f.readlines():
            if line.startswith("+FILE/COMMENT"):
                go = True
            if go:
                print(line.strip())
            if line.startswith("-FILE/COMMENT"):
                go = False


def read_sinex_comments(file):
    """This function reads comments in a SINEX file.

    :param str file: the input SINEX file
    :return: comments block
    :rtype: list of strings
    """
    comments = []
    go = False
    with open(file) as f:
        for line in f.readlines():
            if line.startswith("+FILE/COMMENT"):
                go = True
            if go:
                comments.append(line.strip())
            if line.startswith("-FILE/COMMENT"):
                go = False

        comments.insert(
            -1,
            f"* File created by Geodepy.gnss.py at {datetime.now().strftime('%d-%m-%Y, %H:%M')}",
        )

    return comments


def set_creation_time():
    """This function sets the creation time, in format YY:DDD:SSSSS, for use
    in the SINEX header line

    :return: creation_time
    :rtype: str
    """
    now = datetime.now()
    time_tup = now.timetuple()
    year = str(time_tup.tm_year)[2:]
    doy = time_tup.tm_yday
    doy = "{:03d}".format(doy)
    seconds = (
        now - now.replace(hour=0, minute=0, second=0, microsecond=0)
    ).total_seconds()
    # seconds = '{:.0f}'.format(seconds)
    seconds = "{:05.0f}".format(seconds)
    creation_time = year + ":" + doy + ":" + seconds

    return creation_time


def read_sinex_header_line(file):
    """This function reads the header line of a SINEX file into a string

    :param str file: the input SINEX file
    :return: header_line
    :rtype: str
    """
    with open(file) as f:
        header_line = f.readline()

    return header_line


def read_sinex_custom(fp, start_line, end_line):
    """Read custom line range from SINEX. Useful for
    SINEX with irregular formatting or block that has
    not been accounted for.

    :param str fp: Path to SINEX file.
    :param int start_line:   Beginning line number.
    :param int end_line:     Finishing line number.
    :return list custom:  List of string(s).
    """

    # Setup
    i = start_line
    j = end_line
    f = fp
    custom = []

    with open(fp) as f:
        for line in f.readlines()[i - 1 : j]:
            custom.append(line.strip())

    # Return
    return custom


def read_sinex_estimate(file):
    """This function reads in the SOLUTION/ESTIMATE block of a SINEX file. It
    returns estimate, a list of tuples:

    estimate = [(code, soln, refEpoch, staX, staY, staZ, staX_sd, staY_sd, staZ_sd[, velX, velY, velZ, velX_sd, velY_sd, velZ_sd])...]

    where:
        * code is the stations's 4-character ID
        * soln is the segment of the stations's time series
        * refEpoch is the epoch of the solution in the form YY:DOY:SSSSS (YY
          is the two digit year, DOY is day of year, and SSSSS is the time of
          day in seconds
        * sta[XYZ] is the station coordinates in the Cartesian reference frame
        * sta[XYZ]_sd is the standard deviation of the station coordinates in
          the Cartesian reference frame
        * vel[XYZ] is the station velocity in the Cartesian reference frame
        * vel[XYZ]_sd is the standard deviation of the station velocity in the
          Cartesian reference frame

    Velocities are not included in all SINEX files and so are only returned if
    present.

    :param file: the input SINEX file
    :return: estimate
    """

    # Create data structures and set variables
    lines = []
    estimate = []
    velocities = False
    go = False
    code = ""
    soln = ""
    epoch = ""
    stax = ""
    stay = ""
    staz = ""
    stax_sd = ""
    stay_sd = ""
    staz_sd = ""
    velx = ""
    vely = ""
    velz = ""
    velx_sd = ""
    vely_sd = ""
    velz_sd = ""

    # Read the SOLUTION/ESTIMATE block into a list and determine if there is
    # any velocity information
    with open(file) as f:
        for line in f:
            if line[:18] == "-SOLUTION/ESTIMATE":
                break
            if go and line[:11] == "*INDEX TYPE":
                pass
            elif go:
                if line[7:10] == "VEL":
                    velocities = True
                lines.append(line)
            if line[:18] == "+SOLUTION/ESTIMATE":
                go = True

    for line in lines:
        typ = line[7:11]
        if typ == "STAX":
            code = line[14:18]
            soln = line[23:26].lstrip()
            epoch = line[27:39]
            stax = float(line[47:68])
            stax_sd = float(line[69:80])
        elif typ == "STAY":
            stay = float(line[47:68])
            stay_sd = float(line[69:80])
        elif typ == "STAZ":
            staz = float(line[47:68])
            staz_sd = float(line[69:80])
            if not velocities:
                info = (code, soln, epoch, stax, stay, staz, stax_sd, stay_sd, staz_sd)
                estimate.append(info)
        elif typ == "VELX":
            velx = float(line[47:68])
            velx_sd = float(line[69:80])
        elif typ == "VELY":
            vely = float(line[47:68])
            vely_sd = float(line[69:80])
        elif typ == "VELZ":
            velz = float(line[47:68])
            velz_sd = float(line[69:80])
            info = (
                code,
                soln,
                epoch,
                stax,
                stay,
                staz,
                stax_sd,
                stay_sd,
                staz_sd,
                velx,
                vely,
                velz,
                velx_sd,
                vely_sd,
                velz_sd,
            )
            estimate.append(info)

    return estimate


def read_sinex_matrix(file):
    """This function reads in the SOLUTION/MATRIX_ESTIMATE block of a SINEX
    file. It returns matrix, a list of tuples:

    matrix = [(code, soln, var_x, covar_xy, covar_xz, var_y, covar_yz, var_z[, var_v_x, covar_v_xy, covar_v_xz, var_v_y, covar_v_yz, var_v_z])...]

    where:
        * code is the stations's 4-character ID
        * soln is the segment of the stations's time series
        * var_x is the variance in the X coordinate
        * covar_xy is the covariance between the X and the Y coordinates
        * covar_xz is the covariance between the X and the Z coordinates
        * var_y is the variance in the Y coordinate
        * covar_yz is the covariance between the Y and the Z coordinates
        * var_z is the variance in the Z coordinate
        * var_v_x is the variance in the X velocity
        * covar_v_xy is the covariance between the X and the Y velocities
        * covar_v_xz is the covariance between the X and the Z velocities
        * var_v_y is the variance in the Y velocity
        * covar_v_yz is the covariance between the Y and the Z velocities
        * var_v_z is the variance in the Z velocity

    Velocities are not included in all SINEX files and so their VCV information
    is only returned if they are present.

    :param file: the input SINEX file
    :return: matrix
    """
    '''
    ToDo:
    1. The above order is only valid if the matrix is upper triangle. If it is
       lower triangle, then the covar_xy is actually the var_y. Will need to fix
       this when time permits.
    '''
    # Read in the codes (station names) and solutions, and check for velocities
    data = read_sinex_estimate(file)
    code = []
    soln = []
    velocities = False
    for station in data:
        code.append(station[0])
        soln.append(station[1])
    if len(data[0]) == 15:
        velocities = True

    # Read the SOLUTION/MATRIX_ESTIMATE block into a list and determine if the
    # matrix is upper or lower triangular
    lines = []
    lower_triangular = False
    go = False
    with open(file) as f:
        for line in f:
            if line[:25] == "-SOLUTION/MATRIX_ESTIMATE":
                break
            if go and line[:12] == "*PARA1 PARA2":
                pass
            elif go:
                lines.append(line)
            if line[:25] == "+SOLUTION/MATRIX_ESTIMATE":
                if line[26] == "L":
                    lower_triangular = True
                go = True

    # Create an array containing the matrix elements
    if velocities:
        n = 6 * int(len(code))
    else:
        n = 3 * int(len(code))
    element = zeros((n, n))
    matrix = []
    for line in lines:
        col = line.rstrip().split()
        for i in range(2, len(col)):
            element[int(col[0]) - 1][int(col[1]) + i - 3] = float(col[i])
    if velocities:
        if lower_triangular:
            for i in range(len(code)):
                info = (
                    code[i],
                    soln[i],
                    element[6 * i][6 * i],
                    element[6 * i + 1][6 * i],
                    element[6 * i + 1][6 * i + 1],
                    element[6 * i + 2][6 * i],
                    element[6 * i + 2][6 * i + 1],
                    element[6 * i + 2][6 * i + 2],
                    element[6 * i + 3][6 * i + 3],
                    element[6 * i + 4][6 * i + 3],
                    element[6 * i + 4][6 * i + 4],
                    element[6 * i + 5][6 * i + 3],
                    element[6 * i + 5][6 * i + 4],
                    element[6 * i + 5][6 * i + 5],
                )
                matrix.append(info)
        else:
            for i in range(len(code)):
                info = (
                    code[i],
                    soln[i],
                    element[6 * i][6 * i],
                    element[6 * i][6 * i + 1],
                    element[6 * i][6 * i + 2],
                    element[6 * i + 1][6 * i + 1],
                    element[6 * i + 1][6 * i + 2],
                    element[6 * i + 2][6 * i + 2],
                    element[6 * i + 3][6 * i + 3],
                    element[6 * i + 3][6 * i + 4],
                    element[6 * i + 3][6 * i + 5],
                    element[6 * i + 4][6 * i + 4],
                    element[6 * i + 4][6 * i + 5],
                    element[6 * i + 5][6 * i + 5],
                )
                matrix.append(info)
    else:
        if lower_triangular:
            for i in range(len(code)):
                info = (
                    code[i],
                    soln[i],
                    element[3 * i][3 * i],
                    element[3 * i + 1][3 * i],
                    element[3 * i + 1][3 * i + 1],
                    element[3 * i + 2][3 * i],
                    element[3 * i + 2][3 * i + 1],
                    element[3 * i + 2][3 * i + 2],
                )
                matrix.append(info)
        else:
            for i in range(len(code)):
                info = (
                    code[i],
                    soln[i],
                    element[3 * i][3 * i],
                    element[3 * i][3 * i + 1],
                    element[3 * i][3 * i + 2],
                    element[3 * i + 1][3 * i + 1],
                    element[3 * i + 1][3 * i + 2],
                    element[3 * i + 2][3 * i + 2],
                )
                matrix.append(info)

    return matrix


def read_sinex_sites(file):
    """This function reads in the SITE/ID block of a SINEX file. It returns
    sites, a list of tuples:

    sites = [(site, point, domes, obs, station_description, lon, lat, h)]

    where:
        * site is the site code
        * point is the site's point code
        * domes is the site's dome number
        * obs is the observation technique
        * station_description is a free format desciption of the site
        * lon is the approximate longitude of the site as a DMSAngle object
        * lat is the approximate latitude of the site as a DMSAngle object
        * h is the approximate height of the site

    :param file: the input SINEX file
    :return: sites
    """

    # Read the SITE/ID block into a list
    lines = []
    go = False
    with open(file) as f:
        for line in f:
            if line[:8] == "-SITE/ID":
                break
            if go and line[:8] == "*CODE PT":
                pass
            elif go:
                lines.append(line)
            if line[:8] == "+SITE/ID":
                go = True
    sites = []
    for line in lines:
        site = line[1:5]
        point = line[6:8].lstrip()
        domes = line[9:18]
        obs = line[19:20]
        station_description = line[21:43].lstrip()
        lon = DMSAngle(line[44:55].lstrip())
        lat = DMSAngle(line[56:67].lstrip())
        h = float(line[67:73])
        info = (site, point, domes, obs, station_description, lon, lat, h)
        sites.append(info)

    return sites


def read_disconts(file):
    """This function reads in the SOLUTION/DISCONTINUITY block of a
    SINEX file. It returns disconts , a list of tuples:

    sites = [(site, code1, point, code2, start, end, type)]

    where:
        * site is the site code
        * code1 is unknown
        * point is the site's point code
        * code2 is unknown
        * start is the start time for the point code in YY:DOY:SECOD
        * end is the end time for the point code in YY:DOY:SECOD
        * type is the type of discontinuity; P for position or V for
          velocity

    I could not find the format description for this block.

    :param file: the input discontinuities file
    :return: disconts
    """

    # Read the SOLUTION/DISCONTINUITY block into a list
    lines = []
    go = False
    with open(file) as f:
        for line in f:
            if line[:23] == "-SOLUTION/DISCONTINUITY":
                break
            elif go:
                lines.append(line)
            if line[:23] == "+SOLUTION/DISCONTINUITY":
                go = True
    disconts = []
    for line in lines:
        site = line[1:5]
        code1 = line[5:8].lstrip()
        point = line[8:13].lstrip()
        code2 = line[14:15]
        start = line[16:28]
        end = line[29:41]
        p_or_v = line[42:43]
        info = (site, code1, point, code2, start, end, p_or_v)
        disconts.append(info)

    return disconts


def read_solution_epochs(file):
    """This function reads in the SOLUTION/EPOCHS block of a SINEX file.
    It returns epochs, a list of tuples:

    epochs = [(site, point, sol, obs, start, end, mean)]

    where:
        * site is the site code
        * point is the site's point code
        * sol is the solution number at a site/point
        * obs is the observation technique
        * start is the start time for the solution in YY:DOY:SECOD
        * end is the end time for the solution in YY:DOY:SECOD
        * mean is the mean time for the solution in YY:DOY:SECOD

    :param file: the input SINEX file
    :return: epochs
    """

    # Read the SOLUTION/EPOCHS block into a list
    lines = []
    go = False
    with open(file) as f:
        for line in f:
            if line[:16] == "-SOLUTION/EPOCHS":
                break
            if go and line[:8] == "*Code PT":
                pass
            elif go:
                lines.append(line)
            if line[:16] == "+SOLUTION/EPOCHS":
                go = True
    epochs = []
    # Parse each line, create a tuple and add it to the list
    for line in lines:
        site = line[1:5]
        point = line[6:8].lstrip()
        sol = line[9:13].lstrip()
        obs = line[14:15]
        start = line[16:28]
        end = line[29:41]
        mean = line[42:55].rstrip()
        info = (site, point, sol, obs, start, end, mean)
        epochs.append(info)

    return epochs


def read_sinex_header_block(sinex):
    """
    This function reads in the header block information
    of a SINEX file (All lines before the SITE/ID block).

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    with open(sinex, "r") as f:
        next(f)
        line = f.readline()
        while line:
            block.append(line.rstrip())
            line = f.readline()
            if line.startswith("+SITE/ID"):
                break

    return block


def read_sinex_file_reference_block(sinex):
    """
    This function reads in the +FILE/REFERENCE block
    of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+FILE/REFERENCE"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-FILE/REFERENCE"):
                break
            line = f.readline()

    return block


def read_sinex_input_acknowledgments_block(sinex):
    """
    This function reads in the +INPUT/ACKNOWLEDGMENTS
    block of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+INPUT/ACKNOWLEDGMENTS"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-INPUT/ACKNOWLEDGMENTS"):
                break
            line = f.readline()

    return block


def read_sinex_solution_statistics_block(sinex):
    """
    This function reads in the +SOLUTION/STATISTICS
    block of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SOLUTION/STATISTICS"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SOLUTION/STATISTICS"):
                break
            line = f.readline()

    return block


def read_sinex_site_receiver_block(sinex):
    """
    This function reads in the +SITE/RECEIVER
    block of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SITE/RECEIVER"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SITE/RECEIVER"):
                break
            line = f.readline()

    return block


def read_sinex_site_antenna_block(sinex):
    """
    This function reads in the +SITE/ANTENNA
    block of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SITE/ANTENNA"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SITE/ANTENNA"):
                break
            line = f.readline()

    return block


def read_sinex_site_gps_phase_center_block(sinex):
    """
    This function reads in the +SITE/GPS_PHASE_CENTER
    block of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SITE/GPS_PHASE_CENTER"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SITE/GPS_PHASE_CENTER"):
                break
            line = f.readline()

    return block


def read_sinex_site_eccentricity_block(sinex):
    """
    This function reads in the +SITE/ECCENTRICITY
    block of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SITE/ECCENTRICITY"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SITE/ECCENTRICITY"):
                break
            line = f.readline()

    return block


def read_sinex_site_id_block(sinex):
    """
    This function reads in the SITE/ID block of a SINEX file
    into list of strings.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SITE/ID"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SITE/ID"):
                break
            line = f.readline()

    return block


def read_sinex_solution_epochs_block(sinex):
    """
    This function reads in the SOLUTION/EPOCHS block of a SINEX file
    into list of strings.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SOLUTION/EPOCHS"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SOLUTION/EPOCHS"):
                break
            line = f.readline()
    return block


def read_sinex_solution_estimate_block(sinex):
    """
    This function reads in the SOLUTION/ESTIMATE block of a SINEX
    file into list of strings.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SOLUTION/ESTIMATE"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SOLUTION/ESTIMATE"):
                break
            line = f.readline()
    return block


def read_sinex_solution_apriori_block(sinex):
    """
    This function reads in the +SOLUTION/APRIORI
    block of a SINEX file.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SOLUTION/APRIORI"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SOLUTION/APRIORI"):
                break
            line = f.readline()

    return block


def read_sinex_solution_matrix_estimate_block(sinex):
    """
    This function reads in the SOLUTION/MATRIX_ESTIMATE block of a SINEX
    file into list of strings.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SOLUTION/MATRIX_ESTIMATE"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SOLUTION/MATRIX_ESTIMATE"):
                break
            line = f.readline()

    return block


def read_sinex_solution_matrix_apriori_block(sinex):
    """
    This function reads in the SOLUTION/MATRIX_APRIORI block
    of a SINEX file into list of strings.

    :param str sinex: input SINEX file
    :return: block
    """

    block = []
    go = False
    with open(sinex, "r") as f:
        line = f.readline()
        while line:
            if line.startswith("+SOLUTION/MATRIX_APRIORI"):
                go = True
            if go:
                block.append(line.rstrip())
            if line.startswith("-SOLUTION/MATRIX_APRIORI"):
                break
            line = f.readline()

    return block


def sinex2dataframe_solution_estimate(fp):
    """
    This function reads in a SINEX file and returns
    a dataframe of SOLUTION/ESTIMATE block only.

    :param str sinex: path of input SINEX file
    :return: df
    """

    # Get lines
    lines = read_sinex_solution_estimate_block(fp)

    # Remove non-data lines
    for l in lines[:]:
        if l.startswith("*"):
            lines.remove(l)
        if l.startswith("+"):
            lines.remove(l)
        if l.startswith("-"):
            lines.remove(l)

    # Split by column
    lines = [i.split() for i in lines]

    # Isolate into vectors
    row = np.int_(list(zip(*lines))[0])
    par = list(zip(*lines))[1]
    code = list(zip(*lines))[2]
    pt = list(zip(*lines))[3]
    soln = list(zip(*lines))[4]
    refEpoch = list(zip(*lines))[5]
    unit = list(zip(*lines))[6]
    s = np.int_(list(zip(*lines))[7])
    est = np.float64(list(zip(*lines))[8])
    sigma = np.float64(list(zip(*lines))[9])

    # Organise into DataFrame
    dict_temp = {
        "row": row,
        "par": par,
        "code": code,
        "pt": pt,
        "soln": soln,
        "refEpoch": refEpoch,
        "unit": unit,
        "s": s,
        "est": est,
        "sigma": sigma,
    }
    df = pd.DataFrame(dict_temp)

    # Return
    return df


def sinex2dataframe_solution_apriori(fp):
    """
    This function reads in a SINEX file and returns
    a dataframe of SOLUTION/APRIORI block only.

    :param str sinex: path of input SINEX file
    :return: df
    """

    # Get lines
    lines = read_sinex_solution_apriori_block(fp)

    # Remove non-data lines
    for l in lines[:]:
        if l.startswith("*"):
            lines.remove(l)
        if l.startswith("+"):
            lines.remove(l)
        if l.startswith("-"):
            lines.remove(l)

    # Split by column
    lines = [i.split() for i in lines]

    # Isolate into vectors
    row = np.int_(list(zip(*lines))[0])
    par = list(zip(*lines))[1]
    code = list(zip(*lines))[2]
    pt = list(zip(*lines))[3]
    soln = list(zip(*lines))[4]
    refEpoch = list(zip(*lines))[5]
    unit = list(zip(*lines))[6]
    s = np.int_(list(zip(*lines))[7])
    est = np.float64(list(zip(*lines))[8])
    sigma = np.float64(list(zip(*lines))[9])

    # Organise into DataFrame
    dict_temp = {
        "row": row,
        "par": par,
        "code": code,
        "pt": pt,
        "soln": soln,
        "refEpoch": refEpoch,
        "unit": unit,
        "s": s,
        "est": est,
        "sigma": sigma,
    }
    df = pd.DataFrame(dict_temp)

    # Return
    return df


def sinex2dataframe_solution_matrix_estimate(fp):
    """
    This function reads in a SINEX file and returns
    a dataframe of SOLUTION/MATRIX_ESTIMATE block only.

    :param str sinex: path of input SINEX file
    :return: df
    """

    # Get lines
    lines = read_sinex_solution_matrix_estimate_block(fp)

    # Remove non-data lines
    for l in lines[:]:
        if l.startswith("*"):
            lines.remove(l)
        if l.startswith("+"):
            lines.remove(l)
        if l.startswith("-"):
            lines.remove(l)

    # Split by column
    lines = [i.split() for i in lines]

    # Pad out lines with NaN
    for i in range(len(lines)):
        if len(lines[i]) == 5:
            continue
        if len(lines[i]) == 4:
            lines[i].append(np.nan)
        if len(lines[i]) == 3:
            lines[i].append(np.nan)
            lines[i].append(np.nan)

    # Isolate into vectors
    row = np.int_(list(zip(*lines))[0])
    col = np.int_(list(zip(*lines))[1])
    q1 = np.float64(list(zip(*lines))[2])
    q2 = np.float64(list(zip(*lines))[3])
    q3 = np.float64(list(zip(*lines))[4])

    # Organise into DataFrame
    dict_temp = {
        "row": row,
        "col": col,
        "q1": q1,
        "q2": q2,
        "q3": q3,
    }
    df = pd.DataFrame(dict_temp)

    # Return
    return df


def dataframe2sinex_solution_estimate(df):
    """
    This function reads in a dataframe of the
    SOLUTIONS/ESTIMATE block from a SINEX, then
    converts each row to a string in a list
    ready for writing to SINEX.

    :param dataframe df: dataframe of SOLUTION/ESTIMATE block
    :return: list of strings
    """

    lines = lines = ["+SOLUTION/ESTIMATE"]
    lines.append(
        "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___"
    )

    for i in range(len(df.code)):
        l = " %5i %s   %s  %s    %s %s %-3s  %s %21.14e %.5e" % (
            df.row.iloc[i],
            df.par.iloc[i],
            df.code.iloc[i],
            df.pt.iloc[i],
            df.soln.iloc[i],
            df.refEpoch.iloc[i],
            df.unit.iloc[i],
            df.s.iloc[i],
            df.est.iloc[i],
            df.sigma.iloc[i],
        )
        lines.append(l)

    lines.append("-SOLUTION/ESTIMATE")

    # Return
    return lines


def dataframe2sinex_solution_apriori(df):
    """
    This function reads in a dataframe of the
    SOLUTION/APRIORI block from a SINEX, then
    converts each row to a string in a list
    ready for writing to SINEX.

    :param dataframe df: dataframe of SOLUTION/APRIORI block
    :return: list of strings
    """

    lines = lines = ["+SOLUTION/APRIORI"]
    lines.append(
        "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___"
    )

    for i in range(len(df.code)):
        l = " %5i %s   %s  %s    %s %s %-3s  %s %21.14e %.5e" % (
            df.row.iloc[i],
            df.par.iloc[i],
            df.code.iloc[i],
            df.pt.iloc[i],
            df.soln.iloc[i],
            df.refEpoch.iloc[i],
            df.unit.iloc[i],
            df.s.iloc[i],
            df.est.iloc[i],
            df.sigma.iloc[i],
        )
        lines.append(l)

    lines.append("-SOLUTION/APRIORI")

    # Return
    return lines


def dataframe2matrix_snx_vcv(df, numPar=3):
    """
    This function converts a dataframe created from a
    SINEX VCV based from read_sinex_matrix(), and converts
    it to a full VCV without covariances between sites.

    .. code-block:: text

        i.e.  xx xy xz          0  0  0
              xy yy yz .  .  .  0  0  0
              xz yz zz          0  0  0
                   .     .           .
                   .        .        .
                   .           .     .
              0  0  0           xx xy xz
              0  0  0  .  .  .  xy yy yz
              0  0  0           xz yz zz

    :param DataFrame: dataframe formed from read_sinex_matrix():
    :return: Numpy matrix of VCV with no covariances between sites.
    """

    # Matrix dimensions
    # - To Do: Handle check for velocities (?)
    par = numPar  # (is there a way to automatically determine if velocitie exist?)
    n = len(df.code) * par
    Q = np.zeros((n, n))

    # Matrix elements and rows of dataframe
    i = 0
    j = 0
    r = 0
    while r < len(df.code):

        # variances
        Q[i + 0, j + 0] = df.xx[r]
        Q[i + 1, j + 1] = df.yy[r]
        Q[i + 2, j + 2] = df.zz[r]

        # covariances
        Q[i + 1, j + 0] = df.xy[r]
        Q[i + 0, j + 1] = df.xy[r]
        Q[i + 2, j + 0] = df.xz[r]
        Q[i + 0, j + 2] = df.xz[r]
        Q[i + 2, j + 1] = df.yz[r]
        Q[i + 1, j + 2] = df.yz[r]

        # Counter
        r = r + 1
        i = i + 3
        j = j + 3

    # Return
    return Q


def dataframe2matrix_solution_matrix_estimate(df, tri="L"):
    """
    This function reads in a dataframe of the SINEX
    SOLUTION/MATRIX_ESTIMATE block formed from
    sinex2dataframe_solution_matrix_estimate(), and then
    forms the full VCV matrix from that dataframe.

    :param DataFrame df: dataframe from sinex2dataframe_solution_matrix_estimate().
    :param String tri: String to indicate "upper" or "lower" triagle matrix.
    :return: Numpy matrix of full VCV.
    """

    # Triangularity of matrix
    triangle = tri

    # Matrix size
    n = df.row.iloc[-1]
    Q = np.zeros((n, n))

    if triangle == "L":

        for i in range(len(df.row)):

            # Get matrix indices
            row = int(df.row[i]) - 1
            col = int(df.col[i]) - 1

            # Fill PARA2+0
            q1 = df.q1[i]
            Q[row, col] = q1
            Q[col, row] = q1

            # Fill PARA2+1
            q2 = df.q2[i]
            Q[row, col + 1] = q2
            Q[col + 1, row] = q2

            # Fill PARA2+2
            q3 = df.q3[i]
            Q[row, col + 2] = q3
            Q[col + 2, row] = q3

    if triangle == "U":

        for i in range(len(df.row)):

            # Get matrix indices
            row = int(df.row[i]) - 1
            col = int(df.col[i]) - 1

            if df.col[i] < n - 1:

                # Fill PARA2+0
                q1 = df.q1[i]
                Q[row, col] = q1
                Q[col, row] = q1

                # Fill PARA2+1
                q2 = df.q2[i]
                Q[row, col + 1] = q2
                Q[col + 1, row] = q2

                # Fill PARA2+2
                q3 = df.q3[i]
                Q[row, col + 2] = q3
                Q[col + 2, row] = q3

            if df.col[i] == n - 1:

                # Fill PARA2+0
                q1 = df.q1[i]
                Q[row, col] = q1
                Q[col, row] = q1

                # Fill PARA2+1
                q2 = df.q2[i]
                Q[row, col + 1] = q2
                Q[col + 1, row] = q2

            if df.col[i] == n:

                # Fill PARA2+0
                q1 = df.q1[i]
                Q[row, col] = q1
                Q[col, row] = q1

    # Return
    return Q


def matrix2dataframe_solution_matrix_estimate(m, tri="L"):
    """
    This function reads a VCV in matrix format and writes it
    to dataframe for SINEX SOLUTION/MATRIX_ESTIMATE format. It
    is the format produced from sinex2dataframe_solution_matrix_estimate().

    :param numpy.array() m: A numpy array of the VCV matrix.
    :param str tri: String to indicate "upper" or "lower" triagle matrix.
    :return: A dataframe of SOLUTION/MATRIX_ESTIMATE.
    :rtype: DataFrame
    """

    # Matrix dimensions
    Q = m
    triangle = tri
    n = np.shape(Q)[0]

    if triangle == "L":

        # Make upper triangle NaNs
        Q[np.triu_indices(n, 1)] = np.nan

        row = []
        col = []
        q1 = []
        q2 = []
        q3 = []
        i = 0
        j = 0
        while i < n:
            while j <= i:
                row.append(i + 1)
                col.append(j + 1)
                q1.append(Q[i, j])
                q2.append(Q[i, j + 1])
                q3.append(Q[i, j + 2])
                j += 3
            j = 0
            i += 1

    if triangle == "U":

        # Make lower triangle NaNs
        Q[np.tril_indices(n, -1)] = np.nan

        row = []
        col = []
        q1 = []
        q2 = []
        q3 = []
        i = 0
        j = 0
        while i < n:
            while j < n:
                if j < n - 2:
                    row.append(i + 1)
                    col.append(j + 1)
                    q1.append(Q[i, j])
                    q2.append(Q[i, j + 1])
                    q3.append(Q[i, j + 2])
                    j += 3

                if j == n - 2:
                    row.append(i + 1)
                    col.append(j + 1)
                    q1.append(Q[i, j])
                    q2.append(Q[i, j + 1])
                    q3.append(np.nan)
                    j += 3

                if j == n - 1:
                    row.append(i + 1)
                    col.append(j + 1)
                    q1.append(Q[i, j])
                    q2.append(np.nan)
                    q3.append(np.nan)
                    j += 3

            i += 1
            j = i

    dict_temp = {
        "row": row,
        "col": col,
        "q1": q1,
        "q2": q2,
        "q3": q3,
    }
    df = pd.DataFrame(dict_temp)

    # Return
    return df


def dataframe2sinex_solution_matrix_estimate(df, tri="L"):
    """
    This function reads in a dataframe of the
    SOLUTIONS/MATRIX_ESTIMATE block from a SINEX, the
    converts each row to a string in a list
    ready for writing to SINEX.

    :param dataframe df: dataframe of SOLUTION/ESTIMATE block
    :return: list of strings
    """

    # Upper or lower triangle
    triangle = tri

    lines = [f"+SOLUTION/MATRIX_ESTIMATE {triangle} COVA"]
    lines.append(
        "*PARA1 PARA2 ____PARA2+0__________ ____PARA2+1__________ ____PARA2+2__________"
    )

    for i in range(len(df.row)):
        p1 = df.row.iloc[i]
        p2 = df.col.iloc[i]
        q1 = df.q1.iloc[i]
        q2 = df.q2.iloc[i]
        q3 = df.q3.iloc[i]
        l = " {:5n} {:5n} {:>21.14e} {:>21.14e} {:>21.14e}".format(p1, p2, q1, q2, q3)
        l = l.replace("nan", "")
        lines.append(l)

    lines.append(f"-SOLUTION/MATRIX_ESTIMATE {triangle} COVA")

    # Return
    return lines


def writeSINEX(
    fp,
    header=None,
    comment=None,
    siteID=None,
    solutionEpochs=None,
    solutionEstimate=None,
    solutionMatrixEstimate=None,
    fileReference=None,
    inputAcknowledgments=None,
    solutionStatistics=None,
    siteReceiver=None,
    siteAntenna=None,
    siteGpsPhaseCenter=None,
    siteEccentricity=None,
    solutionApriori=None,
    solutionMatrixApriori=None,
):
    """
    This function writes out SINEX blocks to a new SINEX file. The
    SINEX blocks can be obtained from many of the :literal:`read_sinex_...()`
    functions when writing the same input to output. Or can use the
    :literal:`dataframe2sinex_...()` functions when SINEX manipulations have
    occurred on that specific block. Input arguments are set to 'None'
    by default. This allows user to only write out the required blocks.

    :param str fp: Full filepath to output SINEX.
    :param str header: Single header line. Can get from read_sinex_header_line().
    :param list of str comment:  +FILE/COMMENT block. Can get from read_sinex_comments().
    :param list of str SiteID:  +SITE/ID block. Can get from read_sinex_site_id_block().
    :param list of str SolutionEpochs:  +SOLUTION/EPOCHS block. Can get from read_sinex_solution_epochs_block().
    :param list of str SolutionEstimate:  +SOLUTION/ESTIMATE block. Can get from read_sinex_solution_estimate_block().
    :param list of str SolutionMatricEstimate:  +SOLUTION/MATRIX_ESTIMATE block. Can get from read_sinex_solution_matrix_estimate_block().
    :param list of str fileReference:  +FILE/REFERENCE block. Can get from read_sinex_file_reference_block().
    :param list of str inputAcknowledgments:  +INPUT/ACKNOWLEDGEMENTS block. Can get from read_sinex_input_acknowledgments_block().
    :param list of str siteReceiver:  +SITE/RECEIVER block. Can get from read_sinex_site_receiver_block().
    :param list of str siteAntenna:  +SITE/ANTENNA block. Can get from read_sinex_site_antenna_block().
    :param list of str siteGpsPhaseCenter: +SITE/GPS_PHASE_CENTER block. Can get from read_site_gps_phase_center_block().
    :param list of str siteEccentricity: +SITE/ECCENTRICITY block. Can get from read_sinex_site_eccentricity_block().
    :param list of str solutionApriori: +SOLUTION/APRIORI block. Can get from read_sinex_solution_apriori_block().
    :param list of str solutionMatrixApriori: +SOLUTION/MATRIX_APRIORI block. Can get from read_sinex_solution_matric_apriori_block().

    :return: No return. But a new SINEX file will be written out to the file path (fp).

    """

    # Open File
    with open(fp, "w") as f:

        # Header
        if header == None:
            pass
        else:
            f.write("{}".format(header))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # Comment
        if comment == None:
            pass
        else:
            for i in range(len(comment)):
                f.write("{}\n".format(comment[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # FILE/REFERENCE
        if fileReference == None:
            pass
        else:
            for i in range(len(fileReference)):
                f.write("{}\n".format(fileReference[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # INPUT/ACKNOWLEDGMENTS
        if inputAcknowledgments == None:
            pass
        else:
            for i in range(len(inputAcknowledgments)):
                f.write("{}\n".format(inputAcknowledgments[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SOLUTION/STATISTICS
        if solutionStatistics == None:
            pass
        else:
            for i in range(len(solutionStatistics)):
                f.write("{}\n".format(solutionStatistics[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SITE/ID
        if siteID == None:
            pass
        else:
            for i in range(len(siteID)):
                f.write("{}\n".format(siteID[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SITE/RECEIVER
        if siteReceiver == None:
            pass
        else:
            for i in range(len(siteReceiver)):
                f.write("{}\n".format(siteReceiver[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SITE/ANTENNA
        if siteAntenna == None:
            pass
        else:
            for i in range(len(siteAntenna)):
                f.write("{}\n".format(siteAntenna[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SITE/GPS_PHASE_CENTER
        if siteGpsPhaseCenter == None:
            pass
        else:
            for i in range(len(siteGpsPhaseCenter)):
                f.write("{}\n".format(siteGpsPhaseCenter[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SITE/ECCENTRICITY
        if siteEccentricity == None:
            pass
        else:
            for i in range(len(siteEccentricity)):
                f.write("{}\n".format(siteEccentricity[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SOLUTION/EPOCHS
        if solutionEpochs == None:
            pass
        else:
            for i in range(len(solutionEpochs)):
                f.write("{}\n".format(solutionEpochs[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SOLUTION/ESTIMATE
        if solutionEstimate == None:
            pass
        else:
            for i in range(len(solutionEstimate)):
                f.write("{}\n".format(solutionEstimate[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SOLUTION/APRIORI
        if solutionApriori == None:
            pass
        else:
            for i in range(len(solutionApriori)):
                f.write("{}\n".format(solutionApriori[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SOLUTION/MATRIX_ESTIMATE
        if solutionMatrixEstimate == None:
            pass
        else:
            for i in range(len(solutionMatrixEstimate)):
                f.write("{}\n".format(solutionMatrixEstimate[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # SOLUTION/MATRIX_APRIORI
        if solutionMatrixApriori == None:
            pass
        else:
            for i in range(len(solutionMatrixApriori)):
                f.write("{}\n".format(solutionMatrixApriori[i]))

            f.write(
                "*-------------------------------------------------------------------------------\n"
            )

        # End Line
        f.write("%ENDSNX")

    f.close()


def remove_stns_sinex(sinex, sites):
    """
    This function removes a list sites from a SINEX file

    :param sinex: input SINEX file
    :param sites: list of the sites to be removed
    :return: SINEX file output.snx
    """

    # Open the output file
    with open("output.snx", "w") as out:

        # Get header line and update the creation time and the number of
        # parameter estimates. Write the updated header line to the new file
        header = read_sinex_header_line(sinex)
        old_creation_time = header[15:27]
        creation_time = set_creation_time()
        header = header.replace(old_creation_time, creation_time)
        old_num_params = header[60:65]
        if header[70:71] == "V":
            num_stn_params = 6
        else:
            num_stn_params = 3
        solution_epochs = read_sinex_solution_epochs_block(sinex)
        num_stns_to_remove = 0
        for line in solution_epochs:
            site = line[1:5]
            if site in sites:
                num_stns_to_remove += 1
        del solution_epochs
        num_params = int(old_num_params) - num_stn_params * num_stns_to_remove
        num_params = "{:05d}".format(num_params)
        header = header.replace(str(old_num_params), str(num_params))
        out.write(header)

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +FILE/COMMENT block and write to output file
        comments = read_sinex_comments(sinex)
        for i in comments:
            out.write(f"{i}\n")

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the site ID block and write out the sites not being removed
        site_id = read_sinex_site_id_block(sinex)
        for line in site_id:
            if line.startswith("*") or line.startswith("+") or line.startswith("-"):
                out.write(f"{line}\n")
            else:
                site = line[1:5]
                if site not in sites:
                    out.write(f"{line}\n")
        del site_id

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the solution epochs block and write out the epochs of the
        # sites not being removed
        solution_epochs = read_sinex_solution_epochs_block(sinex)
        for line in solution_epochs:
            if line.startswith("*") or line.startswith("+") or line.startswith("-"):
                out.write(f"{line}\n")
            else:
                site = line[1:5]
                if site not in sites:
                    out.write(f"{line}\n")
        del solution_epochs

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the solution estimate block and write out the estimates of
        # the sites not being removed
        skip = []
        estimate_number = 0
        solution_estimate = read_sinex_solution_estimate_block(sinex)
        for line in solution_estimate:
            if line.startswith("*") or line.startswith("+") or line.startswith("-"):
                out.write(f"{line}\n")
            else:
                site = line[14:18]
                if site in sites:
                    num = int(line[0:6])
                    skip.append(num)
                else:
                    estimate_number += 1
                    number = "{:5d}".format(estimate_number)
                    line = " " + number + line[6:]
                    out.write(f"{line}\n")
        del solution_estimate

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the matrix estimate block and write out minus the sites
        # being removed
        vcv = {}
        solution_matrix_estimate = read_sinex_solution_matrix_estimate_block(sinex)
        if solution_matrix_estimate[0][26:27] == "L":
            matrix = "lower"
        elif solution_matrix_estimate[0][26:27] == "U":
            matrix = "upper"
        out.write(f"{solution_matrix_estimate[0]}\n")
        if solution_matrix_estimate[1].startswith("*"):
            out.write(f"{solution_matrix_estimate[1]}\n")
        for line in solution_matrix_estimate:
            if line.startswith(" "):
                cols = line.split()
                row = cols[0]
                for i in range(2, len(cols)):
                    try:
                        vcv[row].append(cols[i])
                    except KeyError:
                        vcv[row] = []
                        vcv[row].append(cols[i])
        block_end = solution_matrix_estimate[-1]
        del solution_matrix_estimate
        sub_vcv = {}
        sub_row = 0
        for i in range(1, len(vcv) + 1):
            if i not in skip:
                sub_row += 1
                if matrix == "lower":
                    for j in range(i):
                        if j + 1 not in skip:
                            try:
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
                            except KeyError:
                                sub_vcv[str(sub_row)] = []
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
                if matrix == "upper":
                    for j in range(len(vcv) - (i - 1)):
                        if j + i not in skip:
                            try:
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
                            except KeyError:
                                sub_vcv[str(sub_row)] = []
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
        for i in range(1, len(sub_vcv) + 1):
            para1 = "{:5d}".format(i)
            if matrix == "lower":
                j = -2
            elif matrix == "upper":
                j = i - 3
            while sub_vcv[str(i)]:
                j += 3
                para2 = "{:5d}".format(j)
                line = " " + para1 + " " + para2
                n = min([3, len(sub_vcv[str(i)])])
                for k in range(n):
                    val = "{:21.14e}".format(float(sub_vcv[str(i)].pop(0)))
                    line += " " + str(val)
                out.write(line + "\n")
        out.write(block_end)

        # Write out the trailer line
        out.write("%ENDSNX\n")

    return


def remove_velocity_sinex(sinex):
    """
    This function reads in a SINEX file and removes the
    velocity parameters, including the zeros of the SOLUTION/MATRIX_ESTIMATE block,

    :param str sinex: input SINEX file
    :return: SINEX file output.snx
    """

    # From header, confirm that the SINEX has velocity parameters
    header = read_sinex_header_line(sinex)
    header = header.strip()
    if header[-1] == "V":
        pass
    else:
        print(
            "Not removing velocities because SINEX file does not have velocity parameters."
        )
        exit()

    # Open the output file
    with open("output.snx", "w") as out:

        # With header line:
        # - update the creation time
        # - update number of parameter estimates
        # - remove 'V' from parameter list
        # - then write to file
        old_creation_time = header[15:27]
        creation_time = set_creation_time()
        header = header.replace(old_creation_time, creation_time)
        old_num_params = int(header[60:65])
        num_params = int(old_num_params / 2)
        header = header.replace(str(old_num_params), str(num_params))
        header = header.replace("V", "")
        out.write(header)
        out.write("\n")
        del header

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +FILE/COMMENT block and write to output file
        comments = read_sinex_comments(sinex)
        for i in comments:
            out.write(f"{i}\n")

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SITE/ID block and write to file
        site_id = read_sinex_site_id_block(sinex)
        for line in site_id:
            out.write(f"{line}\n")
        del site_id

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SOLUTION/EPOCHS block and write to file
        # - also collect count on number of sites for use later
        numSites = 0
        solution_epochs = read_sinex_solution_epochs_block(sinex)
        for line in solution_epochs:
            out.write(f"{line}\n")
            if line[0] != "+" and line[0] != "*" and line[0] != "-":
                numSites += 1
        del solution_epochs

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SOLUTION/ESTIMATE block:
        # - gather velocity indices
        # - remove velocity rows
        # - change indices to account for removed velocities
        # - write to file
        vel_indices = []
        estimate_number = 0
        solution_estimate = read_sinex_solution_estimate_block(sinex)
        for line in solution_estimate:
            if line[7:10] == "VEL":
                vel_indices.append(int(line[0:6]))
            elif line[0] == "+" or line[0] == "*" or line[0] == "-":
                out.write(f"{line}\n")
            else:
                estimate_number += 1
                number = "{:5d}".format(estimate_number)
                line = " " + number + line[6:]
                out.write(f"{line}\n")
        del solution_estimate

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SOLUTION/MATRIX_ESTIMATE block:
        # - identify if matrix is upper or lower triangle form
        # - form full matrix
        # - remove all rows/columns associated with velocity
        # - write to file with new matrix indices (para1, para2)
        solution_matrix_estimate = read_sinex_solution_matrix_estimate_block(sinex)
        if solution_matrix_estimate[0][26:27] == "L":
            matrix = "lower"
        elif solution_matrix_estimate[0][26:27] == "U":
            matrix = "upper"
        # Size of original matrix
        matrix_dim = numSites * 6
        # Setup matrix of zeros
        Q = zeros((matrix_dim, matrix_dim))
        # Write initial comment line(s), save last comment line, and form matrix
        for line in solution_matrix_estimate:
            if line[0] == "+":
                out.write(f"{line}\n")
                continue
            if line[0] == "*":
                out.write(f"{line}\n")
                continue
            if line[0] == "-":
                block_end = line
            else:
                lineMAT = re.split("\s+", line)
                lineMAT = list(filter(None, lineMAT))
                row = int(lineMAT[0])
                col = int(lineMAT[1])
            if len(lineMAT) >= 3:
                q1 = float(lineMAT[2])
                Q[row - 1, col - 1] = q1
                Q[col - 1, row - 1] = q1
            if len(lineMAT) >= 4:
                q2 = float(lineMAT[3])
                Q[row - 1, col] = q2
                Q[col, row - 1] = q2
            if len(lineMAT) >= 5:
                q3 = float(lineMAT[4])
                Q[row - 1, col + 1] = q3
                Q[col + 1, row - 1] = q3
        # Remove velocity rows/columns from matrix
        num_removed = 0
        for i in vel_indices:
            Q = delete(Q, i - 1 - num_removed, 0)
            Q = delete(Q, i - 1 - num_removed, 1)
            num_removed += 1
        # Write matrix to SINEX (lower triangle)
        if matrix == "lower":
            i = 0
            j = 0
            while i < len(Q):
                while j <= i:
                    out.write(f" {i+1:5d} {j+1:5d} {Q[i,j]:21.14E} ")
                    j += 1
                    if j <= i:
                        out.write(f"{Q[i,j]:21.14E} ")
                        j += 1
                    if j <= i:
                        out.write(f"{Q[i,j]:21.14E}")
                        j += 1
                    out.write(" \n")
                j = 0
                i += 1
        # Write matrix to SINEX (upper triangle)
        if matrix == "upper":
            j = 0
            for i in range(len(Q)):
                j = i
                while j < len(Q):
                    out.write(f" {i+1:5d} {j+1:5d} {Q[i,j]:21.14E} ")
                    j += 1
                    if j < len(Q):
                        out.write(f"{Q[i,j]:21.14E} ")
                        j += 1
                    if j < len(Q):
                        out.write(f"{Q[i,j]:21.14E}")
                        j += 1
                    out.write(" \n")
        # Write out end of block line, and delete large variables
        out.write(block_end)
        del solution_matrix_estimate
        del Q

        # Write out the trailer line
        out.write("%ENDSNX\n")

    return


def remove_matrixzeros_sinex(sinex):
    """
    This function reads in a SINEX file and removes the
    zeros from the SOLUTION/MATRIX_ESTIMATE block only.

    :param str sinex: input SINEX file
    :return: SINEX file output.snx
    """

    # Open the output file
    with open("output.snx", "w") as out:

        # With header line:
        # - update the creation time
        # - then write to file
        header = read_sinex_header_line(sinex)
        old_creation_time = header[15:27]
        creation_time = set_creation_time()
        header = header.replace(old_creation_time, creation_time)
        out.write(header)
        del header

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +FILE/COMMENT block and write to output file
        comments = read_sinex_comments(sinex)
        for i in comments:
            out.write(f"{i}\n")

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SITE/ID block and write to file
        site_id = read_sinex_site_id_block(sinex)
        for line in site_id:
            out.write(f"{line}\n")
        del site_id

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SOLUTION/EPOCHS block and write to file
        solution_epochs = read_sinex_solution_epochs_block(sinex)
        for line in solution_epochs:
            out.write(f"{line}\n")
        del solution_epochs

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SOLUTION/ESTIMATE block
        solution_estimate = read_sinex_solution_estimate_block(sinex)
        for line in solution_estimate:
            out.write(f"{line}\n")
        del solution_estimate

        out.write(
            "*-------------------------------------------------------------------------------\n"
        )

        # Read in the +SOLUTION/MATRIX_ESTIMATE block:
        # - Remove lines that contain only zeros
        solution_matrix_estimate = read_sinex_solution_matrix_estimate_block(sinex)
        for line in solution_matrix_estimate:
            col = line.split()
            numCol = len(col)
            if numCol == 3:
                if col[2] == "0.00000000000000e+00":
                    continue
            if numCol == 4:
                if (
                    col[2] == "0.00000000000000e+00"
                    and col[3] == "0.00000000000000e+00"
                ):
                    continue
            if numCol == 5:
                if (
                    col[2] == "0.00000000000000e+00"
                    and col[3] == "0.00000000000000e+00"
                    and col[4] == "0.00000000000000e+00"
                ):
                    continue
            out.write(line)
        del solution_matrix_estimate

        # Write out the trailer line
        out.write("%ENDSNX\n")

    return

#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
GNSS Module

In Development
"""

from datetime import datetime
from numpy import zeros
from geodepy.angles import DMSAngle


def list_sinex_blocks(file):
    """This script lists the blocks in a SINEX file

    :param str file: the input SINEX file
    """
    blocks = []
    with open(file) as f:
        for line in f.readlines():
            if line.startswith('+'):
                col = line.split(' ')
                block = col[0].replace('+', '')
                blocks.append(block.strip())
    for block in blocks:
        print(block)


def print_sinex_comments(file):
    """This script lists prints comments in a SINEX file

    :param str file: the input SINEX file
    """
    go = False
    with open(file) as f:
        for line in f.readlines():
            if line.startswith('+FILE/COMMENT'):
                go = True
            if go:
                print(line.strip())
            if line.startswith('-FILE/COMMENT'):
                go = False


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
    doy = '{:03d}'.format(doy)
    seconds = (now - now.replace(hour=0, minute=0, second=0, microsecond=0))\
        .total_seconds()
    seconds = '{:.0f}'.format(seconds)
    creation_time = year + ':' + doy + ':' + seconds

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


def read_sinex_estimate(file):
    """This function reads in the SOLUTION/ESTIMATE block of a SINEX file. It
    returns estimate, a list of tuples:

    estimate = [(code, soln, refEpoch, staX, staY, staZ, staX_sd, staY_sd,
                 staZ_sd[, velX, velY, velZ, velX_sd, velY_sd, velZ_sd])...]

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
                if line[7:10] == 'VEL':
                    velocities = True
                lines.append(line)
            if line[:18] == '+SOLUTION/ESTIMATE':
                go = True

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
            if not velocities:
                info = (code, soln, epoch, stax, stay, staz, stax_sd, stay_sd,
                        staz_sd)
                estimate.append(info)
        elif typ == 'VELX':
            velx = float(line[47:68])
            velx_sd = float(line[69:80])
        elif typ == 'VELY':
            vely = float(line[47:68])
            vely_sd = float(line[69:80])
        elif typ == 'VELZ':
            velz = float(line[47:68])
            velz_sd = float(line[69:80])
            info = (code, soln, epoch, stax, stay, staz, stax_sd, stay_sd,
                    staz_sd, velx, vely, velz, velx_sd, vely_sd, velz_sd)
            estimate.append(info)

    return estimate


def read_sinex_matrix(file):
    """This function reads in the SOLUTION/MATRIX_ESTIMATE block of a SINEX
    file. It returns matrix, a list of tuples:

    matrix = [(code, soln, var_x, covar_xy, covar_xz, var_y, covar_yz,
            var_z[, var_v_x, covar_v_xy, covar_v_xz, var_v_y, covar_v_yz,
            var_v_z])...]

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
            if line[:25] == '-SOLUTION/MATRIX_ESTIMATE':
                break
            if go and line[:12] == '*PARA1 PARA2':
                pass
            elif go:
                lines.append(line)
            if line[:25] == '+SOLUTION/MATRIX_ESTIMATE':
                if line[26] == 'L':
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
                info = (code[i], soln[i], element[6 * i][6 * i],
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
                        element[6 * i + 5][6 * i + 5])
                matrix.append(info)
        else:
            for i in range(len(code)):
                info = (code[i], soln[i], element[6 * i][6 * i],
                        element[6 * i][6 * i + 1], element[6 * i][6 * i + 2],
                        element[6 * i + 1][6 * i + 1],
                        element[6 * i + 1][6 * i + 2],
                        element[6 * i + 2][6 * i + 2],
                        element[6 * i + 3][6 * i + 3],
                        element[6 * i + 3][6 * i + 4],
                        element[6 * i + 3][6 * i + 5],
                        element[6 * i + 4][6 * i + 4],
                        element[6 * i + 4][6 * i + 5],
                        element[6 * i + 5][6 * i + 5])
                matrix.append(info)
    else:
        if lower_triangular:
            for i in range(len(code)):
                info = (code[i], soln[i], element[3 * i][3 * i],
                        element[3 * i + 1][3 * i],
                        element[3 * i + 1][3 * i + 1],
                        element[3 * i + 2][3 * i],
                        element[3 * i + 2][3 * i + 1],
                        element[3 * i + 2][3 * i + 2])
                matrix.append(info)
        else:
            for i in range(len(code)):
                info = (code[i], soln[i], element[3 * i][3 * i],
                        element[3 * i][3 * i + 1], element[3 * i][3 * i + 2],
                        element[3 * i + 1][3 * i + 1],
                        element[3 * i + 1][3 * i + 2],
                        element[3 * i + 2][3 * i + 2])
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
            if line[:8] == '-SITE/ID':
                break
            if go and line[:8] == '*CODE PT':
                pass
            elif go:
                lines.append(line)
            if line[:8] == '+SITE/ID':
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
            if line[:23] == '-SOLUTION/DISCONTINUITY':
                break
            elif go:
                lines.append(line)
            if line[:23] == '+SOLUTION/DISCONTINUITY':
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
            if line[:16] == '-SOLUTION/EPOCHS':
                break
            if go and line[:8] == '*Code PT':
                pass
            elif go:
                lines.append(line)
            if line[:16] == '+SOLUTION/EPOCHS':
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
    """This function reads in the header information of a SINEX file

        :param str sinex: input SINEX file
        return: block
        """

    block = []
    with open(sinex, 'r') as f:
        next(f)
        line = f.readline()
        while line:
            block.append(line)
            line = f.readline()
            if line.startswith('+SITE/ID'):
                break

    return block


def read_sinex_site_id_block(sinex):
    """This function reads in the SITE/ID block of a SINEX file

    :param str sinex: input SINEX file
    return: block
    """

    block = []
    go = False
    with open(sinex, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('+SITE/ID'):
                go = True
            if go:
                block.append(line)
            if line.startswith('-SITE/ID'):
                break
            line = f.readline()

    return block


def read_sinex_solution_epochs_block(sinex):
    """This function reads in the SOLUTION/EPOCHS block of a SINEX file

    :param str sinex: input SINEX file
    return: block
    """

    block = []
    go = False
    with open(sinex, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('+SOLUTION/EPOCHS'):
                go = True
            if go:
                block.append(line)
            if line.startswith('-SOLUTION/EPOCHS'):
                break
            line = f.readline()
    return block


def read_sinex_solution_estimate_block(sinex):
    """This function reads in the SOLUTION/ESTIMATE block of a SINEX
    file

    :param str sinex: input SINEX file
    return: block
    """

    block = []
    go = False
    with open(sinex, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('+SOLUTION/ESTIMATE'):
                go = True
            if go:
                block.append(line)
            if line.startswith('-SOLUTION/ESTIMATE'):
                break
            line = f.readline()
    return block


def read_sinex_solution_matrix_estimate_block(sinex):
    """This function reads in the SOLUTION/MATRIX_ESTIMATE block of a SINEX
    file

    :param str sinex: input SINEX file
    return: block
    """

    block = []
    go = False
    with open(sinex, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith('+SOLUTION/MATRIX_ESTIMATE'):
                go = True
            if go:
                block.append(line)
            if line.startswith('-SOLUTION/MATRIX_ESTIMATE'):
                break
            line = f.readline()
    return block


def remove_stns_sinex(sinex, sites):
    """This function removes a list sites from a SINEX file

    :param sinex: input SINEX file
    :param sites: list of the sites to be removed
    :return: SINEX file output.snx
    """

    # Open the output file
    with open('output.snx', 'w') as out:

        # Get header line and update the creation time and the number of
        # parameter estimates. Write the updated header line to the new file
        header = read_sinex_header_line(sinex)
        old_creation_time = header[15:27]
        creation_time = set_creation_time()
        header = header.replace(old_creation_time, creation_time)
        old_num_params = header[60:65]
        if header[70:71] == 'V':
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
        num_params = '{:05d}'.format(num_params)
        header = header.replace(str(old_num_params), str(num_params))
        out.write(header)

        # Read in the site ID block and write out the sites not being removed
        site_id = read_sinex_site_id_block(sinex)
        for line in site_id:
            if line.startswith('*') or line.startswith('+') or \
                    line.startswith('-'):
                out.write(line)
            else:
                site = line[1:5]
                if site not in sites:
                    out.write(line)
        del site_id

        # Read in the solution epochs block and write out the epochs of the
        # sites not being removed
        solution_epochs = read_sinex_solution_epochs_block(sinex)
        for line in solution_epochs:
            if line.startswith('*') or line.startswith('+') or \
                    line.startswith('-'):
                out.write(line)
            else:
                site = line[1:5]
                if site not in sites:
                    out.write(line)
        del solution_epochs

        # Read in the solution estimate block and write out the estimates of
        # the sites not being removed
        skip = []
        estimate_number = 0
        solution_estimate = read_sinex_solution_estimate_block(sinex)
        for line in solution_estimate:
            if line.startswith('*') or line.startswith('+') or \
                    line.startswith('-'):
                out.write(line)
            else:
                site = line[14:18]
                if site in sites:
                    num = int(line[0:6])
                    skip.append(num)
                else:
                    estimate_number += 1
                    number = '{:5d}'.format(estimate_number)
                    line = ' ' + number + line[6:]
                    out.write(line)
        del solution_estimate

        # Read in the matrix estimate block and write out minus the sites
        # being removed
        vcv = {}
        solution_matrix_estimate = \
            read_sinex_solution_matrix_estimate_block(sinex)
        if solution_matrix_estimate[0][26:27] == 'L':
            matrix = 'lower'
        elif solution_matrix_estimate[0][26:27] == 'U':
            matrix = 'upper'
        out.write(solution_matrix_estimate[0])
        if solution_matrix_estimate[1].startswith('*'):
            out.write(solution_matrix_estimate[1])
        for line in solution_matrix_estimate:
            if line.startswith(' '):
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
        for i in range(1, len(vcv)+1):
            if i not in skip:
                sub_row += 1
                if matrix == 'lower':
                    for j in range(i):
                        if j+1 not in skip:
                            try:
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
                            except KeyError:
                                sub_vcv[str(sub_row)] = []
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
                if matrix == 'upper':
                    for j in range(len(vcv)-(i-1)):
                        if j+i not in skip:
                            try:
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
                            except KeyError:
                                sub_vcv[str(sub_row)] = []
                                sub_vcv[str(sub_row)].append(vcv[str(i)][j])
        for i in range(1, len(sub_vcv)+1):
            para1 = '{:5d}'.format(i)
            if matrix == 'lower':
                j = -2
            elif matrix == 'upper':
                j = i - 3
            while sub_vcv[str(i)]:
                j += 3
                para2 = '{:5d}'.format(j)
                line = ' ' + para1 + ' ' + para2
                n = min([3, len(sub_vcv[str(i)])])
                for k in range(n):
                    val = '{:21.14e}'.format(float(sub_vcv[str(i)].pop(0)))
                    line += ' ' + str(val)
                out.write(line + '\n')
        out.write(block_end)

        # Write out the trailer line
        out.write('%ENDSNX\n')

    return

#!/usr/bin/env python3

"""
Geoscience Australia - GeodePy Package
File input/output module

In Development
"""


def rd_snx_est(file):
    """This function reads in the SOLUTION/ESTIMATE block of the given SINEX
    file. It returns solEstimate, a list of tuples:

    solEstimate = [(code, soln, refEpoch, staX, staY, staZ, staX_sd, staY_sd,
                    staZ_sd, velX, velY, velZ, velX_sd, velY_sd, velZ_sd)...]

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

    :param file: the input SINEX file
    :return: solEstimate
    """

    # Open and read in the given SINEX file
    with open(file) as f:
        lines = f.readlines()

    # Read the solution estimates into a list
    solEstimate = []
    go = False
    for line in lines:
        if line[:18] == '-SOLUTION/ESTIMATE':
            go = False
        if go:
            if component == 'BOTH':
                solEstimate.append(tuple(line.split()))
            else:
                if line[7:10] == component:
                    solEstimate.append(tuple(line.split()))
        if line[:11] == '*INDEX TYPE':
            go = True
    return solEstimate
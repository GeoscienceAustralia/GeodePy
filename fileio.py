#!/usr/bin/env python3

"""
Geoscience Australia - GeodePy Package
File input/output module

In Development
"""


def rd_snx_est(file, component):
    """This function reads in the SOLUTION/ESTIMATE block of the given SINEX
    file.

    :param file: the input SINEX file
    :param component: STA for positions only; VEL for velocities only; and BOTH
        for positions and velocities
    :return:
    """

    # Open and read in the given SINEX file
    with open(file) as f:
        lines = f.readlines()

    # Read the solution estimates into a list
    estimate = []
    go = False
    for line in lines:
        if line[:18] == '-SOLUTION/ESTIMATE':
            go = False
        if go:
            if component == 'BOTH':
                estimate.append(tuple(line.split()))
            else:
                if line[7:10] == component:
                    estimate.append(tuple(line.split()))
        if line[:11] == '*INDEX TYPE':
            go = True
    return estimate
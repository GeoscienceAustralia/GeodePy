#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Data Converter - Class Tools Module
"""


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


def first_vel_cfg(cfg_list):
    # Find First Velocity Parameter Inputs in cfg_list
    for group in cfg_list:
        group_header = group[0].lower()
        if group_header.startswith('first vel'):
            wavelength = float(group[1])
            temperature = float(group[2])
            pressure = float(group[3])
            rel_humidity = float(group[4])
            return wavelength, temperature, pressure, rel_humidity
    else:
        return None


def stdev_cfg(cfg_list):
    # Find Survey Apriori Standard Deviations in cfg_list
    for group in cfg_list:
        group_header = group[0].lower()
        if group_header.startswith('standard dev'):
            stdev_angular = float(group[1])
            stdev_pointing = float(group[2])
            stdev_distconst = float(group[3])
            stdev_distppm = float(group[4])
            return stdev_angular, stdev_pointing, stdev_distconst, stdev_distppm
    else:
        return None

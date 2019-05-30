#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Survey Data Converter - Class Tools Module
"""

import itertools
import operator
from statistics import mean
from geodepy.convert import DMSAngle, dec2dms
from geodepy.survey import first_vel_corrn


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
                 face='FL', rounds=0.5, hz_obs=DMSAngle(0, 0, 0), va_obs=DMSAngle(0, 0, 0),
                 sd_obs=0.0, hz_dist=0.0, vert_dist=0.0):
        self.from_id = from_id
        self.to_id = to_id
        self.inst_height = inst_height
        self.target_height = target_height
        self.face = face
        self.rounds = rounds
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
                + '; rounds ' + repr(self.rounds)
                + '; hz_obs ' + repr(self.hz_obs)
                + '; va_obs ' + repr(self.va_obs)
                + '; sd_obs ' + repr(self.sd_obs)
                + '}')

    def changeface(self):
        # Change Horizontal Angle
        if 0 <= self.hz_obs.degree < 180:
            hz_switch = self.hz_obs + DMSAngle(180)
        elif 180 <= self.hz_obs.degree < 360:
            hz_switch = self.hz_obs - DMSAngle(180)
        else:
            raise ValueError('Horizontal Angle out of range (0 to 360 degrees)')
        # Change Vertical Angle
        if 0 <= self.va_obs.degree < 360:
            va_switch = DMSAngle(360) - self.va_obs
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
                           self.rounds,
                           hz_switch,
                           va_switch,
                           self.sd_obs,
                           self.hz_obs,
                           self.vert_dist)


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
        if abs(ob1.hz_obs.dec() - ob2.hz_obs.dec()) < 180:
            meaned_hz = (ob1.hz_obs + ob2.hz_obs) / 2
        elif ob1.hz_obs < ob2.hz_obs:
            ob2_shift = DMSAngle(360) - ob2.hz_obs
            meaned_hz = (ob1.hz_obs + ob2_shift) / 2
            if meaned_hz < DMSAngle(0):
                meaned_hz = meaned_hz + DMSAngle(360)
        elif ob2.hz_obs < ob1.hz_obs:
            ob1_shift = DMSAngle(360) - ob1.hz_obs
            meaned_hz = (ob1_shift + ob2.hz_obs) / 2
            if meaned_hz < DMSAngle(0):
                meaned_hz = meaned_hz + DMSAngle(360)
        meaned_va = (ob1.va_obs + ob2.va_obs) / 2
        meaned_sd = round(((ob1.sd_obs + ob2.sd_obs) / 2), 4)
        return Observation(ob1.from_id,
                           ob1.to_id,
                           ob1.inst_height,
                           ob1.target_height,
                           ob1.face,
                           ob1.rounds + ob2.rounds,
                           meaned_hz,
                           meaned_va,
                           meaned_sd)


def reducesetup(obslist, strict=False, zerodist=False, meanmulti=False):
    """
    Takes a list of Observations from one setup and
    means together FL, FR pairs of Observations.
    :param obslist: List of Observations (i.e. from one InstSetup)
    :param strict: If True, all single-face Obs are ignored. If False, all
    single-face Obs are included and converted to Face Left
    :param zerodist: If True, Obs with Slope Distance of Zero are included.
    If False, these are ignored
    :param meanmulti: If True, multiple rounds of observations are reduced
    to a single FL-sense Observation
    :return: a reduced list of Observations
    """
    # Remove obs with sd_obs == 0
    if not zerodist:
        for ob in obslist:
            if ob.sd_obs == 0:
                obslist.remove(ob)

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
    # Mean multiple repeated measurements into a single FL observation
    if meanmulti:
        multiob = []
        to_ids = set([i.to_id for i in meanedobs])
        for id in to_ids:
            matchedobs = []
            for ob in meanedobs:
                if ob.to_id == id:
                    matchedobs.append(ob)
            ob1 = matchedobs[0]
            repeat_hz = [ob.hz_obs.dec() for ob in matchedobs]
            repeat_hz = sorted(repeat_hz)
            repeat_va = [ob.va_obs.dec() for ob in matchedobs]
            repeat_sd = [ob.sd_obs for ob in matchedobs]
            # Finds where Horizontal Obseravtions are either side of 360, converts those on 360 side to negative
            hz_diff = [j - i for i, j in zip(repeat_hz[:-1], repeat_hz[1:])]
            bigjump = 1e100
            for num, i in enumerate(hz_diff):
                if i > 180:
                    bigjump = num
            for num, ob in enumerate(repeat_hz):
                if num > bigjump:
                    repeat_hz[num] = repeat_hz[num] - 360
            # Mean Repeated Observations
            mean_hz = mean(repeat_hz)
            if mean_hz < 0:
                mean_hz + 360
            mean_va = mean(repeat_va)
            mean_sd = mean(repeat_sd)
            # Compute number of rounds of observations completed
            sum_rounds = 0
            for ob in matchedobs:
                sum_rounds += ob.rounds
            # Output Meaned Observation
            multiob.append(Observation(ob1.from_id,
                                       id,
                                       ob1.inst_height,
                                       ob1.target_height,
                                       ob1.face,
                                       sum_rounds,
                                       dec2dms(mean_hz),
                                       dec2dms(mean_va),
                                       mean_sd))
        meanedobs = multiob
    # Order list of meaned obs
    sorted_meanedobs = sorted(meanedobs, key=operator.attrgetter('hz_obs'))
    return sorted_meanedobs


def first_vel_observations(obslist, params, temp, pressure, rel_humidity):
    """
    Performs a first velocity correction for all observed slope distances in a list of Observations
    :param obslist: List of Observations (i.e. from one InstSetup)
    :param params: Tuple of First Velocity Parameters C and D (see function first_vel_params)
    :param temp: Observed Temperature (degrees Celsius)
    :param pressure: Observed Pressure (hectopascals or millibars)
    :param rel_humidity: Observed Relative Humidity (percentage)
    :return: List of Observations with First Velocity Correction applied
    """
    for obs in obslist:
        obs.sd_obs = obs.sd_obs + first_vel_corrn(obs.sd_obs, params, temp, pressure, rel_humidity)
    return obslist

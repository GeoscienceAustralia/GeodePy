#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Convert Module
"""

from math import sin, cos, atan2, radians, degrees, sqrt


class DMSAngle(object):
    def __init__(self, degree=0, minute=0, second=0.0):
        # Set sign of object based on sign of any variable
        if degree == 0:
            if str(degree)[0] == '-':
                self.sign = -1
            elif minute < 0:
                self.sign = -1
            elif second < 0:
                self.sign = -1
            else:
                self.sign = 1
        elif degree > 0:
            self.sign = 1
        else:  # degree < 0
            self.sign = -1
        self.degree = abs(int(degree))
        self.minute = abs(int(minute))
        self.second = abs(second)

    def __repr__(self):
        if self.sign == -1:
            signsymbol = '-'
        elif self.sign == 1:
            signsymbol = '+'
        else:
            signsymbol = 'error'
        return '{DMSAngle: ' + signsymbol + str(self.degree) + 'd ' +\
               str(self.minute) + 'm ' + str(self.second) + 's}'

    def __add__(self, other):
        try:
            return dec2dms(self.dec() + other.dec())
        except AttributeError:
            raise TypeError('Can only add DMSAngle and/or DDMAngle objects '
                            'together')

    def __radd__(self, other):
        try:
            return dec2dms(other.dec() + self.dec())
        except AttributeError:
            raise TypeError('Can only add DMSAngle and/or DDMAngle objects '
                            'together')

    def __sub__(self, other):
        try:
            return dec2dms(self.dec() - other.dec())
        except AttributeError:
            raise TypeError('Can only subtract DMSAngle and/or DDMAngle '
                            'objects together')

    def __rsub__(self, other):
        try:
            return dec2dms(other.dec() - self.dec())
        except AttributeError:
            raise TypeError('Can only subtract DMSAngle and/or DDMAngle '
                            'objects together')

    def __mul__(self, other):
        try:
            return dec2dms(self.dec() * other)
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __rmul__(self, other):
        try:
            return dec2dms(other * self.dec())
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __truediv__(self, other):
        try:
            return dec2dms(self.dec() / other)
        except TypeError:
            raise TypeError('Division only defined between DMSAngle Object '
                            'and Int or Float')

    def __abs__(self):
        return DMSAngle(self.degree, self.minute, self.second)

    def __neg__(self):
        if self.sign == 1:
            return DMSAngle(-self.degree, -self.minute, -self.second)
        else:  # sign == -1
            return DMSAngle(self.degree, self.minute, self.second)

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def __lt__(self, other):
        return self.dec() < other.dec()

    def __gt__(self, other):
        return self.dec() > other.dec()

    def dec(self):
        return self.sign * (self.degree + (self.minute / 60) +
                            (self.second / 3600))

    def hp(self):
        return self.sign * (self.degree + (self.minute / 100) +
                            (self.second / 10000))

    def ddm(self):
        if self.sign == 1:
            return DDMAngle(self.degree, self.minute + (self.second/60))
        else:  # sign == -1
            return -DDMAngle(self.degree, self.minute + (self.second/60))


class DDMAngle(object):
    def __init__(self, degree=0, minute=0.0):
        # Set sign of object based on sign of any variable
        if degree == 0:
            if str(degree)[0] == '-':
                self.sign = -1
            elif minute < 0:
                self.sign = -1
            else:
                self.sign = 1
        elif degree > 0:
            self.sign = 1
        else:  # degree < 0
            self.sign = -1
        self.degree = abs(int(degree))
        self.minute = abs(minute)

    def __repr__(self):
        if self.sign == -1:
            signsymbol = '-'
        elif self.sign == 1:
            signsymbol = '+'
        else:
            signsymbol = 'error'
        return '{DDMAngle: ' + signsymbol + str(self.degree) + 'd ' + \
               str(self.minute) + 'm}'

    def __add__(self, other):
        try:
            return dec2ddm(self.dec() + other.dec())
        except AttributeError:
            raise TypeError('Can only add DMSAngle and/or DDMAngle objects '
                            'together')

    def __radd__(self, other):
        try:
            return dec2ddm(other.dec() + self.dec())
        except AttributeError:
            raise TypeError('Can only add DMSAngle and/or DDMAngle objects '
                            'together')

    def __sub__(self, other):
        try:
            return dec2ddm(self.dec() - other.dec())
        except AttributeError:
            raise TypeError('Can only add DMSAngle and/or DDMAngle objects '
                            'together')

    def __rsub__(self, other):
        try:
            return dec2ddm(other.dec() - self.dec())
        except AttributeError:
            raise TypeError('Can only add DMSAngle and/or DDMAngle objects '
                            'together')

    def __mul__(self, other):
        try:
            return dec2ddm(self.dec() * other)
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __rmul__(self, other):
        try:
            return dec2ddm(other * self.dec())
        except TypeError:
            raise TypeError('Multiply only defined between DMSAngle Object '
                            'and Int or Float')

    def __truediv__(self, other):
        try:
            return dec2ddm(self.dec() / other)
        except TypeError:
            raise TypeError('Division only defined between DMSAngle Object '
                            'and Int or Float')

    def __abs__(self):
        return DDMAngle(self.degree, self.minute)

    def __neg__(self):
        if self.sign == 1:
            return DDMAngle(-self.degree, -self.minute)
        else:  # sign == -1
            return DDMAngle(self.degree, self.minute)

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def __lt__(self, other):
        return self.dec() < other.dec()

    def __gt__(self, other):
        return self.dec() > other.dec()

    def dec(self):
        return self.sign * (self.degree + (self.minute / 60))

    def hp(self):
        minute_int, second = divmod(self.minute, 1)
        return self.sign * (self.degree + (minute_int / 100) +
                            (second * 0.006))

    def dms(self):
        minute_int, second = divmod(self.minute, 1)
        return self.sign * DMSAngle(self.degree, minute_int, second * 60)


def dec2hp(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    hp = degree + (minute / 100) + (second / 10000)
    return hp if dec >= 0 else -hp


def hp2dec(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    dec = degree + (minute / 60) + (second / 360)
    return dec if hp >= 0 else -dec


def dec2dms(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    return DMSAngle(degree, minute, second) if dec >= 0 \
        else DMSAngle(-degree, minute, second)


def dec2ddm(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    minute = minute + (second / 60)
    return DDMAngle(degree, minute) if dec >= 0 else DDMAngle(-degree, minute)


def hp2dms(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    return DMSAngle(degree, minute, second * 10) if hp >= 0 \
        else DMSAngle(-degree, minute, second * 10)


def hp2ddm(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    minute = minute + (second / 6)
    return DDMAngle(degree, minute) if hp >= 0 else DDMAngle(-degree, minute)


def polar2rect(r, theta):
    x = r * sin(radians(theta))
    y = r * cos(radians(theta))
    return x, y


def rect2polar(x, y):
    r = sqrt(x ** 2 + y ** 2)
    theta = atan2(x, y)
    if theta < 0:
        theta = degrees(theta) + 360
    else:
        theta = degrees(theta)
    return r, theta


def dd2sec(dd):
    minutes, seconds = divmod(abs(dd) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    sec = (degrees * 3600) + (minutes * 60) + seconds
    return sec if dd >= 0 else -sec


def dec2hp_v(dec):
    minutes, seconds = divmod(abs(dec) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    hp = degrees + (minutes / 100) + (seconds / 10000)
    hp[dec <= 0] = -hp[dec <= 0]
    return hp


def hp2dec_v(hp):
    degmin, seconds = divmod(abs(hp) * 1000, 10)
    degrees, minutes = divmod(degmin, 100)
    dec = degrees + (minutes / 60) + (seconds / 360)
    dec[hp <= 0] = -dec[hp <= 0]
    return dec

#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Angles Module
"""

from math import radians


class DECAngle(float):
    """
    Class for working with angles in Decimal Degrees
    Note: GeodePy also supports working with angles in Decimal Degrees as floats
    """

    def __init__(self, dec_angle=0.0):
        """
        :param dec_angle: float Decimal Degrees angle
        """
        super().__init__()
        self.dec_angle = float(dec_angle)

    def __repr__(self):
        if self.dec_angle >= 0:
            return '{DECAngle: +' + str(self.dec_angle) + '}'
        else:  # negative
            return '{DECAngle: ' + str(self.dec_angle) + '}'

    def __add__(self, other):
        try:
            return DECAngle(self.dec() + other.dec())
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __radd__(self, other):
        try:
            return DECAngle(other.dec() + self.dec())
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __sub__(self, other):
        try:
            return DECAngle(self.dec() - other.dec())
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __rsub__(self, other):
        try:
            return DECAngle(other.dec() - self.dec())
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __mul__(self, other):
        try:
            return DECAngle(self.dec() * other)
        except TypeError:
            raise TypeError('Multiply only defined between DECAngle Object '
                            'and Int or Float')

    def __rmul__(self, other):
        try:
            return DECAngle(other * self.dec())
        except TypeError:
            raise TypeError('Multiply only defined between DECAngle Object '
                            'and Int or Float')

    def __truediv__(self, other):
        try:
            return DECAngle(self.dec() / other)
        except TypeError:
            raise TypeError('Division only defined between DECAngle Object '
                            'and Int or Float')

    def __abs__(self):
        return DECAngle(abs(self.dec_angle))

    def __neg__(self):
        return DECAngle(-self.dec())

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def __lt__(self, other):
        return self.dec() < other.dec()

    def __gt__(self, other):
        return self.dec() > other.dec()

    def __str__(self):
        return str(self.dec_angle)

    def __round__(self, n=None):
        return DECAngle(round(self.dec_angle, n))

    def rad(self):
        """
        Convert to radians
        :return: radians
        :rtype: float
        """
        return radians(self.dec_angle)

    def dec(self):
        """
        Convert to Decimal Degrees (float)
        :return: Decimal Degrees
        :rtype: float
        """
        return self.dec_angle

    def hp(self):
        """
        Convert to HP Notation
        :return: HP Notation (DDD.MMSSSS)
        :rtype: float
        """
        return dec2hp(self.dec_angle)

    def hpa(self):
        """
        Convert to HP Notation (class)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: HPAngle
        """
        return HPAngle(self.hp())

    def dms(self):
        """
        Convert to Degrees, Minutes, Seconds Object
        :return: Degrees, Minutes, Seconds Object
        :rtype: DMSAngle
        """
        return dec2dms(self.dec_angle)

    def ddm(self):
        """
        Convert to Degrees, Decimal Minutes Object
        :return: Degrees, Decimal Minutes Object
        :rtype: DDMAngle
        """
        return dec2ddm(self.dec_angle)


class DMSAngle(object):
    """
    Class for working with angles in Degrees, Minutes and Seconds format
    """
    def __init__(self, degree, minute=0, second=0.0):
        """
        :param degree: Angle: whole degrees component (floats truncated)
                       Alt: formatted string '±DDD MM SS.SSS'
        :param minute: Angle: whole minutes component (floats truncated)
        :param second: Angle: seconds component (floats preserved)
        """
        # Convert formatted string 'DDD MM SS.SSS' to DMSAngle
        if type(degree) == str:
            str_pts = degree.split(' ')
            degree = int(str_pts[0])
            minute = int(str_pts[1])
            second = float(str_pts[2])
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
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __radd__(self, other):
        try:
            return dec2dms(other.dec() + self.dec())
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __sub__(self, other):
        try:
            return dec2dms(self.dec() - other.dec())
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __rsub__(self, other):
        try:
            return dec2dms(other.dec() - self.dec())
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

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

    def __str__(self):
        return (str(self.sign * self.degree) + ' ' + str(self.minute) + ' '
                + str(self.second))

    def __round__(self, n=None):
        return DMSAngle(self.sign * self.degree,
                        self.minute,
                        round(self.second, n))

    def rad(self):
        """
        Convert to Radians
        :return: Radians
        :rtype: float
        """
        return radians(self.dec())

    def dec(self):
        """
        Convert to Decimal Degrees (float)
        :return: Decimal Degrees
        :rtype: float
        """
        return self.sign * (self.degree + (self.minute / 60) +
                            (self.second / 3600))

    def deca(self):
        """
        Convert to Decimal Degrees (class)
        :return: Decimal Degrees
        :rtype: DECAngle
        """
        return DECAngle(self.dec())

    def hp(self):
        """
        Convert to HP Notation (float)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: float
        """
        return self.sign * (self.degree + (self.minute / 100) +
                            (self.second / 10000))

    def hpa(self):
        """
        Convert to HP Notation (class)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: HPAngle
        """
        return HPAngle(self.hp())

    def ddm(self):
        """
        Convert to Degrees, Decimal Minutes Object
        :return: Degrees, Decimal Minutes Object
        :rtype: DDMAngle
        """
        if self.sign == 1:
            return DDMAngle(self.degree, self.minute + (self.second/60))
        else:  # sign == -1
            return -DDMAngle(self.degree, self.minute + (self.second/60))


class DDMAngle(object):
    """
    Class for working with angles in Degrees, Decimal Minutes format
    """
    def __init__(self, degree, minute=0.0):
        """
        :param degree: Angle: whole degrees component (floats truncated)
        :param minute: Angle:minutes component (floats preserved)
        """
        # Convert formatted string 'DDD MM.MMMM' to DDMAngle
        if type(degree) == str:
            str_pts = degree.split(' ')
            degree = int(str_pts[0])
            minute = int(str_pts[1])
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
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __radd__(self, other):
        try:
            return dec2ddm(other.dec() + self.dec())
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __sub__(self, other):
        try:
            return dec2ddm(self.dec() - other.dec())
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __rsub__(self, other):
        try:
            return dec2ddm(other.dec() - self.dec())
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

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

    def __str__(self):
        return str(self.sign * self.degree) + ' ' + str(self.minute)

    def __round__(self, n=None):
        return DDMAngle(self.sign * self.degree,
                        round(self.minute, n))

    def rad(self):
        """
        Convert to Radians
        :return: Radians
        :rtype: float
        """
        return radians(self.dec())

    def dec(self):
        """
        Convert to Decimal Degrees (float)
        :return: Decimal Degrees
        :rtype: float
        """
        return self.sign * (self.degree + (self.minute / 60))

    def deca(self):
        """
        Convert to Decimal Degrees (class)
        :return: Decimal Degrees
        :rtype: DECAngle
        """
        return DECAngle(self.dec())

    def hp(self):
        """
        Convert to HP Notation (float)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: float
        """
        minute_int, second = divmod(self.minute, 1)
        return self.sign * (self.degree + (minute_int / 100) +
                            (second * 0.006))

    def hpa(self):
        """
        Convert to HP Notation (class)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: HPAngle
        """
        return HPAngle(self.hp())

    def dms(self):
        """
        Convert to Degrees, Minutes, Seconds Object
        :return: Degrees, Minutes, Seconds Object
        :rtype: DMSAngle
        """
        minute_int, second = divmod(self.minute, 1)
        return self.sign * DMSAngle(self.degree, int(minute_int), second * 60)


class HPAngle(object):
    """
    Class for working with angles in Hewlett-Packard (HP) format
    Note: GeodePy also supports working with angles in HP format as floats
    """
    def __init__(self, hp_angle=0.0):
        """
        :param hp_angle: float HP angle
        """
        self.hp_angle = float(hp_angle)
        hp_dec_str = f'{self.hp_angle:.17f}'.split('.')[1]
        if int(hp_dec_str[0]) > 5:
            raise ValueError(f'Invalid HP Notation: 1st decimal place greater '
                             f'than 5: {self.hp_angle}')
        if len(hp_dec_str) > 2:
            if int(hp_dec_str[2]) > 5:
                raise ValueError(
                    f'Invalid HP Notation: 3st decimal place greater '
                    f'than 5: {self.hp_angle}')

    def __repr__(self):
        if self.hp_angle >= 0:
            return '{HPAngle: +' + str(self.hp_angle) + '}'
        else:  # negative
            return '{HPAngle: ' + str(self.hp_angle) + '}'

    def __add__(self, other):
        try:
            return HPAngle(dec2hp(self.dec() + other.dec()))
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __radd__(self, other):
        try:
            return HPAngle(dec2hp(other.dec() + self.dec()))
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __sub__(self, other):
        try:
            return HPAngle(dec2hp(self.dec() - other.dec()))
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __rsub__(self, other):
        try:
            return HPAngle(dec2hp(other.dec() - self.dec()))
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __mul__(self, other):
        try:
            return HPAngle(dec2hp(self.dec() * other))
        except TypeError:
            raise TypeError('Multiply only defined between Angle objects and '
                            'Int or Float')

    def __rmul__(self, other):
        try:
            return HPAngle(dec2hp(other * self.dec()))
        except TypeError:
            raise TypeError('Multiply only defined between Angle objects and '
                            'Int or Float')

    def __truediv__(self, other):
        try:
            return HPAngle(dec2hp(self.dec() / other))
        except TypeError:
            raise TypeError('Division only defined between HPAngle objects '
                            'and Int or Float')

    def __abs__(self):
        return HPAngle(abs(self.hp_angle))

    def __neg__(self):
        return HPAngle(self.hp_angle.__neg__())

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def __lt__(self, other):
        return self.dec() < other.dec()

    def __gt__(self, other):
        return self.dec() > other.dec()

    def __int__(self):
        return int(self.hp_angle)

    def __float__(self):
        return float(self.hp_angle)

    def __str__(self):
        return str(self.hp_angle)

    def __round__(self, n=None):
        return HPAngle(round(self.hp_angle, n))

    def rad(self):
        """
        Convert to Radians
        :return: Radians
        :rtype: float
        """
        return radians(hp2dec(self.hp_angle))

    def dec(self):
        """
        Convert to Decimal Degrees (float)
        :return: Decimal Degrees
        :rtype: float
        """
        return hp2dec(self.hp_angle)

    def deca(self):
        """
        Convert to Decimal Degrees (class)
        :return: Decimal Degrees
        :rtype: DECAngle
        """
        return DECAngle(self.dec())

    def hp(self):
        """
        Convert to HP Notation (float)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: float
        """
        return float(self.hp_angle)

    def dms(self):
        """
        Convert to Degrees, Minutes, Seconds Object
        :return: Degrees, Minutes, Seconds Object
        :rtype: DMSAngle
        """
        return hp2dms(self.hp_angle)

    def ddm(self):
        """
        Convert to Degrees, Decimal Minutes Object
        :return: Degrees, Decimal Minutes Object
        :rtype: DDMAngle
        """
        return hp2ddm(self.hp_angle)


def dec2hp(dec):
    """
    Converts Decimal Degrees to HP Notation (float)
    :param dec: Decimal Degrees
    :type dec: float
    :return: HP Notation (DDD.MMSSSS)
    :rtype: float
    """
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    hp = degree + (minute / 100) + (second / 10000)
    hp = round(hp, 16)
    return hp if dec >= 0 else -hp


def dec2hpa(dec):
    """
    Converts Decimal Degrees to HP Angle Object
    :param dec: Decimal Degrees
    :type dec: float
    :return: HP Angle Object (DDD.MMSSSS)
    :rtype: HPAngle
    """
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    hp = degree + (minute / 100) + (second / 10000)
    return hp if dec >= 0 else -hp


def hp2dec(hp):
    """
    Converts HP Notation to Decimal Degrees
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: Decimal Degrees
    :rtype: float
    """
    # Check if 1st and 3rd decimal place greater than 5 (invalid HP Notation)
    hp = float(hp)
    hp_dec_str = f'{hp:.17f}'.split('.')[1]
    if int(hp_dec_str[0]) > 5:
        raise ValueError(f'Invalid HP Notation: 1st decimal place greater '
                         f'than 5: {hp}')
    if len(hp_dec_str) > 2:
        if int(hp_dec_str[2]) > 5:
            raise ValueError(f'Invalid HP Notation: 3st decimal place greater '
                             f'than 5: {hp}')
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    dec = degree + (minute / 60) + (second / 360)
    dec = round(dec, 16)
    return dec if hp >= 0 else -dec


def hp2deca(hp):
    """
    Converts HP Notation to DECAngle Object
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: Decimal Degrees Object
    :rtype: DECAngle
    """
    return DECAngle(hp2dec(hp))


def dec2gon(dec):
    """
    Converts Decimal Degrees to Gradians
    :param dec: Decimal Degrees
    :type dec: float
    :return: Gradians
    :rtype: float
    """
    return 10/9 * dec


def dec2dms(dec):
    """
    Converts Decimal Degrees to Degrees, Minutes, Seconds Object
    :param dec: Decimal Degrees
    :type dec: float
    :return: Degrees, Minutes, Seconds Object
    :rtype: DMSAngle
    """
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    return (DMSAngle(degree, minute, second) if dec >= 0
            else DMSAngle(-degree, minute, second))


def dec2ddm(dec):
    """
    Converts Decimal Degrees to Degrees, Decimal Minutes Object
    :param dec: Decimal Degrees
    :type dec: float
    :return: Degrees, Decimal Minutes Object
    :rtype: DDMAngle
    """
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    minute = minute + (second / 60)
    return DDMAngle(degree, minute) if dec >= 0 else DDMAngle(-degree, minute)


def hp2rad(hp):
    """
    Converts HP Notation to radians
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: radians
    :rtype: float
    """
    return radians(hp2dec(hp))


def hp2dms(hp):
    """
    Converts HP Notation to Degrees, Minutes, Seconds Object
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: Degrees, Minutes, Seconds Object
    :rtype: DMSAngle
    """
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    return (DMSAngle(degree, minute, second * 10) if hp >= 0
            else DMSAngle(-degree, minute, second * 10))


def hp2ddm(hp):
    """
    Converts HP Notation to Degrees, Decimal Minutes Object
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: Degrees, Decimal Minutes Object
    :rtype: DDMAngle
    """
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    minute = minute + (second / 6)
    return DDMAngle(degree, minute) if hp >= 0 else DDMAngle(-degree, minute)


def dd2sec(dd):
    """
    Converts angle in decimal degrees to angle in seconds
    :param dd: Decimal Degrees
    :return: Seconds
    """
    minute, second = divmod(abs(dd) * 3600, 60)
    degree, minute = divmod(minute, 60)
    sec = (degree * 3600) + (minute * 60) + second
    return sec if dd >= 0 else -sec


def dec2hp_v(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    hp = degree + (minute / 100) + (second / 10000)
    hp[dec <= 0] = -hp[dec <= 0]
    return hp


def hp2dec_v(hp):
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    dec = degree + (minute / 60) + (second / 360)
    dec[hp <= 0] = -dec[hp <= 0]
    return dec


def angular_typecheck(angle):
    # Converts DMSAngle and DDMAngle Objects to Decimal Degrees
    if type(angle) is DMSAngle or type(angle) is DDMAngle:
        return angle.dec()
    else:
        return angle

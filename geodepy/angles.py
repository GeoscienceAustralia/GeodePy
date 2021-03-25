#!/usr/bin/env python3

"""
Geoscience Australia - Python Geodesy Package
Angles Module

GeodePy supports Angular Notation in 9 different formats

ABRV - FORMAT (type)
--------------------
rad  - Radians (stored as float)
dec  - Decimal Degrees (stored as float)
deca - Decimal Degrees (via DECAngle class)
hp   - Hewlett Packard (HP) Notation (stored as float)
hpa  - Hewlett Packard (HP) Notation (via HPAngle class)
gon  - Gradians (stored as float)
gona - Gradians (via GONAngle class)
dms  - Degrees, Minutes and Seconds Notation (via DMSAngle class)
ddm  - Degrees and Decimal Minutes Notation (via DDMAngle class)

Conversion between all formats is supported as shown below:

Radians to/from Decimal Degrees via builtin math.radians and math.degrees

Formats as floats to all other types via functions in the form abrv2abrv
e.g. gon2hpa()

DECAngle, HPAngle, GONAngle, DMSAngle and DDMAngle class objects via methods in
the form CLASS.abrv()
e.g. HPAngle(value).dec()

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

    def __int__(self):
        return int(self.dec_angle)

    def __float__(self):
        return float(self.dec_angle)

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

    def gon(self):
        """
        Convert to Gradians (float)
        :return: Gradians
        :rtype: float
        """
        return dec2gon(self.dec_angle)

    def gona(self):
        """
        Convert to Gradians (class)
        :return: Gradians
        :rtype: GONAngle
        """
        return GONAngle(dec2gon(self.dec_angle))

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
                    f'Invalid HP Notation: 3rd decimal place greater '
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

    def gon(self):
        """
        Convert to Gradians (float)
        :return: Gradians
        :rtype: float
        """
        return hp2gon(self.hp_angle)

    def gona(self):
        """
        Convert to Gradians (class)
        :return: Gradians
        :rtype: GONAngle
        """
        return GONAngle(self.gon())

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


class GONAngle(object):
    """
    Class for working with angles in Gradians (90 degrees == 100 Gradians)
    Note: GeodePy also supports working with angles in Gradians as floats
    """
    def __init__(self, gon_angle=0.0):
        """
        :param gon_angle: float Gradian angle
        """
        super().__init__()
        self.gon_angle = float(gon_angle)

    def __repr__(self):
        if self.gon_angle >= 0:
            return '{GONAngle: +' + str(self.gon_angle) + '}'
        else:  # negative
            return '{GONAngle: ' + str(self.gon_angle) + '}'

    def __add__(self, other):
        try:
            return GONAngle(dec2gon(self.dec() + other.dec()))
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __radd__(self, other):
        try:
            return GONAngle(dec2gon(other.dec() + self.dec()))
        except AttributeError:
            raise TypeError('Can only add Angle objects with .dec() method '
                            'together')

    def __sub__(self, other):
        try:
            return GONAngle(dec2gon(self.dec() - other.dec()))
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __rsub__(self, other):
        try:
            return GONAngle(dec2gon(other.dec() - self.dec()))
        except AttributeError:
            raise TypeError('Can only subtract Angle objects with .dec() method'
                            ' together')

    def __mul__(self, other):
        try:
            return GONAngle(dec2gon(self.dec() * other))
        except TypeError:
            raise TypeError('Multiply only defined between Angle objects and '
                            'Int or Float')

    def __rmul__(self, other):
        try:
            return GONAngle(dec2gon(other * self.dec()))
        except TypeError:
            raise TypeError('Multiply only defined between Angle objects and '
                            'Int or Float')

    def __truediv__(self, other):
        try:
            return GONAngle(dec2gon(self.dec() / other))
        except TypeError:
            raise TypeError('Division only defined between HPAngle objects '
                            'and Int or Float')

    def __abs__(self):
        return GONAngle(abs(self.gon_angle))

    def __neg__(self):
        return GONAngle(self.gon_angle.__neg__())

    def __eq__(self, other):
        return self.dec() == other.dec()

    def __ne__(self, other):
        return self.dec() != other.dec()

    def __lt__(self, other):
        return self.dec() < other.dec()

    def __gt__(self, other):
        return self.dec() > other.dec()

    def __int__(self):
        return int(self.gon_angle)

    def __float__(self):
        return float(self.gon_angle)

    def __str__(self):
        return str(self.gon_angle)

    def __round__(self, n=None):
        return GONAngle(round(self.gon_angle, n))

    def rad(self):
        """
        Convert to Radians
        :return: Radians
        :rtype: float
        """
        return radians(gon2dec(self.gon_angle))

    def dec(self):
        """
        Convert to Decimal Degrees (float)
        :return: Decimal Degrees
        :rtype: float
        """
        return gon2dec(self.gon_angle)

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
        return gon2hp(self.gon_angle)

    def hpa(self):
        """
        Convert to HP Notation (class)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: HPAngle
        """
        return HPAngle(gon2hp(self.gon_angle))

    def gon(self):
        """
        Convert to Gradians (float)
        :return: Gradians
        :rtype: float
        """
        return float(self.gon_angle)

    def dms(self):
        """
        Convert to Degrees, Minutes, Seconds Object
        :return: Degrees, Minutes, Seconds Object
        :rtype: DMSAngle
        """
        return gon2dms(self.gon_angle)

    def ddm(self):
        """
        Convert to Degrees, Decimal Minutes Object
        :return: Degrees, Decimal Minutes Object
        :rtype: DDMAngle
        """
        return gon2ddm(self.gon_angle)


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
                self.positive = False
            elif minute < 0:
                self.positive = False
            elif second < 0:
                self.positive = False
            else:
                self.positive = True
        elif degree > 0:
            self.positive = True
        else:  # degree < 0
            self.positive = False
        self.degree = abs(int(degree))
        self.minute = abs(int(minute))
        self.second = abs(second)

    def __repr__(self):
        if self.positive:
            signsymbol = '+'
        else:
            signsymbol = '-'
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
        if self.positive:
            return DMSAngle(-self.degree, -self.minute, -self.second)
        else:  # positive == False
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
        if self.positive:
            return (str(self.degree) + ' ' + str(self.minute) + ' '
                    + str(self.second))
        else:
            return ('-' + str(self.degree) + ' ' + str(self.minute) + ' '
                    + str(self.second))

    def __round__(self, n=None):
        if self.positive:
            return DMSAngle(self.degree, self.minute, round(self.second, n))
        else:
            return -DMSAngle(self.degree, self.minute, round(self.second, n))

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
        if self.positive:
            return self.degree + (self.minute / 60) + (self.second / 3600)
        else:
            return -(self.degree + (self.minute / 60) + (self.second / 3600))

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
        if self.positive:
            return self.degree + (self.minute / 100) + (self.second / 10000)
        else:
            return -(self.degree + (self.minute / 100) + (self.second / 10000))

    def hpa(self):
        """
        Convert to HP Notation (class)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: HPAngle
        """
        return HPAngle(self.hp())

    def gon(self):
        """
        Convert to Gradians (float)
        :return: Gradians
        :rtype: float
        """
        return dec2gon(self.dec())

    def gona(self):
        """
        Convert to Gradians (class)
        :return: Gradians
        :rtype: GONAngle
        """
        return GONAngle(self.gon())

    def ddm(self):
        """
        Convert to Degrees, Decimal Minutes Object
        :return: Degrees, Decimal Minutes Object
        :rtype: DDMAngle
        """
        if self.positive:
            return DDMAngle(self.degree, self.minute + (self.second/60))
        else:
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
            minute = float(str_pts[1])
        # Set sign of object based on sign of any variable
        if degree == 0:
            if str(degree)[0] == '-':
                self.positive = False
            elif minute < 0:
                self.positive = False
            else:
                self.positive = True
        elif degree > 0:
            self.positive = True
        else:  # degree < 0
            self.positive = False
        self.degree = abs(int(degree))
        self.minute = abs(minute)

    def __repr__(self):
        if self.positive:
            signsymbol = '+'
        else:
            signsymbol = '-'
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
        if self.positive:
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
        if self.positive:
            return str(self.degree) + ' ' + str(self.minute)
        else:
            return '-' + str(self.degree) + ' ' + str(self.minute)

    def __round__(self, n=None):
        if self.positive:
            return DDMAngle(self.degree, round(self.minute, n))
        else:
            return DDMAngle(-self.degree, -round(self.minute, n))

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
        if self.positive:
            return self.degree + (self.minute / 60)
        else:
            return -(self.degree + (self.minute / 60))

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
        if self.positive:
            return self.degree + (minute_int / 100) + (second * 0.006)
        else:
            return -(self.degree + (minute_int / 100) + (second * 0.006))

    def hpa(self):
        """
        Convert to HP Notation (class)
        :return: HP Notation (DDD.MMSSSS)
        :rtype: HPAngle
        """
        return HPAngle(self.hp())

    def gon(self):
        """
        Convert to Gradians (float)
        :return: Gradians
        :rtype: float
        """
        return dec2gon(self.dec())

    def gona(self):
        """
        Convert to Gradians (class)
        :return: Gradians
        :rtype: GONAngle
        """
        return GONAngle(self.gon())

    def dms(self):
        """
        Convert to Degrees, Minutes, Seconds Object
        :return: Degrees, Minutes, Seconds Object
        :rtype: DMSAngle
        """
        minute_int, second = divmod(self.minute, 1)
        if self.positive:
            return DMSAngle(self.degree, int(minute_int), second * 60)
        else:
            return -DMSAngle(self.degree, int(minute_int), second * 60)


# Functions converting from Decimal Degrees (float) to other formats

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
    return HPAngle(dec2hp(dec))


def dec2gon(dec):
    """
    Converts Decimal Degrees to Gradians
    :param dec: Decimal Degrees
    :type dec: float
    :return: Gradians
    :rtype: float
    """
    return 10/9 * dec


def dec2gona(dec):
    """
    Converts Decimal Degrees to Gradians (class)
    :param dec: Decimal Degrees
    :type dec: float
    :return: Gradians
    :rtype: GONAngle
    """
    return GONAngle(dec2gon(dec))


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


# Functions converting from Hewlett-Packard (HP) format to other formats

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
            raise ValueError(f'Invalid HP Notation: 3rd decimal place greater '
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


def hp2rad(hp):
    """
    Converts HP Notation to radians
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: radians
    :rtype: float
    """
    return radians(hp2dec(hp))


def hp2gon(hp):
    """
    Converts HP Notation to Gradians
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: Gradians
    :rtype: float
    """
    return dec2gon(hp2dec(hp))


def hp2gona(hp):
    """
    Converts HP Notation to Gradians (class)
    :param hp: HP Notation (DDD.MMSSSS)
    :type hp: float
    :return: Gradians
    :rtype: GONAngle
    """
    return GONAngle(hp2gon(hp))


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


# Functions converting from Gradians format to other formats

def gon2dec(gon):
    """
    Converts Gradians to Decimal Degrees
    :param gon: Gradians
    :type gon: float
    :return: Decimal Degrees
    :rtype: float
    """
    return 9/10 * gon


def gon2deca(gon):
    """
    Converts Gradians to DECAngle Object
    :param gon: Gradians
    :type gon: float
    :return: Decimal Degrees Object
    :rtype: DECAngle
    """
    return DECAngle(gon2dec(gon))


def gon2hp(gon):
    """
    Converts Gradians to HP Notation (float)
    :param gon: Gradians
    :type gon: float
    :return: HP Notation (DDD.MMSSSS)
    :rtype: float
    """
    return dec2hp(gon2dec(gon))


def gon2hpa(gon):
    """
    Converts Gradians to HP Angle Object
    :param gon: Gradians
    :type gon: float
    :return: HP Angle Object (DDD.MMSSSS)
    :rtype: HPAngle
    """
    return HPAngle(gon2hp(gon))


def gon2rad(gon):
    """
    Converts Gradians to radians
    :param gon: Gradians
    :type gon: float
    :return: Radians
    :rtype: float
    """
    return radians(gon2dec(gon))


def gon2dms(gon):
    """
    Converts Gradians to Degrees, Minutes, Seconds Object
    :param gon: Gradians
    :type gon: float
    :return: Degrees, Minutes, Seconds Object
    :rtype: DMSAngle
    """
    return dec2dms(gon2dec(gon))


def gon2ddm(gon):
    """
    Converts Gradians to Degrees, Decimal Minutes Object
    :param gon: Gradians
    :type gon: float
    :return: Degrees, Decimal Minutes Object
    :rtype: DDMAngle
    """
    return dec2ddm(gon2dec(gon))


# Miscellaneous other useful functions

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
    # Converts Angle Objects to Decimal Degrees (float) for computations
    supported_types = [DMSAngle, DDMAngle, DECAngle, HPAngle, GONAngle]
    if type(angle) in supported_types:
        return angle.dec()
    else:
        return float(angle)

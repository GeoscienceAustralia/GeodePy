import struct
from datetime import datetime as dt
import numpy as np
import pandas as pd

'''
Created 12/11/2018

@Author: Jaimie Dodd

This script defines classes for reading the AUSGeoid2020 binary NTv2 file, and interpolating the Ellipsoid to AHD 
Separation, Deflection in the Prime Meridian and Deflection in the Prime Vertical of a given latitude and longitude via
bilinear or bicubic interpolation (user specified input or bicubic by default).

Values have been tested against the output of DynAdjust.
'''

class NTv2ReaderBinary(object):
    def __init__(self):
        self.numOrec = 0
        self.numSrec = 0
        self.numFile = 0
        self.gsType = ''
        self.version = ''
        self.systemF = ''
        self.systemT = ''
        self.majorF = 0
        self.minorF = 0
        self.majorT = 0
        self.minorT = 0
        self.subgrids = pd.DataFrame(columns=['CONTAINS', 'SUB_NAME', 'PARENT', 'CREATED', 'UPDATED', 'S_LAT', 'N_LAT',
                                              'E_LONG', 'W_LONG', 'LAT_INC', 'LONG_INC', 'GS_COUNT',
                                              'ELLIPSOID_TO_AHD_SEPARATION', 'DEFLECTION_PRIME_MERIDIAN',
                                              'DEFLECTION_PRIME_VERTICAL'])
        self.ellps2ahdsep = 0
        self.deflectionprimemeridian = 0
        self.deflectionprimevertical = 0

    def ntv2reader(self, file, lat, lon, interpolation_method='bicubic'):
        """
        Function for reading in the binary NTv2 file and

        :param file: NTv2 file (full location as string)
        :param lat: latitude of point of interest (decimal degrees)
        :param lon: longitude of point of interest (decimal degrees)
        :param interpolation_method: 'bilinear' or 'bicubic' ('bicubic' by default)
        :return: tuple containing ellipsoid to AHD separation, deflection in prime meridian and deflection in prime
        vertical.

        """

        # If interpolation method specified is something other than bilinear or bicubic, then raise ValueError
        interpolation_options = ['bilinear', 'bicubic']
        if interpolation_method not in interpolation_options:
            raise ValueError("Invalid Interpolation method. Expected one of: %s" % interpolation_options)

        # Read the NTv2 file as binary
        f = open(file, 'rb')

        # NUM_OREC
        f.seek(8, 1)

        byte = f.read(4)
        self.numOrec = int.from_bytes(byte, byteorder='little')

        f.seek(4, 1)

        # NUM_SREC
        f.seek(8, 1)

        byte = f.read(4)
        self.numSrec = int.from_bytes(byte, byteorder='little')

        f.seek(4, 1)

        # NUM_FILE
        f.seek(8, 1)

        byte = f.read(4)
        self.numFile = int.from_bytes(byte, byteorder='little')

        f.seek(4, 1)

        # GS_TYPE
        f.seek(8, 1)

        byte = f.read(8)
        self.gsType = byte.decode('utf8').strip('\x00').strip()

        # VERSION
        f.seek(8, 1)

        byte = f.read(8)
        self.version = byte.decode('utf8').strip('\x00').strip()

        # SYSTEM_F
        f.seek(8, 1)

        byte = f.read(8)
        self.systemF = byte.decode('utf8').strip('\x00').strip()

        # SYSTEM_T
        f.seek(8, 1)

        byte = f.read(8)
        self.systemT = byte.decode('utf8').strip('\x00').strip()

        # MAJOR_F
        f.seek(8, 1)

        byte = f.read(8)
        self.majorF = struct.unpack('d', byte)[0]

        # MINOR_F
        f.seek(8, 1)

        byte = f.read(8)
        self.minorF = struct.unpack('d', byte)[0]

        # MAJOR_T
        f.seek(8, 1)

        byte = f.read(8)
        self.majorT = struct.unpack('d', byte)[0]

        # MINOR_T
        f.seek(8, 1)

        byte = f.read(8)
        self.minorT = struct.unpack('d', byte)[0]

        # Convert lat and lon to seconds
        if self.gsType == 'SECONDS':
            lat = lat * 3600
            lon = lon * -3600

        # Sub Grids
        for i in range(0, self.numFile):
            self.subgrids = self.subgrids.append(SubGrid().findsubgrid(f, lat, lon, interpolation_method),
                                                 ignore_index=True)

        # Close the NTv2 file
        f.close()

        # Filter subgrids dataframe so that only subgrids with values remain
        self.subgrids = self.subgrids[self.subgrids['CONTAINS']]

        # If more than one subgrid which contains coordinates, then filter dataframe so subgrid with smallest gridwidth
        # remains
        if len(self.subgrids) > 1:
            self.subgrids = self.subgrids.loc[[self.subgrids.loc[:, 'LAT_INC'].idxmin()]]

        # If subgrids dataframe is not empty, then return values for ellipsoid to AHD separation, deflection of prime
        # meridian and deflection of prime vertical
        if not self.subgrids.empty:
            self.ellps2ahdsep = self.subgrids.iloc[0, -3]
            self.deflectionprimemeridian = self.subgrids.iloc[0, -2]
            self.deflectionprimevertical = self.subgrids.iloc[0, -1]
            return self.ellps2ahdsep, self.deflectionprimemeridian, self.deflectionprimevertical
        # If no subgrids contain coordinates, then raise ValueError
        else:
            raise ValueError("The coordinates supplied are outside the extents of the grid.")


class SubGrid(object):

    def __init__(self):
        self.subName = ''
        self.parent = ''
        self.created = ''
        self.updated = ''
        self.sLat = 0
        self.nLat = 0
        self.eLong = 0
        self.wLong = 0
        self.latInc = 0
        self.longInc = 0
        self.gsCount = 0
        self.Nval = 0
        self.Xi = 0
        self.Eta = 0
        self.contains = False

    def findsubgrid(self, f, lat, lon, interpolation_method='bicubic'):
        """
        Function to pull out metadata of a subgrid in the AUSGeoid2020 NTv2 file

        :param f: variable for open file NTv2 being read as binary
        :param lat: latitude of the point of interest (seconds)
        :param lon: longitude of the point of interest (negative seconds)
        :param interpolation_method: interpolation method as specified in ntv2reader() -
        'bilinear' or 'bicubic' ('bicubic' by default)

        :return: subgrid metadata in form of a dictionary or results from bilinear() or bicubic() (depending on
        interpolation method specified) if point lies within subgrid.
        """
        # SUB_NAME
        f.seek(8, 1)

        byte = f.read(8)
        self.subName = byte.decode('utf').strip('\x00').strip()

        # PARENT
        f.seek(8, 1)

        byte = f.read(8)
        self.parent = byte.decode('utf').strip('\x00').strip()

        # CREATED
        f.seek(8, 1)

        byte = f.read(8)
        self.created = dt.strptime(byte.decode('utf').strip('\x00').strip(), '%d%m%Y').strftime('%d/%m/%Y')

        # UPDATED
        f.seek(8, 1)

        byte = f.read(8)
        self.updated = dt.strptime(byte.decode('utf').strip('\x00').strip(), '%d%m%Y').strftime('%d/%m/%Y')

        # S_LAT
        f.seek(8, 1)

        byte = f.read(8)
        self.sLat = round(struct.unpack('d', byte)[0], 3)

        # N_LAT
        f.seek(8, 1)

        byte = f.read(8)
        self.nLat = round(struct.unpack('d', byte)[0], 3)

        # E_LONG
        f.seek(8, 1)

        byte = f.read(8)
        self.eLong = round(struct.unpack('d', byte)[0], 3)

        # W_LONG
        f.seek(8, 1)

        byte = f.read(8)
        self.wLong = round(struct.unpack('d', byte)[0], 3)

        # LAT_INC
        f.seek(8, 1)

        byte = f.read(8)
        self.latInc = round(struct.unpack('d', byte)[0], 6)

        # LONG_INC
        f.seek(8, 1)

        byte = f.read(8)
        self.longInc = round(struct.unpack('d', byte)[0], 6)

        # GS_COUNT
        f.seek(8, 1)

        byte = f.read(4)
        self.gsCount = int.from_bytes(byte, byteorder='little')

        f.seek(4, 1)

        # Check if coordinates fall within subgrid. If yes, return self.contains = True
        if self.sLat <= lat < self.nLat and self.eLong <= lon < self.wLong:
            self.contains = True

        # If self.contains is True, determine number of columns in grid, and row and column of node to bottom right of
        # point of interest, then call relevant interpolation method function
        if self.contains is True:
            # Determine number of columns
            numcols = 1 + int((self.wLong - self.eLong) / self.longInc)

            # Determine row and col numbers of node below right of point
            row = int((lat - self.sLat) / self.latInc)
            col = int((lon - self.eLong) / self.longInc)

            # If interpolation_method == 'bilinear', call bilinear function
            if interpolation_method == 'bilinear':
                return self.bilinear(f, lat, lon, numcols, row, col)
            # If interpolation_method == 'bicubic', call bicubic function
            elif interpolation_method == 'bicubic':
                return self.bicubic(f, lat, lon, numcols, row, col)

        # If point is not in subgrid, skip to end of subgrid and return sub grid metadata in form of dictionary
        else:
            f.seek(16 * self.gsCount, 1)

            return {'CONTAINS': self.contains, 'SUB_NAME': self.subName, 'PARENT': self.parent, 'CREATED': self.created,
                    'UPDATED': self.updated, 'S_LAT': self.sLat, 'N_LAT': self.nLat, 'E_LONG': self.eLong,
                    'W_LONG': self.wLong, 'LAT_INC': self.latInc, 'LONG_INC': self.longInc, 'GS_COUNT': self.gsCount,
                    'ELLIPSOID_TO_AHD_SEPARATION': np.nan, 'DEFLECTION_PRIME_MERIDIAN': np.nan,
                    'DEFLECTION_PRIME_VERTICAL': np.nan}

    def bilinear(self, f, lat, lon, numcols, row, col):

        """
        Function to perform bilinear interpolatoin of the Ellipsoid to AHD Separation, Deflection in Prime Meridian and
        Deflection in Prime Vertical of a point of interest.

        :param f: variable for open file NTv2 being read as binary
        :param lat: latitude of the point of interest (seconds)
        :param lon: longitude of the point of interest (negative seconds)
        :param numcols: number of columns in grid
        :param row: row number of point to the bottom right of the point of interest
        :param col: column number of point to the bottom right of the point of interest

        :return: dictionary of subgrid metadata and values for Ellipsoid to AHD Separation, Deflection in Prime Meridian
        and Deflection in Prime Vertical of point of interest found via bilinear interpolation.
        """

        #   |     |
        # --o-----o--
        #   |4    |3
        #   | *P  |
        # --o-----o--
        #   |2    |1
        #
        # o - Node position.

        # Determine position in the file of the four surrounding nodes
        pos1 = row * numcols + col
        pos2 = pos1 + 1
        pos3 = pos1 + numcols
        pos4 = pos3 + 1

        # Navigate to start of posA node
        f.seek(16 * pos1, 1)

        # Read in values for nodes A and B
        (pos1nval, pos1xi, pos1eta) = self.read_node(f)
        (pos2nval, pos2xi, pos2eta) = self.read_node(f)

        # Navigate to beginning of node C
        f.seek(16*(pos3 - pos2 - 1), 1)

        # Read in values for nodes C and D
        (pos3nval, pos3xi, pos3eta) = self.read_node(f)
        (pos4nval, pos4xi, pos4eta) = self.read_node(f)

        # Determine latitude and longitude of node A
        lat1 = self.sLat + row * self.latInc
        long1 = self.eLong + col * self.longInc

        # Determine interpolation scale factors
        x = round((lon - long1) / self.longInc, 6)
        y = round((lat - lat1) / self.latInc, 6)

        # Call bilinear interpolation function to determine N value, Xi and Eta
        # (Ellipsoid to AHD separation, deflection in prime meridian and deflection in prime vertical).
        self.Nval = round(bilinear_interpolation(pos1nval, pos2nval, pos3nval, pos4nval, x, y), 3)
        self.Xi = round(bilinear_interpolation(pos1xi, pos2xi, pos3xi, pos4xi, x, y), 2)
        self.Eta = round(bilinear_interpolation(pos1eta, pos2eta, pos3eta, pos4eta, x, y), 2)

        # Navigate to the end of the subgrid
        f.seek(16 * (self.gsCount - pos4 - 1), 1)

        # Return dictionary of information
        return {'CONTAINS': self.contains, 'SUB_NAME': self.subName, 'PARENT': self.parent, 'CREATED': self.created,
                'UPDATED': self.updated, 'S_LAT': self.sLat, 'N_LAT': self.nLat, 'E_LONG': self.eLong,
                'W_LONG': self.wLong, 'LAT_INC': self.latInc, 'LONG_INC': self.longInc, 'GS_COUNT': self.gsCount,
                'ELLIPSOID_TO_AHD_SEPARATION': self.Nval, 'DEFLECTION_PRIME_MERIDIAN': self.Xi,
                'DEFLECTION_PRIME_VERTICAL': self.Eta}

    def bicubic(self, f, lat, lon, numcols, row, col):
        """
        Function to perform bicubic interpolatoin of the Ellipsoid to AHD Separation, Deflection in Prime Meridian and
        Deflection in Prime Vertical of a point of interest.

        :param f: variable for open file NTv2 being read as binary
        :param lat: latitude of the point of interest (seconds)
        :param lon: longitude of the point of interest (negative seconds)
        :param numcols: number of columns in grid
        :param row: row number of point to the bottom right of the point of interest
        :param col: column number of point to the bottom right of the point of interest

        :return: dictionary of subgrid metadata and values for Ellipsoid to AHD Separation, Deflection in Prime Meridian
        and Deflection in Prime Vertical of point of interest found via bicubic interpolation.

        """

        #   |     |     |     |
        # --o-----o-----o-----o--
        #   |11   |12   |13   |14
        #   |     |     |     |
        # --o-----o-----o-----o--
        #   |10   |3    |4    |15
        #   |     | *P  |     |
        # --o-----o-----o-----o--
        #   |9    |2    |1    |16
        #   |     |     |     |
        # --o-----o-----o-----o--
        #   |8    |7    |6    |5
        #
        # o - Node position.

        # Determine position in the file of the sixteen surrounding nodes
        pos1 = row * numcols + col
        pos2 = pos1 + 1
        pos3 = pos2 + numcols
        pos4 = pos3 - 1
        pos5 = pos4 - 2 * numcols - 1
        pos6 = pos5 + 1
        pos7 = pos6 + 1
        pos8 = pos7 + 1
        pos9 = pos8 + numcols
        pos10 = pos9 + numcols
        pos11 = pos10 + numcols
        pos12 = pos11 - 1
        pos13 = pos12 - 1
        pos14 = pos13 - 1
        pos15 = pos14 - numcols
        pos16 = pos15 - numcols

        # Navigate to start of posA node
        f.seek(16 * pos5, 1)

        # Read in values for nodes 5-8
        (pos5nval, pos5xi, pos5eta) = self.read_node(f)
        (pos6nval, pos6xi, pos6eta) = self.read_node(f)
        (pos7nval, pos7xi, pos7eta) = self.read_node(f)
        (pos8nval, pos8xi, pos8eta) = self.read_node(f)

        # Navigate to start of pos16 node
        f.seek(16 * (pos16 - pos8 - 1), 1)

        # Read in values for nodes 16, 1, 2, and 9
        (pos16nval, pos16xi, pos16eta) = self.read_node(f)
        (pos1nval, pos1xi, pos1eta) = self.read_node(f)
        (pos2nval, pos2xi, pos2eta) = self.read_node(f)
        (pos9nval, pos9xi, pos9eta) = self.read_node(f)

        # Navigate to start of pos15 node
        f.seek(16 * (pos15 - pos9 - 1), 1)

        # Read in values for nodes 15, 3, 4 and 10
        (pos15nval, pos15xi, pos15eta) = self.read_node(f)
        (pos4nval, pos4xi, pos4eta) = self.read_node(f)
        (pos3nval, pos3xi, pos3eta) = self.read_node(f)
        (pos10nval, pos10xi, pos10eta) = self.read_node(f)

        # Navigate to start of pos14 node
        f.seek(16 * (pos14 - pos10 - 1), 1)

        # Read in values for nodes 11, 12, 13 and 14
        (pos14nval, pos14xi, pos14eta) = self.read_node(f)
        (pos13nval, pos13xi, pos13eta) = self.read_node(f)
        (pos12nval, pos12xi, pos12eta) = self.read_node(f)
        (pos11nval, pos11xi, pos11eta) = self.read_node(f)

        # Determine latitude and longitude of node A
        lat1 = self.sLat + row * self.latInc
        long1 = self.eLong + col * self.longInc

        # Determine interpolation scale factors
        x = round((lon - long1) / self.longInc, 6)
        y = round((lat - lat1) / self.latInc, 6)

        # Call bicubic_interpolation to determine nVal, Xi and Eta
        # (Ellipsoid to AHD separation, deflection in prime meridian and deflection in prime vertical).
        self.Nval = round(bicubic_interpolation(pos1nval, pos2nval, pos3nval, pos4nval, pos5nval, pos6nval, pos7nval,
                                                pos8nval, pos9nval, pos10nval, pos11nval, pos12nval, pos13nval,
                                                pos14nval, pos15nval, pos16nval, x, y), 3)
        self.Xi = round(bicubic_interpolation(pos1xi, pos2xi, pos3xi, pos4xi, pos5xi, pos6xi, pos7xi, pos8xi, pos9xi,
                                              pos10xi, pos11xi, pos12xi, pos13xi, pos14xi, pos15xi, pos16xi, x, y), 2)
        self.Eta = round(bicubic_interpolation(pos1eta, pos2eta, pos3eta, pos4eta, pos5eta, pos6eta, pos7eta, pos8eta,
                                               pos9eta, pos10eta, pos11eta, pos12eta, pos13eta, pos14eta, pos15eta,
                                               pos16eta, x, y), 2)

        # Navigate to the end of the subgrid
        f.seek(16 * (self.gsCount - pos11 - 1), 1)

        # Return dictionary of information
        return {'CONTAINS': self.contains, 'SUB_NAME': self.subName, 'PARENT': self.parent, 'CREATED': self.created,
                'UPDATED': self.updated, 'S_LAT': self.sLat, 'N_LAT': self.nLat, 'E_LONG': self.eLong,
                'W_LONG': self.wLong, 'LAT_INC': self.latInc, 'LONG_INC': self.longInc, 'GS_COUNT': self.gsCount,
                'ELLIPSOID_TO_AHD_SEPARATION': self.Nval, 'DEFLECTION_PRIME_MERIDIAN': self.Xi,
                'DEFLECTION_PRIME_VERTICAL': self.Eta}

    # Define function to read in node values
    @staticmethod
    def read_node(f):
        """
        Function to read in values of nodes

        :param f:  variable for open file NTv2 being read as binary
        :return: Ellipsoid to AHD Separation, Deflection in Prime Meridian and Deflection in Prime Vertical for node.

        """
        # Read ellipsoid to AHD separation value (N)
        byte = f.read(4)
        nval = round(struct.unpack('f', byte)[0], 3)
        # Read deflection in prime meridian value (Xi)
        byte = f.read(4)
        xi = round(struct.unpack('f', byte)[0], 3)
        # Read deflection in prime vertical value (Eta)
        byte = f.read(4)
        eta = round(struct.unpack('f', byte)[0], 3)
        # Skip to beginning of next node
        f.seek(4, 1)
        # Return values
        return nval, xi, eta


# Define function for bilinear interpolation
def bilinear_interpolation(n1, n2, n3, n4, x, y):
    """
    Bilinear interpolation of value for point of interest (P).

    :param n1: value at node 1
    :param n2: value at node 2
    :param n3: value at node 3
    :param n4: value at node 4
    :param x: interpolation scale factor for x axis
    :param y: interpolation scale factor for y axis

    :return: value at node P

    """

    a0 = n1
    a1 = round(n2 - n1, 3)
    a2 = round(n3 - n1, 3)
    a3 = round(n1 + n4 - n2 - n3, 3)
    p = a0 + (a1 * x) + (a2 * y) + (a3 * x * y)
    return p


# Define function for performing bicubic interpolation
def bicubic_interpolation(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, x, y):
    """
    Bicubic interpolation of value for point of interest (P).

    :param n1: value at node 1
    :param n2: value at node 2
    :param n3: value at node 3
    :param n4: value at node 4
    :param n5: value at node 5
    :param n6: value at node 6
    :param n7: value at node 7
    :param n8: value at node 8
    :param n9: value at node 9
    :param n10: value at node 10
    :param n11: value at node 11
    :param n12: value at node 12
    :param n13: value at node 13
    :param n14: value at node 14
    :param n15: value at node 16
    :param n16: value at node 17
    :param x: interpolation scale factor for x axis
    :param y: interpolation scale factor for y axis

    :return: value at node P
    """

    # Define the inverse of the coefficient matrix
    cinv = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                     [-3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0],
                     [2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
                     [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                     [0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1],
                     [0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1],
                     [-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0],
                     [9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2],
                     [-6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2],
                     [2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0],
                     [-6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1],
                     [4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1]])

    # Define x parameters
    # Function values
    x1 = n1
    x2 = n2
    x3 = n3
    x4 = n4
    # X Derivatives
    x5 = round((n2 - n16) / 2, 4)
    x6 = round((n9 - n1) / 2, 4)
    x7 = round((n10 - n4) / 2, 4)
    x8 = round((n3 - n15) / 2, 4)
    # Y Derivatives
    x9 = round((n4 - n6) / 2, 4)
    x10 = round((n3 - n7) / 2, 4)
    x11 = round((n12 - n2) / 2, 4)
    x12 = round((n13 - n1) / 2, 4)
    # Cross Derivatives (XY)
    x13 = round((n3 - n7 - n15 + n5) / 4, 4)
    x14 = round((n10 - n8 - n4 + n6) / 4, 4)
    x15 = round((n11 - n9 - n13 + n1) / 4, 4)
    x16 = round((n12 - n2 - n14 + n16) / 4, 4)

    # Create array from x parameters
    xarr = np.array([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16])

    # Multiply the inverse of the coefficient matrix by the array of x values to give array of alpha values
    alpha = np.matmul(cinv, xarr)

    # Calculate value at the point of interest
    n_p = 0
    for i in range(0, 4):
        for j in range(0, 4):
            n_p = n_p + alpha[i * 4 + j] * x ** i * y ** j

    # Return the value
    return n_p


# TEST OF SCRIPT
# Specify AUSGeoid2020 binary NTv2 file location
ntv2file = "C://Git/Python/NTv2.git/AUSGeoid2020_20180201.gsb"
# Assign class to variable
grids = NTv2ReaderBinary()
# Call ntv2reader to determine values at specific point
values = grids.ntv2reader(ntv2file, -26.948643, 145.6548721, 'bilinear')
# Print values
print(values)

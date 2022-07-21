"""
GeodePy - Python Geodesy Package
ntv2reader Module

Adapted from Jaimie Dodd's ntv2reader.py for AUSGeoid ntv2 files (2018-11-12).
Expanded (2021-03-21) for the more general case to include horizontal transformation grids
"""

import struct
from datetime import datetime as dt
import numpy as np


class NTv2Grid(object):
    def __init__(self, num_orec, num_srec, num_file, gs_type, version, system_f, system_t, major_f, minor_f, major_t,
                 minor_t, file_path):
        self.num_orec = num_orec        # Number of header identifiers
        self.num_srec = num_srec        # Number of sub-header idents
        self.num_file = num_file        # Number of subgrids in file
        self.gs_type = gs_type          # grid shift type
        self.version = version          # grid file version
        self.system_f = system_f        # reference system
        self.system_t = system_t        # reference system
        self.major_f = major_f          # semi major of from system
        self.minor_f = minor_f          # semi minor of from system
        self.major_t = major_t          # semi major of to system
        self.minor_t = minor_t          # semi minor of to system
        self.subgrids = {}
        self.file_path = file_path


class SubGrid(object):
    def __init__(self, sub_name, parent, created, updated, s_lat, n_lat, e_long, w_long, lat_inc, long_inc, gs_count):
        self.sub_name = sub_name        # subgrid name
        self.parent = parent            # parent name
        self.created = created          # date created
        self.updated = updated          # date modified
        self.s_lat = s_lat              # south latitude extents
        self.n_lat = n_lat              # north latitude extents
        self.e_long = e_long            # east longitude extents
        self.w_long = w_long            # west longitude extents
        self.lat_inc = lat_inc          # latitude increment
        self.long_inc = long_inc        # longitude increment
        self.gs_count = gs_count        # total nodes in subgrid

    def ntv2_bilinear(self, lat, lon, num_cols, row, col, f, start_byte):
        """
        Function to perform bicubic interpolation of the four ntv2 fields at a point of interest.

        :param lat: latitude of the point of interest (arc-seconds)
        :param lon: longitude of the point of interest (negative arc-seconds)
        :param num_cols: number of columns in grid
        :param row: row number of point to the bottom right of the point of interest
        :param col: column number of point to the bottom right of the point of interest
        :param f: variable for open file NTv2 being read as binary
        :param start_byte: start index of subgrid

        :return: four field tuple of interpolation results at point of interest.
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
        pos1 = row * num_cols + col
        pos2 = pos1 + 1
        pos3 = pos1 + num_cols
        pos4 = pos3 + 1

        # Navigate to start of subgrid
        f.seek(start_byte, 1)
        # Navigate to start of posA node
        f.seek(16 * pos1, 1)

        # Read in values for nodes 1 and 2
        node_1 = read_node(f)
        node_2 = read_node(f)

        # Navigate to beginning of node 3
        f.seek(16*(pos3 - pos2 - 1), 1)

        # Read in values for nodes 3 and 4
        node_3 = read_node(f)
        node_4 = read_node(f)

        # Determine latitude and longitude of node 1
        lat1 = self.s_lat + row * self.lat_inc
        long1 = self.e_long + col * self.long_inc

        # Determine interpolation scale factors
        x = (lon - long1) / self.long_inc
        y = (lat - lat1) / self.lat_inc

        field_1 = round(bilinear_interpolation(node_1[0], node_2[0], node_3[0], node_4[0], x, y), 6)
        field_2 = round(bilinear_interpolation(node_1[1], node_2[1], node_3[1], node_4[1], x, y), 6)
        field_3 = round(bilinear_interpolation(node_1[2], node_2[2], node_3[2], node_4[2], x, y), 6)
        field_4 = round(bilinear_interpolation(node_1[3], node_2[3], node_3[3], node_4[3], x, y), 6)

        return field_1, field_2, field_3, field_4

    def ntv2_bicubic(self, lat, lon, num_cols, row, col, f, start_byte):
        """
        Function to perform bicubic interpolation of the four ntv2 fields at a point of interest.

        :param lat: latitude of the point of interest (arc-seconds)
        :param lon: longitude of the point of interest (negative arc-seconds)
        :param num_cols: number of columns in grid
        :param row: row number of point to the bottom right of the point of interest
        :param col: column number of point to the bottom right of the point of interest
        :param f: variable for open file NTv2 being read as binary
        :param start_byte: start index of subgrid

        :return: four field tuple of interpolation results at point of interest.

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
        pos1 = row * num_cols + col
        pos2 = pos1 + 1
        pos3 = pos2 + num_cols
        pos4 = pos3 - 1
        pos5 = pos4 - 2 * num_cols - 1
        pos6 = pos5 + 1
        pos7 = pos6 + 1
        pos8 = pos7 + 1
        pos9 = pos8 + num_cols
        pos10 = pos9 + num_cols
        pos11 = pos10 + num_cols
        pos12 = pos11 - 1
        pos13 = pos12 - 1
        pos14 = pos13 - 1
        pos15 = pos14 - num_cols
        pos16 = pos15 - num_cols

        # Navigate to start of subgrid
        f.seek(start_byte, 1)
        # Navigate to start of pos1 node
        f.seek(16 * pos5, 1)

        # Read in values for nodes 5-8
        node_5 = read_node(f)
        node_6 = read_node(f)
        node_7 = read_node(f)
        node_8 = read_node(f)

        # Navigate to start of pos16 node
        f.seek(16 * (pos16 - pos8 - 1), 1)

        # Read in values for nodes 16, 1, 2, and 9
        node_16 = read_node(f)
        node_1 = read_node(f)
        node_2 = read_node(f)
        node_9 = read_node(f)

        # Navigate to start of pos15 node
        f.seek(16 * (pos15 - pos9 - 1), 1)

        # Read in values for nodes 15, 3, 4 and 10
        node_15 = read_node(f)
        node_4 = read_node(f)
        node_3 = read_node(f)
        node_10 = read_node(f)

        # Navigate to start of pos14 node
        f.seek(16 * (pos14 - pos10 - 1), 1)

        # Read in values for nodes 11, 12, 13 and 14
        node_14 = read_node(f)
        node_13 = read_node(f)
        node_12 = read_node(f)
        node_11 = read_node(f)

        # Determine latitude and longitude of node 1
        lat1 = self.s_lat + row * self.lat_inc
        long1 = self.e_long + col * self.long_inc

        # Determine interpolation scale factors
        x = round((lon - long1) / self.long_inc, 6)
        y = round((lat - lat1) / self.lat_inc, 6)

        # Call bicubic interpolation function to interpolate ntv2 fields.
        field_1 = round(bicubic_interpolation(node_1[0], node_2[0], node_3[0], node_4[0],
                                              node_5[0], node_6[0], node_7[0], node_8[0],
                                              node_9[0], node_10[0], node_11[0], node_12[0],
                                              node_13[0], node_14[0], node_15[0], node_16[0],
                                              x, y), 6)
        field_2 = round(bicubic_interpolation(node_1[1], node_2[1], node_3[1], node_4[1],
                                              node_5[1], node_6[1], node_7[1], node_8[1],
                                              node_9[1], node_10[1], node_11[1], node_12[1],
                                              node_13[1], node_14[1], node_15[1], node_16[1],
                                              x, y), 6)
        field_3 = round(bicubic_interpolation(node_1[2], node_2[2], node_3[2], node_4[2],
                                              node_5[2], node_6[2], node_7[2], node_8[2],
                                              node_9[2], node_10[2], node_11[2], node_12[2],
                                              node_13[2], node_14[2], node_15[2], node_16[2],
                                              x, y), 6)
        field_4 = round(bicubic_interpolation(node_1[3], node_2[3], node_3[3], node_4[3],
                                              node_5[3], node_6[3], node_7[3], node_8[3],
                                              node_9[3], node_10[3], node_11[3], node_12[3],
                                              node_13[3], node_14[3], node_15[3], node_16[3],
                                              x, y), 6)

        return field_1, field_2, field_3, field_4


def read_node(f):
    """
    Function to read in values of nodes

    :param f: variable for open file NTv2 being read as binary
    :return: tuple containing the four ntv2 fields at the grid node.
    """

    # field_1: shift lat / geoid separation (N)
    byte = f.read(4)
    field_1 = struct.unpack('f', byte)[0]
    # field 2: shift lon / deflection in prime meridian value (Xi)
    byte = f.read(4)
    field_2 = struct.unpack('f', byte)[0]
    # field 3: reliability lat / deflection in prime vertical value (Eta)
    byte = f.read(4)
    field_3 = struct.unpack('f', byte)[0]
    # field 4: reliability lon / NA
    byte = f.read(4)
    field_4 = struct.unpack('f', byte)[0]

    return field_1, field_2, field_3, field_4


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
    a1 = n2 - n1
    a2 = n3 - n1
    a3 = n1 + n4 - n2 - n3
    p = a0 + (a1 * x) + (a2 * y) + (a3 * x * y)
    return p


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
    x5 = (n2 - n16) / 2
    x6 = (n9 - n1) / 2
    x7 = (n10 - n4) / 2
    x8 = (n3 - n15) / 2
    # Y Derivatives
    x9 = (n4 - n6) / 2
    x10 = (n3 - n7) / 2
    x11 = (n12 - n2) / 2
    x12 = (n13 - n1) / 2
    # Cross Derivatives (XY)
    x13 = (n3 - n7 - n15 + n5) / 4
    x14 = (n10 - n8 - n4 + n6) / 4
    x15 = (n11 - n9 - n13 + n1) / 4
    x16 = (n12 - n2 - n14 + n16) / 4

    # Create array from x parameters
    xarr = np.array([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16])

    # Multiply the inverse of the coefficient matrix by the array of x values to give array of alpha values
    alpha = np.matmul(cinv, xarr)

    # Calculate value at the point of interest
    n_p = 0
    for i in range(0, 4):
        for j in range(0, 4):
            n_p = n_p + alpha[i * 4 + j] * x ** i * y ** j

    return n_p


def read_ntv2_file(ntv2_gsb_file):
    """
    Function to read an ntv2 gsb file and create grid & subgrid objects

    :param ntv2_gsb_file: full path to ntv2 file
    :return: ntv2 grid object
    """

    with open(ntv2_gsb_file, 'rb') as f:
        # NUM_OREC
        f.seek(8, 1)
        byte = f.read(4)
        num_orec = int.from_bytes(byte, byteorder='little')

        f.seek(4, 1)

        # NUM_SREC
        f.seek(8, 1)
        byte = f.read(4)
        num_srec = int.from_bytes(byte, byteorder='little')

        f.seek(4, 1)
        # NUM_FILE
        f.seek(8, 1)
        byte = f.read(4)
        num_file = int.from_bytes(byte, byteorder='little')

        f.seek(4, 1)
        # GS_TYPE
        f.seek(8, 1)
        byte = f.read(8)
        gs_type = byte.decode('utf8').strip('\x00').strip()

        # VERSION
        f.seek(8, 1)
        byte = f.read(8)
        version = byte.decode('utf8').strip('\x00').strip()

        # SYSTEM_F
        f.seek(8, 1)
        byte = f.read(8)
        system_f = byte.decode('utf8').strip('\x00').strip()

        # SYSTEM_T
        f.seek(8, 1)
        byte = f.read(8)
        system_t = byte.decode('utf8').strip('\x00').strip()

        # MAJOR_F
        f.seek(8, 1)
        byte = f.read(8)
        major_f = struct.unpack('d', byte)[0]

        # MINOR_F
        f.seek(8, 1)
        byte = f.read(8)
        minor_f = struct.unpack('d', byte)[0]

        # MAJOR_T
        f.seek(8, 1)
        byte = f.read(8)
        major_t = struct.unpack('d', byte)[0]

        # MINOR_T
        f.seek(8, 1)
        byte = f.read(8)
        minor_t = struct.unpack('d', byte)[0]

        grid = NTv2Grid(
            num_orec=num_orec,
            num_srec=num_srec,
            num_file=num_file,
            gs_type=gs_type,
            version=version,
            system_f=system_f,
            system_t=system_t,
            major_f=major_f,
            minor_f=minor_f,
            major_t=major_t,
            minor_t=minor_t,
            file_path=ntv2_gsb_file
        )

        # read subgrids
        for i in range(0, grid.num_file):
            # SUB_NAME
            f.seek(8, 1)
            byte = f.read(8)
            sub_name = byte.decode('utf').strip('\x00').strip()

            # PARENT
            f.seek(8, 1)
            byte = f.read(8)
            parent = byte.decode('utf').strip('\x00').strip()

            # CREATED
            f.seek(8, 1)
            byte = f.read(8)
            created = dt.strptime(byte.decode('utf').strip('\x00').strip(), '%d%m%Y').strftime('%d/%m/%Y')

            # UPDATED
            f.seek(8, 1)
            byte = f.read(8)
            updated = dt.strptime(byte.decode('utf').strip('\x00').strip(), '%d%m%Y').strftime('%d/%m/%Y')

            # S_LAT
            f.seek(8, 1)
            byte = f.read(8)
            s_lat = round(struct.unpack('d', byte)[0], 3)

            # N_LAT
            f.seek(8, 1)
            byte = f.read(8)
            n_lat = round(struct.unpack('d', byte)[0], 3)

            # E_LONG
            f.seek(8, 1)
            byte = f.read(8)
            e_long = round(struct.unpack('d', byte)[0], 3)

            # W_LONG
            f.seek(8, 1)
            byte = f.read(8)
            w_long = round(struct.unpack('d', byte)[0], 3)

            # LAT_INC
            f.seek(8, 1)
            byte = f.read(8)
            lat_inc = round(struct.unpack('d', byte)[0], 6)

            # LONG_INC
            f.seek(8, 1)
            byte = f.read(8)
            long_inc = round(struct.unpack('d', byte)[0], 6)

            # GS_COUNT
            f.seek(8, 1)
            byte = f.read(4)
            gs_count = int.from_bytes(byte, byteorder='little')

            f.seek(4, 1)
            # skip ahead to next subgrid
            f.seek(16 * gs_count, 1)

            grid.subgrids[sub_name] = SubGrid(
                sub_name=sub_name,
                parent=parent,
                created=created,
                updated=updated,
                s_lat=s_lat,
                n_lat=n_lat,
                e_long=e_long,
                w_long=w_long,
                lat_inc=lat_inc,
                long_inc=long_inc,
                gs_count=gs_count
            )

        return grid


def interpolate_ntv2(grid_object, lat, lon, method='bicubic'):
    """
    Function to interpolate Ntv2Grid objects
    :param grid_object: Ntv2Grid object
    :param lat: latitude (decimal degrees)
    :param lon: longitude (decimal degrees)
    :param method: interpolation strategy, bicubic or bilinear

    :return: tuple of four ntv2 fields
    """

    interpolation_methods = {
        'bicubic',
        'bilinear'
    }
    if method not in interpolation_methods:
        raise ValueError(f'interpolation method "{method}" not supported')

    # convert decimal degrees to arc-seconds
    lat *= 3600
    lon *= -3600

    # determine subgrid for point of interest
    in_subgrids = set()
    for sg in grid_object.subgrids.values():
        if sg.s_lat <= lat < sg.n_lat and sg.e_long <= lon < sg.w_long:
            in_subgrids.add(sg.sub_name)

    # return null fields if no subgrid found, else choose subgrid with finest resolution
    if len(in_subgrids) == 0:
        return None, None, None, None
    else:
        inc = None
        in_grid = None
        for sg in in_subgrids:
            if not inc:
                inc = grid_object.subgrids[sg].lat_inc
                in_grid = grid_object.subgrids[sg]
            else:
                if grid_object.subgrids[sg].lat_inc < inc:
                    inc = grid_object.subgrids[sg].lat_inc
                    in_grid = grid_object.subgrids[sg]

    # Determine number of columns in grid, and row and column of node to bottom right of
    # point of interest, then call relevant interpolation method function

    # determine number of columns
    num_cols = 1 + int((in_grid.w_long - in_grid.e_long) / in_grid.long_inc)

    # determine row and col numbers of node below right of point
    row = int((lat - in_grid.s_lat) / in_grid.lat_inc)
    col = int((lon - in_grid.e_long) / in_grid.long_inc)

    # locate data in gsb_file
    skip_bytes = 176    # grid header length

    with open(grid_object.file_path, 'rb') as f:
        for sg in grid_object.subgrids.values():

            skip_bytes += 176   # subgrid header length
            if sg.sub_name == in_grid.sub_name:
                if method == 'bilinear':
                    results = in_grid.ntv2_bilinear(lat, lon, num_cols, row, col, f, skip_bytes)
                elif method == 'bicubic':
                    results = in_grid.ntv2_bicubic(lat, lon, num_cols, row, col, f, skip_bytes)
            else:
                skip_bytes += sg.gs_count * 16

    return results

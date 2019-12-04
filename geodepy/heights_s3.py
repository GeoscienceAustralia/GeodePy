#___________________________________________________________________________#
# Some notes:
# Written by Jack McCubbine of Geoscience Australia, date: 08/11/2019
# This code contains functions to handle tranformations between GPS and 
# AWVS/AHD and Vice Versa 
# Gridded data used for the varisous reference surfaces are geotiff files
# These allow direct access remotely using "gdal"
#___________________________________________________________________________#
# Import dependencies 
import gdal
import numpy as np
from scipy.interpolate import griddata
import geodepy.Height_filenames as Height_filenames
#___________________________________________________________________________#
# Interpolation functions
def interp_file(Lat,Long,file):
    # Import the DOVPM file
    f = gdal.Open(file)
    # load band (akin to a variable in dataset)
    band = f.GetRasterBand(1)   
    # get the pixel width, height, etc.
    transform = f.GetGeoTransform()
    # Grid resolution (known)
    res=transform[1]
    # convert lat,lon to row,col
    column = (Long - transform[0]) / transform[1]
    row = (Lat - transform[3]) / transform[5]
    # get pixel values surrounding data point
    Surrounding_data=(band.ReadAsArray(np.floor(column-2), np.floor(row-2), 5, 5))
    # convert row,col back to north,east
    Long_c = transform[0] + np.floor(column) * res
    Lat_c = transform[3] - np.floor(row) * res
    # set up matrices for interpolation
    count=-1
    pos=np.zeros((25,2))
    Surrounding_data_v=np.zeros((25,1))
    for k in range(-2,3):
        for j in range(-2,3):
            count=count+1
            pos[count]=(Long_c+j*res,Lat_c-k*res)
            Surrounding_data_v[count]=Surrounding_data[k+2,j+2]          
    interp_val=griddata(pos,Surrounding_data_v,(Long,Lat),method='cubic')
    return interp_val

#___________________________________________________________________________#
# Functions to handle the conversions from one height to another
def GPS_to_AVWS(Lat,Long,GPS_H):
    zeta=interp_file(Lat,Long,Height_filenames.file_AVWS) # AVWS file
    zeta_std=interp_file(Lat,Long,Height_filenames.file_AVWS_STD) # AVWS STD file
    NORMAL_H=GPS_H-zeta
    return [NORMAL_H,zeta_std]

def AVWS_to_GPS(Lat,Long,AVWS_H):
    zeta=interp_file(Lat,Long,Height_filenames.file_AVWS) # AVWS file
    zeta_std=interp_file(Lat,Long,Height_filenames.file_AVWS_STD) # AVWS STD file
    GPS_H=AVWS_H+zeta
    return [GPS_H,zeta_std]

def AHD_to_AVWS(Lat,Long,AHD_H,file_AG2020):
    # Convert to GPS
    GPS_H=AHD_H+interp_file(Lat,Long,Height_filenames.file_AG2020) # AUSGEOID2020 file
    # Convert to AVWS
    Normal_H=GPS_H-interp_file(Lat,Long,Height_filenames.file_AVWS) # AVWS file
    return [Normal_H]

def GPS_to_AHD(Lat,Long,GPS_H):
    N=interp_file(Lat,Long,Height_filenames.file_AG2020) # AUSGEOID2020 file
    N_std=interp_file(Lat,Long,Height_filenames.file_AG2020_STD) # AUSGEOID2020 STD file
    AHD_H=GPS_H-N
    return [AHD_H,N_std]

def AHD_to_GPS(Lat,Long,AHD_H):
    N=interp_file(Lat,Long,Height_filenames.file_AG2020) # AUSGEOID2020 file
    N_std=interp_file(Lat,Long,Height_filenames.file_AG2020_STD) # AUSGEOID2020 STD file
    GPS_H=AHD_H+N
    return [GPS_H,N_std]

def AVWS_to_AHD(Lat,Long,Normal_H):   
    # Convert to GPS    
    GPS_H=Normal_H+interp_file(Lat,Long,Height_filenames.file_AVWS) # AVWS file
    # Convert to AHD
    AHD_H=GPS_H-interp_file(Lat,Long,Height_filenames.file_AG2020) # AUSGEOID2020 file
    return [AHD_H]

def DOV(Lat,Long):   
    # Convert to GPS    
    DOV_PM=interp_file(Lat,Long,Height_filenames.file_DOV_PM) # AVWS file
    # Convert to AHD
    DOV_PV=interp_file(Lat,Long,Height_filenames.file_DOV_PV) # AUSGEOID2020 file
    return [DOV_PM,DOV_PV]

#___________________________________________________________________________#
# Some notes:
# Written by Jack McCubbine of Geoscience Australia, date: 08/11/2019
# This code contains functions to handle tranformations between GPS and 
# AWVS/AHD and Vice Versa 
# Gridded data used for the varisous reference surfaces are geotiff files
# These allow direct access remotely using "gdal"
# File names are hard coded 
# Some exmaples of how the functions can be used are give towards the end 
#___________________________________________________________________________#
# Import dependencies 
import gdal
import numpy as np
from scipy.interpolate import griddata
#___________________________________________________________________________#
# Interpolation functions
def interp_DOVPV(Lat,Long):
    # Import the DOVPM file
    f = gdal.Open('DOV_PV.tif')
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
    interp_DOVPV_val=griddata(pos,Surrounding_data_v,(Long,Lat),method='cubic')
    return interp_DOVPV_val

def interp_DOVPM(Lat,Long):
    # Import the DOVPM file
    f = gdal.Open('DOV_PM.tif')
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
    interp_DOVPM_val=griddata(pos,Surrounding_data_v,(Long,Lat),method='cubic')
    return interp_DOVPM_val

def interp_AUSGeoid2020(Lat,Long):
    # open geotiff file
    f = gdal.Open('AUSGeoid2020_RELEASEV20170908.tif')
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
    interp_AUSGeoid2020_val=griddata(pos,Surrounding_data_v,(Long,Lat),method='cubic')
    return interp_AUSGeoid2020_val

def interp_AUSGeoid2020_STD(Lat,Long):
    # open geotiff file
    f = gdal.Open('AUSGeoid2020_RELEASEV20170908_err.tif')
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
    interp_AUSGeoid2020_STD_val=griddata(pos,Surrounding_data_v,(Long,Lat),method='cubic')
    return interp_AUSGeoid2020_STD_val

def interp_AVWS(Lat,Long):
    # open geotiff file
    f = gdal.Open('AVWS_20191107.tif')
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
    interp_AVWS_val=griddata(pos,Surrounding_data_v,(Long,Lat),method='cubic')
    return interp_AVWS_val

def interp_AVWS_STD(Lat,Long):
    # open geotiff file
    f = gdal.Open('AVWS_STD_20191107.tif')
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
    interp_AVWS_val=griddata(pos,Surrounding_data_v,(Long,Lat),method='cubic')
    return interp_AVWS_val

#___________________________________________________________________________#
# Functions to handle the conversions from one height to another
def GPS_to_AVWS(Lat,Long,GPS_H):
    zeta=interp_AVWS(Lat,Long)
    zeta_std=interp_AVWS_STD(Lat,Long)
    NORMAL_H=GPS_H-zeta
    return [NORMAL_H,zeta_std]

def AVWS_to_GPS(Lat,Long,AVWS_H):
    zeta=interp_AVWS(Lat,Long)
    zeta_std=interp_AVWS_STD(Lat,Long)
    GPS_H=AVWS_H+zeta
    return [GPS_H,zeta_std]

def AHD_to_AVWS(Lat,Long,AHD_H):
    # Convert to GPS
    GPS_H=AHD_H+interp_AUSGeoid2020(Lat,Long) 
    # Convert to AVWS
    Normal_H=GPS_H-interp_AVWS(Lat,Long)
    return [Normal_H]

def GPS_to_AHD(Lat,Long,GPS_H):
    N=interp_AUSGeoid2020(Lat,Long)
    N_std=interp_AUSGeoid2020_STD(Lat,Long)
    AHD_H=GPS_H-N
    return [AHD_H,N_std]

def AHD_to_GPS(Lat,Long,AHD_H):
    N=interp_AUSGeoid2020(Lat,Long)
    N_std=interp_AUSGeoid2020_STD(Lat,Long)
    GPS_H=AHD_H+N
    return [GPS_H,N_std]

def AVWS_to_AHD(Lat,Long,Normal_H):   
    # Convert to GPS    
    GPS_H=Normal_H+interp_AVWS(Lat,Long)
    # Convert to AHD
    AHD_H=GPS_H-interp_AUSGeoid2020(Lat,Long)     
    return [AHD_H]





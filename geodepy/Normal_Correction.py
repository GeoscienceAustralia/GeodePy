# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import gdal
import numpy as np
import geodepy.Height_filenames as Height_filenames
from scipy.interpolate import griddata

def mean_normal_grav(Lat,h):
    # GRS 80 constants
    a=6378137
    b=6356752.3141
    omega=7292115*(10**-11)
    e2=0.00669438002290 
    GM=3986005*10**8
    k=0.001931851353
    # GRS80 normal gravity
    EllGrav=(10**5)*9.7803267715*(1+k*(np.sin(Lat*np.pi/180)**2))/np.sqrt(1-e2*(np.sin(Lat*np.pi/180)**2))
    FA=-((2*(EllGrav/a)*(1+(a-b)/a + omega**2*a**2*b/GM - 2*(a-b)/a*(np.sin(Lat*np.pi/180)**2))*(h**2)/2-3*(EllGrav/a**2)*(h**3)/3)/h)
    mean_normal_g=(EllGrav+FA)*(10**-5)
    return mean_normal_g

def normal_grav(Lat,h):
    # GRS 80 constants
    a=6378137
    b=6356752.3141
    omega=7292115*(10**-11)
    e2=0.00669438002290 
    GM=3986005*10**8
    k=0.001931851353
    # GRS80 normal gravity
    EllGrav=(10**5)*9.7803267715*(1+k*(np.sin(Lat*np.pi/180)**2))/np.sqrt(1-e2*(np.sin(Lat*np.pi/180)**2))
    FA=-(2*EllGrav*h/a)*(1+(a-b)/a+omega**2*a**2*b/GM-2*(a-b)/a*(np.sin(Lat*np.pi/180)**2))+3*(EllGrav*h**2)/(a**2)
    normal_g=(EllGrav+FA)*(10**-5)
    return normal_g

def mean_surface_grav(Lat_A,Long_A,H_A,Lat_B,Long_B,H_B):    
    Surf_Grav_A=interp_grav(Lat_A,Long_A)*(10**-5)+normal_grav(Lat_A,H_A)+0.0419*2.67*H_A*(10**-5)
    Surf_Grav_B=interp_grav(Lat_B,Long_B)*(10**-5)+normal_grav(Lat_B,H_B)+0.0419*2.67*H_B*(10**-5)
    mean_g=(Surf_Grav_A+Surf_Grav_B)/2
    return mean_g

def interp_grav(Lat,Long):
    # Grid resolution (known)
    res=1.0/60 
    # open geotiff file
    f = gdal.Open(Height_filenames.file_GRAV_BA)
    # load band (akin to a variable in dataset)
    band = f.GetRasterBand(1)   
    # get the pixel width, height, etc.
    transform = f.GetGeoTransform()
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
    interp_g=griddata(pos,Surrounding_data_v,(Long,Lat))
    return interp_g

def normal_correction(Lat_A,Long_A,H_A,Lat_B,Long_B,H_B):
    # ellipsoidal gravity at 45 deg. Lat    
    Gamma_0=9.8061992115
    # Normal Gravity at the point
    normal_g_A=mean_normal_grav(Lat_A,H_A)
#    print normal_g_A
    normal_g_B=mean_normal_grav(Lat_B,H_B)
#    print normal_g_B
    dn=H_B-H_A       
    g=mean_surface_grav(Lat_A,Long_A,H_A,Lat_B,Long_B,H_B)     
   # print g
    NC=(dn*(g-Gamma_0)/Gamma_0)+H_A*(normal_g_A-Gamma_0)/Gamma_0-H_B*(normal_g_B-Gamma_0)/Gamma_0
    return NC,g






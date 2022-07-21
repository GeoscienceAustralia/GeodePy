#___________________________________________________________________________#
# Some notes:
# Written by Jack McCubbine of Geoscience Australia, date: 08/11/2019
# This code contains functions to handle tranformations between GPS and 
# AWVS/AHD and Vice Versa 
# Gridded data used for the varisous reference surfaces are geotiff files
# These allow direct access remotely using "gdal"
#___________________________________________________________________________#
# Import dependencies 
import geodepy.constants as cons
import geodepy.geodesy as gg
import gdal
import numpy as np
from scipy.interpolate import griddata
import math as m
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
    zeta=interp_file(Lat, Long, cons.file_AVWS) # AVWS file
    zeta_std=interp_file(Lat, Long, cons.file_AVWS_STD) # AVWS STD file
    NORMAL_H=GPS_H-zeta
    return [NORMAL_H,zeta_std]

def AVWS_to_GPS(Lat,Long,AVWS_H):
    zeta=interp_file(Lat, Long, cons.file_AVWS) # AVWS file
    zeta_std=interp_file(Lat, Long, cons.file_AVWS_STD) # AVWS STD file
    GPS_H=AVWS_H+zeta
    return [GPS_H,zeta_std]

def AHD_to_AVWS(Lat,Long,AHD_H):
    # Convert to GPS
    GPS_H=AHD_H+interp_file(Lat, Long, cons.file_AG2020) # AUSGEOID2020 file
    # Convert to AVWS
    Normal_H=GPS_H-interp_file(Lat, Long, cons.file_AVWS) # AVWS file
    return [Normal_H]

def GPS_to_AHD(Lat,Long,GPS_H):
    N=interp_file(Lat, Long, cons.file_AG2020) # AUSGEOID2020 file
    N_std=interp_file(Lat, Long, cons.file_AG2020_STD) # AUSGEOID2020 STD file
    AHD_H=GPS_H-N
    return [AHD_H,N_std]

def AHD_to_GPS(Lat,Long,AHD_H):
    N=interp_file(Lat, Long, cons.file_AG2020) # AUSGEOID2020 file
    N_std=interp_file(Lat, Long, cons.file_AG2020_STD) # AUSGEOID2020 STD file
    GPS_H=AHD_H+N
    return [GPS_H,N_std]

def AVWS_to_AHD(Lat,Long,Normal_H):   
    # Convert to GPS    
    GPS_H=Normal_H+interp_file(Lat, Long, cons.file_AVWS) # AVWS file
    # Convert to AHD
    AHD_H=GPS_H-interp_file(Lat, Long, cons.file_AG2020) # AUSGEOID2020 file
    return [AHD_H]

def DOV(Lat,Long):   
    # Convert to GPS    
    DOV_PM=interp_file(Lat, Long, cons.file_DOV_PM) # AVWS file
    # Convert to AHD
    DOV_PV=interp_file(Lat, Long, cons.file_DOV_PV) # AUSGEOID2020 file
    return [DOV_PM,DOV_PV]

def GPS_to_AUSGeoid98(Lat,Long,GPS_H):
    N=interp_file(Lat,Long,cons.file_AG98) # AUSGEOID98 file
    AHD_H=GPS_H-N
    return [AHD_H]

def AUSGeoid98_to_GPS(Lat,Long,AHD_H):
    N=interp_file(Lat,Long,cons.file_AG98) # AUSGEOID98 file
    GPS_H=AHD_H+N
    return [GPS_H]

def GPS_to_AUSGeoid09(Lat,Long,GPS_H):
    N=interp_file(Lat,Long,cons.file_AG09) # AUSGEOID09 file
    AHD_H=GPS_H-N
    return [AHD_H]

def AUSGeoid09_to_GPS(Lat,Long,AHD_H):
    N=interp_file(Lat,Long,cons.file_AG09) # AUSGEOID09 file
    GPS_H=AHD_H+N
    return [GPS_H]

def DOV_09(Lat,Long):   
    # Interp PM
    DOV_PM=interp_file(Lat,Long,cons.file_AG09_DOV_PM) # AGQG09 DOV file
    # Interp PV
    DOV_PV=interp_file(Lat,Long,cons.file_AG09_DOV_PV) # AGQG09 DOV file
    return [DOV_PM,DOV_PV]

def DOV_98(Lat,Long):   
    # Interp PM  
    DOV_PM=interp_file(Lat,Long,cons.file_AG98_DOV_PM) # AGQG98 DOV file
    # Interp PV
    DOV_PV=interp_file(Lat,Long,cons.file_AG98_DOV_PV) # AGQG98 DOV file
    return [DOV_PM,DOV_PV]

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
    f = gdal.Open(cons.file_GRAV_BA)
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


def normal_orthometric_correction(lat1, lon1, H1, lat2, lon2, H2):
    """
    Computes the normal-orthometric correction based on Heck (2003).
    See Standard for New Zealand Vertical Datum 2016, Section 3.3

    :param lat1: Latitude at Stn1
    :param lon1: Longitude at Stn1
    :param H1: Physical Height at Stn1
    :param lat2: Latitude at Stn2
    :param lon2: longitude at Stn2
    :param H2: Physical Height at Stn2
    :return: normal-orthometric correction
    """

    f_ng = cons.grs80_ngf
    m_rad = cons.grs80.meanradius

    mid_height = (H1 + H2) / 2
    mid_lat = m.radians((lat1 + lat2) / 2)

    vinc_inv = gg.vincinv(lat1, lon1, lat2, lon2)
    dist = vinc_inv[0]
    az = vinc_inv[1]

    noc = - f_ng / m_rad * mid_height * m.sin(2.0 * mid_lat) * m.cos(m.radians(az)) * dist

    return noc

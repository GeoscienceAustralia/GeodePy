#######################################################
#                GeoDMS2Grid_Redfearn.py              #
#   Josh Batchelor - Geoscience Australia - Nov 2017  #
#                                                     #
#       Script to convert Latitude and Longitude      #
#     coordinates in Degrees, Minutes and Seconds     #
#    into Universal Transverse Mercator Projection    #
#    coordinates. Reads in a comma separated text     #
#    file of the format:                              #
#      Pt,Latitude,Longitude                          #
#    and outputs a comma separated text file of       #
#    the format:                                      #
#      Pt,Zone,Easting,Northing                       #
#                                                     #
#######################################################
from decimal import *
from math import *
import os
import pprint
def GeoDMS2Grid(lat,long):
        # Universal Transverse Mercator Projection Parameters
	Proj = [6378137,Decimal('298.257222101'),500000,10000000,Decimal('0.9996'),6,-177]
	# First Eccentricity
	e2 = (2*(1/Proj[1]))-((1/Proj[1])**2)
	e4 = e2 ** 2
	e6 = e2 * e4
	# Function to convert Degrees Minutes and Seconds into Decimal Degrees
	# Takes DMS as DDD.MMSSSSSSSSSSS where seconds are decimal
	def DMS2dd(DMSraw):
		getcontext().prec = 28
		DMSseg = [int(DMSraw)]
		DMSseg = DMSseg + [Decimal(str(DMSraw - DMSseg[0])).quantize(Decimal('.01'),rounding=ROUND_DOWN)]
		DMSseg = DMSseg + [Decimal(str(DMSraw)) - DMSseg[0] - DMSseg[1]]
		decDeg = DMSseg[0] + (DMSseg[1] / Decimal('0.6')) + (DMSseg[2] / Decimal('0.36'))
		return(decDeg)
	# Convert Lat & Long to Decimal Degrees
	lat = DMS2dd(lat)
	long = DMS2dd(long)
	# Calculate Zone
	zone = int((float(long) - (Proj[6]-(1.5*Proj[5])))/Proj[5])
	CM = float(zone*Proj[5])+(Proj[6]-Proj[5])
	# Calculate Functions
	sin1Lat = sin(radians(lat))
	sin2Lat = sin(2*radians(lat))
	sin4Lat = sin(4*radians(lat))
	sin6Lat = sin(6*radians(lat))
	A0 = 1-(e2/4)-((3*e4)/64)-((5*e6)/256)
	A2 = Decimal(str(float(3)/8))*(e2+(e4/4)+((15*e6)/128))
	A4 = Decimal(str(float(15)/256))*(e4+((3*e6)/4))
	A6 = (35*e6)/3072
	# Meridian Distance
	MerDist1 = (Proj[0]*float(A0)*(radians(lat)))
	MerDist2 = (Proj[0]*float(A2)*sin2Lat)
	MerDist3 = (Proj[0]*float(A4)*sin4Lat)
	MerDist4 = (Proj[0]*float(A6)*sin6Lat)
	MerDist = MerDist1 - MerDist2 + MerDist3 - MerDist4
	# Radii of Curvature
	Rho = float(Proj[0]*(1-e2))/(pow((1-(Decimal(str(e2))*Decimal(str(sin1Lat**2)))),Decimal(str(1.5))))
	Nu = Proj[0]/(pow((1-(Decimal(str(e2))*Decimal(str(sin1Lat**2)))),Decimal(str(0.5))))
	Psi = Nu/Rho
	Psi2 = Psi**2
	Psi3 = Psi*Psi2
	Psi4 = Psi2**2
	# Powers of Latitude
	cos1Lat = cos(radians(lat))
	cos2Lat = cos1Lat**2
	cos3Lat = cos1Lat**3
	cos4Lat = cos1Lat**4
	cos5Lat = cos1Lat**5
	cos6Lat = cos1Lat**6
	cos7Lat = cos1Lat**7
	cos8Lat = cos1Lat**8
	# Powers of Difference in Longitude
	diff1Long = radians(long - Decimal(str(CM)))
	diff2Long = diff1Long**2
	diff3Long = diff1Long**3
	diff4Long = diff1Long**4
	diff5Long = diff1Long**5
	diff6Long = diff1Long**6
	diff7Long = diff1Long**7
	diff8Long = diff1Long**8
	# Powers of Tan Latitude
	tan1Lat = tan(radians(lat))
	tan2Lat = tan1Lat**2
	tan4Lat = tan1Lat**4
	tan6Lat = tan1Lat**6
	# Calculating Easting
	east1Term = Nu*diff1Long*cos1Lat
	east2Term = Nu*diff3Long*cos3Lat*((Psi-tan2Lat)/6)
	east3Term = Nu*diff5Long*cos5Lat*((4*Psi3*(1-6*tan2Lat)+Psi2*(1+8*tan2Lat)-Psi*(2*tan2Lat)+tan4Lat)/120)
	east4Term = Nu*diff7Long*cos7Lat*((61-479*tan2Lat+179*tan4Lat-tan6Lat)/5040)
	East = Proj[2]+Decimal(str(east1Term + east2Term + east3Term + east4Term))*Proj[4]
	# Calculating Northing
	north1Term = Nu*sin1Lat*diff2Long*cos1Lat/2
	north2Term = Nu*sin1Lat*diff4Long*cos3Lat*((4*Psi2+Psi-tan2Lat)/24)
	north3Term = Nu*sin1Lat*diff6Long*cos5Lat*((8*Psi4*(11-24*tan2Lat)-28*Psi3*(1-6*tan2Lat)+Psi2*(1-32*tan2Lat)-Psi*(2*tan2Lat+tan4Lat))/720)
	north4Term = Nu*sin1Lat*diff8Long*cos7Lat*((1385-3111*tan2Lat+543*tan4Lat-tan6Lat)/40320)
	North = Proj[3]+Decimal(str(MerDist + north1Term + north2Term + north3Term + north4Term))*Proj[4]
	####     Output    ####
	return(',' + str(zone) + ',' + str(round(East,3)) + ',' + str(round(North,3)))
####     Enter Filename    ####
print('Enter co-ordinate file:')
fn = input()
####     Open Filename     ####
dataFile = open(fn, 'r')
####   Create Output File  ####
fn_part = (os.path.splitext(fn))
fn_out = fn_part[0] + '_out' + fn_part[1]
outFile = open(fn_out, 'w')
####      Write Output     ####
outFile.write('Pt,Zone,Easting,Northing\n')
for line in dataFile:
    Pt = line.rstrip()
    Pt = Pt.split(",")
    PtNum = Pt[0]
    lat = float(Pt[1])
    long = float(Pt[2])
    output = GeoDMS2Grid(lat, long)
    outFile.write(Pt[0] + output + '\n')
####      Close Files      ####
outFile.close()
dataFile.close()

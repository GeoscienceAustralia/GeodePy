######################################################################################################
#                                          GeoDMS2Grid.py                                            #
#                         Josh Batchelor - Geoscience Australia - Nov 2017                           #
#                                                                                                    #
#       Script to convert Latitude and Longitude coordinates in Degrees, Minutes and Seconds         #
#       into Universal Transverse Mercator Projection coordinates. Reads in a comma separated        #
#       text file of the format:                                                                     #
#        Pt,Latitude,Longitude                                                                       #
#       and outputs a comma separated text file of the format:                                       #
#        Pt,Zone,Easting,Northing                                                                    #
#                                                                                                    #
#       Ref: http://www.icsm.gov.au/gda/tech.html                                                    #
#       Ref: http://www.mygeodesy.id.au/documents/Karney-Krueger%20equations.pdf                     #
#                                                                                                    #
######################################################################################################
from decimal import *
from math import *
import os
import pprint
getcontext().prec = 28
# Universal Transverse Mercator Projection Parameters
Proj = [6378137,Decimal('298.257222101'),500000,10000000,Decimal('0.9996'),6,-177]
# Ellipsoidal Constants
f = 1/Proj[1]
e2 = f*(2-f)
ecc1 = sqrt(e2)
n = f/(2-f)
n = float(n)
n2 = n**2
# Rectifying Radius (Horner Form)
A = Proj[0]/(1+n)*((n2*(n2*(n2*(25*n2+64)+256)+4096)+16384)/16384.)
# Alpha Coefficients (Horner Form)
a2 = (n*(n*(n*(n*(n*(n*((37884525-75900428*n)*n+42422016)-89611200)+46287360)+63504000)-135475200)+101606400))/203212800.
a4 = (n2*(n*(n*(n*(n*(n*(148003883*n+83274912)-178508970)+77690880)+67374720)-104509440)+47174400))/174182400.
a6 = (n**3*(n*(n*(n*(n*(318729724*n-738126169)+294981280)+178924680)-234938880)+81164160))/319334400.
a8 = (n**4*(n*(n*((14967552000-40176129013*n)*n+6971354016)-8165836800)+2355138720))/7664025600.
a10 = (n**5*(n*(n*(10421654396*n+3997835751)-4266773472)+1072709352))/2490808320.
a12 = (n**6*(n*(175214326799*n-171950693600)+38652967262))/58118860800.
a14 = ((13700311101-67039739596*n)*n**7)/12454041600.
a16 = (1424729850961*n**8)/743921418240.
def GeoDMS2Grid(lat,long):
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
	# Conformal Latitude
	sigx = (ecc1*tan(radians(lat)))/sqrt(1+(tan(radians(lat))**2))
	sig = sinh(ecc1*(0.5*log((1+sigx)/(1-sigx))))
	conf_lat = tan(radians(lat))*sqrt(1+sig**2)-sig*sqrt(1+(tan(radians(lat))**2))
	conf_lat = atan(conf_lat)
	# Longitude Difference
	long_diff = radians(Decimal(long) - Decimal(str(CM)))
	# Gauss-Schreiber Ratios
	xi1 = atan(tan(conf_lat)/cos(long_diff))
	eta1x = sin(long_diff)/(sqrt(tan(conf_lat)**2+cos(long_diff)**2))
	eta1 = log(eta1x+sqrt(1+eta1x**2))
	# Transverse Mercator Ratios
	eta2 = a2*cos(2*xi1)*sinh(2*eta1)
	eta4 = a4*cos(4*xi1)*sinh(4*eta1)
	eta6 = a6*cos(6*xi1)*sinh(6*eta1)
	eta8 = a8*cos(8*xi1)*sinh(8*eta1)
	eta10 = a10*cos(10*xi1)*sinh(10*eta1)
	eta12 = a12*cos(12*xi1)*sinh(12*eta1)
	eta14 = a14*cos(14*xi1)*sinh(14*eta1)
	eta16 = a16*cos(16*xi1)*sinh(16*eta1)
	xi2 = a2*sin(2*xi1)*cosh(2*eta1)
	xi4 = a4*sin(4*xi1)*cosh(4*eta1)
	xi6 = a6*sin(6*xi1)*cosh(6*eta1)
	xi8 = a8*sin(8*xi1)*cosh(8*eta1)
	xi10 = a10*sin(10*xi1)*cosh(10*eta1)
	xi12 = a12*sin(12*xi1)*cosh(12*eta1)
	xi14 = a14*sin(14*xi1)*cosh(14*eta1)
	xi16 = a16*sin(16*xi1)*cosh(16*eta1)
	eta = eta1 + eta2 + eta4 + eta6 + eta8 + eta10 + eta12 + eta14 + eta16
	xi = xi1 + xi2 + xi4 + xi6 + xi8 + xi10 + xi12 + xi14 + xi16
	# Transverse Mercator Co-ordinates
	X = A*eta
	Y = A*xi
	# MGA Co-ordinates
	E = Proj[4]*Decimal(str(X))+Proj[2]
	N = Proj[4]*Decimal(str(Y))+Proj[3]
	####     Output    ####
	return(',' + str(zone) + ',' + str(round(E,3)) + ',' + str(round(N,3)))
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
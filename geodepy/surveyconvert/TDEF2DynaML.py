# ----------------------------------------------------------------------
#                          TDEF2DynaML.py
# ----------------------------------------------------------------------
#  Author:  Kent Wheeler
#    Date:  14 May 2020
# Purpose:  Script to convert Trimble data exchange files to DynaML 
#           format, identify and export trivial baselines as a kml
# ----------------------------------------------------------------------
#   Usage:  cmd:\> python TDEF2DynaML.py
# ----------------------------------------------------------------------
#   Notes:  - All files with the extension *.asc in the cwd will be
#             attempted to convert to DynaML
#           - only supports G type measurements, other measurement
#             are ignored
#           - the following baseline details are hardcoded
#                  ControlRec.TimeStatus='E'
#                  ControlRec.EphemerisType='B'
#                  ControlRec.AtReceiver='TRIM'
#                  ControlRec.ToReceiver='TRIM'
#                  ControlRec.FrequencyMode='D'
#                  ControlRec.SurveyTechnique='S'
#                  ControlRec.Solution='FX'
#                  ControlRec.EpochInterval=15
#
import subprocess
import math
import numpy
from numpy import matmul, matrix
from math import sin, cos, radians, sqrt
import datetime, os, sqlite3

class DnaStation:
    def __init__(self):
        self.name = ''
        self.Constraint=''
        self.W_Constraint=''
        self.Type=''
        self.XAxis=''
        self.YAxis=''
        self.Height=''
        self.Description=''
        self.aAxis=0
        self.bAxis=0
        self.ErrAz=0
        self.HorizCoordMethod=''
        self.RelativeHorizAccuracy=''
        self.NonGSNumber=''
        self.SelectPoint='true'
        self.SelectRL='true'
class AdditionalInfoMsr:
    def __init__(self):
        #Additional information is included as a comment in the DynaML file. This can be used for database import
        self.StartDateTime=datetime.datetime(1994, 1, 1, 00, 00,00)
        self.Duration=datetime.datetime(1966, 1, 1, 00, 00,00)
        self.FinishDateTime=datetime.datetime(1994, 1, 1, 00, 00,00)
        self.TimeStatus=''
        self.EphemerisType=''
        self.AtReceiver=''
        self.ToReceiver=''
        self.FrequencyMode=''
        self.SurveyTechnique='SLEV'
        self.Solution=''
        self.EpochInterval=''
        self.Class='LC'
        self.LevelDistance='0.01'
        self.InstrumentModel=''
        self.Derivation='MEAS'
        self.NonGSNumber=''
class DnaMeasurement:
    def __init__(self):
        self.type = ''
        self.vscale='1'
        self.pscale='1'
        self.lscale='1'
        self.hscale='1'
        self.first=''
        self.second=''
        self.stddev=''
        self.total=''
        self.instheight=0
        self.targheight=0
        self.targets=''
        self.value=''
        self.targetstddevs=''
        self.dx=''
        self.dy=''
        self.dz=''
        self.MatrixType=''
        self.Vs=numpy.zeros([3,3])
        self.Ds=numpy.zeros([3,3])
        
class DeviceHeight:
    def __init__(self):
        #Device Height might be the height of instrument at a point or Height of target
        self.StnName=[]
        self.RefHeight=[]
def add_DeviceHeight(self,Stn,Hgt):
    self.StnName.append(Stn)
    self.RefHeight.append(Hgt)
def hms2hp(HMS_Ang):
    #Input: HH MM SS.ssss used by Geolab
    #Output: HH.MMSSSsssss used by DynAdjust
    sign=1
    HMS_Ang=HMS_Ang.upper()
    if HMS_Ang.find('S')!=-1 or HMS_Ang.find('-')!=-1:
        sign=-1
    while HMS_Ang.find('  ')!=-1:
        HMS_Ang=HMS_Ang.replace('  ',' ')
    HMS_Ang=HMS_Ang.replace('S','')
    HMS_Ang=HMS_Ang.replace('E','')
    HMS_Ang=HMS_Ang.replace('.',' ')
    aAng=HMS_Ang.split()
    aAng[0]=str(sign*abs(int(aAng[0])))
    aAng[1]="%02d" % float(aAng[1])
    aAng[2]="%02d" % float(aAng[2])
    return aAng[0] + '.' + aAng[1] + aAng[2]+ aAng[3]

def dec2hp(dec):
    minute, second = divmod(abs(dec) * 3600, 60)
    degree, minute = divmod(minute, 60)
    hp = degree + (minute / 100) + (second / 10000)
    return hp if dec >= 0 else -hp

def hp2dec(hp):
    #Input: HH.MMSSsss
    #Output: dd.dddddd
    degmin, second = divmod(abs(hp) * 1000, 10)
    degree, minute = divmod(degmin, 100)
    dec = degree + (minute / 60) + (second / 360)
    return dec if hp >= 0 else -dec

def FindJobNumber(strg):
    #search a string for 8 consecutive numbers, this is probably the Job Number
    JN=''
    i=0
    while i+7!=len(strg):
        if strg[i:i+8].isnumeric()==True:
            JN=strg[i:i+8]
        i=i+1
    return JN
def kmlFooter():
    return '</Document>\n</kml>'
def kmlHeader(nme):
    strg = '<?xml version="1.0" encoding="UTF-8"?>\n'
    strg = strg + '<kml xmlns="http://earth.google.com/kml/2.0">\n'
    strg = strg + ' <Document>\n'
    strg = strg + '   <name>' + nme + '</name>\n'
    strg = strg + '    <Style id="trivlineStyle">\n'
    strg = strg + '        <LineStyle>\n'
    strg = strg + '            <color>660000ff</color>\n'
    strg = strg + '            <width>6</width>\n'
    strg = strg + '        </LineStyle>\n'
    return strg + '    </Style>'
def MkLine(ln):
    strg = '   <Placemark>\n'
    strg = strg + '      <name>' + str(ln[1])+' - - ' + str(ln[2])+ '</name>\n'
    strg = strg + '    <description>Start: ' + str(ln[7]) + '\nFinish: ' + str(ln[8]) + '</description>\n'
    strg = strg + '    <styleUrl>#trivlineStyle</styleUrl>\n'
    strg = strg + '    <LineString>\n'
    strg = strg + '      <extrude>1</extrude>\n'
    strg = strg + '      <tessellate>1</tessellate>\n'
    strg = strg + '      <altitudeMode>ClampToGround</altitudeMode>\n'
    strg = strg + '      <coordinates>\n'
    strg = strg + '       ' +  str(hp2dec(ln[4])) + ',' +  str(hp2dec(ln[3])) + ',0 ' +  str(hp2dec(ln[6])) + ',' +  str(hp2dec(ln[5])) + ',0 \n'
    strg = strg + '      </coordinates>\n'
    strg = strg + '    </LineString>\n'
    return strg +'  </Placemark>'
def Stn_xml_str(Stn):
    #Output: String for printing to xml that is one complete station
    xml_str='<DnaStation>\n'
    xml_str=xml_str+'<Name>' + Stn.Name + '</Name>\n'
    xml_str=xml_str+'<Constraints>' + Stn.Constraint + '</Constraints>\n'
    xml_str=xml_str+'<Type>' + Stn.Type + '</Type>\n'
    xml_str=xml_str+'<StationCoord>\n'
    xml_str=xml_str+'<Name>' + Stn.Name + '</Name>\n'
    xml_str=xml_str+'<XAxis>' + str(Stn.XAxis) + '</XAxis>\n'
    xml_str=xml_str+'<YAxis>' + str(Stn.YAxis) + '</YAxis>\n'
    xml_str=xml_str+'<Height>' + str(Stn.Height) + '</Height>\n'
    xml_str=xml_str+'</StationCoord>\n'
    xml_str=xml_str+'<Description>'+ Stn.Description+'</Description>\n'
    xml_str=xml_str+'<!--AdditionalInfoStn>\n'
    xml_str=xml_str+'<HorizCoordMethod>' + Stn.HorizCoordMethod + '</HorizCoordMethod>\n'
    xml_str=xml_str+'<RelativeHorizAccuracy>' + Stn.RelativeHorizAccuracy + '</RelativeHorizAccuracy>\n'
    xml_str=xml_str+'<NonGSNumber>' + Stn.NonGSNumber + '</NonGSNumber>\n'
    xml_str=xml_str+'<SelectPoint>' + Stn.SelectPoint + '</SelectPoint>\n'
    xml_str=xml_str+'<SelectRL>' + Stn.SelectRL + '</SelectRL>\n'
    xml_str=xml_str+'</AdditionalInfoStn-->\n'    
    xml_str=xml_str+'</DnaStation>\n'
    return xml_str

def Msr_xml_str(row,re_scale=1):
    #Output: xml string for printing to file. Caters for type G, D, S, B, D, L, H
    Msr=DnaMeasurement()
    Msr.type = row[1]
    Msr.vscale=str(row[2]); Msr.pscale=str(row[3]); Msr.lscale=str(row[4]); Msr.hscale=str(row[5])
    Msr.first=row[6]; Msr.second=row[7]
    Msr.value=row[8]; Msr.stddev=row[9]; Msr.total=row[10]
    Msr.instheight=row[11]; Msr.targheight=row[12]
    Msr.targets=row[13]; Msr.targetstddevs=row[14]
    Msr.dx=str(row[15]); Msr.dy=str(row[16]); Msr.dz=str(row[17])
    Msr.Vs=matrix([[row[18],row[19],row[20]],[row[19],row[21],row[22]],[row[20],row[22],row[23]]])
    
    Msr.value=Msr.value.split(',')
    Msr.targets=Msr.targets.split(',')
    Msr.targetstddevs=Msr.targetstddevs.split(',')
    if Msr.vscale!='': Msr.vscale=str(sqrt(float(Msr.vscale)*sqrt(re_scale))**2)
    if Msr.stddev!='':Msr.stddev=str(float(row[9])*sqrt(re_scale))
    if Msr.targetstddevs!='':
        for i in range(1,len(Msr.targetstddevs)):
            Msr.targetstddevs[i]=str(float(Msr.targetstddevs[i])*sqrt(re_scale))
            
    ControlRec=AdditionalInfoMsr()
    ControlRec.StartDateTime=datetime.datetime.strptime(row[24],"%Y-%m-%d %H:%M:%S")
    if Msr.type=='G':ControlRec.Duration=datetime.datetime.strptime(row[26],"%Y-%m-%d %H:%M:%S")
    ControlRec.TimeStatus=row[26]; ControlRec.EphemerisType=row[28]
    ControlRec.AtReceiver=row[28]; ControlRec.ToReceiver=row[30]
    ControlRec.FrequencyMode=row[30]; ControlRec.SurveyTechnique=row[32]
    ControlRec.Solution=row[32]; ControlRec.EpochInterval=row[34]
    ControlRec.Class=row[34]; ControlRec.LevelDistance=row[36]
    ControlRec.InstrumentModel=row[36]; ControlRec.Derivation=row[38]
    ControlRec.NonGSNumber=row[39]
    
    xml_str='<DnaMeasurement>\n'
    xml_str=xml_str+'<Type>' + Msr.type + '</Type>\n'
    xml_str=xml_str+'<Ignore/>\n'

    if Msr.type == 'G':
        xml_str=xml_str+'<ReferenceFrame>' + GNSSdate2Ref(ControlRec.StartDateTime) + '</ReferenceFrame>\n'
        xml_str=xml_str+'<Epoch>' + ControlRec.StartDateTime.strftime('%d.%m.%Y') + '</Epoch>\n'
        xml_str=xml_str+'<Vscale>' + Msr.vscale + '</Vscale>\n'
        xml_str=xml_str+'<Pscale>' + Msr.pscale + '</Pscale>\n'
        xml_str=xml_str+'<Lscale>' + Msr.lscale + '</Lscale>\n'
        xml_str=xml_str+'<Hscale>' + Msr.hscale + '</Hscale>\n'
    xml_str=xml_str+'<First>' + Msr.first + '</First>\n'
    if Msr.second != '':
        xml_str=xml_str+'<Second>' + Msr.second + '</Second>\n'
    if Msr.type != 'G' and Msr.type != 'D':
        xml_str=xml_str+'<Value>' + Msr.value[1] + '</Value>\n'
        xml_str=xml_str+'<StdDev>' + Msr.stddev + '</StdDev>\n'
    if Msr.type == 'G':
        xml_str=xml_str+'<GPSBaseline>\n'
        xml_str=xml_str+'<X>' + Msr.dx + '</X>\n'
        xml_str=xml_str+'<Y>' + Msr.dy + '</Y>\n'
        xml_str=xml_str+'<Z>' + Msr.dz + '</Z>\n'
        xml_str=xml_str+'<SigmaXX>' + str(Msr.Vs[0,0]) + '</SigmaXX>\n'
        xml_str=xml_str+'<SigmaXY>' + str(Msr.Vs[0,1]) + '</SigmaXY>\n'
        xml_str=xml_str+'<SigmaXZ>' + str(Msr.Vs[0,2]) + '</SigmaXZ>\n'
        xml_str=xml_str+'<SigmaYY>' + str(Msr.Vs[1,1]) + '</SigmaYY>\n'
        xml_str=xml_str+'<SigmaYZ>' + str(Msr.Vs[1,2]) + '</SigmaYZ>\n'
        xml_str=xml_str+'<SigmaZZ>' + str(Msr.Vs[2,2]) + '</SigmaZZ>\n'
        xml_str=xml_str+'</GPSBaseline>\n'
        xml_str=xml_str+'<!--AdditionalInfoMsrG>\n'
        xml_str=xml_str+'<StartDateTime>'+ControlRec.StartDateTime.strftime('%Y-%m-%dT%H:%M:%S')+'</StartDateTime>\n'
        xml_str=xml_str+'<Duration>P'+str(ControlRec.Duration.year-1900)+'Y'+str(ControlRec.Duration.month-1)+'M'+str(ControlRec.Duration.day-1)+'DT'+ControlRec.Duration.strftime('%HH%MM%SS')+'</Duration>\n'
        xml_str=xml_str+'<TimeStatus>'+ControlRec.TimeStatus+'</TimeStatus>\n'
        xml_str=xml_str+'<EphemerisType>'+ControlRec.EphemerisType+'</EphemerisType>\n'
        xml_str=xml_str+'<AtReceiver>'+ControlRec.AtReceiver+'</AtReceiver>\n'
        xml_str=xml_str+'<ToReceiver>'+ControlRec.ToReceiver+'</ToReceiver>\n'
        xml_str=xml_str+'<FrequencyMode>'+ControlRec.FrequencyMode+'</FrequencyMode>\n'
        xml_str=xml_str+'<SurveyTechnique>'+ControlRec.SurveyTechnique+'</SurveyTechnique>\n'
        xml_str=xml_str+'<Solution>'+ControlRec.Solution+'</Solution>\n'
        xml_str=xml_str+'<EpochInterval>'+str(ControlRec.EpochInterval)+'</EpochInterval>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrG-->\n'
    
    if Msr.type == 'L':
        xml_str=xml_str+'<!--AdditionalInfoMsrL>\n'
        xml_str=xml_str+'<SurveyTechnique>'+ControlRec.SurveyTechnique+'</SurveyTechnique>\n'
        xml_str=xml_str+'<LevelDistance>'+ ControlRec.LevelDistance +'</LevelDistance>\n'
        xml_str=xml_str+'<ObsDate>'+ControlRec.StartDateTime.strftime('%Y-%m-%d')+'</ObsDate>\n'
        xml_str=xml_str+'<Derivation>'+ControlRec.Derivation+'</Derivation>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrL-->\n'
    
    if Msr.type == 'S':
        xml_str=xml_str+'<InstHeight>' + str(Msr.instheight) + '</InstHeight>\n'
        xml_str=xml_str+'<TargHeight>' + str(Msr.targheight) + '</TargHeight>\n'
        xml_str=xml_str+'<!--AdditionalInfoMsrS>\n'
        xml_str=xml_str+'<InstrumentModel>'+ControlRec.InstrumentModel+'</InstrumentModel>\n'
        xml_str=xml_str+'<ObsDate>'+ControlRec.StartDateTime.strftime('%Y-%m-%d')+'</ObsDate>\n'
        xml_str=xml_str+'<Derivation>'+ControlRec.Derivation+'</Derivation>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrS-->\n'

    if Msr.type == 'D':
        xml_str=xml_str+'<Value>' + Msr.value[1] + '</Value>\n'
        xml_str=xml_str+'<StdDev>' + Msr.targetstddevs[1] + '</StdDev>\n'
        xml_str=xml_str+'<Total>' + str(Msr.total-1) + '</Total>\n'
        ObsNumber=2
        while ObsNumber<=Msr.total:
            xml_str=xml_str+'<Directions>\n'
            xml_str=xml_str+'<Ignore/>\n'
            xml_str=xml_str+'<Target>' + Msr.targets[ObsNumber] + '</Target>\n'
            xml_str=xml_str+'<Value>' + Msr.value[ObsNumber] + '</Value>\n'
            xml_str=xml_str+'<StdDev>' + Msr.targetstddevs[ObsNumber] + '</StdDev>\n'
            xml_str=xml_str+'</Directions>\n'
            ObsNumber=ObsNumber+1
        xml_str=xml_str+'<!--AdditionalInfoMsrD>\n'
        xml_str=xml_str+'<InstrumentModel>'+ControlRec.InstrumentModel+'</InstrumentModel>\n'
        xml_str=xml_str+'<ObsDate>'+ControlRec.StartDateTime.strftime('%Y-%m-%d')+'</ObsDate>\n'
        xml_str=xml_str+'<Derivation>'+ControlRec.Derivation+'</Derivation>\n'
        xml_str=xml_str+'<Class>'+ControlRec.Class+'</Class>\n'
        xml_str=xml_str+'<NonGSNumber>'+ControlRec.NonGSNumber+'</NonGSNumber>\n'
        xml_str=xml_str+'</AdditionalInfoMsrD-->\n'
    xml_str=xml_str+'<Source></Source>\n'
    xml_str=xml_str+'</DnaMeasurement>\n'
    return xml_str

c_vac = 299792.458
k_0 = 0.9996

# Ellipsoid Constants
class Ellipsoid(object):
    def __init__(self, semimaj, inversef):
        self.semimaj = semimaj
        self.inversef = inversef
        self.f = 1 / self.inversef
        self.semimin = float(self.semimaj * (1 - self.f))
        self.ecc1sq = float(self.f * (2 - self.f))
        self.ecc2sq = float(self.ecc1sq / (1 - self.ecc1sq))
        self.ecc1 = sqrt(self.ecc1sq)
        self.n = float(self.f / (2 - self.f))
        self.n2 = self.n ** 2

# Geodetic Reference System 1980
grs80 = Ellipsoid(6378137, 298.25722210088)

def llh2xyz(lat, lng, ellht, ellipsoid=grs80):
    # Add input for ellipsoid (default: grs80)
    # Convert lat & long to radians
    lat = radians(hp2dec(float(lat)))
    lng = radians(hp2dec(float(lng)))
    ellht=float(ellht)
    # Calculate Ellipsoid Radius of Curvature in the Prime Vertical - nu
    if lat == 0:
        nu = grs80.semimaj
    else:
        nu = ellipsoid.semimaj/(sqrt(1 - ellipsoid.ecc1sq * (sin(lat)**2)))
    # Calculate x, y, z
    x = (nu + ellht) * cos(lat) * cos(lng)
    y = (nu + ellht) * cos(lat) * sin(lng)
    z = ((ellipsoid.semimin**2 / ellipsoid.semimaj**2) * nu + ellht) * sin(lat)
    return x, y, z

def ErrEllip2Ycluster(Stn,w_H):
    #Input: Supply a station with coordinates and error ellipse for coordinate uncertainty
    #Output: xml string for  point cluster (Y-type observation)
    x, y, z = llh2xyz(Stn.XAxis, Stn.YAxis, Stn.Height)

    a=Stn.aAxis/2.44774683068
    b=Stn.bAxis/2.44774683068
    Az=90-Stn.ErrAz
    
    rAz=math.radians(Az)
    rlat=math.radians(float(Stn.XAxis))
    rlng=math.radians(float(Stn.YAxis))

    rl=numpy.zeros([3,3])
    rl[0,0]=-sin(rlng)
    rl[0,1]=-sin(rlat)*cos(rlng)
    rl[0,2]=cos(rlat)*cos(rlng)
    rl[1,0]=cos(rlng)
    rl[1,1]=-sin(rlat)*sin(rlng)
    rl[1,2]=cos(rlat)*sin(rlng)
    rl[2,1]=cos(rlat)
    rl[2,2]=sin(rlat)

    iA=numpy.zeros([3,3])
    iA[0,0]=(cos(rAz)*cos(rAz)*a*a)+(b*b*sin(rAz)*sin(rAz))
    iA[0,1]=(a*a-b*b)*cos(rAz)*sin(rAz)
    iA[1,0]=iA[0,1]
    iA[1,1]=(a*a*sin(rAz)*sin(rAz))+(b*b*cos(rAz)*cos(rAz))
    iA[2,2]=w_H**2
    
    Wt=matmul(matmul(rl,iA),rl.transpose())
    
    xml_str='<DnaMeasurement>\n'
    xml_str=xml_str+'<!--Converted GESMAR 95% Error Ellipse, aAxis: ' + str(Stn.aAxis) + ' bAxis: ' + str(Stn.bAxis) + ' Az: ' + str(Stn.ErrAz) + ' -->\n'
    xml_str=xml_str+'<!--Converted GESMAR 1SD Error Ellipse, aAxis: ' + str(round(a,4)) + ' bAxis: ' + str(round(b,4)) + ' Az: ' + str(Stn.ErrAz) + ' H: ' + str(w_H) + ' -->\n'
    xml_str=xml_str+'<Type>Y</Type>\n'
    xml_str=xml_str+'<Ignore/>\n'
    xml_str=xml_str+'<ReferenceFrame>GDA2020</ReferenceFrame>\n'
    xml_str=xml_str+'<Epoch>01.01.2020</Epoch>\n'
    xml_str=xml_str+'<Vscale>1.000</Vscale>\n'
    xml_str=xml_str+'<Pscale>1.000</Pscale>\n'
    xml_str=xml_str+'<Lscale>1.000</Lscale>\n'
    xml_str=xml_str+'<Hscale>1.000</Hscale>\n'
    xml_str=xml_str+'<Coords>XYZ</Coords>\n'
    xml_str=xml_str+'<Total>1</Total>\n'
    xml_str=xml_str+'<First>' + Stn.Name + '</First>\n'
    xml_str=xml_str+'<Clusterpoint>\n'
    xml_str=xml_str+'<X>'+str(x)+'</X>\n'
    xml_str=xml_str+'<Y>'+str(y)+'</Y>\n'
    xml_str=xml_str+'<Z>'+str(z)+'</Z>\n'
    xml_str=xml_str+'<SigmaXX>'+str(Wt[0,0])+'</SigmaXX>\n'
    xml_str=xml_str+'<SigmaXY>'+str(Wt[0,1])+'</SigmaXY>\n'
    xml_str=xml_str+'<SigmaXZ>'+str(Wt[0,2])+'</SigmaXZ>\n'
    xml_str=xml_str+'<SigmaYY>'+str(Wt[1,1])+'</SigmaYY>\n'
    xml_str=xml_str+'<SigmaYZ>'+str(Wt[1,2])+'</SigmaYZ>\n'
    xml_str=xml_str+'<SigmaZZ>'+str(Wt[2,2])+'</SigmaZZ>\n'
    xml_str=xml_str+'</Clusterpoint>\n'
    xml_str=xml_str+'</DnaMeasurement>\n'
    
    return xml_str

def stn_header():
    xml_str='<?xml version="1.0" encoding="utf-8"?>\n'
    xml_str=xml_str+'<DnaXmlFormat type="Station File" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="DynaML.xsd">\n'
    return xml_str
def msr_header():
    xml_str='<?xml version="1.0" encoding="utf-8"?>\n'
    xml_str=xml_str+'<DnaXmlFormat type="Measurement File" referenceframe="GDA94" epoch="01.01.1994" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="DynaML.xsd">\n'
    return xml_str
def dML_footer():
    xml_str='</DnaXmlFormat>\n'
    return xml_str
def GNSSdate2Ref(obsDate):
    #Use the date of GNSS baseline obsedrvation to determine the reference frame used by broadcast ephemeris
    if obsDate >= datetime.datetime(1900, 1, 1) and obsDate < datetime.datetime(1994, 1, 2):    
        Ref = 'ITRF1991'
    if obsDate >= datetime.datetime(1994, 1, 2) and obsDate < datetime.datetime(1995, 1, 1):    
        Ref = 'ITRF1992'
    if obsDate >= datetime.datetime(1995, 1, 1) and obsDate < datetime.datetime(1996, 6, 30):    
        Ref = 'ITRF1993'
    if obsDate >= datetime.datetime(1996, 6, 30) and obsDate < datetime.datetime(1998, 3, 1):    
        Ref = 'ITRF1994'
    if obsDate >= datetime.datetime(1998, 3, 1) and obsDate < datetime.datetime(1999, 8, 1):    
        Ref = 'ITRF1996'
    if obsDate >= datetime.datetime(1999, 8, 1) and obsDate < datetime.datetime(2001, 12, 2):    
        Ref = 'ITRF1997'
    if obsDate >= datetime.datetime(2001, 12, 2) and obsDate < datetime.datetime(2006, 11, 5):    
        Ref = 'ITRF2000'
    if obsDate >= datetime.datetime(2006, 11, 5) and obsDate < datetime.datetime(2011, 4, 17):    
        Ref = 'ITRF2005'
    if obsDate >= datetime.datetime(2011, 4, 17):    
        Ref = 'ITRF2008'
    return Ref
#####################################################################################
#### Input:Trimble Business Centre *.asc file                                   #####
#### Output: DynaML stn and msr file,                                           #####
####         if the GNSS network contains trivial baslines, it will output      #####
####         these to a google earth kml                                        #####
#####################################################################################

def create_connection(db_file):
    """ create a database connection to a SQLite database """
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except:
        print ("Cannot create a database'")
    finally:
        conn.close()
for f in os.listdir(os.getcwd()):
    if f.endswith('.asc'):
        adjustment_name = f.replace('.asc','')
        # Connect to Sqlite and Open a new database
        dbname = adjustment_name+'.db'
        if not os.path.exists(dbname):
            create_connection(dbname)
            conn = sqlite3.connect(dbname) # or use :memory: to put it in RAM
        else:
            conn = sqlite3.connect(dbname) # or use :memory: to put it in RAM     
        cursor = conn.cursor()
        # Create 4 Tables, Store Stations, Observations, Stations in discrete networks, Stations connected with directions
        cursor.execute('DROP TABLE IF EXISTS STATIONS')
        cursor.execute("""CREATE TABLE IF NOT EXISTS STATIONS (
                      ID integer PRIMARY KEY, STATION_NAME short text, COORD_TYPE text,
                      CONSTRAIN text, W_CONSTRAIN text, LATITUDE double, LONGITUDE double, HEIGHT double, E_HEIGHT double,
                      DESC text, GES_NAME text, GES94_LATITUDE double, GES94_LONGITUDE double, GES_HEIGHT double, HT_ACCURACY text, HT_METHOD text,
                      HZ_ORDER text, HZ_ACCURACY text, HZ_METHOD text, CE double, A_AXIS double, B_AXIS double, ERR_AZ double,
                      GES2020_LATITUDE double, GES2020_LONGITUDE double);""")
        conn.commit()
        cursor.execute('DROP TABLE IF EXISTS OBSERVATIONS')
        cursor.execute("""CREATE TABLE IF NOT EXISTS OBSERVATIONS (
                      ID integer PRIMARY KEY, TYPE text, VSCALE double, PSCALE double, LSCALE double, HSCALE double, FIRST short text, SECOND short text,
                      VALUE text, SDEV text, TOTAL integer, INST_HEIGHT double, TARG_HEIGHT double, TARGETS text, TARGETS_SDEV text, DX double,
                      DY double, DZ double, VS_1_1 double, VS_1_2 double, VS_1_3 double, VS_2_1 double, VS_2_2 double, VS_3_1 double,
                      StartDateTime date,FinishDateTime date, Duration time, TimeStatus text, EphemerisType text, AtReceiver text, ToReceiver text, FrequencyMode text,
                      SurveyTechnique text, Solution text, EpochInterval text, Class text, LevelDistance text, InstrumentModel text,
                      Derivation text, NON_GS text, SESSION integer, TRIVIAL boolean);""")
        conn.commit()
        cursor.execute('DROP TABLE IF EXISTS DIR_TARGETS')
        cursor.execute("""CREATE TABLE IF NOT EXISTS DIR_TARGETS (
                      ID integer PRIMARY KEY, OBSERVATIONS_ID integer, TARGETS short text, VALUE text, TARGETS_SDEV text);""")
        conn.commit()
        cursor.execute('DROP TABLE IF EXISTS NETWORKS')
        cursor.execute("""CREATE TABLE IF NOT EXISTS NETWORKS (
                      ID integer PRIMARY KEY, STATION_NAME text, NETWORK integer);""")
        conn.commit()
        cursor.execute('DROP TABLE IF EXISTS GNSS_SESSIONS')
        cursor.execute("""CREATE TABLE IF NOT EXISTS GNSS_SESSIONS (
                      ID integer PRIMARY KEY, NETWORK integer, SESSION date, STATION_NAME text, StartDateTime text, FinishDateTime text);""")
        conn.commit()
        f = open(adjustment_name+'.asc', 'r')
        # Run through each line of the Trimble Data Exchange file and extract the relevant lines
        lineCount=0; idcount=0; Tgtidcount=0; GNSSmarks='';time_zero = datetime.datetime.strptime('00:00:00', '%H:%M:%S')
        InstHts=DeviceHeight(); TgtHts=DeviceHeight(); CurrentMsr=DnaMeasurement(); ControlRec=AdditionalInfoMsr()
        jobNumber=FindJobNumber(os.getcwd())
        print('Reading the Trimble Data Exchange (.asc) File...')
        for linestr in f.readlines():
            if linestr.startswith('Station='):
                ilinestr=linestr.split(':')
                CurrentStn=DnaStation()
                CurrentStn.Name=ilinestr[2].strip()
                if ilinestr[9].strip()=='0':CurrentStn.Constraint='FF'
                if ilinestr[9].strip()=='1':CurrentStn.Constraint='CC'
                if ilinestr[10].strip()=='1' or ilinestr[11].strip()=='1':CurrentStn.Constraint=CurrentStn.Constraint+'C'
                else:CurrentStn.Constraint=CurrentStn.Constraint+'F'
                CurrentStn.Constraint=CurrentStn.Constraint.strip()
                CurrentStn.Type='LLH'
                CurrentStn.XAxis=dec2hp(-1*float(ilinestr[3].strip()[:-1]))
                CurrentStn.YAxis=dec2hp(float(ilinestr[4].strip()[:-1]))
                CurrentStn.Height=ilinestr[8].strip()
                CurrentStn.NonGSNumber='E'+jobNumber
                stnRec=('||||||||||||||').split('|')
                cursor.execute("INSERT INTO STATIONS (COORD_TYPE, STATION_NAME, CONSTRAIN, W_CONSTRAIN, LATITUDE, LONGITUDE, HEIGHT, DESC, GES_NAME, GES94_LATITUDE, GES94_LONGITUDE, GES_HEIGHT, HT_ACCURACY, HT_METHOD, HZ_ORDER, HZ_ACCURACY, HZ_METHOD, CE, A_AXIS, B_AXIS, ERR_AZ, GES2020_LATITUDE, GES2020_LONGITUDE) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", 
                                [CurrentStn.Type,CurrentStn.Name,CurrentStn.Constraint,CurrentStn.W_Constraint, CurrentStn.XAxis,CurrentStn.YAxis, CurrentStn.Height, CurrentStn.Description, stnRec[0],stnRec[1],stnRec[2],stnRec[3],stnRec[4],stnRec[5],stnRec[6],stnRec[7],stnRec[8],stnRec[9],stnRec[10],stnRec[11],stnRec[12],stnRec[13],stnRec[14]])
                conn.commit()
        
        # Scrape information for Landgate GESMAR control Records
            if linestr.startswith('Vector='):
                ControlRec=AdditionalInfoMsr()
                ilinestr=linestr.split(':')
                stdate=ilinestr[19]
                sttime=ilinestr[20]
                istdate=stdate.split()
                isttime=sttime.split()
                yr=int(istdate[2])
                mth=int(istdate[1])
                ddy=int(istdate[0])
                hr=int(isttime[0])
                mn=int(isttime[1])
                sc=int(float(isttime[2]))
                ControlRec.StartDateTime=datetime.datetime(yr, mth, ddy, hr, mn,sc)
                ControlRec.NonGSNumber='E' + jobNumber
                fhdate=ilinestr[21]
                fhtime=ilinestr[22]
                ifhdate=fhdate.split()
                ifhtime=fhtime.split()
                yr=int(ifhdate[2])
                mth=int(ifhdate[1])
                ddy=int(ifhdate[0])
                hr=int(ifhtime[0])
                mn=int(ifhtime[1])
                sc=int(float(ifhtime[2]))
                ControlRec.FinishDateTime=datetime.datetime(yr, mth, ddy, hr, mn,sc)
                ControlRec.TimeStatus='E'
                ControlRec.EphemerisType='B'
                ControlRec.AtReceiver='TRIM'
                ControlRec.ToReceiver='TRIM'
                ControlRec.FrequencyMode='D'
                ControlRec.SurveyTechnique='S'
                ControlRec.Solution='FX'
                ControlRec.EpochInterval=15
                ControlRec.Class=''
        
                CurrentMsr=DnaMeasurement()
                CurrentMsr.type='G'
                CurrentMsr.first=ilinestr[2]
                CurrentMsr.second=ilinestr[3]
                CurrentMsr.dx=ilinestr[4]
                CurrentMsr.dy=ilinestr[5]
                CurrentMsr.dz=ilinestr[6]
                GNSSmarks=GNSSmarks + ';' + CurrentMsr.first + ';' + CurrentMsr.second
        
                CurrentMsr.MatrixType='COV'
                CurrentMsr.vscale=1
                CurrentMsr.Vs[0,0]=ilinestr[7]; CurrentMsr.Vs[1,0]=ilinestr[8]; CurrentMsr.Vs[2,0]=ilinestr[9]
                CurrentMsr.Vs[0,1]=ilinestr[10]; CurrentMsr.Vs[1,0]=ilinestr[11]
                CurrentMsr.Vs[1,1]=ilinestr[12]
                cursor.execute("INSERT INTO OBSERVATIONS (TYPE, VSCALE, PSCALE, LSCALE, HSCALE, FIRST, SECOND, VALUE, SDEV, TOTAL, INST_HEIGHT, TARG_HEIGHT, TARGETS, TARGETS_SDEV, DX, DY, DZ, VS_1_1, VS_1_2, VS_1_3, VS_2_1, VS_2_2, VS_3_1, StartDateTime, FinishDateTime, Duration, TimeStatus, EphemerisType, AtReceiver, ToReceiver, FrequencyMode, SurveyTechnique, Solution, EpochInterval, Class, LevelDistance, InstrumentModel, Derivation, NON_GS) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", 
                           [CurrentMsr.type, CurrentMsr.vscale, CurrentMsr.pscale, CurrentMsr.lscale, CurrentMsr.hscale, CurrentMsr.first, CurrentMsr.second, CurrentMsr.value,  CurrentMsr.stddev, CurrentMsr.total, CurrentMsr.instheight, CurrentMsr.targheight, CurrentMsr.targets, CurrentMsr.targetstddevs, CurrentMsr.dx, CurrentMsr.dy, CurrentMsr.dz, CurrentMsr.Vs[0,0], CurrentMsr.Vs[0,1], CurrentMsr.Vs[0,2], CurrentMsr.Vs[1,1], CurrentMsr.Vs[1,2], CurrentMsr.Vs[2,2],ControlRec.StartDateTime,ControlRec.FinishDateTime, ControlRec.Duration, ControlRec.TimeStatus, ControlRec.EphemerisType, ControlRec.AtReceiver, ControlRec.ToReceiver, ControlRec.FrequencyMode, ControlRec.SurveyTechnique, ControlRec.Solution, ControlRec.EpochInterval, ControlRec.Class, ControlRec.LevelDistance, ControlRec.InstrumentModel, ControlRec.Derivation, ControlRec.NonGSNumber])
                conn.commit()    
        
            lineCount=lineCount+1
        f.close
        
        ### Write the database to DynaML formats ###
        print('\nWriting DynaML Files...')
        stnout = open(adjustment_name + '.stn.xml', 'w')
        msrout = open(adjustment_name + '.msr.xml', 'w')
        stnout.write(stn_header())
        msrout.write(msr_header())
        sqlstring="SELECT STATIONS.STATION_NAME, STATIONS.CONSTRAIN, STATIONS.COORD_TYPE, STATIONS.LATITUDE, STATIONS.LONGITUDE, STATIONS.HEIGHT, STATIONS.E_HEIGHT, STATIONS.DESC, STATIONS.HZ_METHOD, STATIONS.HZ_ACCURACY, STATIONS.A_AXIS, STATIONS.B_AXIS, STATIONS.ERR_AZ \
                    FROM STATIONS \
                    ORDER BY STATIONS.ID;"
        qry=cursor.execute(sqlstring).fetchall()
        print(' Writing:' + stnout.name)
        for row in qry:
            Stn=DnaStation()
            Stn.Name=row[0]; Stn.Constraint=row[1]; Stn.Type=row[2]
            Stn.XAxis=row[3]; Stn.YAxis=row[4]; Stn.Height=row[5]; Stn.Description=row[7]
            Stn.aAxis=row[10]; Stn.bAxis=row[11]; Stn.ErrAz=row[12]
            Stn.HorizCoordMethod=row[8]; Stn.RelativeHorizAccuracy=row[9]; Stn.NonGSNumber='E'+jobNumber
            stnout.write(Stn_xml_str(Stn))
        sqlstring="SELECT OBSERVATIONS.* \
                    FROM OBSERVATIONS \
                    ORDER BY OBSERVATIONS.ID;"
        qry=cursor.execute(sqlstring).fetchall()
        print(' Writing:' + msrout.name)
        for row in qry:
            msrout.write(Msr_xml_str(row))
        stnout.write(dML_footer())
        msrout.write(dML_footer())
        stnout.close()
        msrout.close()
        
        ##################################################################################
        ### If this is a final adjustment, break the adjustment into discrete networks ###
        ##################################################################################
        if adjustment_name.lower().find('final')!=-1 and cursor.execute("SELECT Count(STATIONS.A_Axis) AS CountOfErrorEllipses FROM STATIONS where STATIONS.A_AXIS<>'';").fetchall()!=0:
            print('\nBreaking adjustment into discrete networks...')
            network_num=1; prevCnt=0
            while prevCnt!=cursor.execute("SELECT Count(STATIONS.STATION_NAME) AS CountOfSTATION_NAME FROM STATIONS;").fetchall():
                n_adjustment_name=adjustment_name +'_' + str(network_num)
                nW_adjustment_name=adjustment_name +'_' + str(network_num) +'W'
                print(' Network: ' + n_adjustment_name)
                sqlstring="INSERT INTO NETWORKS ( STATION_NAME, NETWORK ) \
                    SELECT STATIONS.STATION_NAME, "+ str(network_num) + " AS NETWORK \
                    FROM STATIONS LEFT JOIN NETWORKS ON STATIONS.STATION_NAME = NETWORKS.STATION_NAME \
                    WHERE (((NETWORKS.STATION_NAME) Is Null)) \
                    LIMIT 1;"
                cursor.execute(sqlstring)
                conn.commit()
                while prevCnt!=cursor.execute("SELECT Count(NETWORKS.STATION_NAME) AS CountOfSTATION_NAME FROM NETWORKS;").fetchall():
                    prevCnt = cursor.execute("SELECT Count(NETWORKS.STATION_NAME) AS CountOfSTATION_NAME FROM NETWORKS;").fetchall()
                    sqlstring="INSERT INTO NETWORKS ( STATION_NAME, NETWORK ) \
                            SELECT TMP.STATION_NAME, "+ str(network_num) + " AS NETWORK \
                            FROM (SELECT OBSERVATIONS.SECOND AS STATION_NAME \
                            FROM NETWORKS INNER JOIN OBSERVATIONS ON NETWORKS.STATION_NAME = OBSERVATIONS.FIRST \
                            UNION \
                            SELECT OBSERVATIONS.FIRST AS STATION_NAME \
                            FROM NETWORKS INNER JOIN OBSERVATIONS ON NETWORKS.STATION_NAME = OBSERVATIONS.SECOND \
                            UNION \
                            SELECT DIR_TARGETS.TARGETS AS STATION_NAME \
                            FROM (NETWORKS INNER JOIN OBSERVATIONS ON NETWORKS.STATION_NAME = OBSERVATIONS.FIRST) INNER JOIN DIR_TARGETS ON OBSERVATIONS.ID = DIR_TARGETS.OBSERVATIONS_ID \
                            UNION \
                            SELECT OBSERVATIONS.FIRST AS STATION_NAME \
                            FROM NETWORKS INNER JOIN (OBSERVATIONS INNER JOIN DIR_TARGETS ON OBSERVATIONS.ID = DIR_TARGETS.OBSERVATIONS_ID) ON NETWORKS.STATION_NAME = DIR_TARGETS.TARGETS \
                            UNION \
                            SELECT NETWORKS.STATION_NAME AS STATION_NAME \
                            FROM NETWORKS) AS TMP \
                            LEFT JOIN NETWORKS ON TMP.STATION_NAME = NETWORKS.STATION_NAME \
                            WHERE (((TMP.STATION_NAME)<>\"\") AND ((NETWORKS.STATION_NAME) Is Null));"
                    cursor.execute(sqlstring)
                    conn.commit()
                ### Create non weighted adjustment files for each discrete network ###
                n_stnout = open(adjustment_name +'_' + str(network_num) + '.stn.xml', 'w')
                n_msrout = open(adjustment_name +'_' + str(network_num) + '.msr.xml', 'w')
                n_stnout.write(stn_header())
                n_msrout.write(msr_header())
                sqlstring="SELECT STATIONS.STATION_NAME, STATIONS.CONSTRAIN, STATIONS.COORD_TYPE, STATIONS.LATITUDE, STATIONS.LONGITUDE, STATIONS.HEIGHT, STATIONS.E_HEIGHT, STATIONS.DESC, STATIONS.HZ_METHOD, STATIONS.HZ_ACCURACY, TMP.VALUE, TMP.SDEV, STATIONS.A_AXIS, STATIONS.B_AXIS, STATIONS.ERR_AZ \
                            FROM STATIONS INNER JOIN (NETWORKS LEFT JOIN (SELECT OBSERVATIONS.* \
                            FROM OBSERVATIONS  WHERE (((OBSERVATIONS.TYPE)=\"H\")))  AS TMP ON NETWORKS.STATION_NAME = TMP.FIRST) ON STATIONS.STATION_NAME = NETWORKS.STATION_NAME \
                            WHERE (((NETWORKS.NETWORK)="+ str(network_num) + "));"
                qry=cursor.execute(sqlstring).fetchall()
                for row in qry:
                    Stn=DnaStation()
                    Stn.Name=row[0]; Stn.Constraint=row[1]; Stn.Type=row[2]
                    Stn.XAxis=row[3]; Stn.YAxis=row[4]; Stn.Height=row[5]; Stn.Description=row[7]
                    Stn.aAxis=row[12]; Stn.bAxis=row[13]; Stn.ErrAz=row[14]
                    Stn.HorizCoordMethod=row[8]; Stn.RelativeHorizAccuracy=row[9]; Stn.NonGSNumber='E'+jobNumber
                    n_stnout.write(Stn_xml_str(Stn))
        
                sqlstring="SELECT OBSERVATIONS.* \
                            FROM (NETWORKS INNER JOIN OBSERVATIONS ON NETWORKS.STATION_NAME = OBSERVATIONS.FIRST) \
                            WHERE (((NETWORKS.NETWORK)="+ str(network_num) + ")) \
                            ORDER BY OBSERVATIONS.ID;"
                qry=cursor.execute(sqlstring).fetchall()
                for row in qry:
                    n_msrout.write(Msr_xml_str(row))
            
                n_stnout.write(dML_footer())
                n_msrout.write(dML_footer())
                n_stnout.close()
                n_msrout.close()
                # Run Dynadjust on the network '_n'
                print('  Adjusting: ' + n_adjustment_name)
                subprocess.run("import -n " + n_adjustment_name + " " + n_adjustment_name + ".msr.xml "+ n_adjustment_name + ".stn.xml --flag-unused-stations --remove-ignored-msr")
                subprocess.run("geoid " + n_adjustment_name + " -g \"X:\GA_Processing - AUSGeoid2020\AusGeoid2020_V1.7\AUSGeoid2020_20170908.gsb\" --convert")
                subprocess.run("segment " + n_adjustment_name + " --min-inner-stns 500 --max-block-stns 500")
                subprocess.run("adjust " + n_adjustment_name + " --staged --create-stage-files --output-adj-msr --output-pos-uncertainty --output-adj-gnss-units 1 --max-iterations 20 --free-stn-sd 5 --iteration-threshold 0.0005 --stn-coord-types ENzPLHhXYZ")
                
                # Open the adjusted files and extract the sigma 0 and the Ellipse Height
                n_adj = open(n_adjustment_name + ".phased-stage.adj", 'r')
                lineCount=0
                read_height=False
                for linestr in n_adj.readlines():
                    if linestr[:35].strip() =='Rigorous Sigma Zero':sigma_0=float(linestr[-20:].strip())
                    if linestr.strip() == 'Adjusted Coordinates': lineCount=0; read_height=True
                    if lineCount>=5 and read_height==True and linestr!='\n':
                        sqlstring='UPDATE STATIONS \
                            SET E_HEIGHT = ' + linestr[103:113].strip() +' \
                            WHERE STATION_NAME = \'' + linestr[:20].strip() +'\';'
                        cursor.execute(sqlstring)
                        conn.commit()
                    lineCount=lineCount+1
                n_adj.close
                ### Create weighted adjustment files, scaled by the sigma_0 for each discrete network ###
                w_stnout = open(nW_adjustment_name + '.stn.xml', 'w')
                w_msrout = open(nW_adjustment_name + '.msr.xml', 'w')
                w_stnout.write(stn_header())
                w_msrout.write(msr_header())
                sqlstring="SELECT STATIONS.STATION_NAME, STATIONS.W_CONSTRAIN, STATIONS.COORD_TYPE, STATIONS.LATITUDE, STATIONS.LONGITUDE, STATIONS.HEIGHT, STATIONS.E_HEIGHT, STATIONS.DESC, STATIONS.HZ_METHOD, STATIONS.HZ_ACCURACY, TMP.VALUE, TMP.SDEV, STATIONS.A_AXIS, STATIONS.B_AXIS, STATIONS.ERR_AZ, STATIONS.CONSTRAIN \
                            FROM STATIONS INNER JOIN (NETWORKS LEFT JOIN (SELECT OBSERVATIONS.* \
                            FROM OBSERVATIONS  WHERE (((OBSERVATIONS.TYPE)=\"H\")))  AS TMP ON NETWORKS.STATION_NAME = TMP.FIRST) ON STATIONS.STATION_NAME = NETWORKS.STATION_NAME \
                            WHERE (((NETWORKS.NETWORK)="+ str(network_num) + "));"
                qry=cursor.execute(sqlstring).fetchall()
                for row in qry:
                    Stn=DnaStation()
                    Stn.Name=row[0]; Stn.Constraint=row[1]; Stn.Type=row[2]
                    Stn.XAxis=row[3]; Stn.YAxis=row[4]; Stn.Height=row[5]; Stn.Description=row[7]
                    Stn.aAxis=row[12]; Stn.bAxis=row[13]; Stn.ErrAz=row[14]
                    Stn.HorizCoordMethod=row[8]; Stn.RelativeHorizAccuracy=row[9]; Stn.NonGSNumber='E'+jobNumber
                    if Stn.HorizCoordMethod=='' and GNSSmarks.find(';'+Stn.name+';')!=-1: Stn.HorizCoordMethod='GNSS'
                    else: Stn.HorizCoordMethod='GEOD'
                    w_stnout.write(Stn_xml_str(Stn))
                    if row[15][:2] =='CC':
                        if Stn.Constraint=='FFC':
                            w_msrout.write(ErrEllip2Ycluster(CurrentStn,0.001))
                        else:
                            w_msrout.write(ErrEllip2Ycluster(CurrentStn,100))
                sqlstring="SELECT OBSERVATIONS.* \
                            FROM (NETWORKS INNER JOIN OBSERVATIONS ON NETWORKS.STATION_NAME = OBSERVATIONS.FIRST) \
                            WHERE (((NETWORKS.NETWORK)="+ str(network_num) + ")) \
                            ORDER BY OBSERVATIONS.ID;"
                qry=cursor.execute(sqlstring).fetchall()
                for row in qry:
                    w_msrout.write(Msr_xml_str(row,sigma_0))
             
                w_stnout.write(dML_footer())
                w_msrout.write(dML_footer())
                w_stnout.close()
                w_msrout.close()
                # Run Dynadjust on the network '_n'
                print('  Adjusting: ' + nW_adjustment_name)
                subprocess.run("import -n " + nW_adjustment_name + " " + nW_adjustment_name + ".msr.xml "+ nW_adjustment_name + ".stn.xml --flag-unused-stations --remove-ignored-msr")
                subprocess.run("geoid " + nW_adjustment_name + " -g \"X:\GA_Processing - AUSGeoid2020\AusGeoid2020_V1.7\AUSGeoid2020_20170908.gsb\" --convert")
                subprocess.run("segment " + nW_adjustment_name + " --min-inner-stns 500 --max-block-stns 500")
                subprocess.run("adjust " + nW_adjustment_name + " --staged --create-stage-files --output-adj-msr --output-pos-uncertainty --output-adj-gnss-units 1 --max-iterations 20 --free-stn-sd 5 --iteration-threshold 0.0005 --stn-coord-types ENzPLHhXYZ")
        
                network_num=network_num+1
        
        ####################################
        ### Search for Trivial Baselines ###
        ####################################
        w_kml = open(adjustment_name+'.kml', 'w')
        w_kml.write(kmlHeader(adjustment_name))
        
        #Break into individual networks for each session (event)
        sqlstring="SELECT DISTINCT OBSERVATIONS.StartDateTime AS EventTime \
            FROM OBSERVATIONS WHERE ((OBSERVATIONS.TYPE)='G') \
            ORDER BY EventTime;"
        Events=cursor.execute(sqlstring).fetchall()
        for e in Events:
            network_num=0; prevNtCnt=0
            while prevNtCnt!=cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall():
                prevNtCnt=cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall()
                network_num=network_num+1
                sqlstring="INSERT INTO GNSS_SESSIONS (NETWORK, SESSION, STATION_NAME, StartDateTime, FinishDateTime) \
                SELECT "+ str(network_num) + " AS NETWORK, '"+ str(e[0]) + "' AS SESSION, OBSERVATIONS.FIRST AS STATION_NAME, OBSERVATIONS.StartDateTime, OBSERVATIONS.FinishDateTime \
                FROM OBSERVATIONS LEFT JOIN GNSS_SESSIONS ON OBSERVATIONS.FIRST = GNSS_SESSIONS.STATION_NAME \
                WHERE (((OBSERVATIONS.StartDateTime)<='"+ e[0] + "') AND ((OBSERVATIONS.FinishDateTime)>='"+ e[0] + "') AND ((OBSERVATIONS.TYPE)='G') AND ((GNSS_SESSIONS.STATION_NAME) Is Null)) \
                LIMIT 1"
                cursor.execute(sqlstring)
                conn.commit()
                prevCnt=0
                while prevCnt!=cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall():
                    prevCnt = cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall()
                    sqlstring="INSERT INTO GNSS_SESSIONS (NETWORK, SESSION, STATION_NAME, StartDateTime, FinishDateTime) \
                    SELECT C.NETWORK, C.SESSION, C.STATION_NAME, C.StartDateTime, C.FinishDateTime \
                    FROM (SELECT C_FIRST.NETWORK, C_FIRST.SESSION, C_FIRST.STATION_NAME, C_FIRST.StartDateTime, C_FIRST.FinishDateTime \
                    FROM( \
                    SELECT MARKS_IN_SESS.NETWORK, MARKS_IN_SESS.SESSION, BASE_IN_SESS.FIRST AS STATION_NAME, BASE_IN_SESS.StartDateTime, BASE_IN_SESS.FinishDateTime \
                    FROM (SELECT GNSS_SESSIONS.NETWORK, GNSS_SESSIONS.SESSION, GNSS_SESSIONS.STATION_NAME \
                    FROM GNSS_SESSIONS \
                    WHERE (((GNSS_SESSIONS.NETWORK)="+ str(network_num) + ") AND ((GNSS_SESSIONS.SESSION)='"+ e[0] + "')) \
                    )  AS MARKS_IN_SESS  \
                    INNER JOIN (SELECT OBSERVATIONS.FIRST, OBSERVATIONS.SECOND, OBSERVATIONS.StartDateTime, OBSERVATIONS.FinishDateTime \
                    FROM OBSERVATIONS \
                    WHERE (((OBSERVATIONS.TYPE)='G') AND ((OBSERVATIONS.StartDateTime)<='"+ e[0] + "') AND ((OBSERVATIONS.FinishDateTime)>='"+ e[0] + "')) \
                    )  AS BASE_IN_SESS ON MARKS_IN_SESS.STATION_NAME = BASE_IN_SESS.SECOND) AS C_FIRST \
                    UNION  \
                    SELECT C_SECOND.NETWORK, C_SECOND.SESSION, C_SECOND.STATION_NAME, C_SECOND.StartDateTime, C_SECOND.FinishDateTime \
                    FROM( \
                    SELECT MARKS_IN_SESS.NETWORK, MARKS_IN_SESS.SESSION, BASE_IN_SESS.SECOND AS STATION_NAME, BASE_IN_SESS.StartDateTime, BASE_IN_SESS.FinishDateTime \
                    FROM (SELECT GNSS_SESSIONS.NETWORK, GNSS_SESSIONS.SESSION, GNSS_SESSIONS.STATION_NAME \
                    FROM GNSS_SESSIONS \
                    WHERE (((GNSS_SESSIONS.NETWORK)="+ str(network_num) + ") AND ((GNSS_SESSIONS.SESSION)='"+ e[0] + "')) \
                    )  AS MARKS_IN_SESS INNER JOIN (SELECT OBSERVATIONS.FIRST, OBSERVATIONS.SECOND, OBSERVATIONS.StartDateTime, OBSERVATIONS.FinishDateTime \
                    FROM OBSERVATIONS \
                    WHERE (((OBSERVATIONS.TYPE)='G') AND ((OBSERVATIONS.StartDateTime)<='"+ e[0] + "') AND ((OBSERVATIONS.FinishDateTime)>='"+ e[0] + "')) \
                    )  AS BASE_IN_SESS ON MARKS_IN_SESS.STATION_NAME = BASE_IN_SESS.FIRST) AS C_SECOND)  AS C LEFT JOIN (SELECT GNSS_SESSIONS.STATION_NAME \
                    FROM GNSS_SESSIONS \
                    WHERE (((GNSS_SESSIONS.NETWORK)="+ str(network_num) + ") AND ((GNSS_SESSIONS.SESSION)='"+ e[0] + "')) \
                    )  AS MARKS_IN_SESS ON C.STATION_NAME = MARKS_IN_SESS.STATION_NAME \
                    WHERE (((MARKS_IN_SESS.STATION_NAME) Is Null))"
                    cursor.execute(sqlstring)
                    conn.commit()
        
        ##Print the trivials to kml
        sqlstring="SELECT DISTINCT NETWORK, SESSION \
            FROM GNSS_SESSIONS \
            ORDER BY SESSION, NETWORK"
        session=cursor.execute(sqlstring).fetchall()
        for s in session:
            baseCntSql = "SELECT COUNT(OBSERVATIONS.ID) AS BASE_CNT \
                FROM (GNSS_SESSIONS INNER JOIN OBSERVATIONS ON GNSS_SESSIONS.STATION_NAME = OBSERVATIONS.FIRST) \
                INNER JOIN GNSS_SESSIONS AS GNSS_SESSIONS_1 ON OBSERVATIONS.SECOND = GNSS_SESSIONS_1.STATION_NAME \
                WHERE (((OBSERVATIONS.TYPE)='G') AND ((GNSS_SESSIONS.NETWORK)="+ str(s[0]) + ") AND ((GNSS_SESSIONS.SESSION)='"+ str(s[1]) + "')AND ((GNSS_SESSIONS_1.NETWORK)="+ str(s[0]) + ") AND ((GNSS_SESSIONS_1.SESSION)='"+ str(s[1]) + "'))"
            stationCntSql = "SELECT COUNT(ID) AS STATION_CNT \
                FROM GNSS_SESSIONS \
                WHERE ((GNSS_SESSIONS.NETWORK)="+ str(s[0]) + ") AND ((GNSS_SESSIONS.SESSION)='"+ str(s[1]) + "')"
            stationCnt = cursor.execute(stationCntSql).fetchall()
            baseCnt = cursor.execute(baseCntSql).fetchall()
            if stationCnt[0][0]-1<baseCnt[0][0]:
                kmlsql="SELECT OBSERVATIONS.ID, OBSERVATIONS.FIRST, OBSERVATIONS.SECOND , STATIONS.LATITUDE, STATIONS.LONGITUDE, STATIONS_1.LATITUDE, STATIONS_1.LONGITUDE, GNSS_SESSIONS.StartDateTime, GNSS_SESSIONS.FinishDateTime \
                        FROM (((OBSERVATIONS INNER JOIN GNSS_SESSIONS ON OBSERVATIONS.FIRST = GNSS_SESSIONS.STATION_NAME) INNER JOIN GNSS_SESSIONS AS GNSS_SESSIONS_1 ON OBSERVATIONS.SECOND = GNSS_SESSIONS_1.STATION_NAME) LEFT JOIN STATIONS AS STATIONS_1 ON GNSS_SESSIONS_1.STATION_NAME = STATIONS_1.STATION_NAME) LEFT JOIN STATIONS ON GNSS_SESSIONS.STATION_NAME = STATIONS.STATION_NAME \
                        WHERE (((OBSERVATIONS.TYPE)='G') AND ((GNSS_SESSIONS.SESSION)='"+ str(s[1]) + "') AND ((GNSS_SESSIONS_1.SESSION)='"+ str(s[1]) + "') AND ((GNSS_SESSIONS.NETWORK)="+ str(s[0]) + ") AND ((GNSS_SESSIONS_1.NETWORK)="+ str(s[0]) + "))"
                kml_Trivial=cursor.execute(kmlsql).fetchall()    
                for tv in kml_Trivial:
                    w_kml.write(MkLine(tv))         
        w_kml.write(kmlFooter())
        w_kml.close()
        print('Done :)')
                                
        conn.close()

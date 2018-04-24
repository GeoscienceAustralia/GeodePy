# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 14:08:52 2018

@author: u84157
"""
## Example input
#North Latitude (S): 10
#South Latitude (S): 25
#West Longitude (E): 155
#East Longitude (E): 165
#Directory and file name of AUSGeoid2020 grid: C:\AusGeoidExtract\AUSGeoid2020_20170908_win.dat
#Output Directory and file name of sub-grid: C:\AusGeoidExtract\AUSGeoid2020_20170908_win_SUBGRID.dat
## Get user to input the extents of the file they wish to extract
Latmin=input('North Latitude (S): ')
Latmax=input('South Latitude (S): ')
Longmin=input('West Longitude (E): ')
Longmax=input('East Longitude (E): ')
## Get user to specify input/out files
filename =raw_input('Directory and file name of AUSGeoid2020 grid: ')
f=open(filename)
linestr=f.readline()
filenameout=raw_input('Output Directory and file name of sub-grid: ')
fout=open(filenameout,"w")
fout.writelines(linestr)
## Run through each line of the input file and extract the relevant lines
print("Running")
for linestr in f.readlines():
    Latk=float(linestr[14:16])+float(linestr[17:19])/60+float(linestr[17:19])/3600
    Longk=float(linestr[28:31])+float(linestr[32:35])/60+float(linestr[36:41])/3600
    if Latk>=Latmin and Latk<=Latmax and Longk>=Longmin and Longk<=Longmax   :
        fout.writelines(linestr)
## Close the files
f.close
fout.close
print("Done :)")
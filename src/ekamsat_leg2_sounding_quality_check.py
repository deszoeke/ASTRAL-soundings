# -*- coding: utf-8 -*-
"""
Created on Fri May 10 10:13:20 2024
The script currently checks the sounding dataset (.nc files) for the following errors:
1. non-numeric entry, 
2. quantities exceeding physical limits,
3. height not increasing monotonically.
4. Ascent rate exceeding 20 m/s 

Ref: Loehrer, S. M., Edmands, T. A., & Moore, J. A. (1996). TOGA COARE upper-air sounding data archive: Development and quality control procedures. Bulletin of the American Meteorological Society, 77(11), 2651-2672.

It also plots the soundings on thephigrams and calculates instability parameters: CAPE, CIN,TPW, LCL,LFC, Mixed layer height.
@author: jayes
"""

##############  Load python packages ##################
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')
###############  and functions ##################
import soundings_quality_functions as sqf

# provide path of the sounding dataset
path_in = r"D:\Astral\data\2024\Leg2\Radiosonde\netcdf\\"

file_list = os.listdir(path_in) # list of files in the folder

file_names =[]
num_errors = []
lim_errors = []
ord_errors = []
asc_errors = []
capes = []
cins  = []
pws   = []
els = []
lcls = []
lfcs = []
ts = []
max_heights =[]
min_heights =[]
mix_lay_ds  =[]

for f_name in file_list:
    os.chdir(path_in) # change working dir
    if f_name.endswith('.nc'): # read only nc files
        ncid = Dataset(f_name, 'r')
        
        # Read vriables from the sounding file
        height = ncid["msl"][:]
        P = ncid["P"][:]
        ws = ncid["ws"][:]
        wd = ncid["wd"][:]
        T = ncid["T"][:]
        Td = ncid["Td"][:]
        rh = ncid["rh"][:]
        tm = ncid["time"][:]
        Lat = ncid["lats"][:]
        Lon = ncid["lons"][:]
        asc = ncid["asc"][:]
        
        start_h = 0 # starting index
        data = np.array([height[start_h:], P[start_h:], ws[start_h:], wd[start_h:], T[start_h:], Td[start_h:], rh[start_h:], asc[start_h:]])
        df_data   = pd.DataFrame(np.transpose(data), columns = ['Height (m)','Pressure (hPa)', 'Speed (m/s)', 'Direction (deg)', 'Temperature (C)', 
                                                  'Dewpoint (C)', 'RH (%)', 'Ascent (m/min)'])
        
        # Plot soundings on a tephigram
        sqf.plot_profile(df_data, f_name, path_in, Lat, Lon)        

    
        # get error report and instability paramters
        [num_error, lim_error, ord_error,asc_error,minH, maxH,
         cape,cin,pw,el,lcl,lfc, MLD, df_data_qc] =  sqf.qc_check(df_data, f_name, path_in)
        
        qc_path = r"D:\Astral\data\2024\Leg2\Radiosonde\netcdf\qc_report\qc_soundings\\"
        df_data_qc.to_csv(qc_path + f_name.split('.nc')[0] + '.csv')
                
        num_errors.append(num_error)
        lim_errors.append(num_error)
        ord_errors.append(num_error)
        asc_errors.append(asc_error)
        capes.append(cape)
        cins.append(cin)
        pws.append(pw)
        els.append(el)
        lcls.append(lcl)
        lfcs.append(lfc)
        ts.append(tm)
        max_heights.append(maxH)
        min_heights.append(minH)
        mix_lay_ds.append(MLD)
        file_names.append(f_name)
        
        print("Processed file: "+ f_name)
        
    else:
        print("Not a netCDF file:" + f_name)
        
df_sounding_info = pd.DataFrame(list(zip(file_names,num_errors,lim_errors,ord_errors,asc_errors,
                                         min_heights,max_heights,capes,cins,pws,els,
                                         lcls,lfcs, mix_lay_ds)), 
                  columns = ['File', 'Numerical test', 'Limit test','Order test', 'Ascent test',
                             'Min Height (hPa)','Max Height (hPa)','CPAE (J/kg)', 'CIN (J/kg)','PW (mm)',
                             'EL (hPa)','LCL (hPa)','LFC (hPa)', 'Mixed layer depth (m)'])        
        

path = path_in + "\\qc_report\\"
df_sounding_info.to_csv(path + "qc_report.csv") # Save report

print("Sounding quality report saved in " + path + "qc_report.csv")
#print("Sounding profile plots saved in "  + path + "\\plots\\")
        
        
        
        
        
        
           
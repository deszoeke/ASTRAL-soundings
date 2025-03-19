# -*- coding: utf-8 -*-
"""
Created on Wed May  8 07:28:22 2024

@author: jayesh
Reads the output textfiles from the iMET-4 radiosonde system and stores the data in NetCDF format (L0 files).
Uses functions defined in save_leg2_to_nc.py
"""
import numpy as np
import pandas as pd
import os
from datetime import date
from save_leg2_to_nc import save_sounding_to_nc
from metpy.units import units
import metpy.calc as mpcalc


input_text_file_path = r"C:\ASTRAL\Data\2024\Leg2\Radiosonde\text_files\\"
output_netcdf_file_path = r'C:\ASTRAL\Data\2024\Leg2\Radiosonde\netcdf\L0\\'


os.chdir(input_text_file_path)
# Get file names
fnames = os.listdir(input_text_file_path)
for f_name in fnames:
    if f_name.endswith('.txt'): # read only text files
        try:
            with open(f_name,"r") as f_obj:
                lines = f_obj.readlines()
                
                # get location from the text file
                Lat = float(lines[7].split(":")[1][1:6])
                Lon = float(lines[8].split(":")[1][1:7])
                
                #get date and time from the text file
                a = f_name.split("EKAMSAT")
                y = int(a[1][0:4])
                m = int(a[1][5:6])
                d = int(a[1][6:8])
                h = int(a[1][8:10])
                date_0 = date(1970,1,1)
                date_1 = date(y, m, d)
                no_seconds = 3600*(int( (date_1 - date_0).days * 24 ) + h)
                
                ti = pd.Timestamp(y, m, d, h)
                
                # sonde details from the text file
                sonde_type = lines[11][-7:-1]
                sonde_id = lines[12][-6:-1]
                
                
                # assign file/global attributes
                file_attr = {
                'project'    : 'EKAMSAT',
                'location'  : 'Bay of Bengal - Indian Ocean',
                'research_vessel' : 'Thomas G Thompson',
                'sonde_type': sonde_type, 
                'sonde_id'  : sonde_id,
                'file_title': 'EKAMSAT Leg 2 radiosonde dataset',
                'level'     : 'Level 1',
                'institue'  : 'University of Notre Dame',
                'group'     : 'Environmental Fluid Mechanics Laboratory',
                'group_url'   : 'https://efmlab.nd.edu/',
                'creator_name'   :  'Jayesh Phadtare',
                'creator_email'  :  'jayesh.phadtare@gmail.com',
                'acknowledgements': 'Office of Naval Research'
                }
                
                
                # get the table
                data = np.loadtxt(f_name, skiprows=28)
                
                # Table to pandas DataFrame
                df   = pd.DataFrame(data, columns = ['Time','Height (GP)','Height (MSL)','Lon', 'Lat','Pressure', 
                                                     'Temperature' , 'Dewpoint', 'Mixing ratio', 'RH', 'Density', 
                                                     'Speed', 'Direction', 'vwind', 'uwind', 'Ascent'])
                
                
                p = df["Pressure"]
                T = df["Temperature"]
                Td = df["Dewpoint"]
                q = df["Mixing ratio"]
                
                p = units.Quantity(np.array(p),"hPa")
                T = units.Quantity(np.array(T),"degC")
                Td = units.Quantity(np.array(Td),"degC")
                q = units.Quantity(np.array(q),"")
                
                df["Theta"]  = np.round((mpcalc.potential_temperature(p, T).magnitude),1)
                df["Thetav"] = np.round((mpcalc.virtual_potential_temperature(p, T, 0.001*q).magnitude),1)
                df["Thetae"] = np.round((mpcalc.equivalent_potential_temperature(p, T, Td).magnitude),1)
                
                save_sounding_to_nc(df,no_seconds,f_name,output_netcdf_file_path,file_attr)
    
     
        except FileNotFoundError:
            data = "Not found"
                
                
                
                
                

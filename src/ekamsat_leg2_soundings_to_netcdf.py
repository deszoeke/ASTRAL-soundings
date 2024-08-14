# -*- coding: utf-8 -*-
"""
Created on Wed May  8 07:28:22 2024

@author: jayesh
"""
import numpy as np
import pandas as pd
import os
from datetime import date
from save_leg2_to_nc import save_sounding_to_nc


path = r"D:\Astral\data\2024\Leg2\Radiosonde\\"


os.chdir(path)
# Get file names
fnames = os.listdir(path)

for f_name in fnames:
    if f_name.endswith('.txt'): # read only text files
        try:
            with open(f_name,"r") as f_obj:
                lines = f_obj.readlines()
                Lat = float(lines[7].split(":")[1][1:6])
                Lon = float(lines[8].split(":")[1][1:7])
               
                a = f_name.split("EKAMSAT")
                y = int(a[1][0:4])
                m = int(a[1][5:6])
                d = int(a[1][6:8])
                h = int(a[1][8:10])
                date_0 = date(1900,1,1)
                date_1 = date(y, m, d)
                no_hours = (date_1 - date_0).days * 24 + h
                
                ti = pd.Timestamp(y, m, d, h)
                
                data = f_obj.read()
                
                # get the table
                data = np.loadtxt(f_name, skiprows=28)
                
                # Table to pandas DataFrame
                df   = pd.DataFrame(data, columns = ['Time','Height (GP)','Height (MSL)','Lon', 'Lat','Pressure', 
                                                     'Temperature' , 'Dewpoint', 'Mixing ratio', 'RH', 'Density', 
                                                     'Speed', 'Direction', 'uwind', 'vwind', 'Ascent'])
                
                save_sounding_to_nc(df,no_hours,f_name)
    
     
        except FileNotFoundError:
            data = "Not found"
                
                
                
                
                
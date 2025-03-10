# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:58:00 2025
@author: jphadtar@nd.edu
"""

import xarray as xr;
import os
import glob
import pandas as pd

def L0_to_L1(L0_path, L1_path):
    # file_list = os.listdir(L0_path)
    file_list = [os.path.basename(f) for f in glob.glob(os.path.join(L0_path,"EKAMSAT*.nc"))]
    for file in file_list:
        
        print('Reading ', os.path.join(L0_path, file) + "........")
        
        ds = xr.open_dataset(os.path.join(L0_path, file))
        
        ###########################################################
        ############## Replace non numeric entries by NaN #########
        
        ds_num = ds.map(lambda x: xr.apply_ufunc(pd.to_numeric, x,  kwargs={"errors": "coerce"}))
        for var in ds.data_vars:
            ds[var].values = ds_num[var].values
            
        ###########################################################
        ############## Sort by pressure #########
        
       
        ds = ds.sortby('pressure',  ascending=False)
 
        
        ###########################################################
        ############## Check variable limits and assign flags ######################
    
        #pressure
        pressure_flag = (ds_num.pressure > 1030) & (ds_num.pressure < 0)
        ds["pressure_flag"] = (("pressure"), pressure_flag.data)
    
        #geopotential height
        gp_height_flag = (ds_num.gp_height > 40000) & (ds_num.gp_height < 0)
        ds["gp_height_flag"] = (("pressure"), gp_height_flag.data)

    
        #temperature
        T_flag = (ds_num.T > 55) & (ds_num.T < -100)
        ds["T_flag"] = (("pressure"), T_flag.data)
    
        
        #humidity
        rh_flag = (ds_num.rh > 100) & (ds_num.rh < 0)
        ds["rh_flag"] = (("pressure"), rh_flag.data)

    
        #wind speed
        wspd_flag = (ds_num.ws > 70) & (ds_num.ws < 0)
        ds["wspd_flag"] = (("pressure"), wspd_flag.data)
    
        #wind direction
        wdir_flag = (ds_num.wd > 360) & (ds_num.wd < 0)
        ds["wdir_flag"] = (("pressure"), wdir_flag.data)
    
        #Ascent speed
    
        L1_filename = file.split(".nc")[0] + "_L1.nc"

        ds.to_netcdf(os.path.join(L1_path, L1_filename)) # store as L1 files
        ds.close()
        ds_num.close()
        print('saved to ', os.path.join(L1_path, L1_filename))

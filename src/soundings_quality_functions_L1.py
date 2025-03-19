# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:58:00 2025
@author: jphadtar@nd.edu
"""

import xarray as xr;
import os
import pandas as pd

def L0_to_L1(L0_path, L1_path):
    file_list = os.listdir(L0_path)
    for file in file_list:
        
        print('Reading ', L0_path + file + "........")
        
        ds = xr.open_dataset(L0_path + file)
        
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
        pressure_flag = (ds_num.pressure < 1030) & (ds_num.pressure > 0)
        ds["pressure_flag"] = (("pressure"), pressure_flag.data)
        ds["pressure_flag"].attrs = {   "long_name": "Pressure flag for data quality",
                                        "description": "0: value not within physical limits; 1: value witihin physical limits",
                                    }
    
        #geopotential height
        gp_height_flag = (ds_num.gp_height < 40000) & (ds_num.gp_height > 0)
        ds["gp_height_flag"] = (("pressure"), gp_height_flag.data)
        ds["gp_height_flag"].attrs = {   "long_name": "Geopotential height flag for data quality",
                                         "description": "0: value not within physical limits; 1: value witihin physical limits",
                                     }
    
        #temperature
        T_flag = (ds_num.T < 40) & (ds_num.T > -100)
        ds["T_flag"] = (("pressure"), T_flag.data)
        ds["T_flag"].attrs = {        "long_name": "Temperature flag for data quality",
                                      "description": "0: value not within physical limits; 1: value witihin physical limits",
                                  }
        
        #humidity
        rh_flag = (ds_num.rh <= 100) & (ds_num.rh >= 0)
        ds["rh_flag"] = (("pressure"), rh_flag.data)
        ds["rh_flag"].attrs = {        "long_name": "RH flag for data quality",
                                       "description": "0: value not within physical limits; 1: value witihin physical limits",
                             }
    
        #wind speed
        wspd_flag = (ds_num.ws < 70) & (ds_num.ws >= 0)
        ds["wspd_flag"] = (("pressure"), wspd_flag.data)
        ds["wspd_flag"].attrs = {      "long_name": "wind speed flag for quality",
                                       "description": "0: value not within physical limits; 1: value witihin physical limits",
                                   }
        #wind direction
        wdir_flag = (ds_num.wd < 360) & (ds_num.wd >= 0)
        ds["wdir_flag"] = (("pressure"), wdir_flag.data)
        ds["wdir_flag"].attrs = {        "long_name": "Wind direction flag for quality",
                                         "description": "0: value not within physical limits; 1: value witihin physical limits",
                                  }
        #Ascent speed
        # asc_flag = ds_num.asc < 500
        # ds["asc_flag"] = (("pressure"), asc_flag.data)
        # ds["asc_flag"].attrs = {        "long_name": "radiosonde ascent speed flag for quality",
        #                                  "description": "0: value not within physical limits (> 500 m/min); 1: value witihin physical limits",
        #                           }
        
        # Attributes
        ds.attrs['processing_level']  =  'Level 1'
        
        encoding = {var: {"_FillValue": None} for var in ds.data_vars}

    
        L1_filename = file.split("L0")[0] + "L1.nc"

        ds.to_netcdf(L1_path + L1_filename, encoding=encoding) # store as L1 files
        ds.close()
        ds_num.close()
        print('saved to ', L1_path + L1_filename)
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 08:11:30 2024

@author: jayes
"""
from netCDF4 import Dataset    # Note: python is case-sensitive!
import numpy as np

def save_sounding_to_nc(df,no_hours,fname):
    try: ncfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
     
    
    file_name = fname.split(".txt")[0] + ".nc"

    ncfile = Dataset(r'D:\Astral\data\2024\Leg2\Radiosonde\netcdf\\' + file_name ,mode='w',format='NETCDF4_CLASSIC') 
    
    ################# Create dimension ###############
    height_dim = ncfile.createDimension('msl', len(df))
    ncfile.title='University of Notre Dame: EKAMSAT Leg 2 (R/V TGT)'
    ncfile.subtitle= 'Radiosonde system: InterMet-4'

    ################ Create variables ################
    time = ncfile.createVariable('time', np.int)
    time.units = 'hours since 1900-01-01'
    time.long_name = 'time'
    
    gp_height = ncfile.createVariable('gp_height', np.int, ('msl',))
    gp_height.units = 'meter'
    gp_height.long_name = 'geopotential height'
    
    msl = ncfile.createVariable('msl', np.int, ('msl',))
    msl.units = 'meter'
    msl.long_name = 'mean sea level height'
    
    lats = ncfile.createVariable('lats', np.float32, ('msl',))
    lats.units = 'degrees_north'
    lats.long_name = 'latitude'
    
    lons = ncfile.createVariable('lons', np.float32, ('msl',))
    lons.units = 'degrees_east'
    lons.long_name = 'longitude'
    
    
    P = ncfile.createVariable('P',np.float32,('msl')) # note: unlimited dimension is leftmost
    P.units = 'hPa' # 
    P.standard_name = 'pressure' # this is a CF standard name
    
    T = ncfile.createVariable('T',np.float32,('msl')) # note: unlimited dimension is leftmost
    T.units = 'C' # 
    T.standard_name = 'air temperature' # this is a CF standard name
    
    Td = ncfile.createVariable('Td',np.float32,('msl')) # note: unlimited dimension is leftmost
    Td.units = 'C' # 
    Td.standard_name = 'dewpoint' # this is a CF standard name
    
    
    q = ncfile.createVariable('q',np.float32,('msl')) # note: unlimited dimension is leftmost
    q.units = 'g/kg' # 
    Td.standard_name = 'dewpoint' # this is a CF standard name
    
    
    rh = ncfile.createVariable('rh',np.float32,('msl')) # note: unlimited dimension is leftmost
    rh.units = '%' # 
    rh.standard_name = 'relative humidity' # this is a CF standard name
    
    rho = ncfile.createVariable('rho',np.float32,('msl')) # note: unlimited dimension is leftmost
    rho.units = 'g/m**3' # 
    rho.standard_name = 'density' # this is a CF standard name
    
    ws = ncfile.createVariable('ws',np.float32,('msl')) # note: unlimited dimension is leftmost
    ws.units = 'm/s' # 
    ws.standard_name = 'wind speed' # this is a CF standard name 
    
    wd = ncfile.createVariable('wd',np.float32,('msl')) # note: unlimited dimension is leftmost
    wd.units = 'degrees from the North' # 
    wd.standard_name = 'wind direction' # this is a CF standard name
    
    
    v = ncfile.createVariable('v',np.float32,('msl')) # note: unlimited dimension is leftmost
    v.units = 'm/s' # 
    v.standard_name = 'meridional wind speed' # this is a CF standard name
    
    u = ncfile.createVariable('u',np.float32,('msl')) # note: unlimited dimension is leftmost
    u.units = 'm/s' # 
    u.standard_name = 'zonal wind speed' # this is a CF standard name
    
    asc = ncfile.createVariable('asc',np.float32,('msl')) # note: unlimited dimension is leftmost
    asc.units = 'm/min' # 
    asc.standard_name = 'radiosonde ascent speed' # this is a CF standard name
    
    
    time[:] = no_hours
    gp_height[:] = df['Height (GP)']
    msl[:] = df['Height (MSL)']
    lats[:] = df['Lat'] 
    lons[:] = df['Lon']   
    P[:] = df['Pressure']
    T[:]  = df['Temperature']
    Td[:] = df['Dewpoint']
    q[:]  = df['Mixing ratio']
    rh[:] = df['RH']
    rho[:] = df['Density']
    ws[:] = df['Speed']
    wd[:] = df['Direction']
    u[:]  = df['uwind']
    v[:]  = df['vwind']
    asc[:] = df['Ascent']
       
    
    # close the Dataset.
    print("File created: " + file_name)
    ncfile.close();
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
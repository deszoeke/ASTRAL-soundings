# -*- coding: utf-8 -*-
"""
Created on Wed May  8 08:11:30 2024
@author: jphadtar@nd.edu
Creates a NetCDF file (L0) from the pandas DataFrame (df) passed from ekamsat_leg2_soundings_to_netcdf.py.
"""
from netCDF4 import Dataset    # Note: python is case-sensitive!
from datetime import datetime

def save_sounding_to_nc(df,no_seconds,fname, path,file_attr):
       
    project_name =  fname.split("_HSPOTINT.txt")[0][0:7]
    launch_date  =  fname.split("_HSPOTINT.txt")[0][7:15]
    launch_time  =  fname.split("_HSPOTINT.txt")[0][15:20]
    file_name = project_name + "_RAOB_TGT_" + launch_date + "_" + launch_time + "_L0.nc"

    ncfile = Dataset(path + file_name ,mode='w',format='NETCDF4_CLASSIC')
    
    ################# Create dimension ##############
    
    pressure = ncfile.createDimension('pressure', len(df))  # pressure as vertical coordinate
    time =     ncfile.createDimension('time', None)

    ################# Create attributes ##############
    ncfile.title = file_attr["file_title"]
    ncfile.processing_level = file_attr["level"]
    ncfile.project = file_attr["project"]
    ncfile.location = file_attr["location"]
    ncfile.research_vessel = file_attr["research_vessel"]
    ncfile.sonde_type = file_attr["sonde_type"]
    ncfile.sonde_id = file_attr["sonde_id"]
    ncfile.institue = file_attr['institue']
    ncfile.group = file_attr["group"]
    ncfile.group_url = file_attr["group_url"]
    ncfile.creator_name = file_attr["creator_name"]
    ncfile.creator_email = file_attr["creator_email"]
    ncfile.date_of_creation = datetime.now().strftime("%Y-%m-%d") + " (YYYY-mm-dd)"
    ncfile.acknowledgments = file_attr["acknowledgements"]
    
    
    ################ Create variables ################
    
    time = ncfile.createVariable('time', 'f4',  ('time',))
    time.units = 'seconds since 1970-01-01T00:00:00Z'
    time.standard_name = 'time'
    time.long_name = 'time of radiosonde release in UTC'
        
    
    flight_time = ncfile.createVariable('flight_time', 'f4', ('pressure',))
    flight_time.units = 'seconds'
    flight_time.long_name = 'time since radiosonde release'
    
    
    gp_height = ncfile.createVariable('gp_height', 'f4', ('pressure',))
    gp_height.units = 'meter'
    gp_height.long_name = 'geopotential height'
    gp_height.standard_name = 'geopotential_height'
    gp_height.comment = 'variable output from InterMet system'
    
    msl_height = ncfile.createVariable('msl_height', 'f4', ('pressure',))
    msl_height.units = 'meter'
    msl_height.long_name = 'height above mean sea level'
    msl_height.standard_name = 'height_above_mean_sea_level'
    msl_height.comment = 'variable output from InterMet system'
    
    lats = ncfile.createVariable('lats', 'f4', ('pressure',))
    lats.units = 'degrees_north'
    lats.long_name = 'latitude'
    lats.standard_name = 'latitude'
    lats.comment = 'variable output from InterMet system'
    
    lons = ncfile.createVariable('lons', 'f4', ('pressure',))
    lons.units = 'degrees_east'
    lats.long_name = 'longitude'
    lons.standard_name = 'longitude'
    lons.comment = 'variable output from InterMet system'
    
    
    pressure = ncfile.createVariable('pressure','f4',('pressure',))
    pressure.units = 'hPa'
    pressure.long_name = 'air pressure'
    pressure.standard_name = 'air_pressure'
    pressure.comment = 'variable output from InterMet system'
    
    T = ncfile.createVariable('T','f4',('pressure',))
    T.units = 'C'
    T.long_name = 'air temperature'
    T.standard_name = 'air_temperature'
    T.comment = 'variable output from InterMet system'
    
    Td = ncfile.createVariable('Td','f4',('pressure',)) 
    Td.units = 'C'
    Td.long_name = 'dew point temperature'
    Td.standard_name = 'dew_point_temperature'
    Td.comment = 'variable output from InterMet system'
    
    q = ncfile.createVariable('q','f4',('pressure',))
    q.units = 'g/kg'
    q.long_name = 'humidity mixing ratio'
    q.standard_name = 'humidity_mixing_ratio'
    q.comment = 'variable output from InterMet system'
    
    
    rh = ncfile.createVariable('rh','f4',('pressure',))
    rh.units = '%'
    rh.long_name = 'relative humidity'
    rh.standard_name = 'relative_humidity'
    rh.comment = 'variable output from InterMet system'
    
    rho = ncfile.createVariable('rho','f4',('pressure',))
    rho.units = 'g/m**3'
    rho.long_name = 'air_density'
    rho.standard_name = 'air_density'
    rho.comment = 'variable output from InterMet system'
    
    
    theta = ncfile.createVariable('theta','f4',('pressure',))
    theta.units = 'kelvin'
    theta.long_name = 'potential_temperature'
    theta.standard_name = 'air_potential_temperature'
    theta.comment = 'derived using MetPy 1.6.0'
    
    thetav = ncfile.createVariable('thetav','f4',('pressure',))
    thetav.units = 'kelvin'
    thetav.long_name = 'virtual potential temperature'
    thetav.standard_name = 'air_virtual_potential_temperature'
    thetav.comment = 'derived using MetPy 1.6.0'
    
    thetae = ncfile.createVariable('thetae','f4',('pressure',))
    thetae.units = 'kelvin' 
    thetae.long_name = 'air equivalent potential temperature'
    thetae.standard_name = 'air_equivalent_potential_temperature'
    thetae.comment = 'derived using MetPy 1.6.0'
    
    ws = ncfile.createVariable('ws','f4',('pressure',)) 
    ws.units = 'm/s' 
    ws.standard_name = 'wind_speed'
    ws.long_name = 'wind speed'
    ws.comment = 'variable output from InterMet system'
    
    wd = ncfile.createVariable('wd','f4',('pressure',))
    wd.units = 'degree'
    wd.long_name = 'wind direction'
    wd.standard_name = 'wind_from_direction'
    wd.comment = 'variable output from InterMet system'
    
    
    v = ncfile.createVariable('v','f4',('pressure',))
    v.units = 'm/s'
    v.long_name = 'northward wind velocity'
    v.standard_name = 'northward_wind'
    v.comment = 'variable output from InterMet system'
    
    u = ncfile.createVariable('u','f4',('pressure',))
    u.units = 'm/s' 
    u.long_name = 'eastward wind velocity'
    u.standard_name = 'eastward_wind'
    u.comment = 'variable output from InterMet system'
    
    asc = ncfile.createVariable('asc','f4',('pressure',))
    asc.units = 'm/min' 
    asc.long_name = 'radiosonde ascent speed'
    asc.comment = 'variable output from InterMet system'
    
    
    time[:] = no_seconds
    flight_time[:] = df['Time'].values
    gp_height[:] = df['Height (GP)'].values
    msl_height[:] = df['Height (MSL)'].values
    lats[:] = df['Lat'].values
    lons[:] = df['Lon'].values   
    pressure[:] = df['Pressure'].values
    T[:]  = df['Temperature'].values
    Td[:] = df['Dewpoint'].values
    theta[:] = df["Theta"].values
    thetav[:] = df["Thetav"].values
    thetae[:] = df["Thetae"].values
    q[:]  = df['Mixing ratio'].values
    rh[:] = df['RH'].values
    rho[:] = df['Density'].values
    ws[:] = df['Speed'].values
    wd[:] = df['Direction'].values
    u[:]  = df['uwind'].values
    v[:]  = df['vwind'].values
    asc[:] = df['Ascent'].values
       
    # close the Dataset.
    print("File created: " + file_name)
    ncfile.close();
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
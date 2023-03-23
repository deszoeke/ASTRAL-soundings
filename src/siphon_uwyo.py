# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Python methods and downloads of soundings from UWyo and IGRA2 to NetCDF
#
# NCEI landing page: 
# https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00975

# %% [markdown]
# To make a kernel with a siphon virtualenv from the terminal:
#
# ```
# conda create -n siphon siphon
# conda activate siphon
# conda install ipykernel
# python -m ipykernel install --user --name=siphon
# ```
#
# Then select the `siphon` from the pulldown/window in Jupyter.

# %%
from siphon.simplewebservice.wyoming import WyomingUpperAir
from siphon.simplewebservice.igra2 import IGRAUpperAir
from datetime import datetime
from netCDF4 import Dataset
import os
import pandas

start_time = datetime(2023, 3, 1, 0)
end_time = datetime(2023, 6, 30, 0)

station = "43371" # Trivandrum, Thiruvananthapuram, 'VOTX'?

df = WyomingUpperAir.request_data(start_time, station) # returns Pandas dataframe

# %%
# test stuff

# for column in df:
#     print( df[column].name )
#     print( df[column].values )

# df[column].shape[0]
df[df.columns[0:8]]

# skip column 8

# scalars


df.columns[9]
type(df[df.columns[9]]) is pandas.core.series.Series

df[df.columns[10]][0].strftime('%Y-%m-%d %H:%M:%S') # convert Datestamp to string
type(df[df.columns[10]][0]) is pandas._libs.tslibs.timestamps.Timestamp


# %%
"convert Timestamps into string, otherwise just return input."
def Timestamp2String( t ):
    return ( t.strftime('%Y-%m-%d %H:%M:%S') if type(t) is pandas._libs.tslibs.timestamps.Timestamp else t )

"write a sounding dataframe as a NetCDF4 file."
def sounding2nc( df, filename ):

    # Open a new NetCDF file for writing
    nc_file = Dataset(filename, mode='w')

    # Set the attributes from df cols 9...
    for column in df.columns[9:]:
        #  print( "set attribute "+df[column].name+"=", Timestamp2String(df[column][0]) )
        nc_file.setncattr( df[column].name, Timestamp2String(df[column][0]) )

    # or set all at once from a dict: setncatts(self,attdict)

    # Define dimensions for pressure, temperature, etc. variables
    time_dim = nc_file.createDimension('time', 1)
    level_dim = nc_file.createDimension('level', df[df.columns[0]].shape[0])

    units_list = [ 'hPa', 'm', 'degrees C', 'degrees C', 'degrees', 'm/s', 'm/s', 'm/s' ]

    # Create variables for time, level, latitude, longitude, pressure, and temperature
    for i, column in enumerate( df.columns[0:7] ):
        # create the netcdf variable
        var = nc_file.createVariable(df[column].name, 'f4', ('time', 'level'))
        # write the data to the NetCDF file
        var[:] = df[column][:]
        # supply units attributes
        var.units = units_list[i]

    # Close the NetCDF file
    return nc_file.close()


# %%
# get Trivandrum soundings from 2019, write to netcdf

start_time = datetime(2019, 3, 1, 0)
end_time = datetime(2019, 6, 30, 0)
station = "43371" # Trivandrum, Thiruvananthapuram, 'VOTX'?

for dt in pandas.date_range(start_time, end_time, freq='12H'):
    try:
        df = WyomingUpperAir.request_data(dt, station) # returns Pandas dataframe
        sounding2nc( df, "../data/uwyo/trivandrum/trivandrum"+df.time[0].strftime('%Y%m%d_%H%M')+".nc" )
    except:
        continue
    else:
        continue
        



# %%
start_time = datetime(2019, 3, 1, 0)
end_time = datetime(2019, 6, 30, 0)
station = "43063"  # Pune
for dt in pandas.date_range(start_time, end_time, freq='12H'):
    try:
        df = WyomingUpperAir.request_data(dt, station) # returns Pandas dataframe
        sounding2nc( df, "../data/uwyo/pune/pune"+df.time[0].strftime('%Y%m%d_%H%M')+".nc" )
    except:
        continue
    finally:
        continue

# %%
type(df.columns) # pandas.core.indexes.base.Index
# var_indx = pandas.core.indexes.base.Index(['pressure', 'height', 'temperature', 'dewpoint'])
# header
# df
units_list = [ 'none', 'none', 'seconds?', 'hPa', 'none', 'meters', 'none', 'degrees C', 'none', 'percent', 'degrees', 'm/s', 'datestring', 'm/s', 'm/s', 'degree C']
    
# df.drop(columns='date')
# for i, column in enumerate(df.columns):
#     print( i, column, units_list[i] )

# dfd = df.drop(columns='date')
# for i, column in enumerate(dfd.columns):
#     print( i, column, units_list[i] )

# %%
"write an IGRA2 sounding dataframe as a NetCDF4 file."
def igra2nc( df, header, filename ):

    # Open a new NetCDF file for writing
    nc_file = Dataset(filename, mode='w')

    # Set the attributes from header
    for column in header.columns:
        nc_file.setncattr( header[column].name, Timestamp2String(header[column][0]) )

    # or set all at once from a dict: setncatts(self,attdict)

    # Define dimensions for pressure, temperature, etc. variables
    time_dim = nc_file.createDimension('time', 1)
    level_dim = nc_file.createDimension('level', df[df.columns[0]].shape[0])

    dfd = df.drop(columns='date') # Datestamp redundant and type not allowed in NetCDF
    units_list = [ 'none', 'none', 'seconds?', 'hPa', 'none', 'meters', 'none', 'degrees C', 'none', 'percent', 'degrees', 'm/s', #'datestring', 
                   'm/s', 'm/s', 'degree C']

    # Create variables for time, level, latitude, longitude, pressure, and temperature
    for i, column in enumerate( dfd.columns ):
        # create the netcdf variable
        var = nc_file.createVariable(dfd[column].name, 'f4', ('time', 'level'))
        # write the data to the NetCDF file
        var[:] = dfd[column][:]
        # supply units attributes
        var.units = units_list[i]

    # Close the NetCDF file
    return nc_file.close()

def get_igra2_station( daterange, station, stationname ):
    for dt in daterange: # freq='12H'
        my_path = "../data/igra2/"+stationname+"/"+stationname+dt.strftime('%Y%m%d_%H%M')+".nc"
        if not os.path.exists(my_path) or os.path.getsize(my_path) <= 0:
            try:
                df, header = IGRAUpperAir.request_data(dt.to_pydatetime(), station)
                igra2nc( df, header, my_path )
            except:
                # print exception!
                continue
            finally:
                continue
        else:
            continue


# %%
# test get_igra2_station

stationname = 'pune'
station = 'INM00043063'
start_time = datetime(2019, 3, 1, 0)
end_time = datetime(2019, 6, 30, 0)
dr = pandas.date_range(start_time, end_time, freq='D')
get_igra2_station( dr, station, stationname )

stationname = 'trivandrum'
station = 'INM00043371'
start_time = datetime(2019, 3, 1, 0)
end_time = datetime(2019, 6, 30, 0)
dr = pandas.date_range(start_time, end_time, freq='D')
get_igra2_station( dr, station, stationname )

# %%
# batch rename Trivandrum files
for filename in os.listdir("../data/igra2/thrivandrum/"):
    if filename.startswith("thrivandrum"):
        os.rename("../data/igra2/thrivandrum/"+filename, "../data/igra2/thrivandrum/"+"t"+filename[2:])

# %%
# batch rename Goa -> goa files
for filename in os.listdir("../data/igra2/goa/"):
    if filename.startswith("Goa"):
        os.rename("../data/igra2/goa/"+filename, "../data/igra2/goa/"+"goa"+filename[4:])

# %%
# load all the data for all stations

# list of all stations (lists and tuples) -- try to be pythonic!
X = [ ("INM00043192", "goa"), ("INM00043371", "trivandrum"), ("INM00043063", "pune"), ("CEM00043466", "colombo") ]
stationid = [ x[0] for x in X ]
stationname = [ x[1] for x in X ]
station_dict = dict( zip(stationname, stationid) )
# station_dict.keys()
station_dict.keys()

# %%
# load all the data for all stations
for year in [2018, 2019, 2020, 2021, 2022]:
    dr = pandas.date_range(datetime(year, 3, 1, 0), datetime(year, 6, 30, 0), freq='D')
    for stationname in station_dict.keys():
        get_igra2_station( dr, station_dict[stationname], stationname )

# %%
# get Colombo soundings
for year in [2018, 2019, 2020, 2021, 2022]:
    dr = pandas.date_range(datetime(year, 3, 1, 0), datetime(year, 6, 30, 0), freq='D')
    get_igra2_station( dr, station_dict["colombo"], "colombo" )

# %%

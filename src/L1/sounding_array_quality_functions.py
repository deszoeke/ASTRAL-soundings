# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 14:15:12 2022
This file contains the functions for the following:  
1. Checking the quality of sounding dataset: non-numeric entries, out of range values, sounding descent, incomplete sounding
2. Calculate the thermodynamics instability parameters: CAPE, CIN,TPW, LCL,LFC, Mixed layer height.
3. Plot and save tephigrams.
4. Quality check reports and tephigrams are saved in the 'qc_report' subdirectory in the of datset directory


@author: Jayesh Phadtare (jayesh.phadtare@gmail.com)
Modified by Alex Kinsella to work with the EKAMSAT24 sounding array
"""
import numpy as np
import pandas as pd
import scipy.signal
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import os
from metpy.units import pandas_dataframe_to_unit_arrays
import tephi


# Subroutine for calculating mixed layer height
def get_ml(v,d,h):
    
    ml_found = 0
    
    for i in range(1,len(v),1):
        vm = 0
        rhom = 0
        
        for j in range(0,i,1):
            rhom = np.sum(d[0:j])
            vm   = np.sum(d[0:j]*v[0:j])
                        
        vm = vm/rhom
        
        if ((v[i] - vm)> 0.2):
            ml_found = 1
            break
    if ml_found == 1:  
        return h[i]
    else:
        return -9999


# Subroutine for calculating instability paramters
def get_instability_paramters(df):
    print("Calculating instability_paramters")
    
    df.units =    {'Height (m)' : 'meter',
                         'Pressure (hPa)': 'hPa',
                         'Speed (m/s)': 'knot',
                         'Direction (deg)': 'degrees',
                         'Temperature (C)': 'degC',
                         'Dewpoint (C)' : 'degC',
                         'RH (%)': '',
                         'Ascent (m/min)': ''
                         }

    data = pandas_dataframe_to_unit_arrays(df)

    p_filter = apply_median_filter(data["Pressure (hPa)"])
    
    try:                     
        cape,cin = mpcalc.most_unstable_cape_cin(p_filter, data["Temperature (C)"], data["Dewpoint (C)"])
    except:
        cape = -9999
        cin = -9999
    pw = mpcalc.precipitable_water(p_filter,  data["Dewpoint (C)"])
    eql = mpcalc.el(p_filter, data["Temperature (C)"], data["Dewpoint (C)"])
    lcl = mpcalc.lcl(p_filter, data["Temperature (C)"], data["Dewpoint (C)"])
    
    try:
        lfc = mpcalc.lfc(p_filter, data["Temperature (C)"], data["Dewpoint (C)"])
    except:
        lfc = -9999
    
    q   = mpcalc.specific_humidity_from_dewpoint(p_filter, data["Dewpoint (C)"])
    

    thv = mpcalc.virtual_potential_temperature(p_filter, data["Temperature (C)"], q)
    rho = mpcalc.density(p_filter, data["Temperature (C)"], q)
    
    mld = get_ml(thv.magnitude, rho.magnitude, data["Height (m)"].magnitude)
    
    try:
        return cape.magnitude ,cin.magnitude,pw.magnitude,eql[0].magnitude,lcl[0][0].magnitude,lfc, mld
    except:
        return cape ,cin,pw.magnitude,eql[0].magnitude,lcl[0][0].magnitude,lfc, mld

        

# Subroutine for the sounding quality check

def qc_check(df, filename, path_in):
    
    num_err_found = "Pass" # Error switch for the numerical errors
    lim_err_found = "Pass" # Error switch for the limit errors
    ord_err_found = "Pass" # Error switch for the order errors
    asc_err_found = "Pass" # Error switch for rapid ascent
    the_err_found = 0

    # Initialize all meteorological parameters with default values
    cape = cin = pw = el = lcl = lfc = ml = np.nan
    
    filename = filename.split(".")[0]
    
    ###########################     QC tests: Start          ################################### 
    
    ################## Check if all entries are numerical ###############################
    
    num_test = df.apply(lambda s: pd.to_numeric(s, errors='coerce').notnull().all())
    num_test = num_test.astype("string")
    
    ################## Check if all values are within physical limit ####################
    
    # Height: 0- 30 km
    max_h = df["Height (m)"].max()
    min_h = df["Height (m)"].min()
    
    # Pressure: 0 - 1030 hPa
    max_p = df["Pressure (hPa)"].max()
    min_p = df["Pressure (hPa)"].min()
    
    # Temperature: -99 - 55 C
    max_T = df["Temperature (C)"].max()
    min_T = df["Temperature (C)"].min()
    
    # RH : 0 - 100%
    max_rh = df["RH (%)"].max()
    min_rh = df["RH (%)"].min()
    
    # Speed : 0 - 70 m/s
    max_u = df["Speed (m/s)"].max()
    min_u = df["Speed (m/s)"].min()
    
    # Direction: 0 - 360 deg
    max_d = df["Direction (deg)"].max()
    min_d = df["Direction (deg)"].min()
    
    # Ascent rate: < 20 m/s or 1200 m/min
    max_asc = df["Ascent (m/min)"].max()
    
    ########## Cehck if height/pressure values are in increasing/decreasing order ########
    
    bool_h = pd.to_numeric(df["Height (m)"]).is_monotonic_increasing 
    bool_p = df["Pressure (hPa)"].is_monotonic_decreasing
    ###########################     QC tests: End          ##############################
    
    ########################  Write data to processed file ##############################
    
    # Path for the quality checked files
    path = path_in + "/qc_report/"
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(path)
        
    os.chdir(path)
    
    ################################ Write the QC test report ###########################

    # Numerical test
    with open(filename + "_qc_report.txt", 'w') as f:        
        
        for entry, ind in zip(num_test,num_test.index):
            if entry == "False":
                temp_str = ind + " : Non numeric values detected!!!" + "\n"
                f.writelines(temp_str)
                num_err_found = "Fail"

        f.writelines("Numeric test:" + num_err_found + "\n")
        f.writelines( "\n")

        
        
    # Limit test
    with open(filename + "_qc_report.txt", 'a') as f:        
        f.writelines( "\n") 
        
        # For Altitude
        if (max_h > 40000):
            f.write("GPM AGL (m): One or more values out of range!!! Max detected: " + str(max_h) + " m \n")
            the_err_found = 1
            lim_err_found = "Fail"
            df.loc[df["Height (m)"] > 40000, "Height (m)"] = np.nan
            
        if (min_h < 0):
            f.write("GPM AGL (m): One or more values out of range!!! Min detected: " + str(min_h) + " m \n")
            the_err_found = 1
            lim_err_found = "Fail"
            df.loc[df["Height (m)"] < 0, "Height (m)"] = np.nan
        
        # For Pressure  
        if (max_p > 1030):
            f.write("Pressure (hPa): One or more values out of range!!! Max detected: " + str(max_p) + " hPa \n")
            the_err_found = 1
            lim_err_found = "Fail"
            df.loc[df["Pressure (hPa)"] > 1030, "Pressure (hPa)"] = np.nan
            
        if (min_p < 0):
            f.write("Pressure (hPa): One or more values out of range!!! Min detected: " + str(min_p) + " hPa \n")
            the_err_found = 1
            lim_err_found = "Fail"
            df.loc[df["Pressure (hPa)"] < 0, "Pressure (hPa)"] = np.nan

        
        # For Temperature    
        if (max_T > 55):
            f.write("Temperature (C): One or more values out of range!!! Max detected: " + str(max_T) + " C \n")
            the_err_found = 1
            lim_err_found = "Fail"
            df.loc[df["Temperature (C)"] > 55, "Temperature (C)"] = np.nan

        if (min_T < -99):
            f.write("Temperature (C): One or more values out of range!!! Min detected: " + str(min_T) + " C \n")
            the_err_found = 1
            lim_err_found = "Fail"
            
            df.loc[df["Temperature (C)"] < -99, "Temperature (C)"] = np.nan

            
        # For RH    
        if (max_rh > 100):
            f.write("RH (%): One or more values out of range!!! Max detected: " + str(max_rh) + " % \n")
            the_err_found = 1
            lim_err_found = "Fail"
            
            df.loc[df["RH (%)"] > 100, "RH (%)"] = np.nan

            
        if (min_rh < 0):
            f.write("RH (%): One or more values out of range!!! Min detected: " + str(min_rh) + " % \n")
            the_err_found = 1
            lim_err_found = "Fail"
            
            df.loc[df["RH (%)"] < 0, "RH (%)"] = np.nan

            
        # For wind speed    
        if (max_u > 70):
            f.write("Wind speed (m/s): One or more values out of range!!! Max detected: " + str(max_u) + "  \n")
            lim_err_found = "Fail"
            
            df.loc[df["Wind speed (m/s)"] > 70, "Wind speed (m/s)"] = np.nan

        if (min_u < 0):
            f.write("Wind speed (m/s): One or more values out of range!!! Min detected: " + str(min_u) + "  \n")
            lim_err_found = "Fail"
            
            df.loc[df["Wind speed (m/s)"] < 0, "Wind speed (m/s)"] = np.nan

            
        # For Wind deirection    
        if (max_d > 360):
            f.write("Wind direction (deg): One or more values out of range!!! Max detected: " + str(max_d) + "  \n")
            lim_err_found = "Fail"
            
            df.loc[df["Wind direction (deg)"] > 360, "Wind direction (deg)"] = np.nan

            
        if (min_d < 0):
            f.write("Wind direction (deg): One or more values out of range!!! Min detected: " + str(min_d) + "  \n")
            lim_err_found = "Fail"
            
            df.loc[df["Wind direction (deg)"] < 0, "Wind direction (deg)"] = np.nan

            
            
        # For Ascent rate
        if (max_asc > 1200):
            f.write("Ascent rate (m/min): One or more values out of range (> 1200)!!! Max detected: " + str(max_asc) + " \n")
            asc_err_found = "Fail"           
            
        f.writelines( "\n")  
        
    # Order test
    with open(filename + "_qc_report.txt", 'a') as f:
            
        if bool_h == False:
            f.write("Height (m): Not in monotonically increasing order!!! \n")
            the_err_found = 1
            ord_err_found = "Fail"
        if bool_p == False:
            f.write("Pressure (hPa): Not in monotonically decreasing order!!! \n")
            the_err_found = 1
            ord_err_found = "Fail"

        f.write("Order test:" + ord_err_found + "\n")  
        f.writelines( "\n")
        
        
    ########################### Get meteorological paramters #####################
    #if the_err_found == 0:
    if True: # Calculate for all soundings for now

        [cape,cin,pw,el,lcl,lfc, ml] = get_instability_paramters(df)
            
        with open(os.path.join(path_in, "qc_report", f"{filename}_qc_report.txt"), 'a') as f:
            f.write("Meteorological parameters:\n")
            f.write(f"Max CAPE (J/kg): {round(cape) if np.isfinite(cape) else 'N/A'}\n")
            f.write(f"Max CIN (J/kg): {round(cin) if np.isfinite(cin) else 'N/A'}\n")
            f.write(f"Precipitable water (mm): {round(pw) if np.isfinite(pw) else 'N/A'}\n")
            f.write(f"Equilibrium level (hPa): {el if np.isfinite(el) else 'N/A'}\n")
            f.write(f"Level of liquid condensation (hPa): {lcl if np.isfinite(lcl) else 'N/A'}\n")
            #f.write(f"Level of free convection (hPa): {lfc if np.isfinite(lfc) else 'N/A'}\n")
            f.write(f"Mixed layer depth (m): {ml if np.isfinite(ml) else 'N/A'}\n")
            
            # Handle lfc, which might be an array or contain units
            if lfc is not None: #and hasattr(lfc, 'magnitude'):
                lfc_value = lfc[0].magnitude if not np.isnan(lfc[0].magnitude) else 'N/A'
            else:
                lfc_value = 'N/A'
            
            f.write(f"Level of free convection (hPa): {lfc_value}\n") 
    else:

        with open(os.path.join(path_in, "qc_report", f"{filename}_qc_report.txt"), 'a') as f:
            f.write("Meteorological parameters:\n")
            f.write("Error in thermodynamic profile!!! please see the above report.")
            cape = cin = pw = el = lcl = lfc = np.nan
            ml = np.nan
    
    return  num_err_found,lim_err_found,ord_err_found,asc_err_found,max_p,min_p,cape,cin,pw,el,lcl,lfc_value, ml, df

# Subroutines for plotting profiles

def wind_tuple(s,d,l):
    wt = []
    
    for i in range(0,len(s),1):
        wt.append((s[i],d[i],l[i]))
        
    return wt

def apply_median_filter(p, initial_kernel_size=3, max_kernel_size=301):
    """
    Apply a median filter to the pressure data, starting with the smallest kernel size.
    The kernel size is increased by 2 until the filter is successful (pressure is monotonically decreasing)
    or the maximum kernel size is reached.
    
    Parameters:
    p : array-like
        The pressure data to filter.
    initial_kernel_size : int, optional
        The starting kernel size for the median filter. Default is 3.
    max_kernel_size : int, optional
        The maximum kernel size to try. Default is 29.
    
    Returns:
    p_filtered : array-like
        The filtered pressure data.
    """
    kernel_size = initial_kernel_size
    while kernel_size <= max_kernel_size:
        try:
            print(f"Trying kernel size: {kernel_size}")
            p_filtered = scipy.signal.medfilt(p.magnitude, kernel_size=kernel_size) * units.hPa
            
            # Check if the pressure is monotonically decreasing
            if np.all(np.diff(p_filtered) <= 0):
                print(f"Filtering successful with kernel size: {kernel_size}")
                return p_filtered
            
            print(f"Kernel size {kernel_size} did not produce a monotonically decreasing pressure profile.")
            
        except Exception as e:
            print(f"Kernel size {kernel_size} failed with error: {e}")
        
        # Increment the kernel size by 2
        kernel_size += 2

    # If the loop completes without finding a working kernel size, raise an error
    raise ValueError("Failed to filter pressure data with any kernel size that produces a monotonically decreasing profile.")

def plot_profile(df, filename, path_in, Lat, Lon):
    fs = 24
    fw = "normal"

    a = filename.split("_")
    y = int(a[1][0:4])
    m = int(a[1][5:6])
    d = int(a[1][6:8])
    h = int(a[1][9:10])
    ti = pd.Timestamp(y, m, d, h)
    
    df.units =    {'Height (m)' : 'meter',
                   'Pressure (hPa)': 'hPa',
                   'Speed (m/s)': 'knot',
                   'Direction (deg)': 'degrees',
                   'Temperature (C)': 'degC',
                   'Dewpoint (C)' : 'degC',
                   'RH (%)': '',
                   'Ascent (m/min)': ''
                   }

    data = pandas_dataframe_to_unit_arrays(df)
    
    T = data['Temperature (C)']
    p = data['Pressure (hPa)']
    Td = data['Dewpoint (C)']
    
    # Apply a median filter to the pressure data
    p_filtered = apply_median_filter(p)

    ## change appearance ##
    fig = plt.figure()
    plt.rcParams['font.size'] = fs
    plt.rcParams['font.weight'] = fw
    tpg = tephi.Tephigram(anchor=[(1000, -15), (50, -60)])
    tephi.ISOBAR_TEXT.update({'color': 'k', 'size': 12,'fontweight':'normal'})
    tephi.ISOBAR_LINE.update({'color': 'k','linestyle': 'solid', 'zorder':1})
    tephi.MIXING_RATIO_LINE.update({'color': 'g','zorder':1})
    tephi.MIXING_RATIO_TEXT.update({'color': 'g', 'size': 12,'fontweight':'normal', 'zorder':1})
    tephi.WET_ADIABAT_TEXT.update({'color': 'purple', 'size': 8,'fontweight':'normal'})
    tephi.WET_ADIABAT_LINE.update({'color': 'purple', 'linestyle': 'dashdot'})
    tephi.ISOBAR_SPEC = [(50, None)]

    parcel = mpcalc.parcel_profile(p_filtered, T[0], Td[0]).to('degC')
    tdry =  zip(p_filtered.magnitude, T.magnitude)
    tdew =  zip(p_filtered.magnitude, Td.magnitude)
    tpar =  zip(p_filtered.magnitude, parcel.magnitude)

    tpg.plot(tdry, color = 'red', linewidth = 2, zorder=3, label = "Dry bulb temperature")
    profile = tpg.plot(tdew, color = 'blue', linewidth = 2, linestyle='dashed', zorder=3, label = "Dew point temperature")
    par_profile = tpg.plot(tpar, color = 'magenta', linewidth = 2, linestyle='solid', zorder=3, label = "Adiabatic parcel temperature")

        # Filter out NaN values from wind data before plotting
    valid_indices = ~np.isnan(data["Speed (m/s)"].magnitude) & ~np.isnan(data["Direction (deg)"].magnitude)
    
    if np.any(valid_indices):  # Ensure there is valid data to plot
        ubarbs = wind_tuple(data["Speed (m/s)"].magnitude[valid_indices][::50],
                            data["Direction (deg)"].magnitude[valid_indices][::50],
                            p.magnitude[valid_indices][::50])
        profile.barbs(ubarbs, length=7, color='black', linewidth=2, pivot='middle', gutter=0.07)
    
    tit_tex = str(ti) + " UTC (" + str(Lon[0]) + "$^{\circ}$E," + str(Lat[0]) + "$^{\circ}$N)"

    plt.title(tit_tex, fontsize=14);

    # Save figure
    path = path_in + "/qc_report/plots/"
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(path)
        
    figname = filename.split(".nc")[0] + ".jpeg"
    plt.savefig(path + figname)
    print("Tephigram saved in: " + path + figname)
    plt.close()
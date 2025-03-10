# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 18:02:14 2025
@author: jphadtar@nd.edu
Processes sounding NetCDF files.
"""

import os
import soundings_quality_functions_L1 as sqf1

# Jayesh's path:
#L0_path = r"C:\ASTRAL\Data\2024\Leg2\Radiosonde\netcdf\L0\\"
#L1_path = r"C:\ASTRAL\Data\2024\Leg2\Radiosonde\netcdf\L1\\"

# Simon's laptop path
# /Users/deszoeks/Data/ASTRAL_2024/radiosonde/working/netcdf
# /Users/deszoeks/Data/ASTRAL_2024/radiosonde/tgt/nc
#stempath = r"../data/tgt/nc/" # ./data link in local parent directory

ncpath = os.path.join("..","data","tgt","nc") # ./data link in local parent directory
L0_path = os.path.join(ncpath, "L0")
L1_path = os.path.join(ncpath, "L1")

sqf1.L0_to_L1(L0_path,L1_path)



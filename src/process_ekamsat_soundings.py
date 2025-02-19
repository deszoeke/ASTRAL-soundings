# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 18:02:14 2025
@author: jphadtar@nd.edu
Processes sounding NetCDF files.
"""

import soundings_quality_functions_L1 as sqf1

L0_path = r"C:\ASTRAL\Data\2024\Leg2\Radiosonde\netcdf\L0\\"
L1_path = r"C:\ASTRAL\Data\2024\Leg2\Radiosonde\netcdf\L1\\"

sqf1.L0_to_L1(L0_path,L1_path)



import numpy as np
from datetime import dt

#X = np.loadtxt('rmm1974toRealtime.txt',skiprows=2,comments=('Prelim_','Final_','Missing_')
bsisopath='/Users/sdeszoek/Data/cruises/PISTON_2018/TGT/cruiseshare/indices/BSISO/'
bsisofile=bsisopath+'rmm1974toRealtime.txt'
X = np.genfromtxt(bsisofile,skip_header=2,missing_values=(1.e36, 999))
#RMM values up to "real time". For the last few days, ACCESS analyses are used instead of NCEP
# year, month, day, RMM1, RMM2, phase, amplitude.  Missing Value= 1.E36 or 999

year = (X[:,0]).astype(int)
month = (X[:,1]).astype(int)
day = (X[:,2]).astype(int)
bsiso1 = X[:,3]
bsiso2 = X[:,4]
phase = (X[:,5]).astype(int)
amplitude = X[:,6]

vdt=np.vectorize(dt.datetime)
rmmdate=vdt(year,month,day)
#rmmph=interpolate.interp1d(rmmdate,phase)(time)

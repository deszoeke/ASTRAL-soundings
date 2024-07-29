import numpy as np
from datetime import dt

#X = np.loadtxt('rmm1974toRealtime.txt',skiprows=2,comments=('Prelim_','Final_','Missing_')
rmmpath='/Users/sdeszoek/Data/cruises/PISTON_2018/TGT/cruiseshare/indices/RMM/'
rmmfile=rmmpath+'rmm1974toRealtime.txt'
X = np.genfromtxt(rmmfile,skip_header=2,missing_values=(1.e36, 999))
#RMM values up to "real time". For the last few days, ACCESS analyses are used instead of NCEP
# year, month, day, RMM1, RMM2, phase, amplitude.  Missing Value= 1.E36 or 999

year = (X[:,0]).astype(int)
month = (X[:,1]).astype(int)
day = (X[:,2]).astype(int)
rmm1 = X[:,3]
rmm2 = X[:,4]
phase = (X[:,5]).astype(int)
amplitude = X[:,6]

vdt=np.vectorize(dt.datetime)
rmmdate=vdt(year,month,day)
#rmmph=interpolate.interp1d(rmmdate,phase)(time)

def compsond(rmmph,x):

    nh=np.shape(x)[1]
    y = np.zeros((max(rmmph),nh))
    n = np.zeros((max(rmmph),nh))
    for it in range(np.shape(rmmph)[0]):
        ii=np.isfinite(x[it,:])
        if any(ii):
            y[rmmph[it]-1,ii] += x[it,ii]
            n[rmmph[it]-1,ii] += 1.0
    y = y/n
    return y


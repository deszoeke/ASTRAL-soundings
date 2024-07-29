import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#import string as str
import datetime as dt
from scipy import interpolate

#fileyears=range(2018,2018+1)
#filemonths=range(5,8+1)
fileyears=range(2015,2019) #excludes stop
filemonths=range(5,11)
nrec = 2 * 184 * len(fileyears) # number of sondes, assumes May-Oct inclusive each year, 2 sondes per day
npres=350
ipmax=0

# initialize variables (in global scope)
# housekeeping variables for each sonde
station_num_id=np.zeros(nrec,dtype=int)
station_wban_id=[] # empty list, append with station_wban_id.append("WBAN"). Lists are pythonic.
station_name=[]
time=[]
# sounding time,height variables
pres = np.nan+np.zeros((nrec,npres))
hght = np.nan+np.zeros((nrec,npres))
temp = np.nan+np.zeros((nrec,npres))
dwpt = np.nan+np.zeros((nrec,npres))
relh = np.nan+np.zeros((nrec,npres))
mixr = np.nan+np.zeros((nrec,npres))
drct = np.nan+np.zeros((nrec,npres))
sknt = np.nan+np.zeros((nrec,npres))
thta = np.nan+np.zeros((nrec,npres))
thte = np.nan+np.zeros((nrec,npres))
thtv = np.nan+np.zeros((nrec,npres))

def enfloat(strng,missingval=np.nan):
    """
    Converts a string to a float. Returns missingval if it can't.
    """
    try:
        return float(strng)
    except ValueError:
        return missingval

def get_line(file):
    line = ''
    while(str.strip(line) == ''):
        line = file.readline()
    return line

def read_sounding(line):
    """
    Reads one sounding from the U Wyoming text format and appends data to
    global variables.
    """
    global ipmax
    # the header line could be anything, assume we're at the first line of a record when called:
    #'91408 PTRO Koror, Palau Is Observations at 00Z 01 Aug 2017\n'
    word=str.strip(line).split()
    # find positions of spaces
    p = [ pos for pos, char in enumerate(line) if char==" "]
    # station identifiers
    station_num_id[irec]=int(word[0])
    station_wban_id.append(word[1])
    station_name.append(line[p[1]+1:p[4]])
    # timestamp
    #datestr=line[p[6]+1:p[7]-1]+line[p[7]::]
    datestr=word[-4][0:2]+' '+word[-3]+' '+word[-2]+' '+word[-1]
    #print(datestr)
    time.append(dt.datetime.strptime(datestr, '%H %d %b %Y'))
    
    # data header, always assumed in this form:
    #file.readline() #'\n' # get_line skips blank lines already
    get_line( file ) #'-----------------------------------------------------------------------------\n'
    get_line( file ) #'   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV\n'
    get_line( file ) #'    hPa     m      C      C      %    g/kg    deg   knot     K      K      K \n'
    get_line( file ) #'-----------------------------------------------------------------------------\n'
    # now 10-column positional data table starts
    line = get_line( file )
    ip=0 # loop though vertical coordinate
    while not(str.strip(line).startswith('Station information and sounding indices')):
        #print(line)
        pres[irec,ip] = enfloat(line[0:7])
        hght[irec,ip] = enfloat(line[9:14])
        temp[irec,ip] = enfloat(line[16:21])
        dwpt[irec,ip] = enfloat(line[23:28])
        relh[irec,ip] = enfloat(line[31:35])
        mixr[irec,ip] = enfloat(line[37:42])
        drct[irec,ip] = enfloat(line[46:49])
        sknt[irec,ip] = enfloat(line[52:56])
        thta[irec,ip] = enfloat(line[58:63])
        thte[irec,ip] = enfloat(line[65:70])
        thtv[irec,ip] = enfloat(line[72:77])
        ip += 1
        line = get_line( file )
    ipmax = max(ipmax,ip)
    # skip summary data to end of record
    while not(str.strip(line).startswith('Precipitable water [mm] for entire sounding')):
        line = get_line( file )


# start main loop

irec=0 # initialize global sounding record index, 0-based
# loop over all files
for fy in fileyears:
    for fm in filemonths:
        filename="{:0>4d}{:0>2d}.sed2txt".format(fy,fm)
        print(filename)
        with open(filename,'r') as file:
            # loop over each sounding
            line = get_line( file )
            while str.strip(line).startswith('University of Wyoming'):
                line = get_line( file )
            while not(str.strip(line).startswith(('Description of the data columns or sounding indices.','Description of the'))):
                read_sounding(line) # starts with the present line and reads more lines, appends data to sounding variables
                irec += 1
                line = get_line( file )

# truncate data
station_num_id=station_num_id[:irec]
pres=pres[:irec,:ipmax]
hght=hght[:irec,:ipmax]
temp=temp[:irec,:ipmax]
dwpt=dwpt[:irec,:ipmax]
relh=relh[:irec,:ipmax]
mixr=mixr[:irec,:ipmax]
drct=drct[:irec,:ipmax]
sknt=sknt[:irec,:ipmax]
thta=thta[:irec,:ipmax]
thte=thte[:irec,:ipmax]
thtv=thtv[:irec,:ipmax]

# compute u and v wind
kts2ms=0.51444444444 # converts kts to m/s 463/900?
pio180=np.pi/180
uwnd = -kts2ms*sknt*np.sin(drct*pio180)
vwnd = -kts2ms*sknt*np.cos(drct*pio180)

# time axis for plotting
datn=plt.matplotlib.dates.date2num(time)
adatn=np.tile(plt.matplotlib.dates.date2num(time),(np.shape(hght)[1],1)).transpose() # tile into a 2d grid

# interpolate to standard heights to facilitate plotting
ht=np.array(range(0,20001,50))

def int2ht(hght,xxxx,ht):
    y = np.asarray( [
            interpolate.interp1d( hght[i,:],xxxx[i,:],
                                  bounds_error=False,fill_value=np.nan )(ht)
                for i in range(0,hght.shape[0]) ] )
    return y

th=int2ht(hght,thta,ht)
rh=int2ht(hght,relh,ht)
mr=int2ht(hght,mixr,ht)/1e3 # -> kg/kg
qv=mr/(1+mr)
u=int2ht(hght,uwnd,ht)
v=int2ht(hght,vwnd,ht)

def anom(x,n=0):
    return x-np.nanmean(x,axis=n)               

"""
Read in RMM index
"""
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

def compsond(ph,x):
    nh=np.shape(x)[1]
    y  = np.zeros((max(ph),nh))
    n  = np.zeros((max(ph),nh))
    yy = np.nan+y
    for it in range(np.shape(ph)[0]):
        ii=np.isfinite(x[it,:])
        if any(ii):
            y[ph[it]-1,ii] += x[it,ii]
            n[ph[it]-1,ii] += 1.0
    yy[n>0] = y[n>0]/n[n>0]
    yy[n==0] = np.nan
    return yy

tdsond = ( (np.asarray(time)-time[0]) / dt.timedelta(days=1) ).astype(None) # -> decimal days from time[0]
tdrmm  = ( (         rmmdate-time[0]) / dt.timedelta(days=1) ).astype(None)
iit=np.asarray(time)<=dt.datetime(2018,8,10)
rmmph=( interpolate.interp1d(tdrmm,phase)(tdsond[iit]) ).astype(int)

rmmth  = compsond(rmmph,th)
# composite anomalies
rmmtha = compsond(rmmph,anom(th))
rmmrha = compsond(rmmph,anom(rh))
rmmqva = compsond(rmmph,anom(qv))
rmmua  = compsond(rmmph,anom(u))
rmmva  = compsond(rmmph,anom(v))

fig, ax = plt.subplots(2,2)
ax[0,0].plot(np.transpose(rmmtha),ht/1e3);     ax[0,0].set_title('theta (K)')
ax[0,1].plot(np.transpose(rmmqva)*1e3,ht/1e3); ax[0,1].set_title('qv (kg/kg)')
ax[1,0].plot(np.transpose(rmmua),ht/1e3);      ax[1,0].set_title('u (m/s)')
ax[1,1].plot(np.transpose(rmmva),ht/1e3);      ax[1,1].set_title('v (m/s)')
plt.legend(np.arange(1,9))
fig.set_figheight(6.5)
fig.savefig('Palau_RMMcmpsit.eps')
fig.savefig('Palau_RMMcmpsit.png')

# find gaps between years
start=np.vstack(([[0]],np.argwhere(np.diff(datn)>2)+1))
ends=start[[1,-1,0]]-1

# plot anomalies in 2015
fig, ax = plt.subplots(3,1)
for yi in range(0,3):
    pc=ax[yi].pcolormesh(datn,ht/1e3,anom(th).transpose(),cmap="RdYlBu",vmin=-6,vmax=6)
    ax[yi].set_xlim(datn[start[yi]],datn[ends[yi]])
    ax[yi].xaxis.set_major_locator(mdates.MonthLocator())  # every month
    ax[yi].xaxis.set_major_formatter(mdates.DateFormatter('%m'))
    ax[yi].text(datn[start[yi]]+10,10,'{}'.format(fileyears[yi]))
    plt.colorbar(pc,ax=ax[yi])
    
ax[2].set_xlabel('time (month)')
ax[1].set_ylabel('height (km)')
ax[0].set_title('potential temperature anomaly (K)')
#fig.autofmt_xdate()
plt.show()


np.argwhere(np.isnan(hght[:,0]))[0,0] # finds nan at lowest level
# 20 May 2015

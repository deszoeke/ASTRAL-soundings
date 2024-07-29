import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#import string as str
from datetime import datetime
from scipy import interpolate

fileyears=range(2015,2018) #excludes stop
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

def get_line(file)
    line = ''
    while(line == ''):
        line = str.strip( file.readline() )
    return line

def read_sounding(line):
    """
    Reads one sounding from the U Wyoming text format and appends data to
    global variables.
    """
    global ipmax
    # the header line could be anything, assume we're at the first line of a record when called:
    #'91408 PTRO Koror, Palau Is Observations at 00Z 01 Aug 2017\n'
    word=line.split()
    # find positions of spaces
    p = [ pos for pos, char in enumerate(line) if char==" "]
    # station identifiers
    station_num_id[irec]=int(word[0])
    station_wban_id.append(word[1])
    station_name.append(line[p[1]+1:p[4]])
    # timestamp
    datestr=line[p[6]+1:p[7]-1]+line[p[7]:-1]
    time.append(datetime.strptime(datestr, '%H %d %b %Y'))
    
    # data header, always assumed in this form:
    #file.readline() #'\n' # get_line skips blank lines already
    get_line( file ) #'-----------------------------------------------------------------------------\n'
    get_line( file ) #'   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV\n'
    get_line( file ) #'    hPa     m      C      C      %    g/kg    deg   knot     K      K      K \n'
    get_line( file ) #'-----------------------------------------------------------------------------\n'
    # now 10-column positional data table starts
    line = get_line( file )
    ip=0 # loop though vertical coordinate
    while not(line.startswith('Station information and sounding indices')):
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
    # read to end of record
    while not(line.startswith('Precipitable water [mm] for entire sounding')):
        line = get_line( file )


# start main loop

irec=0 # initialize global sounding record index, 0-based
# loop over all files
for fy in fileyears:
    for fm in filemonths:
        filename="{:0>4d}{:0>2d}.txt".format(fy,fm)
        with open(filename,'r') as file:
            # loop over each sounding
            line = get_line( file )
            while not(line.startswith('Description of the data columns or sounding indices.')):
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


# time axis for plotting
datn=plt.matplotlib.dates.date2num(time)
adatn=np.tile(plt.matplotlib.dates.date2num(time),(np.shape(hght)[1],1)).transpose() # tile into a 2d grid

# interpolate to standard heights to facilitate plotting
ht=np.array(range(0,14001,50))
th = np.asarray( [ interpolate.interp1d(hght[i,:],thta[i,:],bounds_error=False,fill_value=np.nan)(ht) for i in range(0,hght.shape[0]) ] )

def anom(x,n=0):
    return x-np.nanmean(x,axis=n)               

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

"""
Form curl requests to download monthly files from the U. Wyoming web site.
These requests must be run from the command line to get the data

Simon de Szoeke
"""

#import os

stnm = 91408 # station number for Palau
years=range(2016,2018+1)
months=range(5,10+1)
#months=[6,9]
#months=[5,7,8,10]

for year in years:
    for month in months:
        if (month==6 or month==9):
            nday=30
        else:
            nday=31

        cmdstr = ("""curl -s "http://weather.uwyo.edu/cgi-bin/sounding?region=pac&TYPE=TEXT%3ALIST&YEAR="""
                  +str(year)
                  +"&MONTH="
                  +"{:02d}".format(month)
                  +"&FROM=0100&TO="
                  +str(nday)
                  +"23&STNM="
                  +str(stnm)
                  +"""" | sed -e "s/<[a-zA-Z\/][^>]*>//g" > """
                  +str(year)
                  +"{:02d}".format(month)
                  +".txt")
        print(cmdstr)


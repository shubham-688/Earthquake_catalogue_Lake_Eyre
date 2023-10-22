import matplotlib.pyplot as plt
#plt.switch_backend("nbagg")
#plt.style.use('ggplot')
import obspy
from obspy import read, Stream
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import locations2degrees
import datetime
from obspy import UTCDateTime
from obspy.taup import TauPyModel
import sys
import os
from obspy.clients.fdsn import Client
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



##########################################################
client_5g=Client('http://auspass.edu.au:80',user='5g',password='grape71')
inventory = client_5g.get_stations(network="5G",station='*',level='response')
event_info=[]
UTC_times=[]
hour=[]
min=[]

#12 Event(s) in Catalog:
# event_time = "2018-11-27T12:54:14.20 | -29.283, +137.307 | 3.96 mb" #331
# event_time = "2018-12-02T22:29:50.83 | -29.470, +137.472 | 2.34 ml" #336

for line in open("mines_all_1.txt", "r"):
    # line='    1   2020 03 30 08:21:16.188    30076.188   1.4750     -27.5657     136.6053      10.0000     -inf      nan   15   13   28    84.46'
    event_info.append(line.split())
    a=line.split()
    UTC_times.append(UTCDateTime('{},{},{},{}'.format(a[1],a[2],a[3],a[4])))
    UTC_time=UTCDateTime('{},{},{},{}'.format(a[1],a[2],a[3],a[4]))+34200
    hour.append(UTC_time.hour)
    min.append(UTC_time.minute)

for line in open("mines_all_2.txt", "r"):
    # line='    1   2020 03 30 08:21:16.188    30076.188   1.4750     -27.5657     136.6053      10.0000     -inf      nan   15   13   28    84.46'
    event_info.append(line.split())
    a=line.split()
    UTC_times.append(UTCDateTime('{},{},{},{}'.format(a[1],a[2],a[3],a[4])))
    UTC_time=UTCDateTime('{},{},{},{}'.format(a[1],a[2],a[3],a[4]))+34200
    hour.append(UTC_time.hour)
    min.append(UTC_time.minute)

fig, ax=plt.subplots(figsize=(5,3))
ax.scatter(hour,min,c='darksalmon',s=50,alpha=.6,linewidths=.3,edgecolors='maroon')

ax.set_xbound(lower=0,upper=24)
ax.set_ybound(lower=0,upper=60)
ax.xaxis.set_major_locator(MultipleLocator(4))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(5))

ax.grid(alpha=0.4)
ax.set_xlabel('Hour (SA time - UTC+9:30 hours)',fontfamily='helvetica',fontsize=10)
ax.set_ylabel('Minute ',fontfamily='helvetica',fontsize=10)

#fig.show()
plt.savefig('hour_minute_all.png',dpi=400,bbox_inches='tight', pad_inches=0.2)

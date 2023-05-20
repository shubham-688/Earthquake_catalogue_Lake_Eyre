import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from obspy import UTCDateTime
###############

### the file has lat (new and old), long(new and old), depth (new and 10), and resid (new and old)
lat=[]
long=[]
Mlv=[]
error=[]
date=[]
depth=[]
# 8 columns
for line in open('na_mag_p_t.txt','r'):
    line=line.split()
    lat.append(float(line[0]))
    long.append(float(line[1]))
    Mlv.append(float(line[7]))
    date.append(UTCDateTime(line[4]).datetime)


fig, ax = plt.subplots(figsize=(7, 2))
ax.grid(alpha=0.7)
ax.scatter(date, Mlv,c='CadetBlue',alpha=.65)
ax.set_xlabel('Time',fontsize=14,labelpad=5)
ax.set_ylabel('Mlv',fontsize=14,labelpad=5)
plt.savefig('time_ml.png',bbox_inches='tight',dpi=500, pad_inches=0.2)
plt.show()

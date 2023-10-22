import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
# from obspy import UTCDateTime
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
    error.append(float(line[3]))
    depth.append(float(line[2]))

    # date.append(UTCDateTime(line[4]).datetime)

########################### MLv
fig, ax = plt.subplots(figsize=(4,4))
# plt.style.use('seaborn-whitegrid')
# plt.style.use('seaborn-paper')
plt.style.use('seaborn')
# plt.style.use('ggplot')
# plt.style.use('seaborn-paper')
# ax.set_facecolor('whitesmoke')
ax.grid(which='major', axis='both',color='black', linestyle='--',linewidth=.5,alpha=.25,zorder=1)
sns.histplot(data=np.array(Mlv),color='mediumseagreen',alpha=.65,zorder=5)
# ax.set_xlabel('Time',fontsize=14,labelpad=5)
ax.xaxis.set_major_locator(MultipleLocator(.5))
ax.set_xlabel('Mlv',fontsize=12,labelpad=5)
ax.set_ylabel('No. of earthquakes',fontsize=12,labelpad=5)
plt.savefig('ml_hist.png',bbox_inches='tight',dpi=300, pad_inches=0.2)
plt.show()

#################### Residual

fig, ax = plt.subplots(figsize=(4,4))
# plt.style.use('seaborn-whitegrid')
# plt.style.use('seaborn-paper')
plt.style.use('seaborn')
# plt.style.use('ggplot')
# plt.style.use('seaborn-paper')
# ax.set_facecolor('whitesmoke')
ax.grid(which='major', axis='both',color='black', linestyle='--',linewidth=.5,alpha=.25,zorder=1)
sns.histplot(data=np.array(error),color='darksalmon',alpha=.8,zorder=5)
# ax.set_xlabel('Time',fontsize=14,labelpad=5)
ax.xaxis.set_major_locator(MultipleLocator(.1))
ax.set_xlabel('Residual (s)',fontsize=12,labelpad=5)
ax.set_ylabel('No. of earthquakes',fontsize=12,labelpad=5)
plt.savefig('resid_hist.png',bbox_inches='tight',dpi=300, pad_inches=0.2)
plt.show()
#################### Depth Km

fig, ax = plt.subplots(figsize=(4,4))
# plt.style.use('seaborn-whitegrid')
# plt.style.use('seaborn-paper')
plt.style.use('seaborn')
# plt.style.use('ggplot')
# plt.style.use('seaborn-paper')
# ax.set_facecolor('whitesmoke')
ax.grid(which='major', axis='both',color='black', linestyle='--',linewidth=.5,alpha=.25,zorder=1)
sns.histplot(data=np.array(depth),color='slateblue',alpha=.55,zorder=5)
# ax.set_xlabel('Time',fontsize=14,labelpad=5)
ax.xaxis.set_major_locator(MultipleLocator(2.5))
ax.set_xlabel('Depth (km)',fontsize=12,labelpad=5)
ax.set_ylabel('No. of earthquakes',fontsize=12,labelpad=5)
plt.savefig('depth_hist.png',bbox_inches='tight',dpi=300, pad_inches=0.2)
plt.show()

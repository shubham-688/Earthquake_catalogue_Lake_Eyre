import numpy as np
import obspy
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.patches as patches
import sys
#########
def great_circ_calc_dist(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    km = 6371 * c
    return km
def dist_cal(lat1, lon1, lat2, lon2):
    """
    Determine station distance from cross-section pt
    """
    distance=111*np.sqrt((lat1-lat2)**2+(lon1-lon2)**2)
    # st.max_P=max_
    return distance-50
####################
lat=[]
long=[]
Mlv=[]
depth=[]

for line in open('proj.txt','r'):
    line=line.split()
    lat.append(float(line[0]))
    long.append(float(line[1]))
    Mlv.append(75*float(line[3]))
    depth.append(float(line[2]))
#####
cross_lat=-26.56
cross_long=134.5
##
# 135.688004 -27.794399 OOD 162.2
# 138.216202 -30.447201 LCRK 207.0
###########
# azimuth=140.83
##
dist=[]
for i,l in enumerate(lat):
    dist.append(dist_cal(lat[i],long[i],cross_lat,cross_long))

#####
fig, ax=plt.subplots(figsize=(12,5))
# plt.style.use('seaborn-whitegrid')
plt.style.use('seaborn-paper')
ax.set_facecolor('ghostwhite')
ax.grid(which='major', axis='both',color='black', linestyle='--',linewidth=.5,
alpha=.25,zorder=1)

# fig, ax = plt.subplots(figsize=(7, 2))
# ax.grid(alpha=0.7)
ax.scatter(dist, depth,s=Mlv,c='firebrick',alpha=.65)
ax.set_xlabel('Distance (km)',fontsize=14,labelpad=5)
ax.set_ylabel('Depth (km)',fontsize=14,labelpad=5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.scatter(dist_cal(-27.794399,135.688004,cross_lat,cross_long), 4.85,s=110,c='black',clip_on=False,marker='^',alpha=.85)
ax.scatter(dist_cal(-30.447201,138.216202,cross_lat,cross_long), 4.85,s=110,c='black',clip_on=False,marker='^',alpha=.85)
fig.text(.32, .94, 'OOD',fontsize=12, ha='center', va='center',color='black')
fig.text(.895, .94, 'LCRK',fontsize=12, ha='center', va='center',color='black')

plt.ylim(5, 11.75)
plt.xlim(0, 552.5)
plt.gca().invert_yaxis()
plt.savefig('dist_depth.png',bbox_inches='tight',dpi=500, pad_inches=0.2)
plt.show()

sys.exit()
#xy, width, height

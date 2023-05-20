import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys
# from obspy import UTCDateTime
###############
def circular_hist(ax, x, bins=16, density=True, offset=0, gaps=True,color='grey',al=.5):
    """
    Produce a circular histogram of angles on ax.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.PolarAxesSubplot
        axis instance created with subplot_kw=dict(projection='polar').

    x : array
        Angles to plot, expected in units of radians.

    bins : int, optional
        Defines the number of equal-width bins in the range. The default is 16.

    density : bool, optional
        If True plot frequency proportional to area. If False plot frequency
        proportional to radius. The default is True.

    offset : float, optional
        Sets the offset for the location of the 0 direction in units of
        radians. The default is 0.

    gaps : bool, optional
        Whether to allow gaps between bins. When gaps = False the bins are
        forced to partition the entire [-pi, pi] range. The default is True.

    Returns
    -------
    n : array or list of arrays
        The number of values in each bin.

    bins : array
        The edges of the bins.

    patches : `.BarContainer` or list of a single `.Polygon`
        Container of individual artists used to create the histogram
        or list of such containers if there are multiple input datasets.
    """
    # Wrap angles to [-pi, pi)
    x = (x+np.pi) % (2*np.pi) - np.pi

    # Force bins to partition entire circle
    if not gaps:
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    n, bins = np.histogram(x, bins=bins)

    # Compute width of each bin
    widths = np.diff(bins)

    # By default plot frequency proportional to area
    if density:
        # Area to assign each bin
        area = n / x.size
        # Calculate corresponding bin radius
        radius = (area/np.pi) ** .5
    # Otherwise plot frequency proportional to radius
    else:
        radius = n

    # Plot data on ax
    patches = ax.bar(bins[:-1], radius, zorder=1, align='edge', width=widths,
                     edgecolor='grey',color=color, fill=True, linewidth=1,alpha=al)

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels for area plots (they are mostly obstructive)
    if density:
        ax.set_yticks([])

    return n, bins, patches
###
### long lat depth strike dip rake Mlv t axis (azi and plunge) P axis (azi and plunge)
lat=[]
long=[]
depth=[]
strike=[]
dip=[]
rake=[]

# 8 columns
for line in open('fm_Strike_slip_p_t_clean.txt','r'):
    line=line.split()
    lat.append(float(line[1]))
    long.append(float(line[0]))
    depth.append(float(line[2]))
    strike.append(float(line[3]))
    dip.append(float(line[4]))
    rake.append(float(line[5]))
    # date.append(UTCDateTime(line[4]).datetime)

# for i,st in enumerate(strike):
#     if 270 > st > 180:
#         strike[i]=st-180
# ############

fig, ax = plt.subplots(figsize=(4,4),subplot_kw=dict(projection='polar'))
# Visualise by area of bins
circular_hist(ax, np.array(strike),offset=np.pi/2)
ax.set_theta_direction(-1)
ax.axvline(x=2.35619,ls='--',lw='1.75',color='black',alpha=.8)
ax.axvline(x=5.49779,ls='--',lw='1.75',color='black',alpha=.8)
ax.tick_params(which='major',labelsize=12)
plt.text(2.63,.38,'Strike',color='crimson',fontsize=22)
# plt.text(0,.35,'N',color='grey',weight='bold',ha='centre',bbox={'facecolor': 'yellow', 'alpha': 0.95, 'pad': 1},fontsize=13)
plt.text(0,.33,'N',color='grey',weight='bold',bbox={'facecolor': 'white', 'alpha': 0.95, 'pad': 2},ha='center',fontsize=17)
plt.savefig('strike_hist.png',bbox_inches='tight',dpi=300, pad_inches=0.2)

plt.show()

sys.exit()

fig, ax = plt.subplots(figsize=(4,4),subplot_kw=dict(projection='polar'))
# Visualise by area of bins
circular_hist(ax, np.array(dip),offset=np.pi/2,color='rebeccapurple',al=.65)
ax.set_theta_direction(-1)
ax.axvline(x=2.35619,ls='--',lw='1.75',color='black',alpha=.8)
ax.axvline(x=5.49779,ls='--',lw='1.75',color='black',alpha=.8)
ax.tick_params(which='major',labelsize=12)
ax.set_thetamin(0)
ax.set_thetamax(90)
plt.text(1.75,.18,'Dip',color='rebeccapurple',fontsize=22)
# plt.text(0,.35,'N',color='grey',weight='bold',ha='centre',bbox={'facecolor': 'yellow', 'alpha': 0.95, 'pad': 1},fontsize=13)
# plt.text(0,.33,'N',color='grey',weight='bold',bbox={'facecolor': 'white', 'alpha': 0.95, 'pad': 2},ha='center',fontsize=17)
plt.savefig('dip_hist.png',bbox_inches='tight',dpi=300, pad_inches=0.2)

sys.exit()

fig, ax = plt.subplots(figsize=(4,4),subplot_kw=dict(projection='polar'))
# Visualise by area of bins
circular_hist(ax, np.array(rake),offset=np.pi/2,color='darkorange',al=.65)
ax.set_theta_direction(-1)
ax.axvline(x=2.35619,ls='--',lw='1.75',color='black',alpha=.8)
ax.axvline(x=5.49779,ls='--',lw='1.75',color='black',alpha=.8)
ax.tick_params(which='major',labelsize=12)
# ax.set_thetamin(0)
# ax.set_thetamax(90)
plt.text(2.63,.38,'Rake',color='darkorange',fontsize=22)
# plt.text(0,.35,'N',color='grey',weight='bold',ha='centre',bbox={'facecolor': 'yellow', 'alpha': 0.95, 'pad': 1},fontsize=13)
# plt.text(0,.33,'N',color='grey',weight='bold',bbox={'facecolor': 'white', 'alpha': 0.95, 'pad': 2},ha='center',fontsize=17)
plt.savefig('rake_hist.png',bbox_inches='tight',dpi=300, pad_inches=0.2)

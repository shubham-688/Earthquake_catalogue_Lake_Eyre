import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream
import numpy as np
import scipy

import datetime
from collections import deque
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import locations2degrees

from obspy.signal.trigger import classic_sta_lta,z_detect,plot_trigger
from obspy.signal.trigger import recursive_sta_lta,delayed_sta_lta

from eqcorrscan.utils.despike import median_filter

def trigger_onset(charfct, thres1, thres2, max_len=9e99, max_len_delete=False):
    ind1 = np.where(charfct > thres1)[0]
    if len(ind1) == 0:
        return []
    ind2 = np.where(charfct > thres2)[0]
    #
    on = deque([ind1[0]])
    of = deque([-1])
    # determine the indices where charfct falls below off-threshold
    ind2_ = np.empty_like(ind2, dtype=bool)
    ind2_[:-1] = np.diff(ind2) > 1
    # last occurence is missed by the diff, add it manually
    ind2_[-1] = True
    of.extend(ind2[ind2_].tolist())
    on.extend(ind1[np.where(np.diff(ind1) > 1)[0] + 1].tolist())
    # include last pick if trigger is on or drop it
    if max_len_delete:
        # drop it
        of.extend([1e99])
        on.extend([on[-1]])
    else:
        # include it
        of.extend([ind2[-1]])
    #
    pick = []
    while on[-1] > of[0]:
        while on[0] <= of[0]:
            on.popleft()
        while of[0] < on[0]:
            of.popleft()
        if of[0] - on[0] > max_len:
            if max_len_delete:
                on.popleft()
                continue
            of.appendleft(on[0] + max_len)
        pick.append([on[0], of[0]])
    return np.array(pick, dtype=np.int64)

def plot_trigger11(trace, cft,cft1,cft2, thr_on, thr_off, show=False):
    import matplotlib.pyplot as plt
    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts, dtype=np.float32) / df
    fig = plt.figure(figsize=(7,9))
    ax1 = fig.add_subplot(411)
    ax1.plot(t, trace.data, 'k')
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = fig.add_subplot(412, sharex=ax1)
    ax2.plot(t, cft, 'k')
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax3 = fig.add_subplot(413, sharex=ax2)
    ax3.plot(t, cft1/100, 'k')
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax4 = fig.add_subplot(414, sharex=ax3)
    ax4.plot(t, cft2/100, 'k')
    on_off = np.array(trigger_onset(cft, thr_on, thr_off))
    i, j = ax1.get_ylim()
    try:
        ax1.vlines(on_off[:, 0] / df, i, j, color='tomato', lw=2,
                   label="Trigger On")
        ax1.vlines(on_off[:, 1] / df, i, j, color='skyblue', lw=2,
                   label="Trigger Off")
        ax1.legend()
    except IndexError:
        pass
    ax2.axhline(thr_on, color='tomato', lw=1, ls='--')
    ax2.axhline(thr_off, color='skyblue', lw=1, ls='--')
    ax4.set_xlabel("Time after %s [s]" % trace.stats.starttime.isoformat())
    ax2.set_ylabel("STA/LTA")
    ax3.set_ylabel("STA values")
    ax4.set_ylabel("LTA values")
    fig.suptitle(trace.id)
    fig.canvas.draw()
    if show:
        plt.show()

def classic_sta_lta_py(a, nsta, nlta):
    sta = np.cumsum(a ** 2)

    # Convert to float
    sta = np.require(sta, dtype=np.float)

    # Copy for LTA
    lta = sta.copy()

    # Compute the STA and the LTA
    sta[nsta:] = sta[nsta:] - sta[:-nsta]
    sta /= nsta
    lta[nlta:] = lta[nlta:] - lta[:-nlta]
    lta /= nlta

    # Pad zeros
    sta[:nlta - 1] = 0

    # Avoid division by zero by setting zero values to tiny float
    dtiny = np.finfo(0.0).tiny
    idx = lta < dtiny
    lta[idx] = dtiny

    return sta / lta,sta,lta
##############################################################################################################
# AES=['AES20','AES02','AES07', 'AES08', 'AES03', 'AES16','AES11','AES10','AES04','AES15', 'AES12','AES13']
# AEB=['AEB02','AEB20','AEB18', 'AEB11', 'AEB17','AEB16','AEB01','AEB12','AEB13', 'AEB15'] #,'AEB14']
# AES=['AES13','AES12','AES10','AES11']
# AEB=['AEB16','AEB01','AEB12','AEB13', 'AEB15']


#7 Event(s) in Catalog:
event_time = "2018-11-02T21:27:31.54 | -32.100, +138.767 | 2.25 ML" #306
#event_time = "2018-11-16T01:29:46.24 | -31.690, +138.568 | 3.74 mb" #320 #1700
#event_time = "2018-11-27T12:54:14.20 | -29.283, +137.307 | 3.96 mb" #331
#event_time = "2018-12-02T22:29:50.83 | -29.470, +137.472 | 2.34 ml" #336
#event_time = "2018-12-05T21:52:47.69 | -29.833, +141.526 | 3.49 mb"
#event_time = "2018-12-06T12:41:45.77 | -30.562, +138.201 | 3.45 mb" #340
#event_time = "2018-12-11T13:26:35.12 | -32.033, +139.203 | 3.25 mb" #345

#Julian_day
fmt = '%Y-%m-%dT%H:%M:%S'
#s = '2018-12-29T08:33:21'
dt = datetime.datetime.strptime(event_time[:19], fmt)
tt = dt.timetuple()
Jday=tt.tm_yday
year=tt.tm_year
month=tt.tm_mon
day=tt.tm_mday
minu=tt.tm_min
hr=tt.tm_hour
sec=tt.tm_sec
if hr<10:
    hr='%0*d' % (2, hr)
# time_sec=minu*60+sec
# if hr<10:
#     hr='%0*d' % (2, hr)
hr=int(hr)
print('Jday=',Jday, 'Hour=',hr)

###########
x='AES12' #AES12,13 AEB13,15 close...AES20,AEB20,02 far
# AES=['AES20','AES02','AES07', 'AES08', 'AES03', 'AES16','AES11','AES10','AES04','AES15', 'AES12','AES13']
# AEB=['AEB02','AEB20','AEB18', 'AEB11', 'AEB17','AEB16','AEB01','AEB12','AEB13', 'AEB15'] #,'AEB14']
# AES=['AES13','AES12','AES10','AES11']
# hr=12   # 12 for 340/331  ## 331 27 Nov # 340 6 Dec
trace = obspy.read('../../2018/{}/{}*{}0000.EHZ' .format(Jday,x,hr))[0]
df = trace.stats.sampling_rate
trace.filter('bandpass', freqmin=4, freqmax=19)
trace.taper(0.05,type='hann', max_length=4, side='left')
trace=median_filter(trace, multiplier=10, windowlength=4, interp_len=0.5)

###classic_sta_lta # First attempt 13/02: 3.5 - 25; 4 - 1.5 Trigger on/off
# second 4.5;70. 6;1.5

cft = classic_sta_lta_py(trace.data, int(3.5 * df), int(70 * df))
plot_trigger11(trace, cft[0],cft[1],cft[2], 8, 1.5,show=False)
plt.savefig('{}/dsp_Sta_{}.pdf'.format(Jday,x),bbox_inches='tight', pad_inches=0.2)

# z_detect
# cft1 = z_detect(trace.data, int(10 * df))
# plot_trigger(trace, cft1, 1.5, 0.6,show=False)
# plt.savefig('Z_detect.pdf')

## recursive sta/delta
# cft2 = recursive_sta_lta(trace.data, int(5 * df), int(10 * df))
# plot_trigger(trace, cft2, 1.7, 0.5)

client = Client("IRIS")
t_ood= UTCDateTime(2018, month,day , hr, 0,0)
OOD = client.get_waveforms("AU", "OOD", "*", "*Z", t_ood, t_ood + 60 * 60)
trace=OOD[0]
df = trace.stats.sampling_rate
trace.filter('bandpass', freqmin=4, freqmax=19)
trace.taper(0.05,type='hann', max_length=4, side='left')
trace=median_filter(trace, multiplier=10, windowlength=4, interp_len=0.5)

##
cft1 = classic_sta_lta_py(trace.data, int(3.5 * df), int(70 * df))
plot_trigger11(trace, cft1[0],cft1[1],cft1[2], 8, 1.5,show=False)
plt.savefig('{}/cl_Sta_OOD.pdf'.format(Jday),bbox_inches='tight', pad_inches=0.2)

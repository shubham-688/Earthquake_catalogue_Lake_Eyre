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

from obspy.signal.trigger import classic_sta_lta,z_detect,plot_trigger
from obspy.signal.trigger import recursive_sta_lta,delayed_sta_lta

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


x='AEB13' #AEB 13,15 close...AEB20 02 far
hr=12   # 12 for 340/331  ## 331 27 Nov # 340 6 Dec
trace = obspy.read('../../2018/340/{}*{}0000.EHZ' .format(x,hr))[0]
df = trace.stats.sampling_rate
trace.filter('bandpass', freqmin=4, freqmax=19)
trace.taper(0.05,type='hann', max_length=4, side='left')

###classic_sta_lta # First attempt 13/02: 3.5 - 25; 4 - 1.5 Trigger on/off
cft = classic_sta_lta_py(trace.data, int(4.5 * df), int(70 * df))
plot_trigger11(trace, cft[0],cft[1],cft[2], 6, 1.5,show=False)
plt.savefig('classic_SLTA_{}_340.pdf'.format(x),bbox_inches='tight', pad_inches=0.2)

# z_detect
# cft1 = z_detect(trace.data, int(10 * df))
# plot_trigger(trace, cft1, 1.5, 0.6,show=False)
# plt.savefig('Z_detect.pdf')

## recursive sta/delta
# cft2 = recursive_sta_lta(trace.data, int(5 * df), int(10 * df))
# plot_trigger(trace, cft2, 1.7, 0.5)

client = Client("IRIS")
t_ood= UTCDateTime(2018, 12,6 , hr, 0,0)
OOD = client.get_waveforms("AU", "OOD", "*", "*Z", t_ood, t_ood + 60 * 60)
trace=OOD[0]
df = trace.stats.sampling_rate
trace.filter('bandpass', freqmin=4, freqmax=19)
trace.taper(0.05,type='hann', max_length=1, side='left')

##
# cft1 = classic_sta_lta_py(trace.data, int(4.5 * df), int(70 * df))
# plot_trigger11(trace, cft1[0],cft1[1],cft1[2], 6, 1.5,show=False)
# plt.savefig('classic_SLTA_OOD_340.pdf')

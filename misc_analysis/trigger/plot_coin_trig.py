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
from obspy.signal.trigger import recursive_sta_lta,delayed_sta_lta,coincidence_trigger
from obspy.signal.trigger import ar_pick,pk_baer
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

from eqcorrscan.utils.despike import median_filter
import statistics as stat
from pprint import pprint
import time as tmm
import pickle
import os
import sys

def trigger_onset(charfct, thres1, thres2, max_len, max_len_delete=False):
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

# def plot_trigger11(trace, cft,cft1,cft2, thr_on, thr_off, show=False):
#     import matplotlib.pyplot as plt
#     df = trace.stats.sampling_rate
#     npts = trace.stats.npts
#     t = np.arange(npts, dtype=np.float32) / df
#     fig = plt.figure(figsize=(7,9))
#     ax1 = fig.add_subplot(411)
#     ax1.plot(t, trace.data, 'k',lw=.5)
#     ax1.xaxis.set_minor_locator(MultipleLocator(60))
#     plt.setp(ax1.get_xticklabels(), visible=False)
#     ax2 = fig.add_subplot(412, sharex=ax1)
#     ax2.plot(t, cft, 'k',lw=.5)
#     ax2.yaxis.set_minor_locator(MultipleLocator(.5))
#     ax2.xaxis.set_minor_locator(MultipleLocator(60))
#     plt.setp(ax2.get_xticklabels(), visible=False)
#     ax3 = fig.add_subplot(413, sharex=ax2)
#     ax3.plot(t, cft1/100, 'k',lw=.5)
#     plt.setp(ax3.get_xticklabels(), visible=False)
#     ax3.xaxis.set_minor_locator(MultipleLocator(60))
#     ax4 = fig.add_subplot(414, sharex=ax3)
#     ax4.plot(t, cft2/100, 'k',lw=.5)
#     ax4.xaxis.set_major_locator(MultipleLocator(600))
#     ax4.xaxis.set_minor_locator(MultipleLocator(60))
#     on_off = np.array(trigger_onset(cft, thr_on, thr_off))
#     i, j = ax1.get_ylim()
#     try:
#         ax1.vlines(on_off[:, 0] / df, i, j, color='tomato', lw=2,
#                    label="Trigger On")
#         ax1.vlines(on_off[:, 1] / df, i, j, color='skyblue', lw=2,
#                    label="Trigger Off")
#         ax1.legend()
#     except IndexError:
#         pass
#     ax2.axhline(thr_on, color='tomato', lw=1, ls='--')
#     ax2.axhline(thr_off, color='skyblue', lw=1, ls='--')
#     ax4.set_xlabel("Time after %s [s]" % trace.stats.starttime.isoformat())
#     ax2.set_ylabel("STA/LTA")
#     ax3.set_ylabel("STA values")
#     ax4.set_ylabel("LTA values")
#     fig.suptitle(trace.id)
#     fig.canvas.draw()
#     if show:
#         plt.show()

def recursive_sta_lta_py(a, nsta, nlta):
    try:
        a = a.tolist()
    except Exception:
        pass
    ndat = len(a)
    # compute the short time average (STA) and long time average (LTA)
    # given by Evans and Allen
    csta = 1. / nsta
    clta = 1. / nlta
    charfct = [0.0] * len(a)
    sta = [0.0] * len(a)
    lta = [0.0] * len(a)
    sta[0] = 0.
    lta[0] = 1e-99  # avoid zero division
    icsta = 1 - csta
    iclta = 1 - clta
    for i in range(1, ndat):
        sq = a[i] ** 2
        sta[i] = csta * sq + icsta * sta[i-1]
        lta[i] = clta * sq + iclta * lta[i-1]
        charfct[i] = sta[i] / lta[i]
        if i < nlta:
            charfct[i] = 0.
    return np.array(charfct),np.array(sta),np.array(lta)

def plot_trigger_1(trace, cft, thr_on, thr_off, show=False):
    import matplotlib.pyplot as plt
    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts, dtype=np.float32) / df
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, trace.data, 'k',lw=.5)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.plot(t, cft, 'k',lw=.5)
    ax2.yaxis.set_minor_locator(MultipleLocator(.5))
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax1.xaxis.set_major_locator(MultipleLocator(60))
    ax1.xaxis.set_minor_locator(MultipleLocator(30))
    ax2.xaxis.set_major_locator(MultipleLocator(60))
    ax2.xaxis.set_minor_locator(MultipleLocator(30))
    ########
    max_trigger_length=20 # max single trigger length in sec
    kwargs = {'max_len_delete': False}
    kwargs['max_len'] = int(max_trigger_length * df + 0.5)
    ######
    on_off = np.array(trigger_onset(cft, thr_on, thr_off, **kwargs))
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
    ax2.set_xlabel("Time after %s [s]" % trace.stats.starttime.isoformat())
    ax2.set_ylabel("STA/LTA")
    fig.suptitle(trace.id)
    fig.canvas.draw()
    if show:
        plt.show()

def bp_taper_dsp(trace): ## bandpass filter, taper the edges, and despike
    trace.filter('bandpass', freqmin=4, freqmax=19)
    trace.taper(0.05,type='hann', max_length=4, side='left')
    trace=median_filter(trace, multiplier=10, windowlength=4, interp_len=0.5)
    return trace

##############################################################################################################
# This script takes the coincidence_trigger files (in pickel format) detected by
# 'multi_rcr_detect.py' and for each day. Then for each item, it plots the
# seismograms with recursive On/off trigger for the detected stations. I also
# used a max_trigger_length for Recursive by modifying the plot_trigger_1.
# I used 20 sec time window as a Off.

######### taking one day at a time #302 Julian day is 29th Oct when all stations have recordings
###############
#x='AEB14' stopped working
# AES=['AES20','AES02','AES07', 'AES08', 'AES03', 'AES16','AES11','AES10','AES04','AES15', 'AES12','AES13']
ALL=['AES20','AES02','AES07', 'AES08', 'AES03', 'AES16','AES11','AES10','AES04', \
'AES15','AES12','AES13','AEB02','AEB20','AEB18', 'AEB11', 'AEB17','AEB16', \
'AEB01','AEB12','AEB13', 'AEB15','AES09','AEB19']
###########
start = tmm.process_time()
############ Setting params for STA/LTA
STA=3.5
LTA=60
On=5
Off=1
print('STA=',STA, 'LTA=',LTA,'On=',On,'Off=',Off)


######
for day in range(331,332):
    ### getting coincidence trigger from multi station trigger
    print('--> Doing day #',day)
    with open('../trigger/multi_trig/{}_triggers_pickle.txt'.format(day), 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        coinc_trigger = pickle.load(f)
    print('   # triggers =',len(coinc_trigger))
    folder_index=1
    for item in coinc_trigger:
        # item = coinc_trigger[0]
        count=item['coincidence_sum']
        time=item['time']
        Jday= item['time'].julday
        hr= item['time'].hour
        # print('Jday=',Jday)
        # pprint(item)
        if hr<10:
            hr='%0*d' % (2, hr)
        ######
        # print('Total stations in trigger =',count)
        item_index=1
        os.mkdir('{}/trig_{}'.format(day,folder_index))
        for station in item['stations']:
            if station=='OOD':
                client = Client("IRIS")
                OOD=client.get_waveforms("AU", "OOD", "*", "*Z", time-120, time + 300)
                trace=OOD[0]
                df = trace.stats.sampling_rate
                trace=bp_taper_dsp(trace)
                cft1 = recursive_sta_lta_py(trace.data, int(STA * df), int(LTA * df))
                plot_trigger_1(trace, cft1[0], On, Off,show=False)
                plt.savefig('{}/trig_{}/item_{}.pdf'.format(day,folder_index,item_index), \
                bbox_inches='tight', pad_inches=0.2)
            elif station=='MAL49':
                cl = Client('http://auspass.edu.au:80',user='3g', password='orange44')
                MAL=cl.get_waveforms("3G", "MAL49", "*", "*Z", time-120, time + 300)
                trace=MAL[0]
                df = trace.stats.sampling_rate
                trace=bp_taper_dsp(trace)
                cft1 = recursive_sta_lta_py(trace.data, int(STA * df), int(LTA * df))
                plot_trigger_1(trace, cft1[0], On, Off,show=False)
                plt.savefig('{}/trig_{}/item_{}.pdf'.format(day,folder_index,item_index), \
                bbox_inches='tight', pad_inches=0.2)
            else:
                stream= obspy.read('../../2018/{}/{}*0000.EHZ'.format(Jday,station))
                stream=stream.slice(time-120, time+300)
                trace=stream.merge()[0]
                df = trace.stats.sampling_rate
                trace=bp_taper_dsp(trace)
                cft1 = recursive_sta_lta_py(trace.data, int(STA * df), int(LTA * df))
                plot_trigger_1(trace, cft1[0], On, Off,show=False)
                plt.savefig('{}/trig_{}/item_{}.pdf'.format(day,folder_index,item_index), \
                bbox_inches='tight', pad_inches=0.2)
            item_index+=1
            plt.close('all')
        folder_index+=1
        plt.close('all')

### printing run time
runtime = tmm.process_time() - start
print (" Running Time in seconds %s " %runtime)
plt.close('all')


##### The following loop finds the number of times
####   the coincidence_sum values occurs

# trig_count={}
# for x in range(5,26):
#     trig_count["count{}".format(x)]=0
#     var=0
#     for item in coinc_trigger:
#         count=item['coincidence_sum']
#         if count==x:
#             var +=1
#     trig_count["count{}".format(x)]=var
#
# pprint(trig_count)

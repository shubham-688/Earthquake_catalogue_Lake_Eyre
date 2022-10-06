import matplotlib.pyplot as plt
import obspy
import os
from obspy import read,UTCDateTime,Stream
from obspy.signal.trigger import recursive_sta_lta, trigger_onset, classic_sta_lta

import datetime
import numpy as np

from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import locations2degrees
from pprint import pprint
import time as tmm
import sys
#################################################################################

def bp_taper_dsp(trace): ## bandpass filter, taper the edges, and despike
    trace.detrend('linear')
    trace.detrend('demean')
    trace.remove_response(pre_filt=pre_filt,water_level=60)
    trace.filter('bandpass', freqmin=4, freqmax=19)
    trace.taper(0.05,type='hann', max_length=4, side='both')
 #   trace=median_filter(trace, multiplier=10, windowlength=4, interp_len=0.5)

    return trace
##

def recSTALTAPy_h(a, b, nsta, nlta):
    """
    Recursive STA/LTA written in Python.

    .. note::

        There exists a faster version of this trigger wrapped in C
        called :func:`~obspy.signal.trigger.recSTALTA` in this module!

    :type a: NumPy ndarray
    :param a: Seismic Trace
    :type nsta: Int
    :param nsta: Length of short time average window in samples
    :type nlta: Int
    :param nlta: Length of long time average window in samples
    :rtype: NumPy ndarray
    :return: Characteristic function of recursive STA/LTA

    .. seealso:: [Withers1998]_ (p. 98) and [Trnkoczy2012]_
    """
    try:
        a = a.tolist()
    except:
        pass

    try:
        b = b.tolist()
    except:
        pass
    ndat = len(a)
    # compute the short time average (STA) and long time average (LTA)
    csta = 1. / nsta
    clta = 1. / nlta
    sta = 0.
    lta = 1e-99  # avoid zero devision
    charfct = [0.0] * len(a)
    icsta = 1 - csta
    iclta = 1 - clta
    for i in range(1, ndat):
        sq = a[i] ** 2 + b[i] ** 2
        sta = csta * sq + icsta * sta
        lta = clta * sq + iclta * lta
        charfct[i] = sta / lta
        if i < nlta:
            charfct[i] = 0.
    return np.array(charfct)
########################################
############ Setting params for STA/LTA
STA=.2
LTA=15
On=4
Off=2
print('STA=',STA, 'LTA=',LTA,'On=',On,'Off=',Off)
start = tmm.process_time()
pre_filt = [0.005, 0.01, 25, 30] # bandpass filter

##########
client_5g=Client('http://auspass.edu.au:80',user='5g',password='antipasti')
inventory = client_5g.get_stations(network="5G",station='*',level='response')
year=2021
# for Jday in days:
for Jday in range(1,32): # month of Jan
    print('##############################')
    Jday='%0*d' % (3, Jday)
    timeS=UTCDateTime("{},{},0,0,0".format(year,Jday))
    mon=timeS.month
    day=timeS.day
    # if day<10:
    day='%0*d' % (2, day)
    mon='%0*d' % (2, mon)
    date = str(year)+str(mon)+str(day)+'/'
    if not os.path.exists(date):
        os.makedirs(date)
    print('doing day:', timeS, 'Jday=',Jday)

    with open("stations_19.txt", "r") as f:
        for station in f:
            stlo, stla, net, sta, chan, elev = station.split()
            chanz = chan[:2]+"Z"
            chann = chan[:2]+"N"
            chane = chan[:2]+"E"

            try:
                if net=='5G':
                    ste=client_5g.get_waveforms("5G", '{}'.format(sta), "*", "*E", timeS, timeS+24*3600 -1,attach_response=True)
                    stn=client_5g.get_waveforms("5G", '{}'.format(sta), "*", "*N", timeS, timeS+24*3600 -1,attach_response=True)

                    # st = obspy.read('../../../../2018/{}/SAC/{}*.{}.SAC'.format(Jday,sta,chanz))
                    ste.merge(method=1,fill_value=0)
                    stn.merge(method=1,fill_value=0)
                else:
                    client = Client("IRIS")
                    ste=client.get_waveforms("AU", '{}'.format(sta), "*", "*E", timeS, timeS+24*3600 -1,attach_response=True)
                    stn=client.get_waveforms("AU", '{}'.format(sta), "*", "*N", timeS, timeS+24*3600 -1,attach_response=True)

                    ste.merge(method=1,fill_value=0)
                    stn.merge(method=1,fill_value=0)

                    # st.merge()
                tre = ste[0]
                trn = stn[0]


                df = tre.stats.sampling_rate
                tre=bp_taper_dsp(tre)
                trn=bp_taper_dsp(trn)

                tstart = tre.stats.starttime - UTCDateTime(year, int(mon), int(day), 0, 0, 0)
                # tstart = 0

                output = './'+date+net+'.'+sta+'.'+'S.txt'

                # Characteristic function and trigger onsets, see ObsPy website
                cft = recSTALTAPy_h(tre.data,trn.data, int(STA * df), int(LTA * df))
                on_of = trigger_onset(cft, On, Off)


                # Output the triggered
                f = open(output,'w')
                i = 0
                while(i<len(on_of)):
                    trig_on = on_of[i,0]
                    trig_of = on_of[i,1]

                    if max(cft[trig_on:trig_of]) > 8:
                        f.write('{} {} {}\n'.format((tstart+trig_on/df),max(cft[trig_on:trig_of]),0.0))
                        # f.write('{} {} {}\n'.format((tstart+trig_on/df),max(cft[trig_on:trig_of]),amp))
                    i=i+1
                f.close()
                print('done station:',sta,'for day',Jday)
            except:
                # print('no data in/some shit went down for', sta)
                print('-')
    print('########################################\n')


### printing run time
runtime = tmm.process_time() - start
print (" Running Time in seconds %s " %runtime)

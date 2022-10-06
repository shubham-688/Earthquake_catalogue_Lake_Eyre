import matplotlib.pyplot as plt
import obspy
import os
from obspy import read,UTCDateTime,Stream
from obspy.signal.trigger import recursive_sta_lta, trigger_onset, classic_sta_lta
import datetime
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
############ Setting params for STA/LTA
STA=.3
LTA=15
On=5
Off=2
print('STA=',STA, 'LTA=',LTA,'On=',On,'Off=',Off)
start = tmm.process_time()
pre_filt = [0.005, 0.01, 25, 30] # bandpass filter
##########
client_5g=Client('http://auspass.edu.au:80',user='5g',password='salad')
inventory = client_5g.get_stations(network="5G",station='*',level='response')

year=2021 # year 
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
                    st=client_5g.get_waveforms("5G", '{}'.format(sta), "*", "*Z", timeS, timeS+24*3600 -1,attach_response=True)
                    # st = obspy.read('../../../../2018/{}/SAC/{}*.{}.SAC'.format(Jday,sta,chanz))
                    st.merge(method=1,fill_value=0)
                else:
                    client = Client("IRIS")
                    st=client.get_waveforms("AU", '{}'.format(sta), "*", "*Z", timeS, timeS+24*3600 -1,attach_response=True)
                    st.merge(method=1,fill_value=0)
                tr = st[0]

                df = tr.stats.sampling_rate
                tr=bp_taper_dsp(tr)
                tstart = tr.stats.starttime - UTCDateTime(year, int(mon), int(day), 0, 0, 0)
                # tstart = 0

                output = './'+date+net+'.'+sta+'.'+'P.txt'

                # Characteristic function and trigger onsets, see ObsPy website
                cft = recursive_sta_lta(tr.data, int(STA * df), int(LTA * df))
                on_of = trigger_onset(cft, On, Off)

                f = open(output,'w')
                i = 0
                while(i<len(on_of)):
                    trig_on = on_of[i,0]
                    trig_of = on_of[i,1]
                    # trig_off = int(trig_of + (trig_of - trig_on)*4.0)
                    #1000 is from meter to millimeter (mm) see Hutton and Boore (1987)
                    # amp = max(max(abs(datatre[trig_on:trig_off])),max(abs(datatrn[trig_on:trig_off])))*10000
                    if max(cft[trig_on:trig_of]) > 10: #only keeps STA/LTA > 10
                        f.write('{} {} {}\n'.format((tstart+trig_on/df),max(cft[trig_on:trig_of]),0.0))
                        # f.write('{} {} {}\n'.format((tstart+trig_on/df),max(cft[trig_on:trig_of]),amp))
                    i=i+1
                f.close()
                print('done station:',sta,'for day',Jday)
            except:
                # print('no data in/some shit went down for', sta)
                print('-')

### printing run time
runtime = tmm.process_time() - start
print (" Running Time in seconds %s " %runtime)

from math import log10
import numpy as np
import obspy
from obspy.geodetics import gps2dist_azimuth
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import math
import sys
import logging
import copy
import os
from obspy import Trace
from obspy.signal.invsim import simulate_seismometer as seis_sim
from obspy.signal.invsim import estimate_wood_anderson_amplitude_using_response as wa_ampl
from obspy.core.event import (
    Amplitude, Pick, WaveformStreamID, Origin, ResourceIdentifier)
Logger = logging.getLogger(__name__)

###################################
#################
PAZ_WA = {'sensitivity': 2080, 'zeros': [0+0j], 'gain': 1,
          'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j]}
##
PAZ_WA = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
          'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
#paz_WA= {'sensitivity': 2080, 'zeros': [0j], 'gain': 1,
#         'poles':[(-5.49779 - 5.60886j), (-5.49779 + 5.60886j)]}

def _max_p2t(data, delta, return_peak_trough=False):
    """
    Finds the maximum peak-to-trough amplitude and period.
    Originally designed to be used to calculate magnitudes (by
    taking half of the peak-to-trough amplitude as the peak amplitude).
    :type data: numpy.ndarray
    :param data: waveform trace to find the peak-to-trough in.
    :type delta: float
    :param delta: Sampling interval in seconds
    :type return_peak_trough: bool
    :param return_peak_trough:
        Optionally return the peak and trough
    :returns:
        tuple of (amplitude, period, time) with amplitude in the same
        scale as given in the input data, and period in seconds, and time in
        seconds from the start of the data window.
    :rtype: tuple
    """
    turning_points = []  # A list of tuples of (amplitude, sample)
    for i in range(1, len(data) - 1):
        if (data[i] < data[i - 1] and data[i] < data[i + 1]) or\
           (data[i] > data[i - 1] and data[i] > data[i + 1]):
            turning_points.append((data[i], i))
    if len(turning_points) >= 1:
        amplitudes = np.empty([len(turning_points) - 1],)
        half_periods = np.empty([len(turning_points) - 1],)
    else:
        Logger.warning(
            'Turning points has length: ' + str(len(turning_points)) +
            ' data have length: ' + str(len(data)))
        return 0.0, 0.0, 0.0
    for i in range(1, len(turning_points)):
        half_periods[i - 1] = (delta * (turning_points[i][1] -
                                        turning_points[i - 1][1]))
        amplitudes[i - 1] = np.abs(turning_points[i][0] -
                                   turning_points[i - 1][0])
    amplitude = np.max(amplitudes)
    period = 2 * half_periods[np.argmax(amplitudes)]
    delay = delta * turning_points[np.argmax(amplitudes)][1]
    if not return_peak_trough:
        return amplitude, period, delay
    max_position = np.argmax(amplitudes)
    peak = max(
        t[0] for t in turning_points[max_position: max_position + 2])
    trough = min(
        t[0] for t in turning_points[max_position: max_position + 2])
    return amplitude, period, delay, peak, trough

#################
####### this gets mean of a list after removing values outside 2std
def get_stdrem_mean(mags):
    # mags=[i[0] for i in mags_st]
    if len(mags)==0:
        return 'nan'
    elif len(mags)==1:
        mean=np.mean(np.array(mags))
        return round(float(mean),3)
    else:
        mags_ar=np.array(mags)
        mean = np.mean(mags_ar, axis=0)
        sd = np.std(mags_ar, axis=0)

        mags_stdrem = [x for x in mags if (x > mean - 1.5 * sd)]
        mags_stdrem = [x for x in mags_stdrem if (x < mean + 1.5 * sd)]
        mm=np.mean(np.array(mags_stdrem))
        return round(float(mm),3)
###
def unique(list1):
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))
    return unique_list
############
#station='AEB18'
#station='AES12'
model = TauPyModel(model="iasp91")
#
client = Client("IRIS")
# client_aus=Client('http://auspass.edu.au:80')
# inven_S1=client_aus.get_stations(network='S1',station='AUROX',channel='H*',level='response')
client_5g=Client('http://auspass.edu.au:80',user='5g',password='grape71')
inventory = client_5g.get_stations(network="5G",station='*',level='response')
inventory_au = client.get_stations(network="AU",station='MULG,OOD,LCRK',level='response')
##
#paz_WA= {'sensitivity': 2080, 'zeros': [0j], 'gain': 1,'poles':[(-5.49779 - 5.60886j), (-5.49779 + 5.60886j)]}
################
mag_all=[]
Eq_mag_all=[]
# station_name=['AES12','AEB18','AES05','AES04','AES03','AES02','AEB19',\
# 'AES07',' AES08','AEB02','AES09','AES19','AES17','AEB10','AEB08',\
# 'AES11','AES10','AEB15', 'AEB17','AES15','AES14','OOD','LCRK','MULG']

#station_name=['AES03','AEB18','AES09','AEB15']#,'AES13'] # eq 12
#station_name=['AES04','AES10','AES11','AES15'] #eq 4
#event = np.genfromtxt('eq40.txt')
##
eq_num=1
while eq_num < 28: # change here for mine/eq
    print('Doing Earthquake', eq_num)
    eq_file = os.path.join('mine_mag_info/'.format(eq_num), '')#change here for mine/eq
    st_file = os.path.join('mine_mag_info/'.format(eq_num), '')#change here for mine/eq


    station_name=[]
    for stt in open(st_file+'station_'+str(eq_num)+'.txt','r'): # station deets
        station_name.append(stt.split()[0])

    station_name=unique(station_name)
    print('Number of stations',len(station_name))
    ##
    times=[]
    event_lat=[]
    event_long=[]
    event_depth=[]
    for line in open(eq_file+'eq'+str(eq_num)+'.txt','r'): # eq details
        line=line.split()
        line_eq=line.copy()
        event_lat.append(line[0])
        event_long.append(line[1])
        times.append(line[4])
        event_depth.append(line[2])
    #################
    mags_st=[]
    mags_au=[]
    mags_5g_sp=[]
    mags_5g_bb=[]
    for station_5G in station_name:
        # print('Doing', station_5G)
        t=obspy.UTCDateTime(times[0])
    #for station_5G in ['AES12']:
        if station_5G in ['MULG','OOD','LCRK']:

            st_lat=inventory_au.get_coordinates('AU.{}..BHZ'.format(station_5G))['latitude']
            st_long=inventory_au.get_coordinates('AU.{}..BHZ'.format(station_5G))['longitude']
            try:
                response = inventory.get_response('AU.{}..BHZ'.format(station_au),t)
            except:
                pass
            try:
                response = inventory.get_response('AU.{}.00.HHZ'.format(station_au),t)
            except:
                pass
            ### gets waveform around P and S picks for each station
            epi_dist, az, baz = gps2dist_azimuth(float(event_lat[0]),float(event_long[0]), st_lat, st_long)
            epi_dist = epi_dist / 1000
            arrival_st=model.get_travel_times(source_depth_in_km=float(event_depth[0]),\
            distance_in_degree=epi_dist/111.19,phase_list=["P","p"])
            t_arr=arrival_st[0].time
            tsp = 0.73*(epi_dist/6.0) + 8
            st = client.get_waveforms('AU',station_5G,'*','BH*',starttime=t+t_arr-5,\
            endtime=t+t_arr+tsp,attach_response=True)
            ####
            # inventory_1 = client.get_stations(network="AU",station=station_5G,level='response')

            #####
            st1=st.copy()
            st1.detrend("demean")
            st1.detrend("linear")
            st1.filter('bandpass', freqmin=2, freqmax=10)
            st1.taper(0.05,type='hann', max_length=2, side='left')
            tr_z= st1.select(component="Z")[0]
            if max(abs(tr_z.data)) == 0:
                continue
            (amplitude, period, delay)=_max_p2t(tr_z.data,tr_z.stats.delta)
            ampl_wa=wa_ampl(response,amplitude,period)

            #
            # if ampl_wa == 0:
            #     continue

            r=math.sqrt(float(event_depth[0])**2+epi_dist**2)
            mag_v=log10(ampl_wa) + (1.1 * log10(r)) +(0.0013 * r) + 0.7 # SA
            mags_st.append([round(mag_v,2),station_5G,round(r,2)])
            mags_au.append(round(mag_v,2))

        else:
            try:
                st_lat=inventory.get_coordinates('5G.{}..HHZ'.format(station_5G))['latitude']
                st_long=inventory.get_coordinates('5G.{}..HHZ'.format(station_5G))['longitude']
                response = inventory.get_response('5G.{}..HHZ'.format(station_5G),t)

            except:
                pass
            try:
                st_lat=inventory.get_coordinates('5G.{}..EHZ'.format(station_5G))['latitude']
                st_long=inventory.get_coordinates('5G.{}..EHZ'.format(station_5G))['longitude']
                response = inventory.get_response('5G.{}..EHZ'.format(station_5G),t)
            except:
                pass

            epi_dist, az, baz = gps2dist_azimuth(float(event_lat[0]),float(event_long[0]), st_lat, st_long)
            epi_dist = epi_dist / 1000
            arrival_st=model.get_travel_times(source_depth_in_km=float(event_depth[0]),\
            distance_in_degree=epi_dist/111.19,phase_list=["P","p"])
            t_arr=arrival_st[0].time
            tsp = 0.73*(epi_dist/6.0) + 8

            st = client_5g.get_waveforms('5G',station_5G,'','*',starttime=t+t_arr-5,\
            endtime=t+t_arr+tsp,attach_response=True)
            #
            #####
            # inventory_1 = client_5g.get_stations(network="5G",station=station_5G,level='response')

            st1=st.copy()
            st1.detrend("demean")
            st1.detrend("linear")
            st1.filter('bandpass', freqmin=2, freqmax=10)
            st1.taper(0.05,type='hann', max_length=2, side='left')
            tr_z= st1.select(component="Z")[0]
            if max(abs(tr_z.data)) == 0:
                continue
            (amplitude, period, delay)=_max_p2t(tr_z.data,tr_z.stats.delta)
            ampl_wa=wa_ampl(response,amplitude,period)

            # if ampl_wa == 0:
            #     continue
            # print(ampl_z, station_5G)

            r=math.sqrt(float(event_depth[0])**2+epi_dist**2)
            mag_v=log10(ampl_wa) + (1.1 * log10(r)) +(0.0013 * r) + 0.7 # SA
            mags_st.append([round(mag_v,2),station_5G,round(r,2)])
            if st[0].stats.channel[0]=='H':
                mags_5g_bb.append(round(mag_v,2))
            elif st[0].stats.channel[0]=='E':
                mags_5g_sp.append(round(mag_v,2))

    mean_sp=get_stdrem_mean(mags_5g_sp)
    mean_bb=get_stdrem_mean(mags_5g_bb)
    mean_au=get_stdrem_mean(mags_au)


    mags=[i[0] for i in mags_st]
    mean_all=get_stdrem_mean(mags)
    # print('Len of mags_st',len(mags_st))
    print('mean sp, bb, au, all',mean_sp,mean_bb,mean_au,mean_all,'for eq',eq_num)
    eq_num+=1

    eq_mags=[str(mean_sp),str(mean_bb),str(mean_au),str(mean_all)]
    eq_deets_mag=line_eq+eq_mags
    Eq_mag_all.append(eq_deets_mag)
    mag_all.append(eq_mags)
    print('################################')
 # mags_all.append(np.mean(np.array(mags_stdrem)))
with open('mine_na_mag_p_t.txt','w') as fi:
    for event in Eq_mag_all:
        fi.write(event[0]+' '+event[1]+' '+event[2]+' '+event[3]+' '+event[4]+' '+event[5]+' '+event[6]+' '+event[8]+'\n')

with open('mine_just_mags_p_t.txt','w') as fi:
    for event in mag_all:
        fi.write(event[0]+' '+event[1]+' '+event[3]+'\n')
#print(mags_st)

#sys.exit()

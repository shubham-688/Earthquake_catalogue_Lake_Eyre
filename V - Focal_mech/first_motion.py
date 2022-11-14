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
################################################################
model = TauPyModel(model="iasp91")
#
def pre_process(trace, inventory, water_level):
    trace.detrend('simple')
    # Remove response to Velocity
    try:
        trace.remove_response(
            inventory=inventory, output="VEL", water_level=water_level)
    except Exception:
        Logger.error("No response for {trace.id} at {trace.stats.starttime}")
        return None

    # trace
    return trace
##
def get_t_arr(event_lat,event_long,event_depth,st_lat, st_long):
    epi_dist, az, baz = gps2dist_azimuth(float(event_lat),float(event_long), st_lat, st_long)
    epi_dist = epi_dist / 1000
    arrival_st=model.get_travel_times(source_depth_in_km=float(event_depth),\
    distance_in_degree=epi_dist/111.19,phase_list=["P","p"])
    t_arr=arrival_st[0].time
    return t_arr

##
def plot_wave(trace,P_arr,show=True):
    trace.trim(P_arr-10,P_arr+10)
    df = trace.stats.sampling_rate
    npts = trace.stats.npts
    t = np.arange(npts, dtype=np.float32) / df
    p_arr=(P_arr-tr_z.stats.starttime)
    trace_1=trace.copy()
    trace_3=trace.copy()
    trace_6=trace.copy()
    trace_1.filter('bandpass', freqmin=1, freqmax=10)
    trace_1.taper(.5,max_length=.5, side='left')
    trace_3.filter('bandpass', freqmin=3, freqmax=10)
    trace_3.taper(.5,max_length=.5, side='left')
    trace_6.filter('bandpass', freqmin=6, freqmax=10)
    trace_6.taper(.5,max_length=.5, side='left')
    fig = plt.figure(figsize=(15,8))
    ax1 = fig.add_subplot(311)
    ax1.plot(t, trace_1.data, color='navy',lw=1)
    i, j = ax1.get_ylim()
    ax1.vlines(p_arr, i/3, j/3, color='tomato', lw=2,
                   label="P-arrival")
    ax1.axhline(0, color='gray', lw=1, ls='--')
    ax1.set_xticks(np.arange(21))
    ax2 = fig.add_subplot(312)
    ax2.plot(t, trace_3.data, color='navy',lw=1)
    i, j = ax2.get_ylim()
    ax2.vlines(p_arr, i/3, j/3, color='tomato', lw=2,
                   label="P-arrival")
    ax2.axhline(0, color='gray', lw=1, ls='--')
    ax2.set_xticks(np.arange(21))
    ax3 = fig.add_subplot(313)
    ax3.plot(t, trace_6.data, color='navy',lw=1)
    i, j = ax3.get_ylim()
    ax3.vlines(p_arr, i/3, j/3, color='tomato', lw=2,
                   label="P-arrival")
    ax3.axhline(0, color='gray', lw=1, ls='--')
    ax3.set_xticks(np.arange(21))
    ax3.set_xlabel("Time after %s [s]" % trace.stats.starttime.isoformat())
    fig.suptitle(trace.id)
    fig.canvas.draw()
    plt.tight_layout()
    plt.show()
    if show:
        plt.show()
##
###
def unique(list1):
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))
    return unique_list
############
client = Client("IRIS")
client_aus=Client('http://auspass.edu.au:80')
client_5g=Client('http://auspass.edu.au:80',user='5g',password='')
##
inven_S1=client_aus.get_stations(network='S1',station='AUROX',channel='H*',level='response')
inventory_5G = client_5g.get_stations(network="5G",station='*',level='response')
inventory_au = client.get_stations(network="AU",station='MULG,OOD,LCRK',level='response')
########
# starttime=inven_S1.select(station='AUROX')[0][0].start_date
# endtime=inven_S1.select(station='AUROX')[0][0].end_date
##########
# to find fm, do eq's one by one.

eq_num=1 ##

# %reset -f
##################
print('Doing Earthquake', eq_num)
eq_file = os.path.join('../mag_esti/mag_info/', '')#change here for mine/eq
st_file = os.path.join('../mag_esti/mag_info/', '')#change here for mine/eq

station_name_all=[] # all 40 LAke Eyre stations
for i in range(len(inventory_5G.get_contents()['stations'])):
    station_name_all.append(inventory_5G.get_contents()['stations'][i].split()[0].split('.')[1])

station_name=[]
for stt in open(st_file+'station_'+str(eq_num)+'.txt','r'): # station deets
    station_name.append(stt.split()[0])

station_name=unique(station_name)
print('Number of stations',len(station_name))
##
##
for line in open(eq_file+'eq'+str(eq_num)+'.txt','r'): # eq details
    line=line.split()
    event_lat=line[0]
    event_long=line[1]
    times=line[4]
    event_depth=line[2]
    t=obspy.UTCDateTime(times)

print('Station Lat Long:',event_lat,event_long,'\n')
print('Time of event:',t,'\n')

########################

# station_name=[ 'AES07','AES17', 'AES12','AES08','AEB15']
#station_name=[ 'AES10','AES11','AES12', 'AES15','AES08','AES04','AES12','AEB08']
# station_name=[ 'AEB13','AES08','AEB15','AEB16', 'AES04','AEB04','AES12','AES15','AES14']
#station_name=[ 'AEB17', 'AES14','AES02','AEB20','AEB07']

# station_name=[ 'AES13']#,'OOD','LCRK','MULG']

st_all=Stream()
fm_info=[]
for station_5G in station_name:
# for station_5G in ['AEB17']:
    try:
        try:
            st_lat=inventory_5G.get_coordinates('5G.{}..HHZ'.format(station_5G))['latitude']
            st_long=inventory_5G.get_coordinates('5G.{}..HHZ'.format(station_5G))['longitude']
            st_elevation=inventory_5G.get_coordinates('5G.{}..HHZ'.format(station_5G))['elevation']
        except:
            pass
        try:
            st_lat=inventory_5G.get_coordinates('5G.{}..EHZ'.format(station_5G))['latitude']
            st_long=inventory_5G.get_coordinates('5G.{}..EHZ'.format(station_5G))['longitude']
            st_elevation=inventory_5G.get_coordinates('5G.{}..EHZ'.format(station_5G))['elevation']
        except:
            pass
        ##
        t_arr= get_t_arr(event_lat,event_long,event_depth,st_lat, st_long)
        inventory_1 = client_5g.get_stations(network="5G",station=station_5G,level='response')
        st = client_5g.get_waveforms('5G',station_5G,'','*',starttime=t+t_arr-15,\
        endtime=t+t_arr+20,attach_response=True)
        P_arr=t+t_arr ### P arrival in UTCDateTime
        st1=st.copy()
        tr_z= st1.select(component="Z")[0]
        tr_z=pre_process(tr_z,inventory_1,10)
        # tr_z.filter('bandpass', freqmin=1, freqmax=10)
        plot_wave(tr_z,P_arr) # plots waveform usign the defined function
        FM=input("Enter the first motion: ")
        fm_info.append([station_5G,FM,st_lat,st_long,st_elevation])
        st_all.append(tr_z)
        # plt.close()
    except:
        pass
### checking the three AU stations for all eq's
# sys.exit()
for station_au in ['MULG','OOD','LCRK']:
# for station_au in ['OOD']:
    try:
        st_lat=inventory_au.get_coordinates('AU.{}..BHZ'.format(station_au))['latitude']
        st_long=inventory_au.get_coordinates('AU.{}..BHZ'.format(station_au))['longitude']
        st_elevation=inventory_au.get_coordinates('AU.{}..BHZ'.format(station_au))['elevation']
        try:
            response = inventory_au.get_response('AU.{}..BHZ'.format(station_au),t)
        except:
            pass
        try:
            response = inventory_au.get_response('AU.{}.00.HHZ'.format(station_au),t)
        except:
            pass
        ##
        t_arr= get_t_arr(event_lat,event_long,event_depth,st_lat, st_long)
        inventory_1 =  client.get_stations(network="AU",station=station_au,level='response')
        st = client.get_waveforms('AU',station_au,'*','BH*',starttime=t+t_arr-15,endtime=t+t_arr+20,attach_response=True)
        P_arr=t+t_arr ### P arrival in UTCDateTime
        st1=st.copy()
        tr_z= st1.select(component="Z")[0]
        tr_z=pre_process(tr_z,inventory_1,10)
        # tr_z.filter('bandpass', freqmin=1, freqmax=10)
        plot_wave(tr_z,P_arr) # plots waveform usign the defined function
        FM=input("Enter the first motion: ")
        fm_info.append([station_au,FM,st_lat,st_long,st_elevation])
        st_all.append(tr_z)
    # plt.close()
    except:
        pass
# sys.exit()

### saving fm info for each earthquake
## out of window
# AEB13, AES04 eq 29
####
output_f=os.path.join('fm_info/', '')
with open(output_f+'Eq_'+str(eq_num)+'.txt','w') as fi:
# with open(output_f+'Eq_'+str(eq_num)+'_4-10Hz.txt','w') as fi:
    for list in fm_info:
        fi.write(list[0]+' '+list[1]+' '+str(list[2])+' '+str(list[3])+' '+str(list[4])+'\n')

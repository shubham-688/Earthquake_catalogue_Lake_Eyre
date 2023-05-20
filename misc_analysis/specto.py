import matplotlib.pyplot as plt
#plt.switch_backend("nbagg")
#plt.style.use('ggplot')
import obspy
from obspy import read, Stream
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import locations2degrees
import datetime
from obspy import UTCDateTime
from obspy.taup import TauPyModel

#12 Event(s) in Catalog:
#event_time = "2018-11-02T21:27:31.54 | -32.100, +138.767 | 2.25 ML" #306
#event_time = "2018-11-16T01:29:46.24 | -31.690, +138.568 | 3.74 mb" #320 #1700
event_time = "2018-11-27T12:54:14.20 | -29.283, +137.307 | 3.96 mb" #331
event_time = "2018-12-02T22:29:50.83 | -29.470, +137.472 | 2.34 ml" #336
#event_time = "2018-12-05T21:52:47.69 | -29.833, +141.526 | 3.49 mb"
#event_time = "2018-12-06T12:41:45.77 | -30.562, +138.201 | 3.45 mb" #340
#event_time = "2018-12-11T13:26:35.12 | -32.033, +139.203 | 3.25 mb" #345


#st = obspy.read("2018/322/AES16*.EHZ")
#Julian_day
fmt = '%Y-%m-%dT%H:%M:%S'
#s = '2018-12-29T08:33:21'
dt = datetime.datetime.strptime(event_time[:19], fmt)
tt = dt.timetuple()
Jday=tt.tm_yday
year=tt.tm_year
month=tt.tm_mon
day=tt.tm_mday
print(Jday)
min=tt.tm_min
hr=tt.tm_hour
sec=tt.tm_sec

if hr<10:
            hr='%0*d' % (2, hr)
print(int(hr))


utc_t = UTCDateTime(year, month,day , int(hr), min,sec )

t=event_time[:19]
print(t)
#print(t)
###############
AES=['AES20','AES02','AES07', 'AES08', 'AES03', 'AES16','AES11','AES10','AES04','AES15', 'AES12','AES13']
AEB=['AEB02','AEB20','AEB18', 'AEB11', 'AEB17','AEB16','AEB01','AEB12','AEB13', 'AEB15'] #,'AEB14']
list=['AEB01','AEB16','AES12']
########
st = Stream()
for x in list:       ###### EDIT here for changing broadband to short term
	st += obspy.read('../2018/{}/{}*{}0000.EHZ' .format(Jday,x,hr))
########
# Earthquakes' epicenter
lat_long=[float(event_time[25:31]),float(event_time[34:42])]
print(lat_long)
### Reading the station lat_long_name
infile=open('../station/stations.txt',"r")
lat_ori=[]
long_ori=[]
station=[]
for line in infile:
    words = line.split()
    lat_ori.append(float(words[0]))
    long_ori.append(float(words[1]))
    station.append(words[2])
infile.close()
#### Finding the index of station from the station_list
lat=[]
long=[]
for x in list:     ###### EDIT here for changing broadband to shortterm
    for l in range(len(station)):
        if station[l] == x:
            lat.append(lat_ori[l])
            long.append(long_ori[l])
            break
#print(lat)

# Calculating distance from SAC headers lat/lon
i=0
for tr in st:
    tr.stats.distance = gps2dist_azimuth(lat[i], long[i],
                                         lat_long[0], lat_long[1])[0]
    tr.stats.distance_degree = locations2degrees(lat[i], long[i],
                                         lat_long[0], lat_long[1])
    i=i+1
### TAUP
m = TauPyModel(model="ak135")

### FILTERS
#st.filter('bandpass', freqmin=.1, freqmax=20) #AES .1-2
#st.filter('bandpass', freqmin=3, freqmax=20) #AEB 1-10
st.filter('bandpass', freqmin=1, freqmax=40) #AES 3-20 local
#st.filter('highpass', freq=4, corners=4)
print(st)
i=0
for tr in st:
    arrivals = m.get_ray_paths(distance_in_degree=tr.stats.distance_degree,
        source_depth_in_km=10)
    #print(arrivals)
    print(tr.stats.distance_degree)
    first_arrival = utc_t + arrivals[0].time
    tr2 = tr.copy()
    tr2.trim(first_arrival-10, first_arrival+100)
    #print(st2)
    print(first_arrival)
    # tr2.plot(size=(600,300),linewidth=.05,outfile='specto/{}_{}.pdf'.format(Jday,AEB[i]),format='pdf')
    # tr2.spectrogram(log=True,wlen=1,title='5G.{}Z#{} - Record from: {}'.format(AEB[i],Jday,first_arrival-10)
    # ,outfile='specto/Spec_{}_{}.png'.format(Jday,AEB[i]),fmt='png')
    tr2.plot(size=(600,300),linewidth=.05,outfile='specto/{}_{}.pdf'.format(Jday,list[i]),format='pdf')
    tr2.spectrogram(size=(600,300),log=True,wlen=1,per_lap=0.5,title='5G.{}Z- Time after {}'.format(list[i],first_arrival-10)
    ,outfile='specto/Spec_{}_{}.png'.format(Jday,list[i]),fmt='png')
    i=i+1
print(i)
#tr = st[11]

#print(tr.stats)

#st.plot(outfile='../../../Dropbox/{}_AES_hr.pdf'.format(Jday),format='pdf')
#fig = plt.figure()
#ax = plt.subplot(211)
#st2.plot(size=(800,400),linewidth=.3,fig=fig)
#ymin, ymax = ax.get_ylim()
#plt.vlines(t, ymin, ymax, linestyles='dashed',color='r', linewidth=2)
#plt.axvline(x=t,linestyle='dashed',color='r', linewidth=2)
#st2.plot(size=(800,400),linewidth=.4,outfile='../../../Dropbox/{}_AES.pdf'.format(Jday),format='pdf')
#plt.subplot(212, sharex=ax)
#st2.spectrogram(log=True,show=False)
#plt.savefig('specto/{}_AES20.pdf'.format(Jday),format='pdf')

#for tr in st2:
#st2.spectrogram(log=True,wlen=2)#log=True,,outfile='specto/Spec_{}_AES12.pdf'.format(Jday,AES),fmt='pdf')
#

#st.plot(type='dayplot')

#################
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#ax.plot(tr.times("matplotlib"), tr.data, "b-")
#ax.xaxis_date()
#fig.autofmt_xdate()
#plt.show()

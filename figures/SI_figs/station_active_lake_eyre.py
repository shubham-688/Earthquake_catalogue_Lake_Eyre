import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.clients.fdsn import Client
import numpy as np
import sys
#######################
def get_no_of_days(station,year):
    file='list_{}.txt'.format(year)
    count=0
    for line in open(file, "r"):
        line=line.split()[0]

        if line[10:15]==station:
            count=count+1
    return count


client_aus = Client('http://auspass.edu.au:8080',user='5g',password='grape71')
inventory = client_aus.get_stations(network='5G', station='*', level='channel')


station_name=[]
for i in range(len(inventory.get_contents()['stations'])):
    station_name.append(inventory.get_contents()['stations'][i].split()[0].split('.')[1])

####################
list_18=[]
list_19=[]
list_20=[]
list_21=[]
list_22=[]

for station in station_name:
    count=get_no_of_days(station,18)
    list_18.append([station,count])
    ##
    count=get_no_of_days(station,19)
    list_19.append([station,count])
    ##


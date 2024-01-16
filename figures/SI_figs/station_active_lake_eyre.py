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
    count=get_no_of_days(station,20)
    list_20.append([station,count])
    ##
    count=get_no_of_days(station,21)
    list_21.append([station,count])
    ##
    count=get_no_of_days(station,22)
    list_22.append([station,count])
    ##
sys.exit()



# list.sort(key = lambda i: i[1]) #sorts list based on second element which is Jday
name=[i[0] for i in list_18] # gets name back from the list
day_18=[i[1] for i in list_18]
day_19=[i[1] for i in list_19]
day_20=[i[1] for i in list_20]
day_21=[i[1] for i in list_21]
day_22=[i[1] for i in list_22]


# sys.exit()
name_position = np.arange(len(name))  # the label locations
width = 0.3  # the width of the bars

fig, ax = plt.subplots(figsize=(22,8))
rects_18=ax.bar(name_position-.3,day_18,width/2,color='powderblue',alpha=1,label='2018')
rects_19=ax.bar(name_position-.15,day_19,width/2,color='DarkSalmon',alpha=1,label='2019')
rects_20=ax.bar(name_position,day_20,width/2,color='gray',alpha=.55,label='2020')
rects_21=ax.bar(name_position+.15,day_21,width/2,color='seagreen',alpha=.7,label='2021')
rects_22=ax.bar(name_position+.3,day_22,width/2,color='slateblue',alpha=.45,label='2022')

ax.set_ylabel('# of days',fontsize=16)
ax.set_xlabel('Stations',fontsize=16)
ax.set_yticks(np.arange(0, 380, step=20))
ax.set_xticks(name_position)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
ax.set_xticklabels(name)
ax.set_title('# of days stations active in 2018-2022')
ax.legend()
plt.savefig('Station_run_18_22.pdf',bbox_inches='tight', pad_inches=0.2)
plt.show()
# ax.bar_label(rects_18, padding=3)
# ax.bar_label(rects_19, padding=3)
for rect in rects_18:
    height = rect.get_height()
    ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
for rect in rects_19:
    height = rect.get_height()
    ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
for rect in rects_20:
    height = rect.get_height()
    ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
for rect in rects_21:
    height = rect.get_height()
    ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
for rect in rects_22:
    height = rect.get_height()
    ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

fig.tight_layout()
# plt.savefig('Station_run_18_20.pdf',bbox_inches='tight', pad_inches=0.2)
plt.show()

import numpy as np
import os
import datetime
#from obspy import UTCDateTime
###########
fn=1
while fn < 90:

    eq_folder= os.path.join('eq{}/'.format(fn), '')


    temp=[]
    with open(eq_folder+'shnamod.out', 'r') as f:
        for line in f:
            line=line.split()
            if len(line)>0:
                if line[0]=='Latitude':
                    temp.append(line[1])
                if line[0]=='Longitude':
                    temp.append(line[1])
                if line[0]=='Depth':
                    temp.append(line[1])
                if line[0]=='RESIDUAL:':
                    temp.append(line[1])
                    temp.append(line[2])

    with open('relocated_eq_na.txt','a') as out:
        #out.write('\n')
        out.write(' '.join(temp))
        out.write('\n')
    temp=[]
    fn+=1

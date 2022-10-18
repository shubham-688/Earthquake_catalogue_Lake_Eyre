import numpy as np
import obspy
import sys
import datetime
from obspy import UTCDateTime
import os
###################################
#This script gets the stations which recorded a particualr event
##########
events_to_read=50 # 

blk_begin = '  10.0000'#      10.0000'
blk_end = 'end'
eq_num=1
while eq_num < 51:
    eq_folder = os.path.join('eq{}/'.format(eq_num), '')
    data=[]
    write_block = False
    for line in open(eq_folder+'le.evd', 'r'):
   
        line=line.split()
        if line[0][1] in ['S','B']:
            name='A'+line[0]
            data.append(name)
        if line[0] in ['OOD','LCRK','MULG']:
            data.append(line[0])
    mag_fold=os.path.join('mag_info','')
    with open(mag_fold+'station_'+str(eq_num)+'.txt', 'w') as fout:
        fout.write('\n'.join(data))
        #data=[]
    eq_num +=1




    



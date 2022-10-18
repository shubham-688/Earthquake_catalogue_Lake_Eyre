import numpy as np
import obspy
import sys
import datetime
from obspy import UTCDateTime
import os
###################################
#This script gets the lat long depth time of the relocated event
##########
events_to_read=50 # 

blk_begin = '  DATE:'#      10.0000'
blk_end = '  RESIDUAL:'
eq_num=1
while eq_num < 51:
    eq_folder = os.path.join('eq{}/'.format(eq_num), '')
    data=[]
    write_block = False
    for line in open(eq_folder+'shnamod.out', 'r'):
        if line[:7] == blk_begin:
            write_block = True
        if write_block:
            data.append(line)
        if line[:11] == blk_end:
            #data.append(line)
            break
    
    mag_fold=os.path.join('mag_info','')
    with open(mag_fold+'eq'+str(eq_num)+'.txt', 'w') as fout:
        date=data[1].split()[1]+'-'+data[2].split()[1]+'-'+data[3].split()[1]+'T'
        time=data[4].split()[1]+':'+data[5].split()[1]+':'+data[6].split()[1]
        str=data[8].split()[1]+' '+data[9].split()[1]+' '+data[10].split()[1]+' '+data[11].split()[1]
        final=str+' '+date+time
        fout.write(final)
    eq_num +=1

###########

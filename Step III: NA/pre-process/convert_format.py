import obspy
import numpy as np
import datetime
from obspy import UTCDateTime
import os
#############

def myround(x, base=.5):
    return base * round(x/base)
#################################
eq_output="le_all_eq_10k.evd"
#mine_output="le_all_mine.evd"
ff = open(eq_output,'w')
#ff_m = open(mine_output,'w')

with open("phase_all_NA_format.txt", "r") as f:
    for line in f:
        a=line.split()
 #       if float(a[9]) < -29 and float(a[10]) < 136: # mine event...set depth to zero

        if a[0]=='#':
                ff.write('end\n')
                # print(a)
                time=UTCDateTime("{}-{}-{}T{}".format(a[1],a[3],a[4],a[2]))
                ff.write(' EVENT Lake Eyre {} {:0>2} {:0>2}\n'.format(time.year, time.month, time.day))
                ff.write('5G\n')
                ff.write('{}  {:0>2}   {:0>2} {:0>2} {:0>2}  {:.2f}  {:.2f}\n'.format(time.year, time.month, time.day, time.hour, time.minute, round(time.second+time.microsecond*10e-7,2), float(a[8])))
                ff.write('  {:.3f}       {:.4f}\n'.format(float(a[9]),1))
                ff.write('  {:.3f}       {:.4f}\n'.format(float(a[10]),1))
                ff.write('  {:.4f}      {:.4f}\n'.format(15,15)) # original for 10km: a[11],10; and one less space.

        else:
                time_=time+float(a[2])
                if round(abs(float(a[3])),2) == 0 or round(abs(float(a[3])),2) < 0.1:
                    a[3]=0.1
                ff.write('{:<4} {}         {:0>2} {:0>2} {:0>5.2f}  {:.2f}\n'.format(a[0][-4:],a[1],time_.hour,time_.minute,round(time_.second+time_.microsecond*10e-7,2),round(abs(float(a[3])),2)))





ff.write('end')
ff.close()
#ff_m.close()

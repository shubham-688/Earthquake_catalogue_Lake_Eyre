import obspy
import numpy as np
import datetime
from obspy import UTCDateTime
import os
#################################
#scp shubham@es09156:/Users/Shubham/Research/Lake_eyre_data/NA/SHAKENA/pre_process/2018_19/eq* .
#scp -r shubham@es09156:/Users/Shubham/Research/Lake_eyre_data/NA/SHAKENA/pre_process/2018_19/eq* .
# eq_input="le_all_18_19.evd"
# after making a le_all_18_19.evd file, you have to remove the first end from it.
#next step is to get relocated eq info from shmamod.cmd
blocks_to_read = 19 # number of earthquakes..see phase_2018_19_NA_format.txt file
blk_begin = ' EVENT'
blk_end = 'end'

with open('le_all_20.evd', 'r') as f:
    fn = 1
    data = []
    write_block = False
    for line in f:
        if fn > blocks_to_read:
            break
        # line = line.strip()
        if line[:6] == blk_begin:
            write_block = True
        if write_block:
            data.append(line)
        if line[:3] == blk_end:
            write_block = False
            eq_folder= os.path.join('eq{}/'.format(fn), '')
            os.mkdir(eq_folder)
            with open(eq_folder+'le'+'.evd', 'w') as fout:
                fout.write(''.join(data))
                # fout.write(data)

                # fout.write('end')
                data = []
            fn += 1

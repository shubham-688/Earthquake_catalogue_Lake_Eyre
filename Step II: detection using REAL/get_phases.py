import obspy
from obspy import read, Stream
import numpy as np
import datetime
from collections import deque
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import os
import sys
import shutil

# this script picks the relevant events and phases out of phase_all.txt file.
# I have to use this coz i implement gap >180 after REAl, therefore, I need to
# sort the phase_sel.txt as well.
# I do by reading both files, findind common lines and printing P+S pick stations..
##########################
# event_file=open("catalog/final_cat_good/cat_sg200_2018_19.txt","r")
# phase_all=open("catalog/2018_19/phase_all_2018_19.txt","r")

#' 4   2020 01 17 19:15:04.408    69304.408   1.4879     -29.7347     135.4510      10.0000     -inf      nan    8   10   18   177.39'
# now, I do all three years, 2018,19, 20 in one go.
phase_select="catalog/final_cat_good/phase_all_sg200.txt"
phase_select_mine="catalog/final_cat_good/phase_all_sg200_mine.txt"
f = open(phase_select,'w')
fm = open(phase_select_mine,'w')
for line1 in open("catalog/final_cat_good/cat_sg200_ALL.txt","r"):
    i=0
    j=0
    
    #for line2 in open("catalog/phase_all_2020.txt","r"):
    if float(line1.split()[7]) < -29 and float(line1.split()[8]) < 136:
        for line2 in open("catalog/phase_ALL.txt","r"):
            if line1.split()==line2.split():
                i=int(line1.split()[14])+1 # gets the number of P+S picks from line1
                # print(line2,outfile=phase_select)
                # f.write('{}'.format(line2))
            if j<i:
            # print(line2,outfile=phase_select)
                fm.write('{}'.format(line2))
                j+=1
            if j==i and j!=0:
                break
            # i=0
    else:
        for line2 in open("catalog/phase_ALL.txt","r"):
            if line1.split()==line2.split():
                i=int(line1.split()[14])+1 # gets the number of P+S picks from line1
                # print(line2,outfile=phase_select)
                # f.write('{}'.format(line2))
            if j<i:
            # print(line2,outfile=phase_select)
                f.write('{}'.format(line2))
                j+=1
            if j==i and j!=0:
                break
            # i=0

# event_file.close()
# phase_all.close()
f.close()

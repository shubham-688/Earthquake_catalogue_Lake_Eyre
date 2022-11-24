### script for fm usin Hash_py
from hashpy import HashPype, HashError, libhashpy
from hashpy.plotting import FocalMechPlotter
import numpy as np
import obspy
# from obspy.core.event import readEvents
import os
import matplotlib.pyplot as plt
import sys
from obspy.core.event import (Catalog, Event, Origin,Pick,
QuantityError, Arrival, FocalMechanism, MomentTensor, NodalPlanes,
 PrincipalAxes, Axis, NodalPlane, SourceTimeFunction, Tensor,read_events,Magnitude)
import importlib
import full_test_event as fl
from obspy.imaging.beachball import beachball
from obspy.imaging.beachball import beach
importlib.reload(fl)

#
# def plot_fm(event_out,ax):
#     n_fm=len(event_out.focal_mechanisms)
#     for i in range(n_fm):
#         _fm=event_out.focal_mechanisms[i]
#         [strike, dip, rake]=[_fm.nodal_planes.nodal_plane_1.strike,
#         _fm.nodal_planes.nodal_plane_1.dip,_fm.nodal_planes.nodal_plane_1.rake]
#         if rake < 0:
#             rake=rake+360
#         sdr=[strike, dip, rake]
#         # fig.add_subplot(1,n_fm,i+1)
#         beach1=beach(sdr,linewidth=1, facecolor='steelblue',edgecolor='navy',xy=((i)*40,40),width=30)
#         ax.add_collection(beach1)
#
#     return ax
######################################################
# event = read_events('output_.xml').events[0]
# round 1
line='-26.198 135.789 9.567 0.264 2019-08-27T23:42:30.691 1.928 1.996 1.969' #21

print('--',line,'--\n')
# %reset -f

eq_num=21
print('--EQ number:',eq_num,'--\n')
event=fl.full_test_event(line,eq_num)
# Set configuration at creation with a dict...
# ...can from file or interactively, etc
config = { "npolmin" : 7,
           "max_agap": 200,
           "vmodels" : ['ausrem_LE.txt',
                        'ausrem_AU.txt',
                       ]
           }

hp = HashPype(**config)
hp.input(event, format="OBSPY")
hp.load_velocity_models()
hp.generate_trial_data()
hp.calculate_takeoff_angles()

pass1 = hp.check_minimum_polarity()
pass2 = hp.check_maximum_gap()

if pass1 and pass2:
    hp.calculate_hash_focalmech()
    hp.calculate_quality()
    print(hp.output()) # default output is a simple string
else:
    raise HashError("Didn't pass user checks!")

qual=hp.qual[0].decode('UTF-8')
event_out = hp.output(format="OBSPY")
ax=[]
f=plt.figure(figsize=(2.5,1.5))
# plt.clf()
plt.plot([40, 40],[0, 40], "rv", ms=1,alpha=.1)
ax = plt.gca()
n_fm=len(event_out.focal_mechanisms)
for i in range(n_fm):
    _fm=event_out.focal_mechanisms[i]
    [strike, dip, rake]=[_fm.nodal_planes.nodal_plane_1.strike,
    _fm.nodal_planes.nodal_plane_1.dip,_fm.nodal_planes.nodal_plane_1.rake]
    if rake < 0:
        rake=rake+360
    sdr=[strike, dip, rake]
    # fig.add_subplot(1,n_fm,i+1)
    beach1=beach(sdr,linewidth=.75, facecolor='steelblue',edgecolor='navy',xy=((i)*40,40),width=30,alpha=.9)
    ax.add_collection(beach1)

xmax=i*40+20
ax.set_xlim((-20, xmax))
ax.set_ylim((20, 60))

plt.axis('off')
ax.set_title('Earthquake #{}'.format(eq_num))
plt.savefig('#{}.jpg'.format(eq_num),dpi=200,bbox_inches='tight', pad_inches=0.2)
plt.show()

print('\n',event_out.focal_mechanisms)
with open('{}.txt'.format(eq_num), 'w') as f:
    for line in event_out.focal_mechanisms:
        f.write("%s\n" % line)

fmp = FocalMechPlotter(event_out)


# sys.exit()
#
# a=hp.output().split()
# strike=float(a[3])
# dip=float(a[5])
# rake=float(a[7])
# strike2=event_out.preferred_focal_mechanism().nodal_planes.nodal_plane_2.strike
# dip2=event_out.preferred_focal_mechanism().nodal_planes.nodal_plane_2.dip
# rake2=event_out.preferred_focal_mechanism().nodal_planes.nodal_plane_2.rake
# ###
# np1 = [strike, dip, rake]  ## ak
# np2 = [strike2, dip2, rake2]  ## ak
#
# fig = beachball(np1)
# fig2 = beachball(np2)
#
#
# # event_out2=hp.outputOBSPY(hp,event)
# mags=Magnitude(magnitudes=[])
# mags.mag=hp.qmag
# event_out.magnitudes.append(mags)
# print('----- number of solution------',hp.nmult,'\n')
# print('-----Misfit-----',event_out.preferred_focal_mechanism().misfit)


### following code to check take off angle
# for i in range(hp.nmult):
#     print(event_out.focal_mechanisms[i].nodal_planes)

# #
# for ii, arrv in enumerate(event.preferred_origin().arrivals):
#     i=arrv.distance
#     arrivals = model.get_travel_times(source_depth_in_km=10, distance_in_degree=i, phase_list=['p'])
#     arr=arrivals[0]
#     print(i,'deg -----',arr.takeoff_angle,'deg')
#     # print(arrv.distance,'---dist---',arrv.azimuth)

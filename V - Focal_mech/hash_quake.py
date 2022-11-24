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
line='-25.664 137.076 6.350 0.499 2019-10-4T16:35:22.267 3.589 3.765 3.678' #27
line='-28.227 135.887 9.674 0.307 2019-11-5T21:6:10.162 1.17 1.235 1.212'#31
line='-27.882 135.704 10.934 0.418 2020-02-22T9:50:25.164 3.248 3.217 3.222' #40
line='-28.101 135.842 9.931 0.224 2020-03-29T18:54:29.939 1.592 1.792 1.641' #45
line='-27.568 136.662 9.934 0.228 2020-03-30T8:21:15.855 1.286 1.613 1.372' #46

### Didn't pass user checks
# line='-30.008 138.308 9.008 0.256 2021-04-15T9:52:22.269 3.032 3.106 3.064' #24
# line='-30.394 137.061 6.459 0.356 2021-05-29T20:54:21.507 1.786 1.667 1.726' #33
# line='-30.618 137.001 9.443 0.312 2021-07-9T4:30:7.384 1.734 1.755 1.744' #42
# line='-29.752 138.037 6.940 0.229 2021-09-10T20:24:15.993 1.835 1.613 1.653' #53
# line='-27.944 137.821 9.772 0.341 2021-09-25T2:30:50.552 1.598 2.066 1.791' #56
# ##
# line='-28.100 135.703 8.873 0.259 2021-01-1T17:33:57.575 1.822 1.895 1.844' #10
# line='-28.748 135.702 10.450 0.270 2021-01-20T23:33:13.321 2.453 2.478 2.455' #11
# line='-28.695 136.310 9.626 0.205 2021-02-25T6:48:42.411 1.913 2.143 1.998' #18
# line='-28.226 135.941 9.917 0.441 2021-06-20T11:44:23.531 1.785 1.949 1.828' #38
# line='-28.243 136.442 5.166 0.237 2021-07-7T6:43:20.114 1.937 2.084 1.973' #41
# line='-28.489 136.468 9.749 0.353 2021-08-31T12:52:31.042 1.791 1.8 1.779' #51
# line='-28.414 136.051 10.504 0.252 2021-10-18T18:44:41.309 1.9 1.886 1.873' #59
# line='-28.855 135.639 9.115 0.161 2021-12-23T1:29:2.961 1.557 1.694 1.725' #68
# line='-26.511 137.318 7.342 0.578 2022-02-14T22:40:28.263 2.608 2.525 2.542' #75
# line='-28.499 136.193 9.897 0.389 2022-03-6T10:4:13.020 2.777 2.719 2.675' #77
# line='-28.502 136.209 10.320 0.195 2022-03-6T10:35:39.865 2.116 2.121 2.113' #78

print('--',line,'--\n')
# %reset -f

eq_num=46
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

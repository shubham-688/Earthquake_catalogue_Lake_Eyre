import math
import obspy.taup
import numpy as ny
from obspy.taup import TauPyModel
from obspy.taup.taup_create import build_taup_model
import sys

###### build a TauPyModel and cpy it to obspy directory..
# build_taup_model("LA_moho45.nd",'.') # '.' means .npz file is build in same directory
# sudo cp LA_moho45.npz /opt/local/Library/Frameworks/Python.framework/Versions/3.5/lib/python3.5/site-packages/obspy/taup/data/
 #### copying the .npz file to the obspy repository

# model = TauPyModel(model="LA_moho35")
# model = TauPyModel(model="itvel")

model = TauPyModel(model="LA_moho45")
# sys.exit()
with open("ttdb_moho45.txt", "w") as f:
    # f.write("dist dep tp ts tp_slowness ts_slowness tp_hslowness ts_hslowness p_elvecorr s_elvecorr\n")
    # for dep in range(0,11,2): # search 0 - 20 km in depth ori: range(0,21,1)
    for dep in ny.arange(0,11,10):
        for dist in range(5,551,5): # search 0.01 - 1.4 degree in horizontal ori:range(1,141,1)
            dist = dist*0.01
            print(dep,dist)
            arrivals = model.get_travel_times(source_depth_in_km=dep, distance_in_degree=dist, phase_list=["P","p","S","s"])
            #print(arrivals)
            i = 0
            pi = 0
            si = 0
            while(i<len(arrivals)):
                arr = arrivals[i]
                i = i + 1
                if((arr.name == 'P' or arr.name == 'p') and pi == 0):
                    pname = arr.name
                    p_time = arr.time
                    p_ray_param = arr.ray_param*2*ny.pi/360
                    p_hslowness = -1*(p_ray_param/111.19)/math.tan(arr.takeoff_angle*math.pi/180)
                    pi = 1

                if((arr.name == 'S' or arr.name == 's') and si == 0):
                    sname = arr.name
                    s_time = arr.time
                    s_ray_param = arr.ray_param*2*ny.pi/360
                    s_hslowness = -1*(s_ray_param/111.19)/math.tan(arr.takeoff_angle*math.pi/180)
                    si = 1
                if(pi == 1 and si == 1):
                    break

            f.write("{} {} {} {} {} {} {} {} {} {}\n".format(dist, dep, p_time,s_time, p_ray_param, s_ray_param, p_hslowness, s_hslowness, pname, sname))

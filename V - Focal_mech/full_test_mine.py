from obspy import read_events, Catalog, UTCDateTime, read
from obspy.core.event import Pick, WaveformStreamID, Arrival, Amplitude
from obspy.core.event import Event, Origin, Magnitude
from obspy.core.event import EventDescription, CreationInfo, OriginUncertainty
from obspy.core.event import ConfidenceEllipsoid
from obspy.core.util.base import NamedTemporaryFile
import os
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel

########
model = TauPyModel(model="iasp91")

def get_t_arr(event_lat,event_long,event_depth,st_lat, st_long):
    epi_dist, az, baz = gps2dist_azimuth(float(event_lat),float(event_long), st_lat, st_long)
    epi_dist = epi_dist / 1000
    arrival_st=model.get_travel_times(source_depth_in_km=float(event_depth),\
    distance_in_degree=epi_dist/111.19,phase_list=["P","p"])
    t_arr=arrival_st[0].time
    return t_arr,az,epi_dist/111.19

##


def full_test_event(line,number):
    """
    Function to generate a basic, full test event
    """
    line=line.split()
    test_event = Event()
    test_event.origins.append(Origin())
    test_event.origins[0].time = UTCDateTime(line[4])
    test_event.event_descriptions.append(EventDescription())
    test_event.event_descriptions[0].text = 'LE_{}'.format(number)
    test_event.origins[0].latitude = float(line[0])
    test_event.origins[0].longitude = float(line[1])
    test_event.origins[0].depth = float(line[2])
    # test_event.origins[0].depth = 0

    test_event.origins[0].creation_info = CreationInfo(author="SA_le", version=number)
    ouc = OriginUncertainty()
    ce = ConfidenceEllipsoid()
    ouc.confidence_ellipsoid = ce
    test_event.origins[0].origin_uncertainty = ouc
    test_event.origins[0].origin_uncertainty.confidence_ellipsoid.semi_major_axis_length=.5
    test_event.origins[0].origin_uncertainty.confidence_ellipsoid.semi_intermediate_axis_length=.5

    test_event.origins[0].time_errors['Time_Residual_RMS'] = float(line[3])
    test_event.preferred_origin_id = test_event.origins[0].resource_id

    test_event.magnitudes.append(Magnitude())
    test_event.magnitudes[0].mag = float(line[7])
    test_event.magnitudes[0].magnitude_type = 'ML'
    test_event.magnitudes[0].creation_info = CreationInfo(author="SA_le", version=number)
    test_event.magnitudes[0].origin_id = test_event.origins[0].resource_id
    # test_event.preferred_magnitude_id = test_event.magnitudes[0].origin_id

    # Define the test pick
    polarity_file = os.path.join('../mines_fm_info/', '')
    line_num=0
    for st_line in open(polarity_file+'mine_'+str(number)+'.txt','r'): # eq details
        st_line=st_line.split()
        if st_line[0] in ('LCRK','OOD','MULG'):
            network='AU'
        else:
            network='5G'

        _waveform_id=WaveformStreamID(station_code=st_line[0], channel_code='SHZ',
                                          network_code=network)


        t_arr,az,epi_dist= get_t_arr(float(line[0]),float(line[1]),float(line[2]),float(st_line[2]), float(st_line[3]))
        if st_line[1] in ('N'):
            pol_st='undecidable'
        if st_line[1] in ('U','+'):
            pol_st='positive'
        if st_line[1] in ('D','-'):
            pol_st='negative'
        # print(pol_st)
        test_event.picks.append(
            Pick(waveform_id=_waveform_id, onset='impulsive', phase_hint='P',
                 polarity=pol_st, time=UTCDateTime(line[4]) + t_arr,
                 evaluation_mode="manual",creation_info=CreationInfo(author="SA_le", version=number)))
            # Unassociated pick
    ####
        test_event.origins[0].arrivals.append(
            Arrival(time_weight=0, phase=test_event.picks[line_num].phase_hint,
                    pick_id=test_event.picks[line_num].resource_id,distance=epi_dist,
                    azimuth=az))
        line_num+=1
    # test_event.origins[0].arrivals.append(
    #     Arrival(time_weight=2, phase=test_event.picks[1].phase_hint,
    #             pick_id=test_event.picks[1].resource_id,
    #             backazimuth_residual=5, time_residual=0.1, distance=2,
    #             azimuth=28))

    return test_event

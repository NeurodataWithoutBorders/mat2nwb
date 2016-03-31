#!/usr/local/python-2.7.8/bin/python 

# encoding: utf-8
"""
make_nwb.py

Created by Claudia Friedsam on 2015-02-08.
Redesigned by Gennady Denisov on 2016-03-28.

"""

import sys, os
import nwb  
import nwb.nwbts as ts
import nwb.nwbmo as ms
import nwb.nwbep as ep
import nwbco 
import h5py
import datetime
import getpass
import numpy as np
from sets import Set
import re 
import optparse
import h5lib 

# ------------------------------------------------------------------------------

def make_nwb_command_line_parser(parser):
    parser.add_option("-D", "--debug",     action="store_true", dest="debug", help="output debugging info", default=False)
    parser.add_option("-o", "--outfolder", dest="output_folder", help="output folder (default=same as input folder)",metavar="output_folder",default="")
    parser.add_option("-r", "--replace",   action="store_true", dest="replace", help="if the output file already exists, replace it", default=False)
    parser.add_option("-v", "--verbose",   action="store_true", dest="verbose", help="increase the verbosity level of output", default=False)

    return parser

# ------------------------------------------------------------------------------

def check_entry(file_name,obj):
    try:
        return file_name[obj]
    except KeyError:
        print str(obj) +" does not exist"
        return []

# ------------------------------------------------------------------------------

def parse_h5_obj(obj, level = 0, output = [], verbose = 0):
    if level == 0:
        output = []
    try:
        if isinstance(obj, h5py.highlevel.Dataset):
            level = level+1
            if obj.value.any():
                output.append(obj.value)
            else:
                full_name = obj.name.split("/")
                output.append(full_name[-1])
        elif isinstance(obj, h5py.highlevel.Group):
            level = level+1
            if not obj.keys():
                output.append([])
            else:
                for key in obj.keys():
                    parse_h5_obj(obj[key], level, output, verbose)
        else:
            output.append([])
    except KeyError:
        print "Can't find" + str(obj)
        output.append([])
    return output

# ------------------------------------------------------------------------------

# each of nwb_object's hdf5 files have imaging planes and subareas
#   labels consistent within the file, but inconsistent between
#   files. create a map between the h5 plane name and the 
#   identifier used between files
# plane_map = {}
def add_plane_map_entry(plane_map, h5_plane_name, filename):
    toks = filename.split("fov_")
    if len(toks) != 2:
        print "Error parsing %s for imaging plane name" % filename
        sys.exit(1)
    univ_name = "fov_" + toks[1][:5]
    if univ_name not in plane_map:
        #print filename + " -> " + univ_name
        plane_map[h5_plane_name] = univ_name
    return univ_name

# ------------------------------------------------------------------------------

def create_plane_map(orig_h5, plane_map):
    num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    for subarea in range(num_subareas):
        # fetch time array
        grp = orig_h5['timeSeriesArrayHash/value/%d/imagingPlane' %(subarea + 2)]
        grp2 = orig_h5['timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)]
        if grp2.keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(grp2.keys())
        for plane in range(num_planes):
            # if num_planes == 1:
            #     pgrp = grp
            # else:
            #     pgrp = grp["%d"%(plane+1)]
            pgrp = grp["%d"%(plane+1)]
            old_name = "area%d_plane%d" % (subarea+1, plane+1)
            frame_idx = pgrp["sourceFileFrameIdx"]["sourceFileFrameIdx"].value
            # lst = pgrp["sourceFileList"]
            lst = parse_h5_obj(pgrp["sourceFileList"])[0]
            for k in lst:
                # srcfile = str(lst[k][k].value)
                srcfile = str(k)
                add_plane_map_entry(plane_map, old_name, srcfile)
                break

# ------------------------------------------------------------------------------

# fetch start time/date of experiment
def find_exp_time(input_h5):
    # SP data
    if "metaDataHash" in h5lib.get_child_group_names(input_h5):
        d = h5lib.get_value_by_key(input_h5['/metaDataHash'], \
                         "dateOfExperiment")
        t = h5lib.get_value_by_key(input_h5['/metaDataHash'], \
                         "timeOfExperiment")
#       print "d=", d
#       print "t=", t
        dt=datetime.datetime.strptime(d+t, "%Y%m%d%H%M%S")
    # NL data
    elif "dateOfExperiment" in h5lib.get_child_group_names(input_h5) and \
         "timeOfExperiment" in h5lib.get_child_group_names(input_h5):
        try:
            # Python2 version of .h5 file
            d = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["dateOfExperiment", "dateOfExperiment"])).tolist()[0]
            t = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["timeOfExperiment", "timeOfExperiment"])).tolist()[0]
        except:
            # Python3 version of .h5 file
            print "date=", np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["dateOfExperiment"])).tolist()
            print "time=", np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["timeOfExperiment"])).tolist()
            d = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["dateOfExperiment"])).tolist()[0]
            t = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["timeOfExperiment"])).tolist()[0]
#       print "d=", d
#       print "t=", t
        dt=datetime.datetime.strptime(d+t, "%Y%m%d%H%M%S")
    else:
        sys.exit("Cannot extract date and time from input data")
       
    return dt.strftime("%a %b %d %Y %H:%M:%S")

# ------------------------------------------------------------------------------

# create 2-photon time series, pointing to specified filename
# use junk values for 2-photon metadata, for now at least
def create_2p_ts(nwb_object, name, fname, stack_t, plane):
    zero = np.zeros(1)
    twop = ts.TwoPhotonSeries(name, nwb_object, "acquisition")
    twop.set_max_voltage(0)
    twop.set_min_voltage(0)
    twop.set_bits_per_pixel(16)
    twop.set_pmt_gain(1.0)
    twop.set_wavelength(100.0)
    twop.set_indicator("GCaMP6s, 780nm")
    twop.set_imaging_depth(0.150+0.1*plane)
    twop.set_scan_line_rate(16000)
    twop.set_field_of_view([ 600e-6, 600e-6 ])
    twop.set_orientation("rostral is image top")
    twop.set_format("external")
    twop.set_dimension([512, 512, 1])
    twop.set_distance(.15)
    twop.set_external_file(fname)
    twop.set_time(stack_t)
    twop.set_data(zero, "None", 1, 1)
    twop.finalize()

# ------------------------------------------------------------------------------

# create raw data timeseries
# use junk values resolution for now
def create_aq_ts(nwb_object, name, modality, timestamps, rate, data, comments = '', descr = '', link = 0):
    twop = ts.TimeSeries(name, nwb_object, modality)
    if link == 0:
        twop.set_time(timestamps)
    else:
        twop.set_time_as_link(timestamps)
    twop.set_data(data, "Zaber motor steps", 0.1, 1)
    twop.set_comments(comments)
    twop.set_description(descr)
    twop.finalize()

# pull out masterImage arrays and create entries for each in
#   /acquisition/images
# masterImages are store in:
#        tsah::descrHash::[2-7]::value::1::[1-3]::masterImage
# each image has 2 color channels, green and red
reference_image_red = {}
reference_image_green = {}
def create_reference_image(orig_h5, nwb_object, plane_map, area, plane, num_plane = 3):
    area_grp = orig_h5["timeSeriesArrayHash/descrHash"]["%d"%(1+area)]
    if num_plane == 1:
        plane_grp = area_grp["value/1"]
    else:
        plane_grp = area_grp["value/1"]["%d"%(plane)]
    master = plane_grp["masterImage"]["masterImage"].value
    green = np.zeros((512, 512))
    red = np.zeros((512, 512))
    for i in range(512):
        for j in range(512):
            green[i][j] = master[i][j][0]
            red[i][j] = master[i][j][1]
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    #image_plane = "area%d_plane%d" % (area, plane)
    oname = "area%d_plane%d" % (area, plane)
    image_plane = plane_map[oname]
    name = image_plane + "_green"
    fmt = "raw"
    desc = "Master image (green channel), in 512x512, 8bit"
    nwb_object.create_reference_image(green, name, fmt, desc, 'uint8')
    reference_image_green[image_plane] = green
    name = "area%d_plane%d_red" % (area, plane)
    name = image_plane + "_red"
    desc = "Master image (red channel), in 512x512, 8bit"
    nwb_object.create_reference_image(red, name, fmt, desc, 'uint8')
    reference_image_red[image_plane] = red

# ------------------------------------------------------------------------------

# pull out all ROI pixel maps for a particular subarea and imaging plane
#   and store these in the segmentation module
def fetch_rois(orig_h5, plane_map, seg_iface, area, plane, num_planes=3):
    tsah = orig_h5["timeSeriesArrayHash"]
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    #image_plane = "area%d_plane%d" % (area, plane)
    oname = "area%d_plane%d" % (area, plane)
    image_plane = plane_map[oname]           
    # first get the list of ROIs for this subarea and plane
    # if num_planes == 1:
    #     ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["ids"]
    # else:
    #     ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["%d"%plane]["ids"]
    ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["%d"%plane]["ids"]
    roi_ids = ids["ids"].value
    lookup = tsah["value"]["%d"%(area+1)]["ids"]["ids"].value
    for i in range(len(roi_ids)):
        rid = roi_ids[i]
        if num_planes == 1:
            rois = tsah["descrHash"]["%d"%(area+1)]["value"]["1"]
        else:
            rois = tsah["descrHash"]["%d"%(area+1)]["value"]["1"]["%d"%plane]
        # make sure the ROI id is correct
        record = rois["rois"]["%s"%(1+i)]
        x = int(parse_h5_obj(record["id"])[0])
        assert x == int(rid)
        pix = parse_h5_obj(record["indicesWithinImage"])[0]
        # pix = record["indicesWithinImage/indicesWithinImage"].value
        pixmap = []
        for j in range(len(pix)):
            v = pix[j]
            px = int(v / 512)
            py = int(v) % 512
            pixmap.append([py,px])
        weight = np.zeros(len(pixmap)) + 1.0
        print "image_plane=", image_plane
        seg_iface.add_roi_mask_pixels(image_plane, "%d"%x, "ROI %d"%x, pixmap, weight, 512, 512)

# ------------------------------------------------------------------------------

# map between ROI ID in imaging plane and dF/F row valueMatrix:
#        tsah::value::[2-7]::ids::ids
# map of ROI IDs per imaging plane:
#        tsah::value::[2-7]::imagingPlane::[1-3]::ids::ids
# time for dF/F samples
#        tsah::value::[2-7]::valueMatrix
def fetch_dff(orig_h5, nwb_object, plane_map, dff_iface, area, plane, num_planes=3):
    area_grp = orig_h5["timeSeriesArrayHash/value"]["%d"%(area+1)]
    # if num_planes == 1:
    #         plane_ids = area_grp["imagingPlane"]["ids/ids"].value
    #     else:
    #         plane_ids = area_grp["imagingPlane"]["%d"%plane]["ids/ids"].value
    plane_ids = area_grp["imagingPlane"]["%d"%plane]["ids/ids"].value
    area_ids = area_grp["ids/ids"].value
    values = area_grp["valueMatrix/valueMatrix"].value
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    #image_plane = "area%d_plane%d" % (area, plane)
    oname = "area%d_plane%d" % (area, plane)
    image_plane = plane_map[oname]
    t = area_grp["time/time"].value * 0.001
    trial_ids = area_grp["trial/trial"].value
    # for each plane ID, find group idx. df/f is idx'd row in values
    for i in range(len(plane_ids)):
        pid = plane_ids[i]
        row = -1
        for j in range(len(area_ids)):
            if area_ids[j] == pid:
                row = j
                break
        assert row >= 0, "Unable to find row for ROI " + pid
        dff = values[j]
        dff_data = []
        dff_t = []
        dff_trials = []
        for j in range(len(dff)):
            if np.isnan(dff[j]):
                continue
            dff_data.append(dff[j])
            dff_t.append(t[j])
            dff_trials.append(trial_ids[j])
        if len(dff_t) == 0:
            continue
        dff_ts = ts.TimeSeries("%d" % pid, nwb_object, "module", "dff")
        dff_ts.set_description("dF/F for ROI %d" % pid)
        dff_ts.set_comments("For segmentation, see similarly named pixel list under 'segmentation'")
        dff_ts.set_source(oname)
        dff_ts.set_data(dff_data, "dF/F", 1, 1)
        dff_ts.set_dtype('f4')
        dff_ts.set_time(dff_t)
        dff_ts.set_value("roi_segments", "ROIs")
        dff_iface.add_trace(image_plane, dff_ts)
        dff_roi_path = "processing/ROIs/DfOverF/" + image_plane + '/' + str(pid)
        dff_roi_grp = nwb_object.file_pointer[dff_roi_path]
        dff_roi_grp.create_dataset("trial_ids", data=dff_trials)

# ------------------------------------------------------------------------------

def create_behavioral_time_series(nwb_object, hash_folder, keyName, series_type, \
                                  series_name, t, var, description, source, \
                                  options):
    # initialize timeseries
    ts = nwb_object.create_timeseries(series_type, series_name)
    if keyName in ["poleInReach", "rewardCue", "leftReward", "rightReward", \
                   "rewardCue",   "leftLicks", "rightLicks", "whiskerVars", \
                   "touches"]:
        # SP data
        if var == '':
            samples = 0
            for i in range(len(t)):
                if np.isnan(t[i]):
                    continue
                samples += 1
            assert samples == len(t), "NaNs found when reading " + series_name
        else:
            clean_var = []
            clean_t   = []
            lent = len(t)
            for i in range(lent):
                if not np.isnan(t[i]) and not np.isnan(var[i]):
                    clean_var.append(var[i])
                    clean_t.append(t[i])
            t   = np.array(clean_t)
            var = np.array(clean_var)
            if options.verbose:
                print series_name + " has %d nans (removed)" % (lent - len(t))

        if keyName in ["poleInReach", "rewardCue", "leftReward", "rightReward", \
                       "rewardCue", "touches"]:
            # times are stored as 'on' in even intervals, 'off' in odd intervals
            # create list of alternating +1, -1, +1, -1
            ts.set_value("timestamps", t)
            on_off = np.int_(np.zeros(len(t)))
            on_off += -1
            on_off[::2] *= -1
            ts.set_data(on_off)
        elif keyName in ["whiskerVar"]:
            if series_name == "whisker_angle":
# GD unit, conversion and resolution should be read in from input file
                ts.set_data(var, unit="degrees", conversion=1, resolution=0.001)
            else:
# GD what is the unit for curvature?
                ts.set_data(var, "Unknown", 1, 1)
        else: 
            # keyName in ["leftLicks", "rightLicks"]
            # there's no meaningful entries to store in data[], as events are binary
            # store a '1'
            data = np.zeros(len(t))
            data += 1
            ts.set_data(data, "Licks", 1, 1)
        ts.set_time(t)
    else:
        # NL data
        # fill in values
        time = np.array(h5lib.get_value_pointer_by_key(hash_folder, keyName, \
                                                    options.debug)).tolist()
        ts = behavior_helper(ts, time, t, 1)
    ts.set_description(description)
    ts.set_source(source)
    ts.set_comments("Timestamp array stores " + series_name + " times")
    return ts

# ------------------------------------------------------------------------------

def behavior_helper(ts, time, start_times, data_value):
    # times relative to session time
    timestamps = time + start_times
    # delete bad trial points
    timestamps = timestamps[~np.isnan(time)]
    # create dummy data: 1 for 'on', -1 for 'off'
    data = [data_value]*len(timestamps)
    # write data into timeseries
    ts.set_value("timestamps", timestamps)
# GD what are the agguments of 'set_data'?
    ts.set_data(data,'s',0.1,1)

    return ts

# ------------------------------------------------------------------------------

# pole position
# store accessible intervals as BehavioralInterval (pole_accessible)
# store touch times as BehavioralInterval (pole_touch
def process_pole_position(orig_h5, nwb_object, options):
    mod = nwb_object.create_module("Pole")
    
    if "eventSeriesArrayHash" in h5lib.get_child_group_names(orig_h5) and \
       "poleInReach" in h5lib.get_key_list(orig_h5['eventSeriesArrayHash']):
        # SP data
        mod.set_description("Intervals that pole is accessible")
        pole_iface = mod.create_interface("BehavioralTimeSeries")
        pole_iface.set_source("Times and intervals are as reported in Simon's data file")

        keyName = "poleInReach"
        hash_folder = orig_h5['eventSeriesArrayHash']
        grp = h5lib.get_value_pointer_by_key(hash_folder, keyName, \
                                              options.debug)
# GD: should use trialTimeUnit to figure out the multiplier below
        t = grp["eventTimes/eventTimes"].value * 0.001
        source = "Intervals are as reported in somatosensory cortex data file"
        description = h5lib.get_description_by_key(hash_folder, keyName)
        pole_acc = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName, "IntervalSeries", \
                   "pole_accessible", t, '', description, source, options)
        pole_iface.add_timeseries(pole_acc)
        mod.finalize()
    elif "PoleInTime"  in h5lib.get_key_list(orig_h5['trialPropertiesHash']) and \
         "PoleOutTime" in h5lib.get_key_list(orig_h5['trialPropertiesHash']):
        # NL data
        mod.set_description("Intervals that pole is accessible")
        pole_iface = mod.create_interface("BehavioralEvents")
        pole_iface.set_source("Times and intervals are as reported in motor cortex data file, but relative to session start")

        hash_folder = orig_h5["trialPropertiesHash"]
        source = "Times as reported in motor cortex data file, but relative to session start"
        trial_start_times = orig_h5["trialStartTimes/trialStartTimes"].value

        # get relevant pole_in data
        keyName1 = "PoleInTime"
        description = h5lib.get_description_by_key(hash_folder, keyName1)
        pole_in  = create_behavioral_time_series(nwb_object, \
                       hash_folder, keyName1, "IntervalSeries", \
                       "pole_in", trial_start_times, '', description, source, \
                       options)
        pole_iface.add_timeseries(pole_in)

        # same procedure for pole_out
        keyName2 = "PoleOutTime"
        description = h5lib.get_description_by_key(hash_folder, keyName2)
        pole_out = create_behavioral_time_series(nwb_object, \
                       hash_folder, keyName2, "IntervalSeries", \
                       "pole_out", trial_start_times, '', description, source, \
                       options)
        pole_iface.add_timeseries(pole_out)

        mod.finalize()
    else:
        sys.exit("Unable to read pole position")
        
# ------------------------------------------------------------------------------
        
# licks
# BehavioralEvent (lick_left)
# BehavioralEvent (lick_right)
def process_licks(orig_h5, nwb_object):
    mod = nwb_object.create_module("Licks")
    mod.set_description("Lickport contacts, right and left")
    lick_iface = mod.create_interface("BehavioralEvents")
    lick_iface.set_source("Data as reported in somatosensory cortex data file")

    hash_folder = orig_h5['eventSeriesArrayHash']
    source = "Intervals are as reported in somatosensory cortex data file"

    # Handle left licks
    keyName1 = 'leftLicks'
    description = h5lib.get_description_by_key(hash_folder, keyName1)
    grp = h5lib.get_value_by_key(hash_folder, keyName1)
# GD
    t = grp["eventTimes/eventTimes"].value * 0.001
    ts_left = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName1, "IntervalSeries", \
                   "lick_left", t, '', description, source, options)
    lick_iface.add_timeseries(ts_left)

    # Handle right licks      
    keyName2 = 'rightLicks'             
    description = h5lib.get_description_by_key(hash_folder, keyName2) 
    grp = h5lib.get_value_by_key(hash_folder, keyName2)  
# GD
    t = grp["eventTimes/eventTimes"].value * 0.001
    ts_right = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName2, "IntervalSeries", \
                   "lick_right", t, '', description, source, options)
    lick_iface.add_timeseries(ts_right)

    mod.finalize()

# ------------------------------------------------------------------------------
        
# water
# BehavioralInterval (water_left)
# BehavioralInterval (water_right)
def process_water(orig_h5, nwb_object):
    mod = nwb_object.create_module("Water")
    mod.set_description("Water reward intervals, left and right")
    water_iface = mod.create_interface("BehavioralTimeSeries")
    water_iface.set_source("Times and intervals are as reported in Simon's data file")

    hash_folder = orig_h5['eventSeriesArrayHash']
    source = "Intervals are as reported in somatosensory cortex data file"

    # left water
    keyName1 = "leftReward"
    description = h5lib.get_description_by_key(hash_folder, keyName1)
    grp = h5lib.get_value_by_key(hash_folder, keyName1)
# GD
    t = grp["eventTimes/eventTimes"].value * 0.001
    water_left = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName1, "IntervalSeries", \
                   "water_left", t, '', description, source, options)
    water_iface.add_timeseries(water_left)

    # right water 
    keyName2 = "rightReward"
    description = h5lib.get_description_by_key(hash_folder, keyName2)
    grp = h5lib.get_value_by_key(hash_folder, keyName2)
# GD
    t = grp["eventTimes/eventTimes"].value * 0.001
    water_right = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName2, "IntervalSeries", \
                   "water_right", t, '', description, source, options)
    water_iface.add_timeseries(water_right)

    mod.finalize()
        
# ------------------------------------------------------------------------------

# auditory_cue
# BehavioralInterval (auditory_cue)
def process_cue(orig_h5, nwb_object, options):
    mod = nwb_object.create_module("Auditory")

    if "eventSeriesArrayHash" in h5lib.get_child_group_names(orig_h5) and \
                  "rewardCue" in h5lib.get_key_list(orig_h5['eventSeriesArrayHash']):
        # SP data
        mod.set_description("Auditory cue signaling animal to collect reward")
        source = "Intervals are as reported in somatosensory cortex data file"
        cue_iface = mod.create_interface("BehavioralTimeSeries")
        cue_iface.set_source("Intervals are as reported in Simon's data file")

        hash_folder = orig_h5['eventSeriesArrayHash']

        # auditory cue 
        keyName = "rewardCue"
        description = h5lib.get_description_by_key(hash_folder, keyName)
        grp = h5lib.get_value_by_key(hash_folder, keyName)
# GD
        t = grp["eventTimes/eventTimes"].value[0:4] * 0.001
        cue = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName, "IntervalSeries", \
                   "auditory_cue", t, '', description, source, options)
        cue_iface.add_timeseries(cue)      
        mod.finalize()
    elif "CueTime" in h5lib.get_key_list(orig_h5['trialPropertiesHash']):
        # NL data
        mod.set_description("Auditory cue signaling animal to report decision")
        # create interface
        cue_iface = mod.create_interface("BehavioralEvents")
        source = "Times are as reported in Nuo's data file, but relative to session time"
        cue_iface.set_source(source)

        keyName = "CueTime"
        hash_folder = orig_h5['trialPropertiesHash']
        description = h5lib.get_description_by_key(hash_folder, keyName)
        trial_start_times = orig_h5["trialStartTimes/trialStartTimes"].value
        # initialize and process auditory cue time series
        cue = create_behavioral_time_series(nwb_object, \
              orig_h5["trialPropertiesHash"], "CueTime", "IntervalSeries", \
              "auditory_cue", trial_start_times, '', description, source, options)
        cue_iface.add_timeseries(cue)
        mod.finalize()
    else:
        sys.exit("Unable to read cue")

# ------------------------------------------------------------------------------

# time is stored in tsah::value::1::time::time
# angle is stored in tsah::value::1::valueMatrix::valueMatrix[0[
# curvature is stored in tsah::value::1::valueMatrix::valueMatrix[1[
def process_whisker(orig_h5, nwb_object, options):
    # create module
    mod = nwb_object.create_module("Whisker")
    mod.set_description("Whisker angle and curvature (relative) of the whiskers and times when the pole was touched by whiskers")

    # Create interface
    whisker_iface = mod.create_interface("BehavioralEvents")
    whisker_iface.set_source("Data as reported in simon's data file")

    # Create time series
    keyName = 'whiskerVars'
    hash_folder = orig_h5['timeSeriesArrayHash']
    grp   = h5lib.get_value_by_key(hash_folder, keyName)   
# GD scaling must be read in from input file
    t     = grp["time/time"].value * 0.001
    val   = grp["valueMatrix/valueMatrix"].value
    
    # whisker angle time series
    # embedded nans screw things up -- remove them
    # count how many non-nan values, and prepare output array
    angle = val[0]
    source   = "Whisker angle as reported in somatosensory cortex data file"
    description = "Angle of whiskers"
    ts_angle = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName, "IntervalSeries", \
                   "whisker_angle", t, angle, description, source, options)
    whisker_iface.add_timeseries(ts_angle)

    # whisker curvature
    curv  = val[1]
    source = "Whisker curvature as reported in somatosensory cortex data file"
    description = "Curvature (relative) of whiskers"
    ts_curve = create_behavioral_time_series(nwb_object, \
                   hash_folder, keyName, "IntervalSeries", \
                   "whisker_curve", t, curv, description, source, options)
    whisker_iface.add_timeseries(ts_curve)

    # pole touches
    pole_iface = mod.create_interface("BehavioralTimeSeries")

    keyName = "touches"
    hash_folder =  orig_h5["eventSeriesArrayHash"]
    grp = h5lib.get_value_by_key(hash_folder, keyName)

    # protraction touches
# GD
    t = grp["eventTimes/1/1"].value * 0.001
    description = "Intervals that whisker touches pole (protract)"
    source = "Intervals are as reported in somatosensory cortex data file"
    pole_touch_pr = create_behavioral_time_series(nwb_object, hash_folder, \
                        keyName, "IntervalSeries", "pole_touch_protract", \
                        t, '', description, source, options)
    pole_iface.add_timeseries(pole_touch_pr)

    # retraction touches
# GD
    t = grp["eventTimes/2/2"].value * 0.001
    description = "Intervals that whisker touches pole (retract)"
    pole_touch_re = create_behavioral_time_series(nwb_object, hash_folder, \
                        keyName, "IntervalSeries", "pole_touch_retract", \
                        t, '', description, source, options)
    pole_iface.add_timeseries(pole_touch_re)

    mod.finalize()

# ------------------------------------------------------------------------------

def process_pole_touches(orig_h5, nwb_object):
    # add kappaMaxAbsOverTouch to pole_touch_protract
    keyName1 = "touches"
    keyName2 = "kappaMaxAbsOverTouch"
    kappa_ma_pr = h5lib.get_value2_by_key2(orig_h5["eventSeriesArrayHash"], keyName1, \
                                                   "eventPropertiesHash/1", keyName2)
    pole_tp_path = "processing/Whisker/BehavioralTimeSeries/pole_touch_protract"
    pole_tp_grp = nwb_object.file_pointer[pole_tp_path]
    pole_tp_grp.create_dataset("kappa_max_abs_over_touch", data=kappa_ma_pr)
    # add kappaMaxAbsOverTouch to pole_touch_retract
    kappa_ma_re = h5lib.get_value2_by_key2(orig_h5["eventSeriesArrayHash"], keyName1, \
                                                   "eventPropertiesHash/2", keyName2)
    pole_tr_path = "processing/Whisker/BehavioralTimeSeries/pole_touch_retract"
    pole_tr_grp = nwb_object.file_pointer[pole_tr_path]
    pole_tr_grp.create_dataset("kappa_max_abs_over_touch", data=kappa_ma_re)

# ------------------------------------------------------------------------------

def process_stimulus(orig_h5, nwb_object):
    keyName = "StimulusPosition"
    stim_pos = h5lib.get_value_by_key(orig_h5['trialPropertiesHash'], keyName)
    trial_t = orig_h5["trialStartTimes/trialStartTimes"].value * 0.001
    rate = (trial_t[-1] - trial_t[0])/(len(trial_t)-1)
    description = h5lib.get_description_by_key(orig_h5["trialPropertiesHash"], keyName)
    zts = nwb_object.create_timeseries("IntervalSeries","zaber_motor_pos")
    zts.set_time(trial_t)
    zts.set_data(stim_pos, "unknown", 1, 1)
    zts.set_description(description)
    zts.set_comments(keyName)
    zts.set_path("processing/ZMP/BehavioralTimeSeries/zaber_motor_pos")
    zts.finalize()

# ------------------------------------------------------------------------------

# trial start times are stored in: h5::trialStartTimes::trialStartTimes
# trial stop isn't stored. assume that it's twice the duration of other
#   trials -- padding on the high side shouldn't matter
def create_trials(orig_h5, nwb_object):
    trial_id = orig_h5["trialIds/trialIds"].value
    trial_t  = orig_h5["trialStartTimes/trialStartTimes"].value * 0.001
    if "GoodTrials" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
        # NL data
        good_trials = h5lib.get_value_by_key(orig_h5["trialPropertiesHash"],  "GoodTrials")
        ignore_ivals_start = [time for (time, good_trial) in zip(trial_t,good_trials) if good_trial == 0]
        # trial stop isn't stored. assume that it's twice the duration of other
        #   trials -- padding on the high side shouldn't matter
        ival = (trial_t[-1] - trial_t[0]) / (len(trial_t) - 1)
        trial_t = np.append(trial_t, trial_t[-1] + 2*ival)
        ignore_ivals_stop = [time for (time, good_trial) in zip(trial_t[1:],good_trials) if good_trial == 0]
        ignore_intervals = [ignore_ivals_start, ignore_ivals_stop]
    else:
        # SP data
        ival = (trial_t[-1] - trial_t[0]) / (len(trial_t) - 1)
        trial_t = np.append(trial_t, trial_t[-1] + 2*ival)

    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        start = trial_t[i]
        stop  = trial_t[i+1]
        epoch = nwb_object.create_epoch(trial, start, stop)
        # pole_pos_path = "trialPropertiesHash/value/3/3"
        #         pole_pos = str(orig_h5[pole_pos_path].value[i])
        #         epoch.description = ("Stimulus position - in Zaber motor steps (approximately, 10,000 = 1 mm): " + pole_pos)
        if "GoodTrials" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
            # NL data
            raw_path = "descrHash/value/%d" % (trial_id[i])
            raw_file = parse_h5_obj(orig_h5[raw_path])[0]
            if len(raw_file) == 1:
                raw_file = 'na'
            else:
                raw_file = str(raw_file)
            epoch.description = ("Raw Voltage trace data files used to acuqire spike times data: "\
                                 + raw_file \
                                 + "\n\ignore intervals: mark start and stop times of bad trials "\
                                 + " when mice are not performing")
            #epoch.set_ignore_intervals(ignore_intervals)
            # collect behavioral data
            events = {"auditory_cue" :  "Auditory", "pole_in" :  "Pole", "pole_out" :  "Pole"}
            for key in events.keys():
                ts_path = "/processing/" + events[key]  + "/BehavioralEvents/" + key
#               print "ts=", ts
                epoch.add_timeseries(key, ts_path)
            epoch.finalize()
#           print "nwb_object.file_pointer['epochs'].keys()=", nwb_object.file_pointer["epochs"].keys()

#           print "trial=", trial
#           print "nwb_object.file_pointer['epochs'].keys()=", nwb_object.file_pointer["epochs"].keys()
#           print "ignore_intervals=", ignore_intervals
#           print "dir(nwb_object.file_pointer['epochs'])=", dir(nwb_object.file_pointer['epochs'])
# GD        nwb_object.file_pointer["epochs"][trial].create_dataset("ignore_intervals", data=ignore_intervals, dtype='f8')
        else:
            # SP data
            epoch.description = ("Data that belong to " + trial)
            events    = {"lick_right" :  "Licks",   "lick_left"     :  "Licks", \
                      "whisker_angle" :  "Whisker", "whisker_curve" :  "Whisker"}
            intervals = {"water_left" :  "Water",   "water_right"   :  "Water", \
                         "pole_accessible"     :  "Pole",    "auditory_cue"       :  "Auditory", \
                         "pole_touch_protract" :  "Whisker", "pole_touch_retract" :  "Whisker"}
            for key in events:
                ts_path = "/processing/" + events[key] + "/BehavioralEvents/" + key
                epoch.add_timeseries(key, ts_path)
            for key in intervals:
                ts_path = "/processing/" + intervals[key] + "/BehavioralTimeSeries/" + key
                epoch.add_timeseries(key, ts_path)
            epoch.finalize()

# ------------------------------------------------------------------------------

# each subarea has list of trials and ROI ids
# to parse, take each subarea, pull out trials
#   foeach trial, write ROI
#   to get plane, find ROI in imaging plane
def create_trial_roi_map(orig_h5, nwb_object, plane_map):
    fp = nwb_object.file_pointer
    num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    for i in range(num_subareas):
        area = i + 1
        grp_name = "timeSeriesArrayHash/value/%d" % (area + 1)
        block = orig_h5[grp_name]
        trials = block["trial/trial"].value
        ids = block["ids/ids"].value
        # create way to map ROI onto plane
        planemap = {}
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(area + 1)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        if num_planes > 1:
            for j in range(num_planes):
                plane = j + 1
                grp_name = "imagingPlane/%d/ids/ids" % plane
                plane_id = block[grp_name].value
                for k in range(len(plane_id)):
                    planemap["%d"%plane_id[k]] = "%d" % plane
            trial_list = {}
            for j in range(len(trials)):
                tid = trials[j]
                name = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
                #name = "Trial_%d" % trials[j]
                if name not in trial_list:
                    trial_list[name] = j
            for trial_name in trial_list.keys():
                roi_list = []
                plane_list = []
                valid_whisker = get_valid_trials(orig_h5, "whisker")
                valid_Ca = get_valid_trials(orig_h5, "Ca")
                for k in range(len(ids)):
                    roi = "%d" % ids[k]
                    plane = int(planemap[roi])
                    # convert from file-specific area/plane mapping to
                    #   inter-session naming convention
                    #imaging_plane = "area%d_plane%d" % (area, plane)
                    oname = "area%d_plane%d" % (area, plane)
                    imaging_plane = plane_map[oname]
                    s = "processing/ROIs/DfOverF/%s/%s" % (imaging_plane, roi)
                    roi_list.append(ids[k])
                    plane_list.append(imaging_plane)
                grp = nwb_object.file_pointer["epochs/" + trial_name]
                grp.create_dataset("ROIs", data=roi_list)
                grp.create_dataset("ROI_planes", data=plane_list)
                grp.attrs["valid_whisker_data"]
        else:
            print "  Warning: in create_trial_roi_map num_planes == 1"

# ------------------------------------------------------------------------------

# add trial types to epoch for indexing
def get_trial_types(orig_h5, nwb_object):
    fp = nwb_object.file_pointer
    trial_id = orig_h5["trialIds/trialIds"].value
    trial_types_all = []
    trial_type_strings = parse_h5_obj(orig_h5['trialTypeStr'])[0]
    if "PhotostimulationType" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
        # NL data
        photostim_types = h5lib.get_value_by_key(orig_h5["trialPropertiesHash"], "PhotostimulationType")
        print "photostim_types=", photostim_types
        num_trial_types = 8
    else:
        # SP data
        valid_whisker = get_valid_trials(orig_h5, "whisker")
        valid_Ca      = get_valid_trials(orig_h5, "Ca")
        num_trial_types = 6

    # collect all trials (strings)
    for i in range(num_trial_types):
        trial_types_all.append(str(trial_type_strings[i]))
    # write specific entries for the given trial
    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial_name = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        trial_types = []
        trial_type_mat = parse_h5_obj(orig_h5['trialTypeMat'])[0]
        for j in range(num_trial_types):
            if trial_type_mat[j,i] == 1:
                trial_types.append(trial_types_all[j])
        grp = nwb_object.file_pointer["epochs/" + trial_name]
        grp.create_dataset("trial_types", data=trial_types)
        if "PhotostimulationType" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
            # NL data
            ps_type_value = photostim_types[i]
            if ps_type_value == 0:
                photostim_type = "non-stimulation trial"
            elif ps_type_value == 1:
                photostim_type = "PT axonal stimulation"
            elif ps_type_value == 2:
                photostim_type = "IT axonal stimulation"
            else:
                photostim_type = "discard"
            grp.create_dataset("photostimulation_type", data=photostim_type)
        else:
            # SP data
            if i in valid_whisker:
                grp.attrs["valid_whisker_data"] = "True"
            else:
                grp.attrs["valid_whisker_data"] = "False"
            if i in valid_Ca:
                 grp.attrs["valid_Ca_data"] = "True"
            else:
                grp.attrs["valid_Ca_data"] = "False"
        
# ------------------------------------------------------------------------------
        
def get_valid_trials(orig_h5, data):
    ts_path = "timeSeriesArrayHash/descrHash/"
    val = []
    num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    if data == "whisker":
        ids = parse_h5_obj(orig_h5[ts_path + '1/value'])[0]
        # ids = list(orig_h5[ts_path + "1/value/value"].value)
        val = val + list(ids)
        val = list(Set(val))
    if data == "Ca":
        for i in range(2,num_subareas+1):
            ids_path = ts_path + "%d/value/2" % i
            ids = parse_h5_obj(orig_h5[ids_path])[0]
            # ids = list(orig_h5[ids_path].value)
            val = val + list(ids)
        val = list(Set(val))
    return val
        
# ------------------------------------------------------------------------------

def process_metadata(nwb_object, input_h5, options):
    if options.verbose:
        print "Processing metadata"

    genotype = ""
    subject  = ""
    weight   = ""
    animalStrain = ""
    animalSource = ""
    animalGeneModification = ""
    value = ""

    if "metaDataHash" in h5lib.get_child_group_names(input_h5):
        # SP data
        key_list = h5lib.get_key_list(input_h5["metaDataHash"])
        group_pointer = input_h5['/metaDataHash']
    else:
        # NL data
        key_list = h5lib.get_child_group_names(input_h5)
        group_pointer = input_h5
#   print "group_pointer.name=", group_pointer.name
#   print "metadata key_list=", key_list
    for key in key_list:
        try:
            if key in ["behavior", "virus", "fiber", "photostim", "extracellular"]:
                # Note: these are two-level groups
                value = ""
                key1_list = h5lib.get_value_pointer_by_path_items(group_pointer, [key]).keys()
#               print "   key1_list=", key1_list
                for key1 in key1_list:
#                   print "      key1=", key1
                    if key1 in ["siteLocations"]:
                        continue
                    value1 = h5lib.get_value_pointer_by_path_items(group_pointer, [key, key1, key1])[:]          
#                   print "      value1=", value1
                    value2 = [str(v) for v in value1]
                    value += "       " + key1 + ": " + ",".join(value2) + "\n " 
            else:
                if "metaDataHash" in h5lib.get_child_group_names(input_h5):
                    value = h5lib.get_value_by_key(group_pointer,key)
                    if h5lib.item_type(value) == "dataset":
                        value = np.array(value).tolist()
                    elif h5lib.item_type(value) == "group":
                        value = value.name
                else:                          
                    value_list = np.array(h5lib.get_value_pointer_by_path_items(group_pointer, [key, key])).tolist()
                    if len(value_list) == 1:
                        value = value_list[0]
                    elif len(value_list) > 0:
                        value = ",".join(value_list)
                    else:
                        value = ""
        except:
            if options.verbose:
                print "\nWarning: could not extract ", key, " value"

#       print "key=", key, " value=", value

        if key in ["animalGeneCopy", "animalGeneticBackground"]:
            genotype += key + ": " + str(value) + "\n"               
   
        elif re.search("animalStrain", key):
            animalStrain += key + ": " + value + " "

        elif re.search("animalSource", key):
            animalSource += key + ": " + value + " "

        elif re.search("animalGeneModification", key):
            animalGeneModification += key + ": " + value + "; "

        elif key == "animalID":
            nwb_object.set_metadata("subject_id", value)
 
        elif key in ["dateOfBirth"]:
            subject  = key + ": " + str(value) + "\n"

        elif re.search("weight", key):
            weight   += key + ": " + str(value) + "\n"

        elif re.search("citation", key):
            nwb_object.set_metadata("related_publications", value)

        elif re.search("experimentType", key):
            nwb_object.set_metadata("notes", value)

        elif re.search("experimenters", key):
            nwb_object.set_metadata("experimenter", value)

        elif key in ["sex", "species", "age"]:
            nwb_object.set_metadata(key, value)

        elif re.search("referenceAtlas", key):
            nwb_object.set_metadata("reference_atlas", value)

        elif re.search("whiskerConfig", key):
            nwb_object.set_metadata("whisker_configuration", value)

        elif key in ["behavior", "virus", "fiber", "photostim", "extracellular", \
                     "surgicalManipulation"]:
            nwb_object.set_metadata(key, value)

    nwb_object.set_metadata("genotype", animalGeneModification + "\n" + genotype + "\n")
    nwb_object.set_metadata("description", subject + " " + animalStrain + "  " + animalSource)
    nwb_object.set_metadata("weight",   weight)

    nwb_object.set_metadata_from_file("surgery",
                             os.path.join(os.environ['NWB_DATA'],"surgery.txt"))
    nwb_object.set_metadata_from_file("data_collection",
                             os.path.join(os.environ['NWB_DATA'], "data_collection.txt"))
    nwb_object.set_metadata_from_file("experiment_description",
                             os.path.join(os.environ['NWB_DATA'],"experiment_description.txt"))

# ------------------------------------------------------------------------------

def process_behavioral_data(orig_h5, nwb_object, options):
    if options.verbose:
        print "Processing behavioral time series data"

    # Both NL and SP data:
    process_pole_position(orig_h5, nwb_object, options)
    process_cue(orig_h5, nwb_object, options)

    # Only SP data:
    if "whiskerVars" in h5lib.get_key_list(orig_h5['timeSeriesArrayHash']):
        process_whisker(orig_h5, nwb_object, options)
        process_licks(orig_h5, nwb_object)
        process_water(orig_h5, nwb_object)
        process_pole_touches(orig_h5, nwb_object)
        process_stimulus(orig_h5, nwb_object)

# ------------------------------------------------------------------------------

# initialize units, so that the field is there, no matter what
def set_default_units(orig_h5, nwb_object):
    fp = nwb_object.file_pointer
    trial_id = orig_h5["trialIds/trialIds"].value
    units = []
    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial_name = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        trial_grp = nwb_object.file_pointer["epochs/" + trial_name]
        trial_grp.create_dataset("units", data=units)

# ------------------------------------------------------------------------------

# collect unit information for a given trial
def get_trial_units(orig_h5, nwb_object, unit_num):
    fp = nwb_object.file_pointer
    for i in range(unit_num):
        i = i+1
        unit = "unit_%d%d" % (int(i/10), i%10)
        grp_name = "eventSeriesHash/value/%d" % i
        grp_top_folder = orig_h5[grp_name]
        trial_ids = grp_top_folder["eventTrials/eventTrials"].value
        trial_ids = Set(trial_ids)
        for trial_num in trial_ids:
            tid = trial_num
            trial_name = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
#           print "\ni=", i, " trial_num=", trial_num, " trial_name=", trial_name
            # check if unit fiels is there and if get values
            try:
                trial_grp = nwb_object.file_pointer["epochs/" + trial_name]
            except:
                trial_grp = nwb_object.file_pointer.create_group("epochs/" + trial_name)

            try:
                units = list(trial_grp["units"].value)
            except KeyError:
                units = []
                trial_grp.create_dataset("units", data=units)
            # update unit information
            units.append(unit)
            del trial_grp["units"]
            trial_grp["units"] = units

# ------------------------------------------------------------------------------

def create_epochs(orig_h5, nwb_object, options):
    if options.verbose:
        print "Creating epochs"

    create_trials(orig_h5,   nwb_object)
    get_trial_types(orig_h5, nwb_object)

    if "GoodTrials" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
        # NL data  only
        set_default_units(orig_h5, nwb_object)
        unit_num = len(orig_h5['eventSeriesHash/value'].keys())
        get_trial_units(orig_h5, nwb_object, unit_num)

# ------------------------------------------------------------------------------

def process_image_data(orig_h5, nwb_object, plane_map):
    # store master images
    print "Creating reference images"
    num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            sys.stdout.write('_')
    print ""
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            create_reference_image(orig_h5, nwb_object, plane_map, subarea+1, plane+1, num_planes)
            sys.stdout.write('.')
            sys.stdout.flush()
    print ""

    mod =  nwb_object.create_module("ROIs")
    mod.set_description("Segmentation (pixel-lists) and dF/F (dffTSA) for all ROIs")
    dff_iface = mod.create_interface("DfOverF")
    dff_iface.set_source("2Photon time series under acquisition are the bases for these ROIs. Those time series are named similarly to the image planes")
    seg_iface = mod.create_interface("ImageSegmentation")
    # pull out image segmentation data. do it by subarea and imaging plane,
    #   as that's how data is stored in the source file
    print "Reading ROI and dF/F"
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            sys.stdout.write('_')
    print ""
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            fetch_rois(orig_h5,             plane_map, seg_iface, subarea+1, plane+1, num_planes)
            fetch_dff( orig_h5, nwb_object, plane_map, dff_iface, subarea+1, plane+1, num_planes)

            sys.stdout.write('.')
            sys.stdout.flush()
    print ""
    mod.finalize()

    # add reference images to image segmentation
    # TODO
    for k in plane_map.keys():
        plane = plane_map[k]
        img = reference_image_red[plane]
        seg_iface.add_reference_image(plane, "%s_red_channel"%plane, img)
        img = reference_image_green[plane]
        seg_iface.add_reference_image(plane, "%s_green_channel"%plane, img)

    print "Creating map between trials and ROIs"
    create_trial_roi_map(orig_h5, nwb_object, plane_map)


    # generate time series for /acquisition
    # data is categorized into 7 chunks -- 6 are 2photon subareas and 1 is
    #   whisker data
    # process 2photon here
    # image stack file links stored in:
    #   timeSeriesArrayHash
    #        value
    #            2-7
    #                imagingPlane
    #                    1-3
    #                        sourceFileList
    print "Creating entries for source .tif"
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            sys.stdout.write('_')
    print ""
    for subarea in range(num_subareas):
        # fetch time array
        grp = orig_h5["timeSeriesArrayHash"]["value"]["%d"%(subarea+2)]
        t = 0.001 * grp["time"]["time"].value
        # now move into imaging plane group, to generate time series for
        #   2photon image stacks
        grp = grp["imagingPlane"]
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            # if num_planes == 1:
            #                 pgrp = grp
            #             else:
            #                 pgrp = grp["%d"%(plane+1)]
            pgrp = grp["%d"%(plane+1)]
            frame_idx = pgrp["sourceFileFrameIdx"]["sourceFileFrameIdx"].value
            lst = parse_h5_obj(pgrp["sourceFileList"])[0]
            cnt = 0
            srcfile = {}
            for k in range(len(lst)):
                srcfile[str(k+1)] = lst[k]
                cnt += 1
            filemap = {}
            lastfile =-1
            lastframe = 1
            stack_t = []
            fname = None
            zero = np.zeros(1)
            assert len(t) == len(frame_idx[0])
            for i in range(len(frame_idx[0])):
                filenum = frame_idx[0][i]
                if lastfile < 0:
                    lastfile = filenum
                framenum = frame_idx[1][i]
                stack_t.append(t[i])
                # check for embedded NaNs
                if np.isnan(filenum):
                    continue
                if np.isnan(framenum):
                    continue
                # use fname as a flag. if it's not None then there's data
                #   to write
                if fname is None:
                    # convert from file-specific area/plane mapping to
                    #   inter-session naming convention
                    oname = "area%d_plane%d" % (subarea+1, plane+1)
                    nname = plane_map[oname]
                    name = "%s_%d" % (nname, filenum)
                    fname = srcfile["%d"%filenum]
                # make sure frames and file numbers are sequential
                assert (lastfile == filenum and framenum == lastframe+1) or framenum == 1
                if lastfile != filenum:
                    if i>0:
                        if not np.isnan(frame_idx[0][i-1] ) and not np.isnan(frame_idx[1][i-1]):
                            create_2p_ts(nwb_object, name, fname, stack_t, plane)
                            stack_t = []
                            fname = None
                lastframe = framenum
                lastfile = filenum
            # make sure we write out the last entry
            if fname is not None:
                create_2p_ts(nwb_object, name, fname, stack_t, plane)
            sys.stdout.write('.')
            sys.stdout.flush()

# ------------------------------------------------------------------------------

def read_probe_locations_matrix(meta_h5, options):
    num_locs = len(h5lib.get_value_pointer_by_path_items(meta_h5, \
                   ["extracellular", "siteLocations"]).keys())
    M = np.zeros([num_locs, 3])
    for i in range(num_locs):
        probe_id = i + 1
        coords = np.array(h5lib.get_value_pointer_by_path_items(meta_h5, \
                   ["extracellular", "siteLocations", \
                    str(probe_id), str(probe_id)])).tolist()
        M[i][:] = coords[:]
    return M

# ------------------------------------------------------------------------------

def detect_shanks(M):
    P = np.zeros([M.shape[0], 1]) # vector of already processed row ids
#   print "\nP=", P
    x1 = M[0,0]
    y1 = M[0,1]
    cols1 = np.where(M[:, 0] == x1)
    cols2 = np.where(M[:, 1] == y1)
    cols = np.intersect1d(cols1, cols2)
#   print "\ncols0=", cols
    shank_size = np.prod(cols.shape)
#   print "shank_size =", shank_size
    curr_shank_id = 0
    num_shanks = 0
    shank_coords = []
    for i in range(M.shape[0]):
        if P[i] == 1:
            continue
        shank_coords.append([M[i, 0], M[i, 1]])
        cols1 = np.where(M[:,0] == M[i, 0])
        cols2 = np.where(M[:,1] == M[i, 1])
        cols = np.intersect1d(cols1, cols2)
        P[cols] = 1
        num_shanks = num_shanks + 1
    return (num_shanks, shank_size, shank_coords)    

# ------------------------------------------------------------------------------

def get_description(meta_h5, options):
    extra_grp = h5lib.get_value_pointer_by_path_items(meta_h5, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["siteLocations"]:
            continue
        key1_list = h5lib.get_value_pointer_by_path_items(extra_grp, [key]).keys()
        for key1 in key1_list:
#           print "      key=", key, " key1=", key1
            value1 = h5lib.get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
            value2 = [str(v) for v in value1]
            value += "       " + key1 + ": " + ",".join(value2) + "\n "
#           print "key=", key, " value=", value
    return value

# ------------------------------------------------------------------------------

def add_shank_group(nwb_object, meta_h5, name, location, \
                    description_data, options):
    group_path = nwbco.EXTRA_CUSTOM(name)
    grp = nwb_object.file_pointer.create_group(group_path)
#   print "description_data=", description_data
    grp.create_dataset("description", data=description_data)

    ################
    # patch-clamp metadata
    # TODO add error checking

# ------------------------------------------------------------------------------

def process_ephys_electrode_map(nwb_object, meta_h5, options):
    M = read_probe_locations_matrix(meta_h5, options)

    num_shanks, shank_size, shank_coords = detect_shanks(M)

    shank = []
    for i in range(1, (num_shanks+1)):
        for j in range(shank_size):
            shank.append("shank" + str(i))

    fp = nwb_object.file_pointer

    # Creating dataset "electrode_map"
    probe = M.tolist()
    map_ephys = fp.create_dataset(nwbco.EXTRA_ELECTRODE_MAP, data=probe) # dataset

    # Creating dataset "electrode_group"
    sz = 0
    for i in range(len(shank)):
        sz = max(sz, len(shank[i]))
    stype = "S%d" % (sz + 1)
    dset = fp.create_dataset(nwbco.EXTRA_ELECTRODE_GROUP, (len(shank),), dtype=stype)
    for i in range(len(shank)):
        dset[i] = shank[i]

    # Creating the shank groups 
    probe_type = h5lib.get_value_pointer_by_path_items(meta_h5, \
                     ["extracellular", "probeType", "probeType"]).value
#   print "\nprobe_type=", probe_type
    description_data = get_description(meta_h5, options)
    rloc = h5lib.get_value_pointer_by_path_items(meta_h5, \
               ["extracellular", "recordingLocation", "recordingLocation"]).value
#   print "\nrloc=", rloc
    for i in range(num_shanks):
        loc = str(rloc[0])
        P = str(shank_coords[i][0])
        Lat = str(shank_coords[i][1])
        add_shank_group(nwb_object, meta_h5, "shank" + str(i+1), \
            "loc: " + loc + ", P: " + P + ", Lat: " + Lat + ", recordingLocation=" + rloc, \
            description_data, options)

# ------------------------------------------------------------------------------

def process_raw_data(orig_h5, nwb_object, options):
    # raw data section
    # lick trace is stored in acquisition
    # photostimulation wave forms is stored in stimulus/processing
    if options.verbose:
        print "Reading raw data"
    # get times
    grp_name = "timeSeriesArrayHash/value/time/time"
    timestamps = orig_h5[grp_name].value
    # calculate sampling rate
    rate = (timestamps[-1] - timestamps[0])/(len(timestamps)-1)
    # get data
    grp_name = "timeSeriesArrayHash/value/valueMatrix/valueMatrix"
    lick_trace = orig_h5[grp_name][:,0]
    aom_input_trace = orig_h5[grp_name][:,1]
    laser_power = orig_h5[grp_name][:,2]
    # get descriptions
    comment1 = parse_h5_obj(orig_h5["timeSeriesArrayHash/keyNames"])[0]
    comment2 = parse_h5_obj(orig_h5["timeSeriesArrayHash/descr"])[0]
    comments = comment1 + ": " + comment2
    grp_name = "timeSeriesArrayHash/value/idStrDetailed"
    descr = parse_h5_obj(orig_h5[grp_name])[0]
    # create timeseries
    #create_aq_ts(nwb_object, "lick_trace", "acquisition",timestamps, rate, lick_trace, comments, descr[0])
    lick_ts = nwb_object.create_timeseries("IntervalSeries","lick_trace")
    lick_ts.set_path("acquisition/timeseries/lick_trace")
    lick_ts.set_data(lick_trace, "unknown", 1, 0)
    lick_ts.set_comments(comments[0])
    lick_ts.set_description(descr[0])
    lick_ts.set_time(timestamps)
    lick_ts.finalize()

    #create_aq_ts(nwb_object, "laser_power", "acquisition", "acquisition/timeseries/lick_trace/",\
    # rate, laser_power, comments, descr[2], 1)
    laser_ts = nwb_object.create_timeseries("IntervalSeries","laser_power")
    laser_ts.set_data(lick_trace, "Watts", 1000.0, 0)
    laser_ts.set_comments(comments[0])
    laser_ts.set_description(descr[2])
    laser_ts.set_time_as_link("acquisition/timeseries/lick_trace/")
    laser_ts.set_time(timestamps)
    laser_ts.set_path("/acquisition/timeseries")
    laser_ts.finalize()

    #create_aq_ts(nwb_object, "aom_input_trace", "stimulus", "acquisition/timeseries/lick_trace/",\
    # rate, aom_input_trace,comments, descr[1], 1)
    aom_ts = nwb_object.create_timeseries("IntervalSeries","aom_input_trace")
    aom_ts.set_path("stimulus/presentation/aom_input_trace")
    aom_ts.set_data(lick_trace, "Volts", 1.0, 0)
    aom_ts.set_comments(comments[0])
    aom_ts.set_description(descr[1])
    aom_ts.set_time_as_link("acquisition/timeseries/lick_trace/")
    aom_ts.set_time(timestamps)
    aom_ts.finalize()

# ------------------------------------------------------------------------------

def process_extracellular_spike_time(orig_h5, nwb_object, options):
    # Create module 'Units' for ephys data
    # Interface 'UnitTimes' contains spike times for the individual units
    # Interface 'EventWaveform' contains waveform data and electrode information
    # Electrode depths and cell types are collected in string arrays at the top level
    # of the module
    if options.verbose:
        print "Reading Event Series Data"
    # create module unit
    mod = nwb_object.create_module("Units")
    mod.set_description("Spike times and waveforms")
    # create interfaces
    spk_waves_iface = mod.create_interface("EventWaveform")
    spk_waves_iface.set_source("Waveform data as reported in Nuo's data file")
    spk_times_iface = mod.create_interface("UnitTimes")
    spk_times_iface.set_source("UnitTimes data as reported in Nuo's data file")
    # top level folder
    grp_name = "eventSeriesHash/value"
    # determine number of units
    unit_num = len(orig_h5[grp_name].keys())
    # initialize cell_types and electrode_depth arrays with default values
    cell_types = ['unclassified']*unit_num
    electrode_depths = [float('NaN')]*unit_num
    unit_descr = parse_h5_obj(orig_h5['eventSeriesHash/descr'])[0]
    # process units
    for i in range(unit_num):
        i = i+1
        unit = "unit_%d%d" % (int(i/10), i%10)
        # initialize timeseries
        spk = nwb_object.create_timeseries("SpikeEventSeries", unit)
        # get data
        grp_name = "eventSeriesHash/value/%d" % i
        grp_top_folder = orig_h5[grp_name]
        timestamps = grp_top_folder["eventTimes/eventTimes"]
        trial_ids = grp_top_folder["eventTrials/eventTrials"]
        # calculate sampling rate
        rate = (timestamps[-1] - timestamps[0])/(len(timestamps)-1)
        waveforms = grp_top_folder["waveforms/waveforms"]
        sample_length = waveforms.shape[1]
        channel = grp_top_folder["channel/channel"].value
        # read in cell types and update cell_type array
        cell_type = parse_h5_obj(grp_top_folder["cellType"])[0]
        if  'numpy' in str(type(cell_type)):
            cells_conc = ' and '.join(map(str, cell_type))
            cell_type = cells_conc
        else:
            cell_type = str(cell_type)
        # try:
        #         cell_type = grp_top_folder["cellType/cellType"][0,0]
        #     except KeyError:
        #         cell_type_1 = grp_top_folder["cellType/1/1"][0,0]
        #         cell_type_2 = grp_top_folder["cellType/2/2"][0,0]
        #         cell_type = cell_type_1 + " and " + cell_type_2
        cell_types[i-1] = cell_type
        try:
            # read in electrode depths and update electrode_depths array
            depth = parse_h5_obj(grp_top_folder["depth"])[0]
            # depth = grp_top_folder["depth/depth"][0,0]
            electrode_depths[i-1] = depth
        except:
            if options.verbose:
                print "Could not extract electrode_depths"
        # fill in values for the timeseries
        spk.set_value("source", "Spike data as reported in Nuo's data file")
        spk.set_path("/processing/Units/EventWaveform")
#       print "\n\nspk.spec=", spk.spec["_attributes"]["source"]
        spk.set_value("sample_length", sample_length)
#       print "\n\nspk.value[sample_length]=", spk.value["sample_length"]
        spk.set_time(timestamps)
        spk.set_data(waveforms, "Volts", 0.1,1)
        spk.set_value("electrode_idx", [channel[0]])
#       print "\n\nelectrode_idx=", [channel[0]]
        spk.set_description("single unit %d with cell-type information and approximate depth, waveforms sampled at 19531.25Hz" %i)
        spk_waves_iface.add_timeseries(spk)
        # add spk to interface
        description = unit_descr[i-1]
        spk_times_iface.add_unit(unit, timestamps, description, "Data from Nuo's file")
        spk_times_iface.append_unit_data(unit, trial_ids, "trial_ids")
    mod.finalize()
    grp = nwb_object.file_pointer["processing/Units"]
    grp.create_dataset("CellTypes", data=cell_types)
    grp.create_dataset("ElectrodeDepths", data=electrode_depths)

# ------------------------------------------------------------------------------

def collect_analysis_information(orig_h5, nwb_object, options):
    if options.verbose:
        print "Collecting analysis information"

    trial_start_times = orig_h5["trialStartTimes/trialStartTimes"].value
    trial_types_all = []
    trial_type_strings = parse_h5_obj(orig_h5['trialTypeStr'])[0]
    # collect all trials (strings)
    for i in range(8):
        trial_types_all.append(str(trial_type_strings[i]))
    # trial_type_strings = orig_h5["trialTypeStr"]
    # trial_types_all = []
    # for i in range(8):
    #     trial_types_all.append(trial_type_strings["%d/%d" %(i+1,i+1)].value[0,0])
    trial_type_mat = orig_h5['trialTypeMat/trialTypeMat'].value
    good_trials = orig_h5['trialPropertiesHash/value/4/4'].value
    grp = nwb_object.file_pointer["analysis"]
    grp.create_dataset("trial_start_times", data = trial_start_times)
    grp.create_dataset("trial_type_string", data = trial_types_all)
    grp.create_dataset("trial_type_mat", data = trial_type_mat)
    grp.create_dataset("good_trials", data = good_trials)

# ------------------------------------------------------------------------------

#ep.ts_self_check()
#sys.exit(0)
def produce_nwb(data_path, metadata_path, output_nwb, options):
    orig_h5 = h5py.File(data_path, "r")
    meta_h5 = ""
    if len(metadata_path) > 0:
        meta_h5 = h5py.File(metadata_path, "r")

    print "Input:", data_path, " ", metadata_path
    print "Output:", output_nwb, "\n"
    if options.replace and os.path.exists(output_nwb):
        os.remove(output_nwb)

    # each of simon's hdf5 files have imaging planes and subareas
    #   labels consistent within the file, but inconsistent between
    #   files. create a map between the h5 plane name and the 
    #   identifier used between files

    vargs = {}
    if len(metadata_path) > 0:                        
        # NL data
        vargs["start_time"] = find_exp_time(meta_h5)
    else:                                             
        # SP data
        vargs["start_time"] = find_exp_time(orig_h5)

    vargs["filename"]       = output_nwb
    vargs["identifier"]     = nwb.create_identifier("Neurodata testing (Nwb_Object)")
    vargs["description"]    = os.path.join(os.environ['NWB_DATA'],"experiment_description.txt")

#   print "vargs=", vargs
    nwb_object = nwb.NWB(**vargs)

    # Process metadat
    if "metaDataHash" in h5lib.get_child_group_names(orig_h5):
        # SP data
        process_metadata(nwb_object, orig_h5, options)
    elif not len(meta_h5) == 0:
        # NL data
        process_metadata(nwb_object, meta_h5, options)
    else:
        sys.exit("\nCannot process metadata. Check your input.")

    # Process behavioral time series
    process_behavioral_data(orig_h5, nwb_object, options)

    # Create epochs
    create_epochs(orig_h5, nwb_object, options)

    if len(metadata_path) == 0:
        # SP data
        plane_map = {}
        create_plane_map(orig_h5, plane_map)
        process_image_data(orig_h5, nwb_object, plane_map)
    else:
        # NL data
        process_ephys_electrode_map(nwb_object, meta_h5, options)
        process_raw_data(orig_h5, nwb_object, options)
        process_extracellular_spike_time(orig_h5, nwb_object, options)
        collect_analysis_information(orig_h5, nwb_object, options)
            
    if options.verbose:
        print "Closing file"

    nwb_object.close()
    print "Done"

# ------------------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog data_h5 [meta_data_h5] [options (-h to list)]"

    parser = optparse.OptionParser(usage=usage)
    parser = make_nwb_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.verbose:
        print "len(args)=", len(args)
    if len(args) in [1, 2]:
        data_path     = args[0]
        if not re.search(".h5", data_path):
            sys.exit("\nScript make_mwb.py accepts as input .h5 data file")
        metadata_path = ""
        if len(args) == 2:
            metadata_path = args[1]
            if not re.search(".h5", metadata_path):
                sys.exit("\nScript make_mwb.py accepts as input .h5 metadata file")
        
        data_basename = os.path.basename(data_path)
        if len(options.output_folder) == 0:
            options.output_folder = os.path.dirname(data_path)
        output_path = os.path.join(options.output_folder, data_basename.split(".")[0] + ".nwb")

        produce_nwb(data_path, metadata_path, output_path, options)
    else:
        parser.print_usage()
        sys.exit(2)


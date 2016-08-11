#!/usr/local/python-2.7.8/bin/python
#
# make_nwb.py
#
# Created by Claudia Friedsam on 2015-02-08.
# Redesigned by Gennady Denisov on 2016-03-28.

import sys
import os

import nwb
from nwb import nwb_file
from nwb import nwb_utils as ut

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
    parser.add_option("-D", "--debug",   action="store_true",  dest="debug", help="output debugging info", default=False)
    parser.add_option("-e", "--no_error_handling", action="store_false", dest="handle_errors", help="handle_errors", default=True)
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
def add_plane_map_entry(plane_map, h5_plane_name, filename, options):
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

def create_plane_map(orig_h5, plane_map, options):
    if options.handle_errors:
        num_subareas = 1
        if '1' in orig_h5['timeSeriesArrayHash/descrHash'].keys():
            num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    else:
        num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1

    for subarea in range(num_subareas):
        # fetch time array
#       disp(['path10=' 'timeSeriesArrayHash/value/%d/imagingPlane' %(subarea + 2)]);
#       disp(['path20=' 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)]);
        if options.handle_errors:
            try:
                grp = orig_h5['timeSeriesArrayHash/value/%d/imagingPlane' %(subarea + 2)]
                grp2 = orig_h5['timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)]
            except:
                try:
                    grp = orig_h5['timeSeriesArrayHash/value/imagingPlane']
                    grp2 = orig_h5['timeSeriesArrayHash/descrHash/value']
                except:
                    print "Cannot create plane map"
                    break
        else:
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
#           print "\nsubarea=", subarea, " plane=", plane
            if options.handle_errors:
                try:
                    pgrp = grp["%d"%(plane+1)]
                except:
                    # Warning: only one imaging plane available (instead of 3)
                    pgrp = grp
            else:
                pgrp = grp["%d"%(plane+1)]

            old_name = "area%d_plane%d" % (subarea+1, plane+1)
            frame_idx = pgrp["sourceFileFrameIdx"]["sourceFileFrameIdx"].value
            lst = parse_h5_obj(pgrp["sourceFileList"])[0]
            for k in lst:
                # srcfile = str(lst[k][k].value)
                srcfile = str(k)
                add_plane_map_entry(plane_map, old_name, srcfile, options)
                break

# ------------------------------------------------------------------------------

# fetch start time/date of experiment
def find_exp_time(input_h5):
    # SP data
    child_groups = h5lib.get_child_group_names(input_h5)
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
        print "Nuo"
        try:
            # Python2 version of .h5 file
            d = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["dateOfExperiment", "dateOfExperiment"])).tolist()[0]
            t = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["timeOfExperiment", "timeOfExperiment"])).tolist()[0]
        except:
            # Python3 version of .h5 file
#           print "date=", np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
#                        ["dateOfExperiment"])).tolist()
#           print "time=", np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
#                        ["timeOfExperiment"])).tolist()
            d = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["dateOfExperiment"])).tolist()[0]
            t = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["timeOfExperiment"])).tolist()[0]
        print "d=", d
        print "t=", t
        dt=datetime.datetime.strptime(d+t, "%Y%m%d%H%M%S")
    # JY data
    elif "date" in h5lib.get_child_group_names(input_h5):
        d = np.array(h5lib.get_value_pointer_by_path_items(input_h5, \
                         ["date", "date"])).tolist()[0]
#       print "date=", "20"+d
        dt = datetime.datetime.strptime("20"+d, "%Y%m%d")
    else:
        sys.exit("Cannot extract date and time from input data")
       
    return dt.strftime("%a %b %d %Y %H:%M:%S")

# ------------------------------------------------------------------------------

# create 2-photon time series, pointing to specified filename
# use junk values for 2-photon metadata, for now at least
def create_2p_tsa(nwb_object, img_plane, ext_file, starting_frame, \
                  timestamps, name):
#def create_2p_tsa(simon, plane_name, external_file, starting_frame, timestamps):
    twop = nwb_object.make_group("<TwoPhotonSeries>", img_plane, \
               path='/acquisition/timeseries', attrs={ \
               "source": "Device 'two-photon microscope'",\
                "description": \
               "2P image stack, one of many sequential stacks for this field of view"})
    twop.set_dataset("format", "external")
    # need to convert name to utf8, otherwise error generated:
    # TypeError: No conversion path for dtype: dtype('<U91').  Added by jt.
    # fname_utf =  fname.encode('utf8')
    twop.set_dataset("external_file", ext_file, attrs={"starting_frame": starting_frame})
    twop.set_dataset("dimension", [512, 512])
    #twop.set_value("pmt_gain", 1.0)
    twop.set_dataset("scan_line_rate", 16000.0)
    twop.set_custom_dataset("channel_name", name+"_red and "+name+"_green")
    twop.set_dataset("field_of_view", [ 600e-6, 600e-6 ])
    twop.set_dataset("imaging_plane", img_plane)
    twop.set_dataset("timestamps", timestamps)

# ------------------------------------------------------------------------------

# save frames for 2-photon series.  These are used in routine create_2p_tsa
def save_2p_frames(external_file, starting_frame, timestamps, fname, stack_t):
    starting_frame.append(len(timestamps))
    timestamps.extend(stack_t)
    external_file.append(fname.encode('utf8'))

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

# pull out masterImage arrays and create entries for each in
#   /acquisition/images
# masterImages are store in:
#        tsah::descrHash::[2-7]::value::1::[1-3]::masterImage
# each image has 2 color channels, green and red
reference_image_red = {}
reference_image_green = {}
def create_reference_image(orig_h5, nwb_object, master_shape, plane_map, area, \
                           plane, options, num_plane = 3):
    area_grp = orig_h5["timeSeriesArrayHash/descrHash"]["%d"%(1+area)]
    if num_plane == 1:
        plane_grp = area_grp["value/1"]
    else:
        plane_grp = area_grp["value/1"]["%d"%(plane)]
    master = plane_grp["masterImage"]["masterImage"].value
    master1_shape = np.array(master).shape
    green = np.zeros([master1_shape[0],master1_shape[1]])
    red   = np.zeros([master1_shape[0],master1_shape[1]])
#   print "\nmaster1_shape=", master1_shape
    for i in range(master1_shape[0]):
        for j in range(master1_shape[1]):
            if len(master1_shape) == 3 or not options.handle_errors:
                green[i][j] = master[i][j][0]
                red[i][j]   = master[i][j][1]
            else:
                # Warning: only one master image is available, so the green and red signals will be identical
                green[i][j] = master[i][j]
                red[i][j]   = master[i][j]
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    oname = "area%d_plane%d" % (area, plane)
    master_shape[oname] = master1_shape
    image_plane = plane_map[oname]
    name = image_plane + "_green"
    fmt = "raw"
    desc = "Master image (green channel), in 512x512, 8bit"
    nwb_object.set_dataset("<image_X>", green, name=name, dtype='uint8', attrs={
                           "description":desc, "format": fmt})
    reference_image_green[image_plane] = green
    name = "area%d_plane%d_red" % (area, plane)
    name = image_plane + "_red"
    desc = "Master image (red channel), in 512x512, 8bit"
    nwb_object.set_dataset("<image_X>", red,   name=name, dtype='uint8', attrs={
                           "description":desc, "format": fmt})
    reference_image_red[image_plane] = red
    return master_shape

# ------------------------------------------------------------------------------

# pull out all ROI pixel maps for a particular subarea and imaging plane
#   and store these in the segmentation module
def fetch_rois(orig_h5, master_shape, plane_map, seg_iface, area, plane, \
               options, num_planes=3):
    tsah = orig_h5["timeSeriesArrayHash"]
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    #image_plane = "area%d_plane%d" % (area, plane)
    oname = "area%d_plane%d" % (area, plane)
    master1_shape = master_shape[oname]
    image_plane = plane_map[oname]           

    # first get the list of ROIs for this subarea and plane
    if options.handle_errors:
        try:
            ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["%d"%plane]["ids"]
        except:
            # Warning: only one imaging plane is available (instead of 3)
            ids = tsah["value"]["%d"%(area+1)]["imagingPlane"]["ids"]
    else:
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
        try:
            record = rois["rois"]["%s"%(1+i)]
            x = int(parse_h5_obj(record["id"])[0])
            assert x == int(rid)
        except:
            print "Missing ROI for area=", area, " plane=", plane, " id=", i
            continue
        pix = parse_h5_obj(record["indicesWithinImage"])[0]
        # pix = record["indicesWithinImage/indicesWithinImage"].value
        pixmap = []
        for j in range(len(pix)):
            v = pix[j]
            px = int(v  / master1_shape[1])
            py = int(v) % master1_shape[0]
            pixmap.append([py,px])
        weight = np.zeros(len(pixmap)) + 1.0
#       print "image_plane=", image_plane, " oname=", oname
        ut.add_roi_mask_pixels(seg_iface, image_plane, "%d"%x, "ROI %d"%x, pixmap, \
            weight, master1_shape[1], master1_shape[0])

# ------------------------------------------------------------------------------

def fetch_dff(orig_h5, dff_iface, seg_iface, plane_map, area, plane, \
              options, num_planes=3):
    area_grp = orig_h5["timeSeriesArrayHash/value"]["%d"%(area+1)]
    if options.handle_errors:
        try:
            plane_ids = area_grp["imagingPlane"]["%d"%plane]["ids/ids"].value
        except:
            # Warning: only one imaging plane is available (instead of 3)
            plane_ids = area_grp["imagingPlane"]["ids/ids"].value
    else:
        plane_ids = area_grp["imagingPlane"]["%d"%plane]["ids/ids"].value

    area_ids = area_grp["ids/ids"].value
    # store dff in matrix and supply that to time series
    dff_data = area_grp["valueMatrix/valueMatrix"].value
    # convert from file-specific area/plane mapping to
    #   inter-session naming convention
    oname = "area%d_plane%d" % (area, plane)
    image_plane = plane_map[oname]
    t = area_grp["time/time"].value * 0.001
    # create array of ROI names for each matrix row
    roi_names = []
    trial_ids = area_grp["trial/trial"].value
    # for each plane ID, find group idx. df/f is idx'd row in values
    for i in range(len(plane_ids)):
        roi_names.append("ROI%d"%plane_ids[i])
    dff_ts = dff_iface.make_group("<RoiResponseSeries>", image_plane)
    dff_ts.set_dataset("data", dff_data, attrs={"unit":"dF/F",
        "conversion": 1.0, "resolution":0.0})
    dff_ts.set_dataset("timestamps", t)
    dff_ts.set_dataset("roi_names", roi_names)
    #- dff_ts.set_value_as_link("segmentation_interface", seg_iface)
    dff_ts.make_group("segmentation_interface", seg_iface)
    trial_ids = area_grp["trial/trial"].value
    dff_ts.set_custom_dataset("trial_ids", trial_ids)

# ------------------------------------------------------------------------------

# Function: create_time_series
# Returns: pointer to the created group
# Arguments: series_name, series_type, series_path, nwb_object,
#            var, t, description, source, comments,
#            hash_group , keyName, options
# Action:
# - creates a group <series_name>
#   of type <series_type>
#   at the path <spath>
#   relative to <nwb_group>
# - creates datasets "data", "num_samples" and "timestamps"
#   inside of the series group
# - populates the datasets from the arrays var and t             
# - populates dataset "timestamps" with the data from input array t
#   using <keyName>
# - sets the attributes "description", "source" and "comments"
#   to the time series group using the input strings
#   <description>, <source> and <comments>
# - perform additional data handling using <h5_hash_group> and <keyName>

# NOTE: argument t has different meaning depending on the keyName
#       (it may be either start time or actual time)

def create_time_series(series_name, series_type, series_path, \
                       orig_h5, nwb_group, group_attrs, var, t, data_attrs, \
                       hash_group_pointer, keyName, options):
    # initialize timeseries
    if options.verbose:
        print "    Creating group ", series_name, " type=", series_type, " path=", series_path
    if len(series_path) > 0:
        ts = nwb_group.make_group(series_type, series_name, path = series_path, \
                                  attrs = group_attrs)
    else:
        ts = nwb_group.make_group(series_type, series_name, attrs = group_attrs)
#   if keyName in ["poleInReach", "rewardCue", "leftReward", "rightReward", \
#                  "leftLicks", "rightLicks", "whiskerVars", "touches"]:

    # Get rid of NaNs   
#   if len(var) == 0:
#       t_not_nan = t[np.isfinite(np.array(t))]
#       if options.verbose:
#           print series_name + " has %d nans (removed)" % (len(t) - len(t_not_nan))
#       t = t_not_nan
#   else:
#       print "np.array(t).shape=", np.array(t).shape, " np.array(var).shape=", np.array(var).shape
#       t_not_nan = t[np.isfinite(np.array(t)) & np.isfinite(np.array(var))]
#       v_not_nan =     var[np.isfinite(np.array(t)) & np.isfinite(np.array(var))]
#       if options.verbose:
#           print series_name + " has %d nans (removed)" % (len(t) - len(t_not_nan))
#       t   = t_not_nan
#       var = v_not_nan
#   if options.verbose:
#           print series_name + " has %d nans (removed)" % (len(t) - len(t_not_nan))

    # 1)  Data = on_off 
    if series_name in ["pole_accessible",                \
                       "water_left",      "water_right", \
                       "pole_touch_protract","pole_touch_retract"] \
       or (series_name == "auditory_cue" and keyName == "rewardCue"):

        timestamps = t
        on_off = np.int_(np.zeros(len(t)))
        on_off += -1
        on_off[::2] *= -1
        data = on_off

        # SP data only
        if series_name == "pole_touch_protract":
            kappa_ma_path = "eventSeriesArrayHash/value/2/eventPropertiesHash/1/value/2/2"
            kappa_ma = orig_h5[kappa_ma_path].value
            ts.set_custom_dataset("kappa_max_abs_over_touch", kappa_ma)
        elif series_name == "pole_touch_retract":
            kappa_ma_path = "eventSeriesArrayHash/value/2/eventPropertiesHash/2/value/2/2"
            kappa_ma = orig_h5[kappa_ma_path].value
            ts.set_custom_dataset("kappa_max_abs_over_touch", kappa_ma)

    # 2) Data = [1]*len(timestamps)
    elif series_name in ["lick_left", "lick_right",\
                         "pole_in",   "pole_out",  \
                         "lick_time"]              \
       or (series_name == "auditory_cue" and keyName == "CueTime"):

        data = [1] * len(t)
        timestamps = t

    # 3) Data = value
    elif series_name in ["whisker_angle", "whisker_curve", \
                         "lick_trace", "aom_input_trace",\
                         "simple_optogentic_stimuli"] \
       or keyName in ["whiskerVars", "Ephys"]:

        data = var
        timestamps = t
        
        if series_name == "simple_optogentic_stimuli":
            ts.set_dataset("site", "site 1")

    else:
       sys.exit("Unknown key "     + keyName     + \
                " or series_name " + series_name + \
                " in create_time_series")

    data_attrs['keyName'] = keyName
    ts.set_dataset("data", data, attrs=data_attrs)
    ts.set_dataset("timestamps", timestamps)
    ts.set_dataset("num_samples", len(timestamps))
    return ts

# ------------------------------------------------------------------------------

# pole position
# store accessible intervals as BehavioralInterval (pole_accessible)
# store touch times as BehavioralInterval (pole_touch
def process_pole_position(orig_h5, nwb_object, options):
#   mod = nwb_object["processing"].make_group("<Module>", "Pole")
    mod = nwb_object.make_group("<Module>", "Pole", path="/processing", \
                                attrs={"description":"Intervals that pole is accessible"})
    pole_iface = mod.make_group("BehavioralEpochs", attrs={
            "source": "Pole intervals"})

    if "eventSeriesArrayHash" in h5lib.get_child_group_names(orig_h5) and \
       "poleInReach" in h5lib.get_key_list(orig_h5['eventSeriesArrayHash']):
        # SP data
        keyName = "poleInReach"
        hash_group_pointer = orig_h5['eventSeriesArrayHash']
        grp = h5lib.get_value_pointer_by_key(hash_group_pointer , keyName, \
                                              options.debug)
# GD: should use trialTimeUnit to figure out the multiplier below
        t = grp["eventTimes/eventTimes"].value * 0.001
#       source = "Intervals are as reported in somatosensory cortex data file"
#       series_path="/stimulus/presentation"
        series_path = "/stimulus/presentation"
        description = h5lib.get_description_by_key(hash_group_pointer , keyName)
        group_attrs = {"source" : "Intervals are as reported in Simon's data file"}
        data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0, \
                       "description" : description}
        pole_acc = create_time_series("pole_accessible", "<IntervalSeries>", series_path,\
                       orig_h5, nwb_object, group_attrs, '', t, data_attrs, \
                       hash_group_pointer, keyName, options)
    elif "PoleInTime"  in h5lib.get_key_list(orig_h5['trialPropertiesHash']) and \
         "PoleOutTime" in h5lib.get_key_list(orig_h5['trialPropertiesHash']):
        # NL or JY data
#       mod.set_description("Intervals that pole is accessible")

        hash_group_pointer = orig_h5["trialPropertiesHash"]
        source = "Times as reported in motor cortex data file, but relative to session start"
        trial_start_times = orig_h5["trialStartTimes/trialStartTimes"].value
        grp = orig_h5["trialPropertiesHash/value/"]
        time = grp["1/1"].value
        series_path = "/stimulus/presentation"

        # get relevant pole_in data
        keyName1 = "PoleInTime"
        time = grp["1/1"].value
        t = time + trial_start_times
        description = h5lib.get_description_by_key(hash_group_pointer, keyName1)
        group_attrs = {"source": source, "description" : description}
        data_attrs  = {"resolution":0.1,"conversion":1.0} 
        pole_in  = create_time_series("pole_in", "<IntervalSeries>", series_path,\
                       orig_h5, nwb_object, group_attrs, '', t, data_attrs, \
                       hash_group_pointer, keyName1, options)

        # same procedure for pole_out
        keyName2 = "PoleOutTime"
        time = grp["2/2"].value
        t = time + trial_start_times
        description = h5lib.get_description_by_key(hash_group_pointer, keyName2)
        group_attrs = {"source": source, "description" : description}
        data_attrs  = {"resolution":float('nan'),"unit":"unknown","conversion":1.0}    
        pole_out = create_time_series("pole_out", "<IntervalSeries>", series_path,\
                       orig_h5, nwb_object, group_attrs, '', t, data_attrs, \
                       hash_group_pointer, keyName2, options)
    else:
        sys.exit("Unable to read pole position")
        
# ------------------------------------------------------------------------------
        
# licks
# BehavioralEvent (lick_left)
# BehavioralEvent (lick_right)
def process_licks(orig_h5, nwb_object):
    mod = nwb_object.make_group("<Module>", "Licks", attrs = \
                                {"description" : "Lick port contact times"}) 

    if not 'LickTime'  in h5lib.get_key_list(orig_h5["trialPropertiesHash"]) and\
       not "EphusVars" in h5lib.get_key_list(orig_h5["timeSeriesArrayHash"]):
        # SP data
        lick_iface = mod.make_group("BehavioralEvents", attrs={
                         "source": "Lick Times as reported in Simon's data file"})
        mod.set_custom_dataset("description", "Lickport contacts, right and left")
        lick_iface.set_attr("source", "Data as reported in somatosensory cortex data file")
        hash_group_pointer  = orig_h5['eventSeriesArrayHash']
        source = "Intervals are as reported in somatosensory cortex data file"

        group_attrs = {"description": "Left lickport contact times (beam breaks left)",
                       "source": "Times as reported in Simon's data file",
                       "comments": "Timestamp array stores lick times"}

        # Handle left licks
        keyName1 = 'leftLicks'
        description = h5lib.get_description_by_key(hash_group_pointer, keyName1)
        grp = h5lib.get_value_by_key(hash_group_pointer, keyName1)
        t = grp["eventTimes/eventTimes"].value * 0.001
        data_attrs = {"description": description, \
                      "source": "Times as reported in Simon's data file",
                      "unit":"Licks", "conversion": 1.0, "resolution": 1.0}
        ts_left = create_time_series("lick_left", "<TimeSeries>", "",\
                       orig_h5, lick_iface, group_attrs, '', t, data_attrs, \
                       hash_group_pointer, keyName1, options)

        # Handle right licks      
        keyName2 = 'rightLicks'             
        description = h5lib.get_description_by_key(hash_group_pointer, keyName2) 
        grp = h5lib.get_value_by_key(hash_group_pointer, keyName2)  
        t = grp["eventTimes/eventTimes"].value * 0.001
        data_attrs = {"description": description, \
                      "source": "Times as reported in Simon's data file",
                      "unit":"Licks", "conversion": 1.0, "resolution": 1.0}
        ts_right = create_time_series("lick_right", "<TimeSeries>", "",\
                       orig_h5, lick_iface, group_attrs, '', t, data_attrs, \
                       hash_group_pointer, keyName2, options)
    elif "EphusVars" in h5lib.get_key_list(orig_h5["timeSeriesArrayHash"]):
        # NL data
        keyName = "EphusVars"
        hash_group_pointer = orig_h5["timeSeriesArrayHash"]
        lick_trace = hash_group_pointer["value/valueMatrix/valueMatrix"][:,0]
        grp_name   = "timeSeriesArrayHash/value/time/time"
        timestamps = orig_h5[grp_name].value
        series_path = "/acquisition/timeseries"
        timestamps = hash_group_pointer["value/time/time"].value
        description = parse_h5_obj(hash_group_pointer["value/idStrDetailed/idStrDetailed"])[0][0]
        comment1 = keyName
        comment2 = h5lib.get_description_by_key(hash_group_pointer, keyName)     
        comments = comment1 + ": " + comment2
        print "comment1=", comment1, " comment2=", comment2, " comments=", comments
        data_attrs={"conversion":1.0, "unit":"unknown", "resolution":float('nan')}  
        group_attrs={"description" : description, "comments" : comments, \
                    "source": "Times as reported in Nuo's data file"}
        lick_ts = create_time_series("lick_trace", "<TimeSeries>", series_path,\
                      orig_h5, nwb_object, group_attrs, lick_trace, timestamps, data_attrs, \
                      hash_group_pointer, keyName, options)
    else:    
        # JY data
        lick_iface = mod.make_group("BehavioralEvents", attrs={
                         "source": "Lick Times as reported in Jianing's data file"})
        hash_group_pointer = orig_h5["trialPropertiesHash"]

        # Handle lick times
        keyName = "LickTime"
        description = h5lib.get_description_by_key(hash_group_pointer , keyName)
        data_attrs = {"resolution":0.1,"conversion":1.0,\
                      "description" : description, \
                      "source" : "Lick times as reported in Jianing's data file"}
        trial_start_times = orig_h5["trialStartTimes/trialStartTimes"].value
        grp = h5lib.get_value_by_key(hash_group_pointer, keyName)
        start_t = []
        num_trials = len(grp.keys())
        for k in sorted([int(k) for k in grp.keys()]):  
            dt = grp[str(k) + "/" + str(k)].value * 0.001
            for i in range(0, len(dt)):
                start_t.append(trial_start_times[int(k-1)])
        ts_time = create_time_series("lick_time", "<TimeSeries>", "",\
                       orig_h5, lick_iface, {}, '', start_t, data_attrs, \
                       hash_group_pointer, keyName, options)

# ------------------------------------------------------------------------------
        
# water
# BehavioralInterval (water_left)
# BehavioralInterval (water_right)
def process_water(orig_h5, nwb_object):
    hash_group_pointer  = orig_h5['eventSeriesArrayHash']
    source = "Intervals are as reported in somatosensory cortex data file"
    series_path="/stimulus/presentation"

    # left water
    keyName1 = "leftReward"
    description1 = h5lib.get_description_by_key(hash_group_pointer, keyName1)
    grp = h5lib.get_value_by_key(hash_group_pointer , keyName1)
    t = grp["eventTimes/eventTimes"].value * 0.001
    group_attrs = {"source" : source}
    data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0,\
                   "description" : description1}
    water_left = create_time_series("water_left", "<IntervalSeries>", \
                     series_path, orig_h5, nwb_object, group_attrs, '', t, \
                     data_attrs, hash_group_pointer, keyName1, options)

    # right water 
    keyName2 = "rightReward"
    description2 = h5lib.get_description_by_key(hash_group_pointer , keyName2)
    grp = h5lib.get_value_by_key(hash_group_pointer , keyName2)
    t = grp["eventTimes/eventTimes"].value * 0.001
    group_attrs = {"source" : source}
    data_attrs  = {"unit": "None", "conversion": 1.0, "resolution":0.0,\
                   "description" : description2}   
    water_right = create_time_series("water_right", "<IntervalSeries>", \
                     series_path, orig_h5, nwb_object, {}, '', t, \
                     data_attrs, hash_group_pointer, keyName2, options)
        
# ------------------------------------------------------------------------------

# auditory_cue
# BehavioralInterval (auditory_cue)
def process_cue(orig_h5, nwb_object, options):

    series_path = "/stimulus/presentation"

    if "eventSeriesArrayHash" in h5lib.get_child_group_names(orig_h5) and \
                  "rewardCue" in h5lib.get_key_list(orig_h5['eventSeriesArrayHash']):
        # SP data
        # auditory cue 
        hash_group_pointer  = orig_h5['eventSeriesArrayHash']
        keyName = "rewardCue"
        description = h5lib.get_description_by_key(hash_group_pointer, keyName)
        grp = h5lib.get_value_by_key(hash_group_pointer, keyName)
        t = grp["eventTimes/eventTimes"].value * 0.001
        group_attrs = {"description" : "Intervals when auditory cue presented",\
                       "source": "Intervals are as reported in Simon's data file"}
        data_attrs = {"unit": "None", "conversion": 1.0, "resolution": 0.0}
        cue_ts = create_time_series("auditory_cue", "<IntervalSeries>", series_path,
                     orig_h5, nwb_object, group_attrs, '', t, data_attrs, \
                     hash_group_pointer, keyName, options)
    elif "CueTime" in h5lib.get_key_list(orig_h5['trialPropertiesHash']):
        # NL data
        keyName = "CueTime"
        trial_start_times = orig_h5["trialStartTimes/trialStartTimes"].value
        grp = orig_h5["trialPropertiesHash/value/"]
        time = grp["3/3"].value
        t = time + trial_start_times
        hash_group_pointer  = orig_h5['trialPropertiesHash']
        description = h5lib.get_description_by_key(hash_group_pointer, keyName)
        group_attrs = {"comments" : description,\
                       "description" : keyName, \
                       "source" : "Times are as reported in Nuo's data file, but relative to session time"}
        data_attrs  = {"resolution":float('nan'),"unit":"unknown","conversion":1.0, \
                       "description" : description, \
                       "source": "Times are as reported in Nuo's data file, but relative to session time"}
        trial_start_times = orig_h5["trialStartTimes/trialStartTimes"].value
        cue_ts = create_time_series("auditory_cue", "<IntervalSeries>", series_path,\
                     orig_h5, nwb_object, group_attrs, '', t, data_attrs, \
                     hash_group_pointer, keyName, options)

# ------------------------------------------------------------------------------

def process_whisker(orig_h5, nwb_object, options):
    # create module
    mod = nwb_object.make_group("<Module>", "Whisker")
    mod.set_custom_dataset("description", \
        "Whisker angle and curvature (relative) of the whiskers and times when the pole was touched by whiskers")

    # Create interface
    whisker_iface = mod.make_group("BehavioralTimeSeries", attrs={
        "source": "Whisker data"})

    # Create time series
    keyName = 'whiskerVars'
    hash_group_pointer  = orig_h5['timeSeriesArrayHash']
    grp   = h5lib.get_value_by_key(hash_group_pointer , keyName)   
# GD scaling must be read in from input file
    t     = grp["time/time"].value * 0.001
    
    if 'eventSeriesArrayHash' in orig_h5.keys() and \
       'touches' in h5lib.get_key_list(orig_h5['eventSeriesArrayHash']):
        # SP data

        var   = grp["valueMatrix/valueMatrix"].value 
        if options.handle_errors:
            try:
                grp = orig_h5["timeSeriesArrayHash/value/1"]
            except:
                grp = orig_h5["timeSeriesArrayHash/value"]
        else:
            grp = orig_h5["timeSeriesArrayHash/value/1"]
        descr = parse_h5_obj(grp["idStrs"])[0]
 
        # whisker angle time series
        # embedded nans screw things up -- remove them
        # count how many non-nan values, and prepare output array
        angle = var[0]
        group_attrs={"description": descr[0],
                     "source": "Whisker angle as reported in Simon's data file"}
        data_attrs ={"unit": "degrees", "conversion":1.0, "resolution":0.001} 
        ts_angle = create_time_series("whisker_angle", "<TimeSeries>", "",\
                       orig_h5, whisker_iface, group_attrs, angle, t, data_attrs, \
                       hash_group_pointer, keyName, options)
        ts_angle.set_attr("description", "Angle of whiskers")
        ts_angle.set_attr("source", "Whisker angle as reported in Simon's data file")

        # whisker curvature
        curv  = var[1]
        group_attrs={"description": descr[1],
                     "source": "Curvature (relative) of whiskers as reported in Simon's data file"}
        data_attrs={"unit":"Unknown", "conversion": 1.0, "resolution": 1.0}
        ts_curve = create_time_series("whisker_curve", "<TimeSeries>", "",\
                       orig_h5, whisker_iface, group_attrs, curv, t, data_attrs, \
                       hash_group_pointer , keyName, options)

        # pole touches
        pole_iface = mod.make_group("BehavioralEpochs", attrs={
                          "source": "Pole intervals as reported in Simon's data file"})
        keyName = "touches"
        hash_group_pointer  =  orig_h5["eventSeriesArrayHash"]
        grp = h5lib.get_value_by_key(hash_group_pointer , keyName)

        # protraction touches
        t = grp["eventTimes/1/1"].value * 0.001
        group_attrs = {"description" : "Intervals that whisker touches pole (protract)",\
                       "source" : "Intervals are as reported in Simon's data file"}
        if len(t) > 0 or not options.handle_errors:
            pole_touch_pr = create_time_series("pole_touch_protract", "<IntervalSeries>",\
                                "", orig_h5, pole_iface, group_attrs, '', t, {},\
                                hash_group_pointer, keyName, options)


        # retraction touches
        t = grp["eventTimes/2/2"].value * 0.001
        group_attrs = {"description" : "Intervals that whisker touches pole (retract)",\
                       "source" : "Intervals are as reported in Simon's data file"}
        if len(t) > 0 or not options.handle_errors:
            pole_touch_re = create_time_series("pole_touch_retract", "<IntervalSeries>",\
                                "", orig_h5, pole_iface, group_attrs, '', t, {},\
                                hash_group_pointer, keyName, options)
    else:   
        # JY data
        group_attrs = {"source" : "Times as reported in Jianing's data file",
                       "description" : "Time moment that whisker touches pole"}
        num_vars = len(np.array(grp["id/id"]).tolist())
        if len(t) > 0:
            valueMatrix = np.array(grp["valueMatrix/valueMatrix"])
            idStr         = np.array(grp["idStr/idStr"])
            try:
                idStrDetailed = np.array(grp["idStrDetailed/idStrDetailed"])
            except:
                idStrDetailed = idStr
            for i in range(0, num_vars):
                var         = valueMatrix[i]
                series_name = idStr[i]
                description = idStrDetailed[i]
                whisker_var = create_time_series(series_name, "<TimeSeries>", "", \
                                  orig_h5, whisker_iface, group_attrs, var, t, {},\
                                  hash_group_pointer, keyName, options)

# ------------------------------------------------------------------------------

def process_pole_touches(orig_h5, nwb_object, options):
    keyName2 = "kappaMaxAbsOverTouch"
    kappa_ma_pr = h5lib.get_value2_by_key2(orig_h5["eventSeriesArrayHash"], keyName1, \
                                                   "eventPropertiesHash/1", keyName2)
    pole_tp_path = "processing/Whisker/BehavioralEpochs/pole_touch_protract"
    if options.handle_errors:
        try:
            pole_tp_grp = nwb_object.file_pointer[pole_tp_path]
            pole_tp_grp.create_dataset("kappa_max_abs_over_touch", data=kappa_ma_pr)
        except:
            print "Cannot create dataset processing/Whisker/BehavioralTimeSeries/pole_touch_protract/kappa_max_abs_over_touch"
    else:
        pole_tp_grp = nwb_object.file_pointer[pole_tp_path]
        pole_tp_grp.create_dataset("kappa_max_abs_over_touch", data=kappa_ma_pr)

    # add kappaMaxAbsOverTouch to pole_touch_retract
    kappa_ma_re = h5lib.get_value2_by_key2(orig_h5["eventSeriesArrayHash"], keyName1, \
                                                   "eventPropertiesHash/2", keyName2)
    pole_tr_path = "processing/Whisker/BehavioralEpochs/pole_touch_retract"
    try:
        pole_tr_grp = nwb_object.file_pointer[pole_tr_path]
        pole_tr_grp.create_dataset("kappa_max_abs_over_touch", data=kappa_ma_re)
    except: 
        print "Cannot create dataset processing/Whisker/BehavioralTimeSeries/pole_touch_retract/kappa_max_abs_over_touch"

# ------------------------------------------------------------------------------

def process_stimulus(orig_h5, nwb_object):
    keyName = "StimulusPosition"
    stim_pos = h5lib.get_value_by_key(orig_h5['trialPropertiesHash'], keyName)
    trial_t = orig_h5["trialStartTimes/trialStartTimes"].value * 0.001
    rate = (trial_t[-1] - trial_t[0])/(len(trial_t)-1)
    description = h5lib.get_description_by_key(orig_h5["trialPropertiesHash"], keyName)
    zts = nwb_object.make_group("<TimeSeries>", "zaber_motor_pos", path="/stimulus/presentation",\
                           attrs={"description": description})
    zts.set_dataset("timestamps", trial_t)
    zts.set_dataset("data", stim_pos, attrs={"unit":"unknown",\
                    "conversion": 1.0, "resolution":1.0})

# ------------------------------------------------------------------------------

def process_intracellular_ephys_data(orig_h5, meta_h5, nwb_object, options):
    mod = nwb_object.make_group("<Module>", "Intracellular_ephys")
    keyName = 'Ephys'
    hash_group_pointer  = orig_h5['timeSeriesArrayHash']
    mod_descr = h5lib.get_description_by_key(hash_group_pointer, keyName)
    mod.set_attr("description", mod_descr)

    # Create interface
    ephys_iface = mod.make_group("FilteredEphys", attrs = \
                     {"source" : "Intracellular ephys data as reported in Jianing's data file"})

    # Create time series
    grp   = h5lib.get_value_by_key(hash_group_pointer , keyName)
    time = np.transpose(np.array(grp["time/time"]))
    t = []
    for i in range(time.shape[0]):
        t = t + time[i,:].tolist()
    source = "Times as reported in Jianing's data file"
    num_ts = 0
    if len(t) > 0:
        valueMatrix = np.array(grp["valueMatrix/valueMatrix"])
        idStr       = np.array(grp["idStr/idStr"])
        try:
            idStrDetailed = np.array(grp["idStrDetailed/idStrDetailed"])
        except:
            idStrDetailed = idStr
        num_var = len(idStr)                   
#       print "num_var=", num_var 
        for i in range(num_var):
            var = valueMatrix[i]
            series_name = idStr[i]
            if options.handle_errors:
                try:
                    electrode_idx = meta_h5["intracellular/recording_coord_location/recording_coord_location"]
                except:
                    electrode_idx = meta_h5["intracellular/recordingLocation/recordingLocation"]
            else:
                electrode_idx = meta_h5["intracellular/recording_coord_location/recording_coord_location"]
            group_attrs = {"description" : idStr[i], "source" : source,\
                           "electrode_idx" : electrode_idx}
            ephys_var = create_time_series(series_name, "<ElectricalSeries>", "", \
                        orig_h5, ephys_iface, group_attrs, var, np.array(t), {},\
                        hash_group_pointer, keyName, options)

# ------------------------------------------------------------------------------

# trial start times are stored in: h5::trialStartTimes::trialStartTimes
# trial stop isn't stored. assume that it's twice the duration of other
#   trials -- padding on the high side shouldn't matter
def create_trials(orig_h5, nwb_object, options):
    trial_id = orig_h5["trialIds/trialIds"].value
    if options.verbose:
        print "\nCreating trials with ids: ", trial_id
    trial_t  = orig_h5["trialStartTimes/trialStartTimes"].value * 0.001
    # trial stop isn't stored. assume that it's twice the duration of other
    #   trials -- padding on the high side shouldn't matter
    ival     = (trial_t[-1] - trial_t[0]) / (len(trial_t) - 1)
    trial_t  = np.append(trial_t, trial_t[-1] + 2*ival)

    if "GoodTrials" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
        # NL data only
        good_trials = h5lib.get_value_by_key(orig_h5["trialPropertiesHash"],  "GoodTrials")
        ignore_ivals_start = [time for (time, good_trial) in zip(trial_t,good_trials) if good_trial == 0]
        # trial stop isn't stored. assume that it's twice the duration of other
        #   trials -- padding on the high side shouldn't matter
        ignore_ivals_stop = [time for (time, good_trial) in zip(trial_t[1:],good_trials) if good_trial == 0]
        ignore_intervals  = [ignore_ivals_start, ignore_ivals_stop]

    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        start = trial_t[i]
        stop  = trial_t[i+1]
        # pole_pos_path = "trialPropertiesHash/value/3/3"
        #         pole_pos = str(orig_h5[pole_pos_path].value[i])
        #         epoch.description = ("Stimulus position - in Zaber motor steps (approximately, 10,000 = 1 mm): " + pole_pos)
        if "GoodTrials" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
            # NL data
            epoch = ut.create_epoch(nwb_object, trial, start, stop)
            tags = []
            if good_trials[i] == 1:
                tags.append("Good trial")
            else:
                tags.append("Non-performing")
            for j in range(len(epoch_tags[trial])):
                tags.append(epoch_tags[trial][j])
            epoch.set_dataset("tags", tags)

            # keep with tradition and create a units field, even if it's empty
            if trial not in epoch_units:
                units = []
            else:
                units = epoch_units[trial]
            epoch.set_custom_dataset("units", units)

            raw_path = "descrHash/value/%d" % (trial_id[i])
            raw_file = parse_h5_obj(orig_h5[raw_path])[0]
            if len(raw_file) == 1:
                raw_file = 'na'
            else:
                raw_file = str(raw_file)

            epoch.set_dataset("description", "Raw Voltage trace data files used to acuqire spike times data: "\
                                 + raw_file \
                                 + "\n\ignore intervals: mark start and stop times of bad trials "\
                                 + " when mice are not performing")
            #epoch.set_ignore_intervals(ignore_intervals)
            # collect behavioral data
            events = {"auditory_cue"    : "/stimulus/presentation/auditory_cue",\
                      "pole_in"         : "/stimulus/presentation/pole_in",\
                      "pole_out"        : "/stimulus/presentation/pole_out",\
                      "lick_trace"      : "/acquisition/timeseries/lick_trace",\
                      "aom_input_trace" : "/stimulus/presentation/aom_input_trace",\
                      "simple_optogentic_stimuli" \
                                        : "/stimulus/presentation/simple_optogentic_stimuli"}
            for key in events.keys():
                if options.verbose:
                    print "Creating epoch ts ", key, " ... "
                ut.add_epoch_ts(epoch, start, stop, key, events[key])
            
        elif "TrialName" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
            # JY data
#           print "    ...trial=", trial
            epoch = nwb_object.make_group("<epoch_X>", trial, attrs = \
                        {"description" : "Data that belong to " + trial})
            ts = "/processing/Licks/BehavioralEvents/lick_time"              
            ut.add_epoch_ts(epoch, start, stop, "lick_time", ts)
            # loop through whiskerVars
            keyName = 'whiskerVars'
            hash_group_pointer  = orig_h5['timeSeriesArrayHash']
            grp   = h5lib.get_value_by_key(hash_group_pointer , keyName)

            valueMatrix = np.array(grp["valueMatrix/valueMatrix"])
            idStr         = np.array(grp["idStr/idStr"])
            try:
                idStrDetailed = np.array(grp["idStrDetailed/idStrDetailed"])
            except:
                idStrDetailed = idStr
            num_vars = len(np.array(grp["id/id"]).tolist())
            for i in range(0, num_vars):
                var         = valueMatrix[i]
                series_name = str(idStr[i])
                ts = "/processing/Whisker/BehavioralTimeSeries/" + series_name
#               print "    ... Adding epoch ts for series_name=", series_name
                ut.add_epoch_ts(epoch, start, stop, series_name, ts)

        else:
            # SP data
            epoch = nwb_object.make_group("<epoch_X>", trial)
            epoch.set_dataset("start_time", start)
            epoch.set_dataset("stop_time", stop)

            if trial in epoch_roi_list:
                epoch.set_custom_dataset("ROIs", epoch_roi_list[trial])
                epoch.set_custom_dataset("ROI_planes", epoch_roi_planes[trial])
            tags = []
            if trial in epoch_trial_types:
                for j in range(len(epoch_trial_types[trial])):
                    tags.append(epoch_trial_types[trial][j])
            epoch.set_dataset("tags", tags)
            epoch.set_dataset("description", "Data that belong to " + trial)

            ts = "/processing/Licks/BehavioralEvents/lick_left"
            ut.add_epoch_ts(epoch, start, stop, "lick_left", ts)
            ts = "/processing/Licks/BehavioralEvents/lick_right"
            ut.add_epoch_ts(epoch, start, stop, "lick_right", ts)
            ts = "/stimulus/presentation/water_left"
            ut.add_epoch_ts(epoch, start, stop, "water_left", ts)
            ts = "/stimulus/presentation/water_right"
            ut.add_epoch_ts(epoch, start, stop, "water_right", ts)
            ts = "/stimulus/presentation/pole_accessible"
            ut.add_epoch_ts(epoch, start, stop, "pole_accessible", ts)
            ts = "/processing/Whisker/BehavioralEpochs/pole_touch_protract"
            if options.handle_errors:
                try:
                    ut.add_epoch_ts(epoch, start, stop, "pole_touch_protract", ts)
                except:
                    print "cannot add epoch for time series=", ts
            else:
                ut.add_epoch_ts(epoch, start, stop, "pole_touch_protract", ts)
            ts = "/processing/Whisker/BehavioralEpochs/pole_touch_retract"
            if options.handle_errors:
                try:
                    ut.add_epoch_ts(epoch, start, stop, "pole_touch_retract", ts)
                except:
                    print "cannot add epoch for time series=", ts
            else:
                ut.add_epoch_ts(epoch, start, stop, "pole_touch_protract", ts)
            ts = "/stimulus/presentation/auditory_cue"
            ut.add_epoch_ts(epoch, start, stop, "auditory_cue", ts)
            ts = "/processing/Whisker/BehavioralTimeSeries/whisker_angle"
            ut.add_epoch_ts(epoch, start, stop, "whisker_angle", ts)
            ts = "/processing/Whisker/BehavioralTimeSeries/whisker_curve"
            ut.add_epoch_ts(epoch, start, stop, "whisker_curve", ts)

# ------------------------------------------------------------------------------

# each subarea has list of trials and ROI ids
# to parse, take each subarea, pull out trials
#   foeach trial, write ROI
#   to get plane, find ROI in imaging plane
epoch_roi_list = {}
epoch_roi_planes = {}
def create_trial_roi_map(orig_h5, nwb_object, plane_map, options):
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
            print "  Warning: in create_trial_roi_map num_planes == 1"

        for j in range(num_planes):
            plane = j + 1
            try:
                grp_name = "imagingPlane/%d/ids/ids" % plane
                plane_id = block[grp_name].value
            except:
                print "Warning: only one imaging plane detected for area ", area
                grp_name = "imagingPlane/ids/ids"
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
            valid_whisker = get_valid_trials(orig_h5, "whisker", options)
            valid_Ca = get_valid_trials(orig_h5, "Ca", options)
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
            epoch_roi_list[trial_name] = roi_list
            epoch_roi_planes[trial_name] = plane_list

# ------------------------------------------------------------------------------

# add trial types to epoch for indexing
epoch_tags = {}
epoch_trial_types = {}
def get_trial_types(orig_h5, nwb_object, options):
#   fp = nwb_object.file_pointer
    trial_id = orig_h5["trialIds/trialIds"].value
    trial_types_all = []
    trial_type_strings = parse_h5_obj(orig_h5['trialTypeStr'])[0]
    num_trial_types    = len(trial_type_strings)
    if "PhotostimulationType" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
        # NL data
        photostim_types = h5lib.get_value_by_key(orig_h5["trialPropertiesHash"], "PhotostimulationType")
    elif "StimulusPosition" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
        # SP data
        valid_whisker = get_valid_trials(orig_h5, "whisker", options)
        valid_Ca      = get_valid_trials(orig_h5, "Ca", options)
    else:
        # JY data
        trial_types = h5lib.get_value_by_key(orig_h5["trialPropertiesHash"], 'TrialName')

    # collect all trials (strings)
    for i in range(num_trial_types):
        trial_types_all.append(str(trial_type_strings[i]))
    # write specific entries for the given trial
    for i in range(len(trial_id)):
        tid = trial_id[i]
        trial_name = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
        epoch_tags[trial_name] = []
        trial_types = []
        trial_type_mat = parse_h5_obj(orig_h5['trialTypeMat'])[0]
        for j in range(num_trial_types):
            if trial_type_mat[j,i] == 1:
                epoch_tags[trial_name].append(trial_types_all[j])
                trial_types.append(trial_types_all[j])
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
            epoch_tags[trial_name].append(photostim_type)                              
        elif "StimulusPosition" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
            # SP data
            if i in valid_whisker:
                trial_types.append("Valid whisker data")
            else:
                trial_types.append("Invalid whisker data")
            if i in valid_Ca:
                trial_types.append("Valid Ca data")
            else:
                trial_types.append("Invalid Ca data")
            epoch_trial_types[trial_name] = trial_types
#       else:
#           # JY data
#           nwb_object.create_dataset("trial_name", data=trial_types)
        
# ------------------------------------------------------------------------------
        
def get_valid_trials(orig_h5, data, options):
    ts_path = "timeSeriesArrayHash/descrHash/"
    val = []
    if options.handle_errors:
        num_subareas = 1
        if '1' in orig_h5['timeSeriesArrayHash/descrHash'].keys():
            num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1
    else:
        num_subareas = len(orig_h5['timeSeriesArrayHash/descrHash'].keys()) - 1

    if data == "whisker":
        if options.handle_errors:
            try:                      
                ids = parse_h5_obj(orig_h5[ts_path + '1/value'])[0]
            except:
                ids = parse_h5_obj(orig_h5[ts_path + 'value'])[0]
        else:
            ids = parse_h5_obj(orig_h5[ts_path + '1/value'])[0]
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

def set_metadata(group, keyname, value):
    print "    set_metadata: keyname=", keyname, " value=", value
    if keyname in ["extracellular", "intracellular"]:
        group.set_dataset("description", value)
    else:
        print "keyname=", keyname, " value=", value
        group.set_custom_dataset(keyname, value)   

# ------------------------------------------------------------------------------

def set_metadata_from_file(group, keyname, file_name):
    group.set_dataset(keyname, ut.load_file(os.path.join(os.environ['NWB_DATA'],file_name)))

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
    general_group = nwb_object.make_group("general", abort=False)
    subject_group = general_group.make_group("subject", abort=False)

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
        if key in ["extracellular", "intracellular"]:
            ephys_group = general_group.make_group(key + "_ephys", abort=False)
        try:
            if key in ["behavior", "virus", "fiber", "photostim", "extracellular", "intracellular"]:
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

        if key in ["animalGeneCopy", "animalGeneticBackground", "animalGeneModification"]:
            genotype += key + ": " + str(value) + "\n"               
   
        elif re.search("animalStrain", key):
            animalStrain += key + ": " + value + " "

        elif re.search("animalSource", key):
            animalSource += key + ": " + value + " "

        elif key == "animalID":
            set_metadata(subject_group, "subject_id", value)
 
        elif key in ["dateOfBirth"]:
            subject  = key + ": " + str(value) + "\n"

        elif re.search("citation", key):
            set_metadata(general_group, "related_publications", value)

        elif re.search("experimentType", key):
            set_metadata(general_group, "notes", value)

        elif re.search("experimenters", key):
            set_metadata(general_group, "experimenter", value)

        elif key in ["sex", "species", "age", "cell", "weight"]:
            set_metadata(subject_group, key, value)

        elif re.search("referenceAtlas", key):
            general_group.set_custom_dataset("reference_atlas", value)

        elif re.search("whiskerConfig", key):
            general_group.set_custom_dataset("whisker_configuration", value)

        elif key in ["virus", "fiber", "photostim", "surgicalManipulation"]:
            general_group.set_custom_dataset(key, value)

        elif key == "behavior":
            task_kw = map(str,parse_h5_obj(input_h5["behavior/task_keyword"])[0])
            nwb_object.set_custom_dataset("task_keyword", task_kw)

        elif key in ["extracellular", "intracellular"]:
            ephys_group.set_custom_dataset(key, value)

    set_metadata(subject_group, "genotype", genotype + "\n")
    set_metadata(subject_group, "description", animalStrain + "  " + animalSource + "  " + subject)

    set_metadata_from_file(general_group, "surgery", "surgery.txt")
    set_metadata_from_file(general_group, "data_collection", "data_collection.txt")
    set_metadata_from_file(general_group, "experiment_description", "experiment_description.txt")

# ------------------------------------------------------------------------------

def process_behavioral_data(orig_h5, nwb_object, options):

    # NL, SP and JY data:
    if options.verbose:
        print "    Processing pole position ..."
    process_pole_position(orig_h5, nwb_object, options)

    if options.verbose:
        print "    Processing licks data ..."
    process_licks(orig_h5, nwb_object)

    if ("eventSeriesArrayHash" in h5lib.get_child_group_names(orig_h5) and \
                   "rewardCue" in h5lib.get_key_list(orig_h5['eventSeriesArrayHash'])) or \
                     "CueTime" in h5lib.get_key_list(orig_h5['trialPropertiesHash']):
        if options.verbose:
            print "    Processing cue ..."
        process_cue(orig_h5, nwb_object, options)

    if "whiskerVars" in h5lib.get_key_list(orig_h5['timeSeriesArrayHash']):
        # SP and JY data:
        print "    Processing whisker data ..."
        process_whisker(orig_h5, nwb_object, options)
        
        if not 'LickTime' in h5lib.get_key_list(orig_h5['trialPropertiesHash']):
            # SP data only
            print "    Processing water data ..."
            process_water(orig_h5, nwb_object)
#           print "    Processing touches data ..."
#           process_pole_touches(orig_h5, nwb_object, options)
            print "    Processing stimulus data ..."
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
epoch_units = {}
def get_trial_units(orig_h5, nwb_object, unit_num):
#   fp = nwb_object.file_pointer
    for i in range(unit_num):
        i = i+1
        unit = "unit_%d%d" % (int(i/10), i%10)
        grp_name = "eventSeriesHash/value/%d" % i
        grp_top_folder = orig_h5[grp_name]
        trial_ids = grp_top_folder["eventTrials/eventTrials"].value
        trial_ids = Set(trial_ids)
        print "trial_ids=", trial_ids
        for trial_num in trial_ids:
            tid = trial_num
            trial_name = "Trial_%d%d%d" % (int(tid/100), int(tid/10)%10, tid%10)
            if trial_name not in epoch_units:
                epoch_units[trial_name] = []
            epoch_units[trial_name].append(unit)

# ------------------------------------------------------------------------------

def create_epochs(orig_h5, nwb_object, options):
    if options.verbose:
        print "Creating epochs"

    print "    Getting trial types ..."
    get_trial_types(orig_h5, nwb_object, options)

    if "GoodTrials" in h5lib.get_key_list(orig_h5["trialPropertiesHash"]):
        # NL data  only
        print "    Getting trial units ..."
        unit_num = len(orig_h5['eventSeriesHash/value'].keys())
        get_trial_units(orig_h5, nwb_object, unit_num)

    print "    Creating trials ..."    
    create_trials(orig_h5, nwb_object, options)

# ------------------------------------------------------------------------------

def process_ROIs_and_dFoverF(orig_h5, dff_iface, seg_iface, num_subareas, \
                             plane_map, master_shape, options):
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
            fetch_rois(orig_h5, master_shape, plane_map, seg_iface, subarea+1, \
                       plane+1, options, num_planes)
            fetch_dff( orig_h5, dff_iface, seg_iface, plane_map, subarea+1, \
                       plane+1, options, num_planes)

            sys.stdout.write('.')
            sys.stdout.flush()
    print ""
    print "Writing ROI and dF/F"

# ------------------------------------------------------------------------------

def create_master_shape(plane_map, num_subareas, orig_h5, nwb_object, options):
    master_shape = {}
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            master_shape = create_reference_image(orig_h5, nwb_object, \
                               master_shape, plane_map,subarea+1, \
                               plane+1, options, num_planes)
            sys.stdout.write('.')
            sys.stdout.flush()
    return master_shape

# ------------------------------------------------------------------------------

def add_reference_images_to_image_segmentation(seg_iface, plane_map, \
        reference_image_red, reference_image_green):
    for k in plane_map.keys():
        plane = plane_map[k]
        try:
            img = reference_image_red[plane]
            ut.add_reference_image(seg_iface, plane, "%s_red_channel"%plane, img)
        except:
            print "Cannot store red reference image"
        try:
            img = reference_image_green[plane]
            ut.add_reference_image(seg_iface, plane, "%s_green_channel"%plane, img)
        except:
            print "Cannot store green reference image"

# ------------------------------------------------------------------------------

def process_image_data(orig_h5, nwb_object, plane_map, options):
    # store master images
    num_subareas = len(h5lib.get_key_list(orig_h5['timeSeriesArrayHash'])) - 1
    master_shape = create_master_shape(plane_map, num_subareas, orig_h5, \
                                       nwb_object, options)

    print "Creating reference images"
    num_subareas = len(h5lib.get_key_list(orig_h5['timeSeriesArrayHash'])) - 1
    for subarea in range(num_subareas):
        plane_path = 'timeSeriesArrayHash/descrHash/%d/value/1' %(subarea + 2)
        if orig_h5[plane_path].keys()[0] == 'masterImage':
            num_planes = 1
        else:
            num_planes = len(orig_h5[plane_path].keys())
        for plane in range(num_planes):
            sys.stdout.write('_')
# GD: strange: the two loops are almost exactly repeated
#     the role of the 1st set of loops is just to wrire symbol '_' ti screen ?!!
    print ""

 # create empty time series for whisker video
    whisker_vid = nwb_object.make_group("<ImageSeries>", "whisker_video", path="/acquisition/timeseries")
    whisker_vid.set_dataset("format", "unknown");
    whisker_vid.set_dataset("external_file", ["whisker_data.tar.gz"], attrs={"starting_frame": [0]})
    whisker_vid.set_dataset("dimension", [1024, 1024])
    whisker_vid.set_attr("description", "Place-holder for whisker tracking video. Video information and timestamps not available")
    whisker_vid.set_dataset("timestamps", "extlink:external_timestamps_link,/path/to/timestamps")
    whisker_vid.set_dataset("num_samples", 0)

    mod = nwb_object.make_group("<Module>", "ROIs")
    mod.set_custom_dataset("description", "Segmentation (pixel-lists) and dF/F (dffTSA) for all ROIs")

    dff_iface = mod.make_group("DfOverF")
    dff_iface.set_attr("source", "This module's ImageSegmentation interface")

    seg_iface = mod.make_group("ImageSegmentation", attrs={"source": "Simon's datafile"})

    process_ROIs_and_dFoverF(orig_h5, dff_iface, seg_iface, num_subareas, \
                             plane_map, master_shape, options)

    # add reference images to image segmentation
    # TODO
    add_reference_images_to_image_segmentation(seg_iface, plane_map, reference_image_red, \
                                               reference_image_green) 

    print "Creating map between trials and ROIs"
    create_trial_roi_map(orig_h5, nwb_object, plane_map, options)


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
            if options.handle_errors:
                try:
                    pgrp = grp["%d"%(plane+1)]
                except:
                    # Warning: only one imaging plane is available (instead of 3)
                    pgrp = grp
            else:
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
            nname = None
            fname = None
            zero = np.zeros(1)
            assert len(t) == len(frame_idx[0])
            # following arrays used to make external_file as an array, reducing number of image_series
            external_file = []
            starting_frame = []
            timestamps = []
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
                if not (lastfile == filenum and framenum == lastframe+1) and \
                   not                          framenum == 1:
                    # Warning: framenum or filenum does not start from 1
                    continue
#               assert (lastfile == filenum and framenum == lastframe+1) or framenum == 1
                if lastfile != filenum:
                    if i>0:
                        if not np.isnan(frame_idx[0][i-1] ) and \
                           not np.isnan(frame_idx[1][i-1]):
                            if options.handle_errors:
                                try:
                                     save_2p_frames(external_file, starting_frame, \
                                                    timestamps, fname, stack_t)
                                except:
                                    print "Warning: unable to create_2p_ts for name=", name
                            else:
                                save_2p_frames(external_file, starting_frame, \
                                               timestamps, fname, stack_t)
                            stack_t = []
                            fname = None
                lastframe = framenum
                lastfile = filenum
            # make sure we write out the last entry
            if fname is not None:
                if options.handle_errors:
                    try:
                        save_2p_frames(external_file, starting_frame, \
                                       timestamps, fname, stack_t)
                    except:
                        print "Warning: unable to create_2p_ts for name=", name
                else:
                    save_2p_frames(external_file, starting_frame, \
                                   timestamps, fname, stack_t)
            create_2p_tsa(nwb_object, nname, external_file, starting_frame, \
                          timestamps, name)
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
        if key in ["recordingType", "penetrationN", "groundCoordinates","referenceCoordinates", \
                   "extracellularDataType", "cellType", "identificationMethod","spikeSorting"]:
            key1_list = h5lib.get_value_pointer_by_path_items(extra_grp, [key]).keys()
            for key1 in key1_list:
#               print "      key=", key, " key1=", key1
                value1 = h5lib.get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                value2 = [str(v) for v in value1]
                value += "       " + key1 + ": " + ",".join(value2) + "\n "
#               print "key=", key, " value=", value
    return value

# ------------------------------------------------------------------------------

def get_device(meta_h5, options):
    extra_grp = h5lib.get_value_pointer_by_path_items(meta_h5, \
                                ["extracellular"])
    value = ""
    for key in extra_grp.keys():
        if key in ["probeSource", "probeType", "ADunit", "amplifierRolloff"]:
            key1_list = h5lib.get_value_pointer_by_path_items(extra_grp, [key]).keys()
            for key1 in key1_list:
#               print "      key=", key, " key1=", key1
                value1 = h5lib.get_value_pointer_by_path_items(extra_grp, [key, key1])[:]
                value2 = [str(v) for v in value1]
                value += "       " + key1 + ": " + ",".join(value2) + "\n "
    return value

# ------------------------------------------------------------------------------

#add empty ElectricalSeries to acquisition
def create_empty_acquisition_series(name, num):
    #- vs = nuo.create_timeseries("ElectricalSeries", name, "acquisition")
    vs = nuo.make_group("<ElectricalSeries>", name, path="/acquisition/timeseries")
    data = [0]
    vs.set_attr("comments","Acquired at 19531.25Hz")
    vs.set_attr("source", "Device 'ephys-acquisition'")
    vs.set_dataset("data", data, attrs={"unit": "none", "conversion": 1.0, "resolution": float('nan')})
    timestamps = [0]
    vs.set_dataset("timestamps", timestamps)
    el_idx = 8 * num + np.arange(8)
    vs.set_dataset("electrode_idx", el_idx)
    vs.set_attr("description", "Place-holder time series to represent ephys recording. Raw ephys data not stored in file")

# ------------------------------------------------------------------------------

def process_ephys_electrode_map(nwb_object, meta_h5, options):
    M = read_probe_locations_matrix(meta_h5, options)

    # probe = M.tolist()
    probe = []
    sites = parse_h5_obj(check_entry(meta_h5, "extracellular/siteLocations"))
    assert len(sites) == 32, "Expected 32 electrode locations, found %d"%len(sites)
    for i in range(len(sites)):
        probe.append(sites[i])
        probe[-1] = probe[-1] * 1.0e-6
    probe = np.asarray(probe)

    num_shanks, shank_size, shank_coords = detect_shanks(M)
    shank = []
    for i in range(1, (num_shanks+1)):
        for j in range(shank_size):
            shank.append("shank" + str(i))
    gg = nwb_object.make_group("general", abort=False)
    ee = gg.make_group("extracellular_ephys", abort=False)
    ee.set_dataset('electrode_map', probe)
    ee.set_dataset('electrode_group', shank)
    ee.set_dataset('filtering', "Bandpass filtered 300-6K Hz")

    ephys_device_txt = "32-electrode NeuroNexus silicon probes recorded on a PCI6133 National Instrimunts board. See 'general/experiment_description' for more information"
    nwb_object.set_dataset("<device_X>", ephys_device_txt, name="ephys-acquisition")
    nwb_object.set_dataset("<device_X>", "Stimulating laser at 473 nm", name="optogenetic-laser")

    # Creating the shank groups 
    probe_type = h5lib.get_value_pointer_by_path_items(meta_h5, \
                     ["extracellular", "probeType", "probeType"]).value
#   print "\nprobe_type=", probe_type
    rloc = h5lib.get_value_pointer_by_path_items(meta_h5, \
               ["extracellular", "recordingLocation", "recordingLocation"]).value
    description = get_description(meta_h5, options)
    device      = get_device(     meta_h5, options)
    if options.verbose:
        print "description=", description
        print "device=", device
        print "rloc=", rloc
    for i in range(num_shanks):
        loc = str(rloc[0])
        P = str(shank_coords[i][0])
        Lat = str(shank_coords[i][1])
        location = "loc: " + loc + ", P: " + P + ", Lat: " + Lat + ", recordingLocation=" + rloc
        eg = ee.make_group("<electrode_group_X>", "shank" + str(i))
        eg.set_dataset("location",    location)
        eg.set_dataset("description", description)
        eg.set_dataset("device",      device)
    
# ------------------------------------------------------------------------------

def process_laser_and_aom_input_data(orig_h5, nwb_object, options):
    # raw data section
    # lick trace is stored in acquisition
    # photostimulation wave forms is stored in stimulus/processing
    if options.verbose:
        print "Processing laser and aom data ..."
    # get times
    grp_name = "timeSeriesArrayHash/value/time/time"
    timestamps = orig_h5[grp_name].value
    # calculate sampling rate
    rate = (timestamps[-1] - timestamps[0])/(len(timestamps)-1)
    # get descriptions
    comment1 = parse_h5_obj(orig_h5["timeSeriesArrayHash/keyNames"])[0][0]
    comment2 = parse_h5_obj(orig_h5["timeSeriesArrayHash/descr"])[0][0]
    comments = comment1 + ": " + comment2
    grp_name = "timeSeriesArrayHash/value/idStrDetailed"
    description = parse_h5_obj(orig_h5[grp_name])[0]

    # laser data
    keyName1 = "EphusVars"
    keyName2 = "laser_power" 
    hash_group_pointer = orig_h5["timeSeriesArrayHash"]               
    laser_power = h5lib.get_value2_by_key2(orig_h5["timeSeriesArrayHash"], \
                      keyName1, "",keyName2)
    group_attrs = {"description" : description[2], "comments" : comments, \
                   "source" : "Nuo's data file"}
    data_attrs = {"resolution":float('nan'), "unit":"Watts", "conversion":1000.0}
    laser_ts = create_time_series("simple_optogentic_stimuli", \
                       "<OptogeneticSeries>", "/stimulus/presentation", \
                       orig_h5, nwb_object, group_attrs, laser_power, timestamps, \
                       data_attrs, hash_group_pointer, keyName1, options);

    # aom input data
    keyName1 = "EphusVars"
    keyName2 = "aom_input_trace"    
    descr = description[1]
    aom_input_trace= h5lib.get_value2_by_key2(orig_h5["timeSeriesArrayHash"], keyName1, \
                                                      "", keyName2)
    group_attrs = {"description" : description[1], "comments" : comments, \
                   "source" : "Nuo's data file"}
    data_attrs = {"resolution":float('nan'), "unit":"Volts", "conversion":1.0}
    aom_ts = create_time_series("aom_input_trace", "<TimeSeries>", \
                       "/stimulus/presentation", orig_h5, nwb_object, group_attrs, \
                       aom_input_trace, timestamps, data_attrs, \
                       hash_group_pointer, keyName1, options);

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
    mod = nwb_object.make_group("<Module>", "Units")
    mod.set_custom_dataset('description', 'Spike times and waveforms')

    # create interfaces
    spk_waves_iface = mod.make_group("EventWaveform")                
    spk_waves_iface.set_attr("source", "Data as reported in Nuo's file")
    spk_times_iface = mod.make_group("UnitTimes")  
    spk_times_iface.set_attr("source", "EventWaveform in this module")
          
    # top level folder
    grp_name = "eventSeriesHash/value"
    # determine number of units
    unit_num = len(orig_h5[grp_name].keys())
    # initialize cell_types and electrode_depth arrays with default values
    cell_types = ['unclassified']*unit_num
    n = max(range(unit_num)) + 1
    electrode_depths = np.zeros(n)
    unit_descr = parse_h5_obj(orig_h5['eventSeriesHash/descr'])[0]
    # process units
    for i in range(unit_num):
        i = i+1
        unit = "unit_%d%d" % (int(i/10), i%10)
        # initialize timeseries
        spk = spk_waves_iface.make_group("<SpikeEventSeries>", unit)
        # get data
        grp_name = "eventSeriesHash/value/%d" % i
        grp_top_folder = orig_h5[grp_name]
        timestamps = grp_top_folder["eventTimes/eventTimes"]
        trial_ids = grp_top_folder["eventTrials/eventTrials"]
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
        cell_types[i-1] = unit + " - " + cell_type
        try:
            # read in electrode depths and update electrode_depths array
            depth = parse_h5_obj(grp_top_folder["depth"])[0]
            # depth = grp_top_folder["depth/depth"][0,0]
            electrode_depths[i-1] = depth
            electrode_depths[i-1] = 0.001 * depth
        except:
            if options.verbose:
                print "Could not extract electrode_depths"
        # fill in values for the timeseries
        spk.set_custom_dataset("sample_length", sample_length)
        spk.set_attr("source", "---")
        spk.set_dataset("timestamps", timestamps)
        spk.set_dataset("data", waveforms, attrs={"resolution":float('nan'), "unit":"Volts", "conversion":0.1})
        spk.set_dataset("electrode_idx", [channel[0]])
        spk.set_attr("description", cell_type)
        # spk_waves_iface.add_timeseries(spk)
        # add spk to interface
        #description = unit_descr[i-1] + " -- " + cell_type
        # spk_times_iface.add_unit(unit, timestamps, cell_type, "Data from processed matlab file")
        ug = spk_times_iface.make_group("<unit_N>", unit)
        ug.set_dataset("times", timestamps)
        ug.set_dataset("source", "Data from processed matlab file")
        ug.set_dataset("unit_description", cell_type)
        # spk_times_iface.append_unit_data(unit, "trial_ids", trial_ids)
        ug.set_custom_dataset("trial_ids", trial_ids)
    # spk_times_iface.set_value("CellTypes", cell_types)
    spk_times_iface.set_custom_dataset("CellTypes", cell_types)
    spk_times_iface.set_custom_dataset("ElectrodeDepths", electrode_depths)

# ------------------------------------------------------------------------------

# Extract all keys from input data and check if there are keys not in the pre-defined list
def check_keys(orig_h5, meta_h5, options):              
    known_keys=[ 'pole_in_reach', 'touches', 'left_licks', 'right_licks', 'left_reward', 'right_reward', 'reward_cue', 'species', 'animal_strain1', 'animal_source1', 'animal_gene_modification1', 'animal_strain2', 'animal_source2', 'animal_gene_modification2', 'sex', 'date_of_experiment', 'time_of_experiment', 'experimenters', 'whisker_vars', 'dffTSA', 'LWaterValveTime', 'RWaterValveTime', 'StimulusPosition', 'EphusVars' ,  'PoleInTime', 'PoleOutTime', 'CueTime', 'GoodTrials', 'PhotostimulationType', 'animal_gene_copy', 'animal_gene_modification', 'animal_genetic_background', 'animal_ID', 'animal_source', 'animal_strain', 'behavior', 'citation', 'date_of_birth', 'experiment_type', 'extracellular', 'photostim', 'reference_atlas', 'surgical_manipulation', 'unit', 'virus', 'weight', 'weight_before', 'whisker', 'fiber', 'Ephys', 'cell', 'PolePos', 'LickTime', 'TrialName', 'PassiveTouch', 'intracellular']
    known_keys_camel=[ 'poleInReach', 'touches', 'leftLicks', 'rightLicks', 'leftReward', 'rightReward', 'rewardCue', 'species', 'animalStrain1', 'animalSource1', 'animalGeneModification1', 'animalStrain2', 'animalSource2', 'animalGeneModification2', 'sex', 'dateOfExperiment', 'timeOfExperiment', 'experimenters', 'whiskerVars', 'dffTSA', 'LWaterValveTime', 'RWaterValveTime', 'StimulusPosition', 'EphusVars' ,  'PoleInTime', 'PoleOutTime', 'CueTime', 'GoodTrials', 'PhotostimulationType', 'animalGeneCopy', 'animalGeneModification', 'animalGeneticBackground', 'animalID', 'animalSource', 'animalStrain', 'behavior', 'citation', 'dateOfBirth', 'experimentType', 'extracellular', 'photostim', 'referenceAtlas', 'surgicalManipulation', 'unit', 'virus', 'weight', 'weightBefore', 'whisker', 'fiber', 'Ephys', 'cell', 'PolePos', 'LickTime', 'TrialName', 'PassiveTouch', 'intracellular']
    Jaining_keys = ['weightAfter', 'whiskerConfig']
    known_keys = known_keys + Jaining_keys + known_keys_camel
    if options.verbose:
        print "known_keys=", sorted(known_keys)
    unknown_keys = []
    all_keys = h5lib.get_all_keys(orig_h5, meta_h5)
    if options.verbose:
        print "\nAll keys=", sorted(all_keys)
    for k in all_keys:
        if not k in known_keys and \
           not re.search("unit", k) and \
           not re.search("dffTSA", k):
            unknown_keys.append(k)
#   if len(unknown_keys) > 0:
#       sys.exit("\nUnknown keys found: " + str(unknown_keys))
    return 

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
    
    trial_type_mat = orig_h5['trialTypeMat/trialTypeMat'].value
    good_trials = orig_h5['trialPropertiesHash/value/4/4'].value
    grp = nwb_object.make_group("analysis", abort=False)
    grp.set_custom_dataset("trial_start_times", trial_start_times)
    grp.set_custom_dataset("trial_type_string", trial_types_all)
    grp.set_custom_dataset("trial_type_mat",    trial_type_mat)
    grp.set_custom_dataset("good_trials",       good_trials)

# ------------------------------------------------------------------------------

#ep.ts_self_check()
#sys.exit(0)
def produce_nwb(data_path, metadata_path, output_nwb, options):
    orig_h5 = h5py.File(data_path, "r")
    meta_h5 = ""
    if len(metadata_path) > 0:
        meta_h5 = h5py.File(metadata_path, "r")

    check_keys(orig_h5, meta_h5, options)              

    print "Input:", data_path, " ", metadata_path
    print "Output:", output_nwb, "\n"
    if options.replace and os.path.exists(output_nwb):
        os.remove(output_nwb)

    # each of simon's hdf5 files have imaging planes and subareas
    #   labels consistent within the file, but inconsistent between
    #   files. create a map between the h5 plane name and the 
    #   identifier used between files

    vargs = {}
    session_id = os.path.basename(output_nwb)[0:-4]
    if len(metadata_path) > 0: 
        if "date" in h5lib.get_child_group_names(orig_h5):
            # JY data                        
            vargs["start_time"] = find_exp_time(orig_h5)
        else:
            # NL data
            vargs["start_time"] = find_exp_time(meta_h5)
    else:                                             
        # SP data
        vargs["start_time"] = find_exp_time(orig_h5)

    vargs["file_name"]      = output_nwb
    vargs["identifier"]     = ut.create_identifier(session_id)                              
    vargs["description"]    = os.path.join(os.environ['NWB_DATA'],"experiment_description.txt")
    vargs["mode"] = "w"

#   print "vargs=", vargs
    nwb_object = nwb_file.open(**vargs)

    # Process metadata
    print "Processing metadata ..."
    if "metaDataHash" in h5lib.get_child_group_names(orig_h5):
        # SP data
        process_metadata(nwb_object, orig_h5, options)
    elif not len(meta_h5) == 0:
        # NL data
        process_metadata(nwb_object, meta_h5, options)
    else:
        sys.exit("\nCannot process metadata. Check your input.")

    # Process behavioral time series
    print "Processing behavioral data ..."
    process_behavioral_data(orig_h5, nwb_object, options)

    if len(metadata_path) == 0:
        # SP data
        plane_map = {}
        create_plane_map(orig_h5, plane_map, options)
        print "After create_plane_map: plane_map=", plane_map
        if len(plane_map.keys()) > 0:
            process_image_data(orig_h5, nwb_object, plane_map, options)
    else:
        # NL or JY data
        if "extracellular" in meta_h5.keys():
            # NL data only
            print "Processing extracellular ephys data ..."
            process_ephys_electrode_map(nwb_object, meta_h5, options)
            process_laser_and_aom_input_data(orig_h5, nwb_object, options)
            process_extracellular_spike_time(orig_h5, nwb_object, options)
            collect_analysis_information(orig_h5, nwb_object, options)
        if "intracellular" in meta_h5.keys():
            # JY data only
            print "Processing intracellular ephys data ..."
            process_intracellular_ephys_data(orig_h5, meta_h5, nwb_object, options)

    # Create epochs
    print "Creating epochs ..."
    create_epochs(orig_h5, nwb_object, options)

    if options.verbose:
        print "Closing file"

    nwb_object.close()

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


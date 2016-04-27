#! /usr/local/python-2.7.6/bin/python     

import sys
import os 
import h5py
import re        
import numpy as np
import h5lib

# ------------------------------------------------------------------------------

# If the top hash group is specified, return its keys; 
# otherwise, return the keys of all top hash groups
if __name__ == "__main__":
    if len(sys.argv) in [2,3,4]:
        h5_data    = sys.argv[1]
        hash_group = ''
        if len(sys.argv) >= 3:
            hash_group = sys.argv[2]
        match_string = ''
        if len(sys.argv) == 4:
            match_string = sys.argv[3]
    else:
        sys.exit("\nUsage: h5_get_keys.py <h5_file_or_dir> [ <hash_group> [ match_string ] ] \n")

    file_list = h5lib.get_file_list(h5_data, match_string)
          
    for h5_file in file_list: 
        orig_h5 = h5py.File(h5_file, "r")
        top_groups = h5lib.get_child_group_names(orig_h5)
        if len(hash_group) > 0:
            path_items = [hash_group, "keyNames", "keyNames"]
            key_list = h5lib.get_value_pointer_by_path_items(orig_h5, path_items)
            print h5_file, ",", hash_group, ": keys=", np.array(key_list).tolist()
        elif re.search("meta", h5_file):
            key_list = orig_h5.keys()
            print h5_file, ": keys=", np.array(key_list).tolist()
            for k in key_list:
                key_list2 = orig_h5[k].keys()
                if len(key_list2) > 1:
                    print h5_file, ", ", k, ": keys2=", np.array(key_list2).tolist()
        else:
            for group in top_groups:
                if not re.search("Hash", group):
                    continue
                group_keys = orig_h5[group].keys()
                if len(group_keys) < 2:
                    print h5_file, ",", group, ": keys=[]"
                    continue
                path_items = [group, "keyNames", "keyNames"]
                print "path_items=", path_items
                key_list = h5lib.get_value_pointer_by_path_items(orig_h5, path_items)
                print h5_file, ",", group, ": keys=", np.array(key_list).tolist()

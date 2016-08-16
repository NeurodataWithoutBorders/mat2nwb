#!/usr/local/python-2.7.8/bin/python

import sys,os 
import h5py
import re        
import numpy as np
import h5lib

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) >= 1:
        h5_data      = sys.argv[1]
        path_items = []
        if len(sys.argv) >= 2:
            path_items   = sys.argv[2:]
    else:
        sys.exit("\nUsage: h5_get_value_by_path_items.py <h5_file_or_dir> [ <group(s)> ] \n")

    verbose = 0

    # Compile a file list
    file_list = h5lib.get_file_list(h5_data, "meta")
    print "file_list=", file_list

    # Handle files in the file list 
    for h5_file in file_list: 
        orig_h5 = h5py.File(h5_file, "r")

        # Extract path in a given file
        value_pointer = ""
        try:
            value_pointer = h5lib.get_value_pointer_by_path_items(orig_h5, path_items)
            print h5_file + ": value_pointer.name=", value_pointer.name
            keys = value_pointer.keys()
            print "group: num_items=", len(keys), " items=" + str(keys)
        except:
            try:
                data = value_pointer
                print "dataset: data=", data, " type=", type(data), " shape=", np.array(data).shape, " dtype=", np.array(data).dtype
                print "data=", np.array(data).tolist()
            except:
                try:
                    data_ref = data[0][0]
                    data = orig_h5[data_ref]
                    print "dataset: data=", data, " dtype=", np.array(data).dtype
                except:
                    print "value_pointer=", value_pointer
                    sys.exit("Path is not valid")


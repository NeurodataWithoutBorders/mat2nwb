#! /usr/local/python-2.7.6/bin/python     

import sys,os 
import h5py
import re
import numpy as np
import h5lib

# ------------------------------------------------------------------------------

def parse_h5_item(h5_root, item_path, verbose):
        item_pointer = h5_root[item_path]
        try:
            # item is a group
            keys = item_pointer.keys()
            if verbose:
#               print "group: path=", item_path, " members=" + str(keys)
                print "group: path=", item_path, " num_members=", len(keys)
            for k in keys:
                if len(item_path) == 1:                # item_path == '/'
                    item_path1 = item_path + k
                else:
                    item_path1 = item_path + '/' + k
                parse_h5_item(h5_root, item_path1, verbose)
        except:
            # item is a dataset
            try:
                data = np.array(item_pointer)
                if verbose:
                    print "dataset: path=", item_path, " , shape=", data.shape, \
                          " , dtype=", data.dtype, " data=", data.tolist()
            except:
                sys.exit("Path " + path + " is not valid")

# ------------------------------------------------------------------------------

if __name__ == "__main__":
#   print "len(sys.argv)=", len(sys.argv)
    if len(sys.argv) in [2,3]:
        h5_data  = sys.argv[1]
        verbose = 1
        if len(sys.argv) == 3:
            verbose = int(sys.argv[2])
#           print "verbose=", verbose
    else:
        sys.exit("\nUsage: h5_dump_all.py <h5_file_or_dir> \n")

    file_list = h5lib.get_file_list(h5_data, "")

    # Parse the contexts of entire file, and dump its summary as text          
    for h5_file in file_list: 
        print "\nParsing file ", h5_file, "\n"
        base_name = os.path.basename(h5_file)
        h5_root  = h5py.File(h5_file, "r")
        path = '/'
        parse_h5_item(h5_root, path, verbose)


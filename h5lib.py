import numpy as np
import re
import os, sys

verbose = 0    

# ------------------------------------------------------------------------------

def get_child_group_names(parent_group_pointer):
    return parent_group_pointer.keys()

# ------------------------------------------------------------------------------

def get_key_list(hash_group_pointer):
    keyNames_dataset = hash_group_pointer['keyNames/keyNames']
    key_list = np.array(keyNames_dataset).tolist()
    return key_list                                         

# ------------------------------------------------------------------------------

def get_value_by_key(hash_group_pointer, key):
    if verbose:
        print "In get_value_by_key: group_pointer.name=", hash_group_pointer.name
        print "                                    key=", key
    key_list = get_key_list(hash_group_pointer)
    if verbose:
        print "In get_value_by_key: partial_key=", key, " key_list=", key_list
    key_name = ""
    for k in key_list:
        if re.search(key, k):
            key_name = k
    if verbose:
        print "In get_value_by_key: key_name=", key_name, "\n"
    ind = key_list.index(key_name) + 1
    if verbose:
        print "ind=", ind, "\n"
    value_group = hash_group_pointer['value']
    if verbose:
        print "In get_value_by_key: value_group.name=",  value_group.name
        print "                     value_group.keys()=", value_group.keys()
    if unicode(ind) in value_group.keys():
        value_pointer = value_group[str(ind) + '/' + str(ind)]
        if verbose:
            print "item_type(value_pointer)=", item_type(value_pointer)
            print "In get_value_by_key: return value name=", value_pointer.name
        if item_type(value_pointer) == "dataset":
            value = np.array(value_pointer).tolist()
    else:
        value = value_group['value'][ind]
        if verbose:
            print "In get_value_by_key: return value name=", value
    return value

# ------------------------------------------------------------------------------

def get_value_path_by_partial_key(hash_group_pointer, partial_key):
    key_list = get_key_list(hash_group_pointer)
    key_name = ""
    for k in key_list:
        if re.search(partial_key, k):
            key_name = k
    ind = key_list.index(key_name)
    value_group = hash_group_pointer['value/value']
    value_name = value_group[ind]
    return value_name

# ------------------------------------------------------------------------------

def get_value_pointer_by_path_items(orig_h5, path_items):  
    path = ""
    if len(path_items) > 0 and not path_items == ['']:
        for i in range(0, len(path_items)):
            if i == 0:
                path += path_items[i]
            else:
                path += '/' + path_items[i]
        value_pointer = orig_h5[path]
#       print "path=", path, " value=", value_pointer
    else:
        value_pointer = orig_h5
    if verbose:
        print "    In get_value_pointer_by_path_items: path=", path
        if item_type(value_pointer) == "group":
            print "                                        key_list=", value_pointer.keys()
    return value_pointer

# ------------------------------------------------------------------------------

def get_key_index(key_list, partial_key):
    key_name = ""
    for k in key_list:
        if re.search(partial_key, k):
            key_name = k
    ind = key_list.index(key_name)
    return ind

# ------------------------------------------------------------------------------
# Returns a pointer to the value 
def get_value_pointer_by_key(hash_group_pointer, partial_key, verbose):
    key_list = get_key_list(hash_group_pointer)             
    ind = get_key_index(key_list, partial_key) + 1
    value_group_pointer = hash_group_pointer['value/' + str(ind)]
    if verbose:
        print "   In get_value_pointer_by_key: value_group_items=", value_group_pointer.keys(), " ind=", ind
        print "                                item_type(value_group_pointer)=", item_type(value_group_pointer)
    if str(ind) in value_group_pointer.keys():
        if verbose:
            print "    In get_value_pointer_by_key: case 1"
        value_pointer1 = value_group_pointer[str(ind)]     
        if item_type(value_pointer1) == "group" and str(ind) in value_pointer1.keys():
            if verbose:
                print "    In get_value_pointer_by_key: case 11"
            value_pointer = value_pointer1[str(ind)]
        else:
            if verbose:
                print "    In get_value_pointer_by_key: case 12"
                print "                                 item_type(value_pointer1)=", item_type(value_pointer1)
                print "                                 value_pointer1.name=", value_pointer1.name
            value_pointer = value_pointer1 
    else:
        if verbose:
            print "    In get_value_pointer_by_key: case 2"
        value_pointer = value_group_pointer
    if item_type(value_pointer) == "dataset":
        value_pointer = np.array(value_pointer).tolist()              
    return value_pointer  

# ------------------------------------------------------------------------------

def get_description_by_key(hash_group_pointer, partial_key):
    key_list = get_key_list(hash_group_pointer)
    ind = get_key_index(key_list, partial_key) + 1
    descr_data = np.array(hash_group_pointer['descr/descr']).tolist()
    if verbose:
        print "   In get_description_by_key: descr_data_items=", descr_data, " ind=", ind
        print "Function get_description_by_key returns: ", descr_data[ind-1]
    return descr_data[ind]

# ------------------------------------------------------------------------------

def get_file_list(data, match_string):
    # Compile a file list
    file_list = []
    if os.path.isfile(data) and (re.search(".h5", data) or re.search(".nwb", data) or re.search(".borg", data)):
        file_list = [ data ]
    elif os.path.isdir(data):
        file_list = []
        if verbose:
            print "num_files=", len(os.listdir(data))
        for file in os.listdir(data):
            if len(match_string) > 0 and not re.search(match_string, file):
                continue
            if (re.search(".h5", file) or re.search(".nwb", data)) \
                and os.path.isfile(os.path.join(data, file)):
                file_path = os.path.join(data, file)
                file_list.append(file_path)
    else:
        sys.exit("Cannot process " + data)
    return file_list
   
# ------------------------------------------------------------------------------

def item_type(item_pointer):
    item_type = ""
    try:
        keys = item_pointer.keys()
        item_type = "group"
    except:
        try:
            keys = np.array(item_pointer).tolist() 
            item_type = "dataset"
        except:
            item_type = "data_item"
    return item_type

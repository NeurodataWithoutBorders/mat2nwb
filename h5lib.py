import numpy as np
import re
import os, sys

verbose = 0    

# ------------------------------------------------------------------------------

def get_child_group_names(parent_group_pointer):
    return parent_group_pointer.keys()

# ------------------------------------------------------------------------------

def get_group_keys(group_pointer):
    key_list = []
    for k in group_pointer.keys():
        if not re.search('#', k):
            key_list.append(k)
    return sorted(key_list)

# ------------------------------------------------------------------------------

def get_data_from_refs_dataset(orig_h5, dataset_pointer):
    path = dataset_pointer.name
    if dataset_pointer.dtype.kind == 'O':
        ref = orig_h5[path][0][0]
        obj = orig_h5[ref]
        value = ''.join(i for i in obj[:])
    else:
        obj = orig_h5[path]
        value = ''.join(i.encode(ascii) for i in obj[:])
    return value

# ------------------------------------------------------------------------------

def get_key_list(hash_group_pointer):
    try:
        keyNames_dataset = hash_group_pointer['keyNames/keyNames']
        key_list = np.array(keyNames_dataset).tolist()
    except:
        key_list = []
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
    ind = key_list.index(key_name)
    if verbose:
        print "ind=", ind, "\n"
    value_group = hash_group_pointer['value']
    if verbose:
        print "\nIn get_value_by_key: value_group.name=",  value_group.name
        print "\n                   value_group.keys()=", value_group.keys()
        print "\nlen(value_group.keys())=", len(value_group.keys()), " ind=", ind  
    if len(value_group.keys()) == 1 and "value" in value_group.keys():
        value_pointer = value_group['value']
        if item_type(value_pointer) == "dataset":
            if verbose:
                print "Case1: returning an element of a dataset"
            value = value_pointer[ind]
        else:
            if verbose:
                print "Case2: value/value type is ", item_type(value_pointer)
            sys.exit("Unsapported case in get_value_by_key()")
    else:
        if verbose:
            print "level2 key_list=", value_group.keys(), ";  ind+1=", ind+1
        if str(ind+1) in value_group.keys():
            value_pointer = value_group[str(ind+1)]
            if item_type(value_pointer) == "dataset":
                if verbose:
                    print "Case3: value/" + str(ind+1) + " is a dataset"
                value = np.array(value_pointer).tolist()
            else:
                if str(ind+1) in value_pointer.keys() and \
                   item_type(value_pointer[str(ind+1)]) == "dataset":
                    if verbose:
                        print "Case4: value/" + str(ind+1) + "/" + str(ind+1) + " is ", item_type(value_pointer)
                    value = np.array(value_pointer[str(ind+1)]).tolist()
                else:
                    if verbose:
                        print "Case5: value/" + str(ind+1) + " is ", item_type(value_pointer)
                    value =  value_pointer
        else:
            value = value_group
    return value

# ------------------------------------------------------------------------------

def get_value2_by_key2(hash_group1_pointer, key1, hash_group2, key2):
    if verbose:
        print "\n\nEntering get_value2_by_key2 ..."
    key_list = get_key_list(hash_group1_pointer)
    key_name = ""
    for k in key_list:
        if re.search(key1, k):
            key_name = k
    if verbose:
        print "In get_value_by_key: key_name=", key_name, "\n"
    ind = key_list.index(key_name)
    if verbose:
        print "ind=", ind, "\n"
    value_group = hash_group1_pointer['value']
    if verbose:
        print "\nIn get_value_by_key: value_group.name=",  value_group.name
        print "\n                   value_group.keys()=", value_group.keys()
        print "\nlen(value_group.keys())=", len(value_group.keys()), " ind=", ind
        print "In get_value2_by_key2: level2 key_list=", value_group.keys()
    print "str(ind+1)=", str(ind+1), " value_group.keys()=", value_group.keys()
    if str(ind+1) in value_group.keys():
        value_pointer = value_group[str(ind+1)]
    else:
        value_pointer = value_group
    if verbose:
        print "value_pointer.keys()=", value_pointer.keys(), "; hash_group2=", hash_group2
    if 'valueMatrix' in value_pointer.keys() and \
       'idStr'       in value_pointer.keys():
       # NL data
       print " NL data"
       key_list2 = np.array(value_pointer['idStr/idStr']).tolist()
#      try:
       ind2 = key_list2.index(key2)
       if verbose:
           print "ind2=", ind2, "\n"
       value2 = value_pointer['valueMatrix/valueMatrix'][:,ind2]
#      except:
#           sys.exit("\nCannot determine value2 in get_value2_by_key2")
    else: 
        # SP data
        try:
            value_pointer2 = value_pointer[hash_group2]
            value2 = get_value_by_key(value_pointer2, key2)
        except:
            sys.exit("\nCannot determine value2 in get_value2_by_key2")
    return value2

# ------------------------------------------------------------------------------

def get_value_pointer_by_path_items(orig_h5, path_items):  
    path = ""
    if len(path_items) > 0 and not path_items == ['']:
        for i in range(0, len(path_items)):
            if i == 0:
                path += path_items[i]
            else:
                path += '/' + path_items[i]
        if verbose:
            print "path=", path
        value_pointer = orig_h5[path]
#       print "value=", value_pointer
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
    if len(key_list) > 1: 
        ind = get_key_index(key_list, partial_key) + 1
    else:
        ind = 1
    value_group_pointer = hash_group_pointer['value']
    if str(ind) in hash_group_pointer['value/'].keys():
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

def get_all_keys(orig_h5, meta_h5):
    all_keys = []
    # Extract data keys
    top_groups = get_child_group_names(orig_h5)
    for group in top_groups:
        if not re.search("Hash", group):
            continue
        group_keys = orig_h5[group].keys()
        if len(group_keys) < 2:
            continue
        path_items = [group, "keyNames", "keyNames"]
        print "path_items=", path_items
        key_list = get_value_pointer_by_path_items(orig_h5, path_items)
        for k in key_list:
            if not k in all_keys:
                all_keys.append(k)
    # Extract metadata keys
    if len(meta_h5) > 0:
        top_groups = get_child_group_names(meta_h5)
        for k in top_groups:
            if not k in key_list:
                all_keys.append(k)
    return all_keys

# ------------------------------------------------------------------------------

def get_description_by_key(hash_group_pointer, partial_key):
    key_list = get_key_list(hash_group_pointer)
    ind = get_key_index(key_list, partial_key) + 1
    descr_data = np.array(hash_group_pointer['descr/descr']).tolist()
    if verbose:
        print "\nIn get_description_by_key: len(descr_data)=", len(descr_data), \
              "\n     descr_data_items=", descr_data, " ind=", ind
        print "\nFunction get_description_by_key returns: ", descr_data[ind-1]
    return descr_data[ind-1]

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


# ------------------------------------------------------------------------------

def get_data_by_key(nwb_root, item_path, verbose):
        item_pointer = h5_root[item_path]
        try:
            # item is group
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
            # item is dataset
            try:
                data = np.array(item_pointer)
                if verbose:
                    print "dataset: path=", item_path, " , shape=", data.shape, \
                          " , dtype=", data.dtype, " data=", data.tolist()
            except:
                sys.exit("Path " + path + " is not valid")


#!/usr/local/python-3.5.1/bin/python3

#
# Copyright (C) 2016 by Howard Hughes Medical Institute.
#               Author: Gennady Denisov
#

import scipy.io
import numpy
import h5py
import os
import time
import shutil, numpy
import sys, re, optparse
import scipy.sparse as sps

__version__ = "0.1"

def mat2hdf5_command_line_parser(parser):
    parser.add_option("-o", "--outfolder", dest="output_folder", help="parent directory for the output tracker folder",metavar="output_folder",default=".")
    parser.add_option("-i", "--inputfolder",dest="input_folder", help="folder containing choreography results", metavar="input_folder", default=".")
    parser.add_option("-m", "--metadata", dest="metadata", help="(comma-separated list of) metadata file(s) to be combined with data file(s)",metavar="metadata",default="")
    parser.add_option("-v", "--verbose",   action="store_true", dest="verbose",               help="increase the verbosity level of output", default=False)
    
    return parser

# ----------------------------------------------------------------------

def get_array_type(array):
    len_shape = len(array.shape)
    if len_shape == 0:
        array_type = type(array)
    elif array.shape[0] == 0:
        array_type = 'None'
    elif len_shape == 1:
        array_type = type(array[0])
    else:
        array_type = type(array[0][0]) 
    return array_type

# ----------------------------------------------------------------------

def item_type(item_value):
    item_type = ""
#   if sps.issparse(sps.coo_matrix(item_value)):
    if type(item_value).__name__ in ["csc_matrix", "csr_matrix"]:
        return "sparse_matrix"
    if hasattr(item_value , '__dict__'):
        return "dictionary"
    if hasattr(item_value , '__array__') and len(item_value.shape) > 0:
        return "array"
    try:
        if isinstance(item_value, basestring):
            return "string"
        if isinstance(item_value, unicode):
            return "unicode"
    except:
        if isinstance(item_value, bytes):
            return "string"
        if isinstance(item_value, str):
            return "unicode"
    if isinstance(item_value, int):
        return "int"
    if isinstance(item_value, float):
        return "float"
    if item_value.dtype in ["uint8", "uint16"]:
        return "uint"
    if isinstance(item_value, complex):
        return "complex"
    if item_value.dtype == "object":
        return "object"
    return item_type

# ----------------------------------------------------------------------

def parse_item(item_name, item_value, hdf5_group, level, upper_folder, options):
    if item_type(item_value) == "sparse_matrix":
        if options.verbose:
            print ("...Handling1 sparse_matrix ", item_name, \
                   " of dtype=", item_value.dtype, " and shape=", item_value.shape)
        array = numpy.zeros(item_value.shape, dtype = item_value.dtype)
        array[:] = item_value.todense()
        dset = hdf5_group.create_dataset(item_name, item_value.shape, \
                   dtype=item_value.dtype, data=array)
    elif item_type(item_value) == "dictionary":  # is folder
#       if re.search("scipy.sparse", type(item_value)):
#       if sps.issparse(sps.coo_matrix(item_value)):
#           print "type(", item_name, ")=", type(item_value).__name__
#           try:
#               print "class_name=", item_value._class_.__name
#           except:
#               print item_name, " has no class attribute"
        keys = item_value.__dict__.keys()   
        if options.verbose:
            print ("...Handling2 item ", item_name, \
                   " of type 'dictionary' at level ", level, \
                   " in parent folder ", upper_folder, " keys =", keys) 
        if len(keys) > 0:
            for k in sorted(keys):
                if not re.search('__', k) and not k[0]=='_':
                    hdf5_subfolder = hdf5_group.create_group(k)
                    parse_item(k, item_value.__dict__[k], hdf5_subfolder, \
                               level+1, item_name, options)
    elif item_type(item_value) == "array": 
        if options.verbose:
             print ("...Handling3 item ", item_name, " in parent folder ", upper_folder, \
                    " item type 'array' of shape ", item_value.shape, \
                    "at level ", level, " item_value=", item_value)
             if item_value.shape[0] > 0:
                 print ("   item type=", item_type(item_value[0]))
        len_shape = len(item_value.shape)
        itype = get_array_type(item_value )
        idtype = item_value.dtype
#       if not item_type(item_value[0]) == "string":
        ishape = item_value.shape
#       else:
#           ishape = tuple([len(i) for i in item_value])
#       print "ishape=", ishape

#       if not itype == 'None':
        if options.verbose:
            print ("   array len=", item_name.__len__(), \
                   " shape=", ishape, " type=", itype, " dtype=", idtype)
        if  not str(idtype) == 'object':
            if options.verbose:
                print ("item_name=", item_name, " ishape=", ishape, " idtype=", idtype.type)
            if len_shape > 0:
#               dset = hdf5_group.create_dataset(item_name, ishape, dtype=idtype, data=item_value )
                try:
                    dset = hdf5_group.create_dataset(item_name, ishape, data=item_value )
                except:
                    data = numpy.array([a.encode('utf8') for a in item_value]) # explicitly encoding
                    dset = hdf5_group.create_dataset(item_name, ishape, dtype=data.dtype, data=data)
            else: # variable-length string
                dt = h5py.special_dtype(vlen=str)
                dset = hdf5_group.create_dataset(item_name, (1,1), dtype=dt, data=str(item_value ))
        elif len_shape ==1 and item_type(item_value[0]) == "array"  and item_value[0].dtype.type is 'numpy.unicode_':
            if options.verbose:
                print ("    ...Handling4 cellular array ", item_name, \
                       " of strings, shape=", ishape, \
                       " type(item_value[0])=", item_value[0].dtype.type, \
                       " item_value[0]=", item_value[0], " at level ", level)
            dset = hdf5_group.create_dataset(item_name, ishape[0], dtype='numpy.unicode_', data=item_value ) 
        elif len_shape ==1 and item_type(item_value[0]) == "array":
            if options.verbose:
                print ("    ...Handling5 cellular array ", item_name, \
                       " of objects of type array", " at level ", level)
            for d in range(0, int(ishape[0])):
                hdf5_subfolder = hdf5_group.create_group(str(d+1)) 
                parse_item(str(d+1), item_value[d], hdf5_subfolder, level+1, \
                           item_name, options)   
        elif len_shape ==1 and item_type(item_value[0]) == "dictionary":
            if options.verbose:
                print ("    ...Handling6 cellular array of structures ", item_name, \
                       " ishape=", ishape, " at level ", level)
            for d in range(0, int(ishape[0])):
                hdf5_subfolder = hdf5_group.create_group(str(d+1))
                if item_type(item_value[d]) == "dictionary":
                    keys = item_value[d].__dict__.keys()
                    if options.verbose:
                        print ("       d=", str(d), " keys=", keys)
                    if len(keys) > 0:
                        for k in sorted(keys):
                            if not re.search('__', k) and not k[0]=='_':
                                hdf5_subfolder2 = hdf5_subfolder.create_group(k)
                                value = item_value[d].__dict__[k]
                                value_str = str(value)
                                value_len = len(value_str) 
                                if value_len > 0 and not value_str=='[]' \
                                                 and not value_str=='[[] []]':
                                    if options.verbose:
                                        print ("      item with name=", k, \
                                               " value=", value , " value_len=", value_len)
                                    parse_item(k, item_value[d].__dict__[k], \
                                               hdf5_subfolder2, level+1,item_name, options)        
                elif item_type(item_value[d]) == "array":
                    ishape_d = item_value[d].shape
                    itype_d  = item_type(item_value[d][0])
                    if options.verbose:
                        print ("    ...Handling7 tricky array ", d, " of shape ", ishape_d, " and dtype", item_value[d].dtype,\
                               " in parent folder ", upper_folder + "/" + item_name, " at level ", level)
                    dset = hdf5_subfolder.create_dataset(str(d), ishape_d, dtype=item_value[d].dtype, data=item_value[d])
                else:
                    print ("    ...Cannot1 handle item ", item_name, "of type ", item_type(item_value[d]), \
                           " in upper_folder", upper_folder)
        elif item_type(item_value[0]) in ["string", "unicode"]:
#           ishape = tuple([len(i) for i in item_value ])
            if options.verbose:
                print ("...Handling8 array of strings ", item_name, "item_value.dtype=", item_value.dtype, " ishape=", ishape)
                print ("item_value=", item_value, " at level ", level)
#           data = numpy.chararray(ishape, itemsize=100)
#           data[:] = item_value
            data = numpy.array(item_value, dtype=object)
            if options.verbose:
               print ("data=", data)
            if item_type(item_value[0]) == "unicode":
                dt = h5py.special_dtype(vlen=str)
            else:
                dt = h5py.special_dtype(vlen=bytes)
            dset = hdf5_group.create_dataset(item_name, data.shape, dtype=dt, data=data)
        elif  item_type(item_value[0]) == "int":    
            if options.verbose:
                print ("...Handling9 array of type int and size", len(item_value))
            data = numpy.array(numpy.int32(item_value))
            if options.verbose:
                print ("data=", data, " ishape=", ishape)
            dset = hdf5_group.create_dataset(item_name, ishape, dtype = numpy.int32, data=data)
        else:
            print ("    ...Cannot2 handle array ",item_name, "of idtype ",idtype,\
                   " and type ", item_type(item_value[0]), "in upper_folder", upper_folder)
    elif item_type(item_value) in ["string", "unicode"]:
        if options.verbose:
            print ("...Handling10 string2 item ", item_name, " in parent folder ", upper_folder, \
                   " item type ", item_type(item_value), \
                   "at level ", level, " item_value=", item_value)
        data = numpy.chararray((1,), itemsize=100)
        data[:] = item_value
        if options.verbose:
           print ("data=", data)
        try:
            dt = h5py.special_dtype(vlen=str)
        except:
            dt = h5py.special_dtype(vlen=unicode)
        dset = hdf5_group.create_dataset(item_name, (1,), dtype=dt, data=data)
    elif item_type(item_value) == "float":
        if options.verbose:
            print ("...Handling11 item ", item_name, " in parent folder ", upper_folder, \
                   " item type ", item_type(item_value), \
                   "at level ", level, " item_value=", item_value)
        dset = hdf5_group.create_dataset(item_name, dtype=numpy.float_, \
                                          data=numpy.float_(item_value ), \
                                          shape=(1,))
    elif item_type(item_value) == "int":
        if options.verbose:
            print ("...Handling12 item ", item_name, " in parent folder ", upper_folder, \
                   " item type ", item_type(item_value), \
                   "at level ", level, " item_value=", item_value)
        dset = hdf5_group.create_dataset(item_name, dtype=numpy.int_, \
                                          data=numpy.int_(item_value ), \
                                          shape=(1,))
    else:  
        # Create annotation
        print ("...Cannot3 handle item ", item_name, \
               " of type ", item_type(item_value), " at level=", str(level), \
               " in upper_folder=", upper_folder, " item_value=", item_value )

# ----------------------------------------------------------------------

def create_hdf5_file(mat_file_name, hdf5_file_name, options):

    # Read the input file and initialize hdf5 object
#   print "opening mat file ..."
    mat = scipy.io.loadmat(mat_file_name,squeeze_me=True,chars_as_strings=True,struct_as_record=False)
#   print "opening hdf5 file ..."
    f  = h5py.File(hdf5_file_name, 'w')  # 'w" stands for truncating if file exists

    # Parse mat file and store its contents in hdf5 object
    print ("parsing ...")
    header  = mat['__header__']
    version = mat['__version__']
#   print "Attributes=", dir(mat)
    keys   = mat.keys()
    if options.verbose:
        print ("keys=", keys, " groups:", f.keys())
    if len(keys) > 0:
        for k in sorted(keys):
            if not re.search('__', k):
                parse_item(k, mat[k], f, 1, 'top', options)

    # Add informations from metadata file
    if len(options.metadata) > 0:
        mat2 = scipy.io.loadmat(options.metadata,squeeze_me=True,struct_as_record=False)
        metadata = f.create_group('metadata')
        keys2 = mat2.keys()
        if options.verbose:
            print ("metadata keys=", keys, " groups:", metadata.keys())
        for k2 in sorted(keys2):
            if not re.search('__', k2):
                parse_item(k2, mat2[k2], metadata, 1, 'metadata', options)

    print ("...Done")

# ----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog mat_file_name [options (-h to list)]"   

    parser = optparse.OptionParser(usage=usage, version="%%prog %s" % __version__)
    parser = mat2hdf5_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.verbose:
        print ("len(args)=", len(args))
    if len(args) >= 1:
        # Extract the name of an input file
        mat_file_name = str(args[0:1][0])
        if not mat_file_name.split(".")[1] == "mat":
            sys.exit("\nNot a MAT file")
        hdf5_file_name = mat_file_name[0:len(mat_file_name)-3] + "h5"
        print ("Input:  ", mat_file_name)
        print ("Output: ", hdf5_file_name)

        create_hdf5_file(mat_file_name, hdf5_file_name, options)
    else:
        parser.print_usage()
        sys.exit(2)

"""
mat2nwb.py     

Purpose: convert all the MAT files in the specified input directory to the NWB format

Copyright (c) 2015 HHMI. All rights reserved.

"""

import os, sys 
import h5py
import getpass
import numpy as np
import re, optparse
import shutil
try:
    import commands, Set # python-2  
except:
    import subprocess    # python-3

python_version = str(sys.version_info.major) + "." + \
                 str(sys.version_info.minor)

print("python_version=" + str(python_version))

try:
    nwb_home = os.environ['NWB_HOME']
except:
    nwb_home = ""
try: 
    nwb_data = os.environ['NWB_DATA']
except:
    nwb_data = ""

# ------------------------------------------------------------------------------

def mat2nwb_command_line_parser(parser):
    parser.add_option("-A", "--project",dest="project_code",help="code to be used with qsub",metavar="project_code", default="nwb")
    parser.add_option("-D", "--debug",dest="debug",help="don't delete intermediate outputs", action="store_true", default=False)
    parser.add_option("-o", "--outfolder", dest="output_folder", help="output folder (default=same as input folder)",metavar="output_folder",default="")
    parser.add_option("-M", "--motor_cortex", action="store_true", dest="mc", help="perform conversion oof data from motor cortex", default=False)
    parser.add_option("-n", "--node",   dest="node", help="id of the cluster node to be used", metavar="node",  default=0)
    parser.add_option("-p", "--processing_start",dest="processing_start",help="start processing from mat2h5 (=1) or make_*_nwb (=2) ",metavar="processing_start",default=1)
    parser.add_option("-P", "--processing_end",dest="processing_end",help="complete processing at step mat2h5 (=1) or make_*_nwb (=2)",metavar="processing_end",default=2)
    parser.add_option("-r", "--replace",   action="store_true", dest="replace", help="if the output file already exists, replace/overwrite it", default=False)
    parser.add_option("-S", "--somatosensory_cortex", action="store_true", dest="ssc", help="perform conversion of data from somatosensory cortex", default=False)
    parser.add_option("-v", "--verbose",   action="store_true", dest="verbose",  help="increase the verbosity level of output", default=False)
    parser.add_option("-w", "--overwrite",   action="store_true", dest="overwrite", help="perform processing even if output already exists", default=False)
    parser.add_option("-z", "--zmin",dest="zmin",help="min file # to be processed", metavar="zmin", default=0)
    try:
        parser.add_option("-Z", "--zmax",dest="zmax",help="max file # to be processed", metavar="zmax", default=sys.maxint)  # python-2
    except:
        parser.add_option("-Z", "--zmax",dest="zmax",help="max file # to be processed", metavar="zmax", default=sys.maxsize) # python-3
    return parser

# ------------------------------------------------------------------------------

def get_input_type(input_data, options):
    input_type = ""
    print("input_data=" + str(input_data)+ " os.path.isdir(input_data)=" + str(os.path.isdir(input_data)))
    if   os.path.isfile(input_data) or \
         os.path.isfile(os.path.join(nwb_data,input_data)):
        if mimetypes.guess_type(input_data)[0] == 'text/plain':
            input_type = "fof"
    elif os.path.isdir(input_data) or \
         os.path.isdir(os.path.join(nwb_data,input_data)):     
        input_type = "dir"
    return input_type

# ------------------------------------------------------------------------------

def compile_file_list(input_data, input_type, options):
    file_list = []
#   print "input_data=", input_data, " input_type=", input_type
    num_files = 0
    if input_type == "dir":
        if options.mc:
            file_list = [[],[]]
            for file in os.listdir(input_data):
                if re.search(".mat", file) and \
                  ( re.search("data_structure_ANM", file)  or \
                   (re.search("data_structure_JY",  file) and not re.search("meta", file)) or \
                   (re.search("_data.mat", file) and not re.search("voltage", file))):
                   data_path = os.path.join(input_data, file)
                   if re.search("data_structure_ANM", file):     # NL data
                       meta_data_path = os.path.join(input_data, "meta_data_ANM" + file[18:])
                   elif re.search("data_structure_JY",  file):   # JY data
                       meta_data_path = os.path.join(input_data, "meta" + file)
                   else:                                         # DG data
                       meta_data_path = os.path.join(input_data, file[0:-8] + "meta.mat")
                   print ("data_path=" + data_path)
                   print ("meta_data_path=" + meta_data_path)
                   print ("meta_data_exists=" + str(os.path.exists(meta_data_path)))
                   if os.path.exists(meta_data_path):
                       print (" ... num_files=" + str(num_files) + " zmin=" + str(options.zmin) + " zmax=" +str(options.zmax))
                       if num_files >= int(options.zmin) and num_files <= int(options.zmax):
                           file_list[0].append(data_path)
                           file_list[1].append(meta_data_path)
                           num_files += 1
            print ("    len(file_list[0])=" + str(len(file_list[0])))
            print ("    len(file_list[1])=" + str(len(file_list[1])))
        elif options.ssc:
#           print "num_input files=", len(os.listdir(input_data))
            for file in os.listdir(input_data):
                if re.search("data_struct.mat", file):
                   if num_files >= int(options.zmin) and num_files <= int(options.zmax):
                       data_path = os.path.join(input_data, file)
                       file_list.append(data_path)
                   num_files += 1
#           print "    len(file_list)=", len(file_list)
#   print "    num_files=", num_files
    if num_files == 0:
        sys.exit("Incorrectly specified input data: num_files == 0")
    return file_list

# ------------------------------------------------------------------------------

def create_array_job(outfolderpath, input_data, file_list, options):
    if options.mc:
        conversion_shell_script_path = \
            os.path.join(nwb_data, "Conversion_script.mc.sh")
    elif options.ssc:
        conversion_shell_script_path = \
            os.path.join(outfolderpath, "Conversion_script.ssc.sh")

    scr = open(conversion_shell_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    if options.mc:
        scr.write("#$ -t 1-" + str(len(file_list[0])) + "\n")
    else:
        scr.write("#$ -t 1-" + str(len(file_list)) + "\n")
    command1 = ""
    command2 = ""
    if not int(options.processing_start) == 2:
        command1 = "python3 " + os.path.join(nwb_data, "mat2nwb.py") 
    if not int(options.processing_end) == 1:
        command2 = "python  " + os.path.join(nwb_data, "mat2nwb.py")
    opts = " " + os.path.join(nwb_data, str(input_data)) + \
        " -n $SGE_TASK_ID " + \
        " -p " + str(options.processing_start) + \
        " -P " + str(options.processing_end)   + \
        " -o " + outfolderpath 

    if len(command1) > 0:
        command1 += opts 
        if options.verbose:
            command1 += " -v "
        if options.mc:
            command1 += " -M "
        if options.ssc:
            command1 += " -S "
        if options.debug:
            command1 += " -D "
        command1 += "\n"
        scr.write(command1)

    if len(command2) > 0:
        command2 += opts
        if options.verbose:
            command2 += " -v "
        if options.mc:
            command2 += " -M "
        if options.ssc:
            command2 += " -S "
        if options.debug:
            command2 += " -D "
        command2 += "\n"
        scr.write(command2)

    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", \
                  conversion_shell_script_path]))
    scr.write("\n")
    scr.close()

    return conversion_shell_script_path

# ------------------------------------------------------------------------------

def submit_job(conversion_shell_script_path, options):
    base_command = os.path.join(os.environ['SGE_ROOT'], "bin", "lx-amd64", "qsub")
    prog_name = "mc.conv"
    if options.ssc:
        prog_name = "ssc.conv"
    command = ""
    jobid = 0
    command = base_command + " -V -N " + prog_name
    if not options.debug:
        command += " -o /dev/null -e /dev/null "
    if len(options.project_code) > 0:
        command += "  -A " + options.project_code
    command += " -pe batch 2"
    command += " " + conversion_shell_script_path
    try:
        res = commands.getstatusoutput(command)   # python-2
    except:
        res = subprocess.getstatusoutput(command) # python-3
    jobid = (res[1].split()[2]).split(".")[0]
    print (str(res[1])+ "\n")
    if options.verbose:
        print ("Submit conversion job command=" + command)
    return jobid

# ------------------------------------------------------------------------------

def high_level_processing(input_data, input_type, options):
    # Create a list of files to be processed 
    file_list = compile_file_list(input_data, input_type, options)
    if options.verbose:
        print ("file_list=" + str(file_list))

    print ("input_data=" + str(input_data))
    if len(options.output_folder) == 0:
        if input_type == "dir":
            if  os.path.isabs(input_data):
                outfolderpath = input_data
            else:
                outfolderpath = os.path.join(nwb_data, input_data)
        else:
            if os.path.isabs(input_data):
                outfolderpath = os.path.dirname(input_data)
            else:
                outfolderpath = nwb_data 
    else:
        if os.path.isabs(options.output_folder):
            # Absolute path is specified 
            outfolderpath = options.output_folder
        else:
            outfolderpath = os.path.join(nwb_data, options.output_folder)

    if options.ssc and len(file_list) == 0:
        sys.exit("No valid input data found. Please, check your input")
   
    if options.mc and len(file_list[0]) == 0:
        sys.exit("No valid input data found. Please, check your input")

    conversion_shell_script_path = \
        create_array_job(outfolderpath, input_data, file_list, options)
    if options.verbose:
        print ("conversion_shell_script_path=" + str(conversion_shell_script_path))
  
    submit_job(conversion_shell_script_path, options)

# ------------------------------------------------------------------------------

def low_level_processing(input_data, input_type, options):
    print ("... Enter low_level_processing; options.mc=" + str(options.mc) + " options.ssc=" + str(options.ssc)+ " options.processing_start=" + 
          str(options.processing_start))
    file_list = compile_file_list(input_data, input_type, options)
    node = int(options.node)
    if options.verbose:
        print ("node=", options.node)
        print ("file_list=", str(file_list))
    if options.mc:
        data_file      = file_list[0][node - 1]
        meta_data_file = file_list[1][node - 1]
        if int(options.processing_start) == 1:
            executable_mat2h5 = os.path.join(nwb_home, "mat2h5.py")
            command_data = executable_mat2h5 + " " + data_file 
            if options.replace:
                command_data += " -r "
            if options.verbose:
                command_data += " -v "
                print ("Running command: "+ command_data)
            os.system(command_data)
            command_meta_data = executable_mat2h5 + " " + meta_data_file
            if options.verbose:
                command_meta_data += " -v "
                print ("Running command: "+ command_meta_data)
            os.system(command_meta_data)
        if int(options.processing_end) == 2:
            executable_h52nwb = sys.executable + " " + os.path.join(nwb_home, "make_nwb.py")
            command = executable_h52nwb + " " + data_file.split(".")[0] + ".h5" + " " +\
                                           meta_data_file.split(".")[0] + ".h5"
            if options.replace:
                command += " -r "
            if options.verbose:
                command += " -v "
                print ("Running command: "+ command)
            os.system(command)
    elif options.ssc:
        data_file      = file_list[node - 1]
        if int(options.processing_start) == 1:
            executable_mat2h5 = os.path.join(nwb_home, "mat2h5.py")
            command_data = executable_mat2h5 + " " + data_file
            if options.replace:
                command_data += " -r "
            if options.verbose:
                command_data += " -v "
                print ("Running command: "+ command_data)
            os.system(command_data)
        if int(options.processing_end) == 2:
            executable_h52nwb = sys.executable + " " + os.path.join(nwb_home, "make_nwb.py")
            command = executable_h52nwb + " " + data_file.split(".")[0] + ".h5"
            if options.replace:
                command += " -r "
            if options.verbose:
                command += " -v "
                print ("Running command: "+ command)
            os.system(command)

# ------------------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog input_data [options (-h to list)] \nwhere\n\
    input_data=\n\
    - path to a directory of input MAT files, or \n\
    - path to a text file of files (FOF), which contains a list of paths to the MAT files to be converted"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = mat2nwb_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.verbose:
        print ("len(args)=" + str(len(args)))
    if len(args) == 1:
        if options.mc and options.ssc:
            sys.exit("Please, specify only one of options '-M' or '-S'")
        if not options.mc and not options.ssc:
            sys.exit("Please, specify one of options '-M' or '-S'")
        input_data = args[0]
        input_type = get_input_type(input_data, options)
        if len(options.output_folder) == 0:
            if input_type == "dir":
                options.output_folder = input_data
            else:
                options.output_folder = os.path.dirname(input_data)
        if options.verbose:
            print ("options.node=" + str(options.node))
        if int(options.node) == 0:
            high_level_processing(input_data, input_type, options)
        else:
            low_level_processing(input_data, input_type, options)
    else:
        parser.print_usage()
        sys.exit(2)


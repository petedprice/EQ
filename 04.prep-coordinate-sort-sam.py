#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' Prepares COORDINATE SORT
Takes an infolder containing folders of sam files and prepares .sh scripts to run
samtools sort. Makes a new folder for the scripts (called wdpath/scripts).
'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder containing folders of sam files which to sort")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if is_number(f)]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if name.endswith("HISAT2.sam"):
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    infolders = list_folder(args.infolder)
    print "Number of infolders =", len(infolders)

    #Make new folder for scripts
    for infolder in infolders:
        script_folder_path = os.path.dirname(infolder)+"/coord_sort_scripts"
        print "Making folder for scripts in", script_folder_path
        os.makedirs(script_folder_path)
        break

    for infolder in infolders:
        infile = list_files(infolder)
        if len(infile) == 1:
            infile = infile[0]
            wd = os.path.dirname(infile)
            sample_id = os.path.basename(infile).split("_")[0]
            scriptname = script_folder_path+"/"+sample_id+"_coord_sort.sh"
            bamoutfile = wd+"/"+os.path.basename(infile).split(".sam")[0]+"_coord_sorted"

        with open(scriptname, "w") as outfile:
            outfile.write("#!/bin/bash -l")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# Batch script to run an OpenMP threaded job on Legion with the upgraded")
            outfile.write("\n")
            outfile.write("# software stack under SGE.")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 1. Force bash as the executing shell.")
            outfile.write("\n")
            outfile.write("#$ -S /bin/bash")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).")
            outfile.write("\n")
            outfile.write("#$ -l h_rt=5:0:0")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 3. Request 1 gigabyte of RAM for the entire job (independent of thread number)")
            outfile.write("\n")
            outfile.write("#$ -l mem=10G")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 4. Set the name of the job.")
            outfile.write("\n")
            outfile.write("#$ -N sort_"+sample_id)
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 5. Select 12 threads (the most possible on Legion).")
            outfile.write("\n")
            outfile.write("#$ -pe smp 1")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 7. Set the working directory to somewhere in your scratch space.  This is")
            outfile.write("\n")
            outfile.write("# a necessary step with the upgraded software stack as compute nodes cannot")
            outfile.write("\n")
            outfile.write("# write to $HOME.")
            outfile.write("\n")
            outfile.write("# Replace <your_UCL_id> with your UCL user ID :)")
            outfile.write("\n")
            outfile.write("#$ -wd "+wd)
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 8. Run the application.")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("echo started")
            outfile.write("\n")
            outfile.write("date")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("source /usr/local/extras/Genomics/.bashrc")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("samtools view -Su ")
            outfile.write(str(infile))
            outfile.write(" | samtools sort - ")
            outfile.write("-o ")
            outfile.write(bamoutfile)
            outfile.write("\n")
            outfile.write("\n")

            outfile.write("echo finished")
            outfile.write("\n")
            outfile.write("date")

if __name__ == "__main__":
    main()

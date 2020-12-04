#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' Prepares HISAT2
Takes a folder containing paired-end fastq files and prepares .sh scripts to run
HISAT2. Makes a new folder for the scripts (called wdpath/scripts), and a new folder for 
each sample (called wdpath/sample).Assumes naming WTCHG_208_1_forward_paired.fastq.
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
                    help="Infolder of fastq files on which to run HISAT2")
parser.add_argument("wdpath", type=str,
                    help="Working directory path")
parser.add_argument("hisatindex", type=str,
                    help="Path and name of HISAT2 index")
parser.add_argument("qualityscores", type=str,
                    help="Type of quality scoring eg phred33")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("fastq")]
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    infiles = list_folder(args.infolder)
    print "Number of infiles (fastq) =", len(infiles)

    #Creates dictionary of forward and reverse pairs
    filedictionary = defaultdict(list)
    for infile in infiles:
        name = os.path.basename(infile)
        if name.split("_")[-2] == "forward":
            uniqid = name.split("_1_")[0]
        elif name.split("_")[-2] == "reverse":
            uniqid = name.split("_2_")[0]
        filedictionary[uniqid].append(infile)
    print "Number of pairs of infiles (fastq) =", len(filedictionary)

    #check each one has forward and reverse
    for sample in filedictionary:
        check = []
        for s in filedictionary[sample]:
            check.append(s.split("_paired.fastq")[0].split("_")[-1])
        if len(check) == 2 and str("forward") in check and str("reverse") in check:
                pass
        else:
            print "ERROR-pairs incorrectly assigned"
            print filedictionary
            break

    #Make new folder for scripts
    if args.wdpath.endswith("/"):
        wdpath = args.wdpath[:-1]
    else:
        wdpath = args.wdpath

    script_folder_path = wdpath+"/HISAT2_scripts"
    print "Making folder for scripts in", script_folder_path
    os.makedirs(script_folder_path)

    #Make new folders and .sh scripts
    for sample in filedictionary:
        for s in filedictionary[sample]:
            if s.split("_paired.fastq")[0].split("_")[-2] == str("1"):
                forward = s
            elif s.split("_paired.fastq")[0].split("_")[-2] == str("2"):
                reverse = s

        sample_id = sample.split("_")[-1]
        wd = wdpath+"/"+sample_id
        print "Making folder for sample", wd
        os.makedirs(wd)

        scriptname = script_folder_path+"/"+sample_id+"_HISAT2.sh"
        samoutput = wd +"/"+sample_id+"_HISAT2.sam"
        metfile = wd +"/"+sample_id+"_HISAT2.alignstats"

        with open(scriptname, "w") as outfile:
            outfile.write("#!/bin/bash")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("#$ -P ressexcon")
            outfile.write("\n")
            outfile.write("#$ -q ressexcon.q")
            outfile.write("\n")            
            outfile.write("\n")            
            outfile.write("# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).")
            outfile.write("\n")
            outfile.write("#$ -l h_rt=10:0:0")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 3. Request 1 gigabyte of RAM for the entire job (independent of thread number)")
            outfile.write("\n")
            outfile.write("\n")          
            outfile.write("#$ -l mem=5G")
            outfile.write("\n")
            outfile.write("\n")
            outfile.write("# 4. Set the name of the job.")
            outfile.write("\n")
            outfile.write("\n")   
            outfile.write("#$ -N HISAT2_"+sample_id)
            outfile.write("\n")
            outfile.write("\n")            
            outfile.write("# 5. Select 12 threads (the most possible on Legion).")
            outfile.write("\n")
            outfile.write("\n")            
            outfile.write("#$ -pe smp 12")
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
            outfile.write("hisat2")
            outfile.write(" -x ")
            outfile.write(args.hisatindex)
            outfile.write(" -1 ")
            outfile.write(forward)
            outfile.write(" -2 ")
            outfile.write(reverse)
            outfile.write(" --")
            outfile.write(args.qualityscores)
            outfile.write(" -q -p 12 --no-discordant --no-mixed --no-unal --dta")
            outfile.write(" -S ")
            outfile.write(samoutput)
            outfile.write(" --met-file ")
            outfile.write(metfile)
            outfile.write("\n")
            outfile.write("\n")

            outfile.write("echo finished")
            outfile.write("\n")
            outfile.write("date")

if __name__ == "__main__":
    main()

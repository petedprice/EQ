#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' QUALITY TRIM
Takes a folder containing paired-end fastq.gz files, finds forward and reverse pairs, and 
quality trims with Trimmomatic. 
'''
#==============================================================================
import argparse
import sys
import os
import time
from subprocess import Popen, list2cmdline
from itertools import product
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder of fastq.gz files to trim")
parser.add_argument("pathtotrimmomatic", type=str,
                    help="Path to trimmomatic eg /Users/alison/Programs/Trimmomatic-0.35/trimmomatic-0.35.jar")
parser.add_argument("pathtoadapters", type=str,
                    help="Path to adapter sequences")
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
def cpu_count():
    ''' Returns the number of CPUs in the system
    '''
    num = 1
    if sys.platform == 'win32':
        try:
            num = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (ValueError, KeyError):
            pass
    elif sys.platform == 'darwin':
        try:
            num = int(os.popen('sysctl -n hw.ncpu').read())
        except ValueError:
            pass
    else:
        try:
            num = os.sysconf('SC_NPROCESSORS_ONLN')
        except (ValueError, OSError, AttributeError):
            pass

    return num

def exec_in_row(cmds):
    ''' Exec commands one after the other until finished. This is helpful
    if a program is already parallelized, but we have to submit many jobs
    after each other.'''
    if not cmds:
        return  # empty list

    def done(p):
        return p.poll() is not None

    def success(p):
        return p.returncode == 0

    def fail():
        sys.exit(1)

    for task in cmds:
        print task
        p = Popen(task, shell=True)
        p.wait()

    if done(p):
            if success(p):
                print "done!"
            else:
                fail()

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("fastq.gz")]
#==============================================================================
#Main==========================================================================
#==============================================================================

def main():
    infiles = list_folder(args.infolder)
    print "Number of infiles (fastq.gz) =", len(infiles)

    #Creates dictionary of forward and reverse pairs
    filedictionary = defaultdict(list)
    for infile in infiles:
        name = os.path.basename(infile)
        uniqid = name.split(".fastq")[0][:-2]
        filedictionary[uniqid].append(infile)
    print "Number of pairs of infiles (fastq.gz) =", len(filedictionary)

    #check each one has forward and reverse
    for sample in filedictionary:
        check = []
        for s in filedictionary[sample]:
            check.append(s.split(".fastq")[0].split("_")[-1])
        if len(check) == 2 and str(1) in check and str(2) in check:
                pass
        else:
            print "ERROR-pairs incorrectly assigned"
            print filedictionary
            break

    #make trimmomatic commands
    print "Adapter file =", args.pathtoadapters
    print "Path to trimmomatic =", args.pathtotrimmomatic
    maketrimmomaticrun = []
    for sample in filedictionary:
        for s in filedictionary[sample]:
            if s.split(".fastq")[0].split("_")[-1] == str("1"):
                forward = s
            elif s.split(".fastq")[0].split("_")[-1] == str("2"):
                reverse = s
        output_forward_paired = forward.split(".fastq.gz")[0] + "_forward_paired.fastq.gz"
        output_reverse_paired = reverse.split(".fastq.gz")[0] + "_reverse_paired.fastq.gz"
        output_forward_unpaired = forward.split(".fastq.gz")[0] + "_forward_unpaired.fastq.gz"
        output_reverse_unpaired = reverse.split(".fastq.gz")[0] + "_reverse_unpaired.fastq.gz"
        trimmomaticcommand = ["java -jar "+args.pathtotrimmomatic+" PE -"+args.qualityscores+" "+
                                forward+" "+
                                reverse+" "+
                                output_forward_paired+" "+
                                output_forward_unpaired+" "+
                                output_reverse_paired+" "+
                                output_reverse_unpaired+" "+
                                "ILLUMINACLIP:"+args.pathtoadapters+
                                ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:95"]
        maketrimmomaticrun.append(trimmomaticcommand)

    #run trimmomatic
    print "Number of commands to run =", len(maketrimmomaticrun)
    exec_in_row(maketrimmomaticrun)

if __name__ == "__main__":
    main()

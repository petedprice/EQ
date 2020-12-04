#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' CONCATENATE SAMPLES
Takes a folder containing paired-end fastq files, finds forward and reverse pairs for each sample, 
and concatenates them. Produces one forward and one reverse file per sample.
Assumes naming is Illumina format eg WTCHG_117425_208_1.fastq where 208 is the sample number.
Important that the order of files concatenated is preserved between forward and reserve files for each sample.
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
                    help="Infolder of fastq files to concatenate")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("paired.fastq")]
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

    #Creates dictionary of forward and reverse pairs for each sample - ensures file order is preserved. Very important.
    samplefiledictionary = defaultdict(list)
    for file in filedictionary:
        sample = file.split("_")[-1]
        for f in filedictionary[file]:
            splitter = "_"+sample+"_"
            uniqid = sample+"_"+f.split(splitter)[1].split("_")[0]
            samplefiledictionary[uniqid].append(f)
    print "Number of samples =", len(samplefiledictionary)/2

    #check order is preserved
    checksamplefiledictionary = defaultdict(list)
    for sample in samplefiledictionary:
        sample_id = sample.split("_")[0]
        checksamplefiledictionary[sample_id].append(samplefiledictionary[sample])
    for sample_id in checksamplefiledictionary:
        splitter = "_"+sample_id+"_"
        checkforward = []
        checkreverse = []
        for c in checksamplefiledictionary[sample_id]:
            for file in c:
                lane = file.split(splitter)[0].split("_")[-1]
                direction = file.split(splitter)[1].split("_")[0]
                if direction == 1:
                    checkforward.append(lane)
                elif direction == 2:
                    checkreverse.append(lane)
        if checkforward == checkreverse:
            pass
        else:
            print "ERROR-lane order not preserved"
            print samplefiledictionary

    # #make concatenate commands
    makeconcatenaterun = []
    for sample in samplefiledictionary:
        files = samplefiledictionary[sample]
        direction = sample.split("_")[-1]
        if direction == "1":
            outfile = os.path.dirname(files[0])+"/"+os.path.basename(files[0]).split("_")[0]+"_"+sample+"_forward_paired.fastq"
        elif direction == "2":
            outfile = os.path.dirname(files[0])+"/"+os.path.basename(files[0]).split("_")[0]+"_"+sample+"_reverse_paired.fastq"

        infilenames = files[0]
        for file in samplefiledictionary[sample][1:]:
            infilenames = infilenames+" "+file
        
        concatenatecommand = ["cat "+str(infilenames)+" > "+outfile]
        makeconcatenaterun.append(concatenatecommand)

    #run concatenate
    print "Number of commands to run =", len(makeconcatenaterun)
    exec_in_row(makeconcatenaterun)

if __name__ == "__main__":
    main()

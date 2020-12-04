#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' EXTRACT expression for reciprocal orthologs
This script takes a file containing read counts and extracts expression for reciprocal orthologs. 
Prints read counts into one file.
'''
#==============================================================================
import argparse
import sys
import os
import cPickle as pickle
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("pklfile", type=str,
                    help="A pkl file of reciprocal orthologs")
parser.add_argument("infile", type=str,
                    help="File containing read counts")
parser.add_argument("focalspecies", type=str,
                    help="Species name for which expression has been calculated eg Poeciliareticulata")
parser.add_argument("outfile", type=str,
                    help="An outfile")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def read_pkl(pkl_file):
    '''Reads in a pickled file and returns the stored object'''
    with open(pkl_file, "rb") as infile:
        return pickle.load(infile)

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    
    #extract reciprocal orthologs from pickle file
    recip_besthit_all = read_pkl(args.pklfile)
    reciprocal_orthologs = {}
    print "Reciprocal besthits =", len(recip_besthit_all)   
    for recip in recip_besthit_all:
        for r in recip:
            species = r.split("_")[0]
            splitter = species+"_"
            transcript = r.split(splitter)[1]
            if species == args.focalspecies:
                gene = transcript.split(".")[0]+"."+transcript.split(".")[1]
                focal = gene
            else:
                ortholog = transcript
        reciprocal_orthologs[gene] = ortholog
    print "Reciprocal besthits =", len(reciprocal_orthologs)  
    
    #extract genes with expression
    counttotal = 0
    countorthologs = 0
    with open(args.outfile, "w") as outfile:
        with open(args.infile, "r") as infile:
            for line in infile:
                if line.startswith("MSTRG"):
                    counttotal += 1
                    gene = line.split("\t")[0]
                    if gene in reciprocal_orthologs:
                        countorthologs += 1
                        ortholog = reciprocal_orthologs[gene]
                        outfile.write(line)
                else:
                    outfile.write(line)
    print "Genes expressed =", counttotal
    print "Reciprocal orthologs expressed =", countorthologs
     
if __name__ == '__main__':
    main()

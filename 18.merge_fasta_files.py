#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' Merge fasta files
Takes fasta files and removes replicate sequences. Prints merged fasta file into new file.
(!Note: Assumes that if a gene is present in more than one fasta file the gene's
sequence is the same in all those files)
'''
#==============================================================================
import argparse
import sys
import os
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("fasta_folder", type=str,
                    help="Folder of fasta files")
parser.add_argument("merge_fasta", type=str,
                    help="A fasta file of merged fasta files")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]

def merge_fasta(fastafile, fastadict):
    ''' Read fasta file and outputs a dictionary 
    key [seq name] values [sequence] [header]
    '''
    try:
        with open(fastafile, "r") as infile:
            count = 0
            for line in infile.readlines():
                line = line.rstrip()
                if line[0] == ">":
                    name = line[0:]
                    count += 1
                    fastadict[name] = ["", name]
                else:
                    fastadict[name][0] += line
            print "Number of sequences in fasta file =", count
            return fastadict
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % fastafile
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    fasta_folder = list_folder(args.fasta_folder)
    print "Number of files =", len(fasta_folder)

    # Creates dictionary of gene and sequence
    fastadict = {}
    for fastafile in fasta_folder:
        print fastafile
        fasta_dict = merge_fasta(fastafile, fastadict)
    print "Number of sequences in merged file =", len(fastadict)

    # Write merged fasta file
    with open(args.merge_fasta, 'w') as outfile:
        for key, value in fasta_dict.iteritems():
            outfile.write(value[1])
            outfile.write("\n")
            outfile.write(value[0])
            outfile.write("\n")

    

if __name__ == '__main__':
    main()
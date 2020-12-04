#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' Merged exon sequence
Takes a GTF file and bedtools fasta file containing exon sequences. Extracts
exon coordinates for each transcript and merges exons to produce one sequence for
each transcript. 
'''
#==============================================================================
import argparse
import sys
from collections import defaultdict
from collections import OrderedDict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gtf", type=str,
                    help="Reference GTF file")
parser.add_argument("fasta", type=str,
                    help="Fasta file of sequences")
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
def read_scaffolds(scaff):
    ''' Read fasta file and outputs a dictionary 
    key [seq name] values [sequence] [header]
    '''
    SG={}
    try:
        with open(scaff, "r") as infile:
            for line in infile.readlines():
                line = line.rstrip()
                if line[0] == ">":
                    name = line[1:]
                    SG[name] = ["", name]
                else:
                    SG[name][0] += line
            return SG
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % scaff
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #get fasta sequences
    fasta = read_scaffolds(args.fasta)
    print "Number of unique fasta sequences =", len(fasta)

    #get exon coordinates for each isoform
    count = 0
    coordinate_dict = defaultdict(list)
    scaffold_dict = {}
    gene_list = []
    with open(args.gtf,"r") as infile:
        for line in infile:
            count += 1
            line = line.rstrip()
            if not line.startswith("#"):
                if not line.split("\t")[2] == "transcript":
                    #want to ignore transcript as gives start-end of entire gene
                    scaffold = line.split("\t")[0]
                    gene = line.split("\t")[8].split(";")[0].split('"')[1]
                    gene_list.append(gene)
                    transcript = line.split("\t")[8].split(";")[1].split('"')[1]
                    coordinate1 = line.split("\t")[3]
                    coordinate2 = line.split("\t")[4]
                    coord = [coordinate1, coordinate2]
                    #GTF file is ordered so that
                    #transcript
                    #first exon (smallest start positon)
                    #last exon (largest start positon)
                    #therefore exon order is preserved
                    coordinate_dict[transcript].append(coord)
                    scaffold_dict[transcript] = scaffold
    print "Number of lines in GTF file =", count
    print "Number of genes =", len(set(gene_list))
    print "Number of transcripts =", len(coordinate_dict)

    #check exon order is preserved
    for transcript in coordinate_dict:
        start = None
        for coord in coordinate_dict[transcript]:
            if start is None:
                start = float(coord[0])
            elif float(coord[0]) < start:
                print "ERROR - exon order not preserved"
                print coordinate_dict[transcript]

    #merge transcripts
    with open(args.outfile, "w") as outfile:
        for transcript in coordinate_dict:
            transcript_sequence = []
            scaffold = scaffold_dict[transcript]
            currentisoform = None
            #extract sequence for each exon
            for coord in coordinate_dict[transcript]:
                #bedtools fasta file has substracted one from the start coordinate in the GTF file
                start = str(float(coord[0])-1).split(".")[0]
                end = coord[1]
                if float(start) > float(end):
                    print "ERROR - coordinates reversed"
                    print coord
                else:
                    name = scaffold+":"+start+"-"+end
                    seq = fasta[name][0]
                    transcript_sequence.append(seq)
            #print sequence for whole transcript
            outfile.write(">"+transcript)
            outfile.write("\n")
            for s in transcript_sequence:
                outfile.write(s)
            outfile.write("\n")                 

if __name__ == '__main__':
    main()
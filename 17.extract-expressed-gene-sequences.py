#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import csv
import os
''' Extract transcript sequences for expressed genes
Takes a fasta file with the longest transcript sequence for each gene and extracts 
genes in the read count file. Assumes naming follows StringTie format eg 
transcript = MSTRG.22287.1 gene = MSTRG.22287.
Assumes one transcript sequences per gene.
'''
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("fasta", type=str,
                    help="A fasta file of longest transcripts")
parser.add_argument("read_counts", type=str,
                    help="A txt file of read counts across samples")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def read_fasta(source):
    ''' Read fasta file and outputs a dictionary 
    key [seq name] values [sequence] [header]
    '''
    SG = {}
    gene_transcript_id = {}
    try:
        with open(source, "r") as file:
            for line in file.readlines():
                if line[0] == ">":
                    name = line[1:].rstrip().split()[0]
                    gene = name.split(".")[0]+"."+name.split(".")[1]
                    header = line[1:].rstrip()
                    SG[name] = ["", header]
                    gene_transcript_id[gene] = name
                else:
                    SG[name][0] += line.rstrip()
            return SG, gene_transcript_id
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % source
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

def read_expr(expr):
    ''' Extract gene names in read counts file and 
    append to a list '''
    SG=[]
    try:
        with open(expr, "r") as infile:
            for line in infile.readlines():
                if line.startswith("MSTRG"):
                    geneid = line.split("\t")[0]
                    SG.append(geneid)
            return SG
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % expr
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    fasta_dict, gene_transcript_id = read_fasta(args.fasta)
    print "Number of sequences in fasta file = ", len(fasta_dict)

    expr_list = read_expr(args.read_counts)
    print "Number of genes in expression file = ", len(expr_list)

    outfile = args.read_counts+".fasta"
    count = 0
    print "Printing sequences to ...", outfile
    with open(outfile, "w") as out:
        for sequence in fasta_dict:
            gene = sequence.split(".")[0]+"."+sequence.split(".")[1]
            if gene in expr_list:
                transcript_id = gene_transcript_id[gene]
                count += 1
                out.write(">"+transcript_id)
                out.write("\n")
                out.write(fasta_dict[transcript_id][0])
                out.write("\n")
    print "Number of genes printed =", count

if __name__ == '__main__':
    main()

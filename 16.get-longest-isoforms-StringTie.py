#!/usr/bin/python
# -*- coding: utf-8 -*-

''' Get LONGEST ISOFORM
This script processes an fasta file and picks the longest isoform for each gene.
Outputs a new fasta file with ending in _longest.fasta. Assumes naming Assumes naming follows 
StringTie format eg 
transcript = MSTRG.22287.1 gene = MSTRG.22287
'''
#==============================================================================
import argparse
import sys, os, time
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("genome", type=str,
                    help="An input genome in fasta format")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================

def read_fasta(source):
    ''' Read fasta file and outputs a dictionary 
    key [seq name] values [sequence] [header]
    '''
    SG = {}
    try:
        with open(source, "r") as file:
            for line in file.readlines():
                if line[0] == ">":
                    name = line[1:].rstrip().split()[0]
                    header = line[1:].rstrip()
                    SG[name] = ["", header]
                else:
                    SG[name][0] += line.rstrip()
            return SG
    except IOError:
        print("File does not exit!")


def find_isoforms(genome):
    ''' Reads dictionary and outputs another dictionary 
    key [gene name] values [isoform name, length of isoform]
    '''
    isoforms = {}
    for seq in genome:
        gname = genome[seq][1].split(".")[0]+"."+genome[seq][1].split(".")[1]
        if gname not in isoforms:
            isoforms[gname] = [(seq, len(genome[seq][0]))]
        else:
            isoforms[gname].append((seq, len(genome[seq][0])))
    return isoforms


def write_longest_isoform(isoforms, genome, fname):
    ''' Reads dictionary, picks longest isoform and prints to fasta
    '''
    count = 0
    with open(fname, "w") as outfile:
        print(fname)
        for seq in isoforms:
            if len(isoforms[seq]) > 1:
                m = max(zip(*isoforms[seq])[1])
                longest = [p[0] for p in isoforms[seq] if p[1] is m][0]
                header = ">"+ genome[longest][1]+"\n"
                count += 1
                outfile.write(header)
                outfile.write(genome[longest][0]+"\n")
            else:
                count += 1
                longest = isoforms[seq][0][0]
                header = ">"+ genome[longest][1]+"\n"
                outfile.write(header)
                outfile.write(genome[longest][0]+"\n")
    print("Number of longest isoforms =", count)

#==============================================================================

def main():

    fname = args.genome.split("/")[-1].split(".")[0]+"_longest.fasta"
    genome = read_fasta(args.genome)
    print("Number of transcripts =", len(genome))
    isoforms = find_isoforms(genome)
    write_longest_isoform(isoforms, genome, args.genome[:-3]+"_longest.fasta")

if __name__ == "__main__":
    main()

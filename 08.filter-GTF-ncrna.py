#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' FILTER GTF file to remove ncRNA
Takes a folder of blast tophit files and filters the GTF to remove these transcripts.
In addition, it removes all other transcripts of a given gene from the GTF file if any one 
transcript maps to ncRNA.
'''
#==============================================================================
import argparse
import sys
import csv
import os
from itertools import combinations
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("blastfolder", type=str,
                    help="Infolder containing blast output files")
parser.add_argument("gtf", type=str,
                    help="GTF file")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".tophits")]

def get_tophits(file):
    '''Returns a list of all query transcripts in blast output file'''
    transcriptlist = []
    with open(file, "r") as infile:
        for line in infile:
            line = line.rstrip()
            transcript = line.split(",")[0]
            transcriptlist.append(transcript)
    return transcriptlist

def get_gene_transcript_id(file):
    gene_isoform_dict = {}
    with open(file,"r") as infile:
        for line in infile:
            line = line.rstrip()
            if not line.startswith("#"):
                if line.split("\t")[2] == "transcript":
                    gene = line.split("\t")[8].split(";")[0].split('"')[1]
                    transcript = line.split("\t")[8].split(";")[1].split('"')[1]

                    gene_isoform_dict[transcript] = gene
    return gene_isoform_dict

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    #read blast files
    blast_tophits = list_folder(args.blastfolder)
    print "Number of blast tophit files =", len(blast_tophits)

    #get list of isoforms with blast hits to ncrna
    removealltranscriptlist = []
    for blast in blast_tophits:
        print blast
        removetranscriptlist = get_tophits(blast)
        removealltranscriptlist = removealltranscriptlist+removetranscriptlist
        print "No. of transcript with hits to ncrna =", len(removetranscriptlist)
    removealltranscriptlist = set(removealltranscriptlist)
    print "Total no. of transcript to remove from GTF =", len(removealltranscriptlist)

    #get gene to transcript id
    gene_transcript_dict = get_gene_transcript_id(args.gtf)
    print "Total no. of transcripts in GTF =", len(gene_transcript_dict)

    #get genes to remove
    removeallgenelist = []
    for transcript in gene_transcript_dict:
        if transcript in removealltranscriptlist:
            gene = gene_transcript_dict[transcript]
            removeallgenelist.append(gene)
    removeallgenelist = set(removeallgenelist)
    print "Total no. of genes to remove from GTF =", len(removeallgenelist)

    #filter and print gtf file
    count = 0
    original = 0
    outfilename = args.gtf[:-4]+"_ncrnafiltered.gtf"
    print "Writing to ... ", outfilename
    with open(outfilename,"w") as outfile:
        with open(args.gtf, "r") as infile:
            for line in infile.readlines():
                original += 1
                if line.startswith("#"):
                    outfile.write(line)
                else:
                    gene = line.split("gene_id")[1].split('"')[1]
                    transcript = line.split("transcript_id")[1].split('"')[1]
                    if gene not in removeallgenelist:
                        outfile.write(line)
                        count += 1
                        if transcript in removealltranscriptlist:
                            print "ERROR - discrepancy between transcript and gene list"
                            print gene, transcript

    print "Original number of lines in GTF =", original
    print "Number of filtered lines in GTF =", count
    
if __name__ == '__main__':
    main()
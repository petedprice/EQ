#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' EXTRACT counts annotated
Takes a file of read counts across samples extracted from HTseq-count. Extracts read counts for genes on scaffolds
assigned to chromosomes and prints read counts into one file. Extracts positional information for genes on scaffolds
assigned to chromosomes and prints positional information into another file.
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
parser.add_argument("read_counts", type=str,
                    help="A txt file of read counts across samples")
parser.add_argument("assigned_scaffolds", type=str,
                    help="A txt file of scaffolds assigned to chromosomes"
                    "FORMAT: Scaffold, Chromosome, LG, Start"
                    "FORMAT: scaffold80274,NC_024331.1,Chr1,15649593.0")
parser.add_argument("gtf", type=str,
                    help="GTF file")
parser.add_argument("read_counts_outfile", type=str,
                    help="An outfile of read counts")
parser.add_argument("gene_position_outfile", type=str,
                    help="An outfile of the position of genes within the genome")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    #get coordinates for each gene from GTF file
    gene_index = {}
    gene_position = {}
    with open(args.gtf, "r") as infile:
            for line in infile.readlines():
                if not line.startswith("#"):
                    if line.split("\t")[2] == "transcript":
                        gene = line.split("gene_id")[1].split('"')[1]
                        scaffold = line.split("\t")[0]
                        start = float(line.split("\t")[3])
                        gene_index[gene] = scaffold
                        if gene in gene_position:
                            if start < gene_position[gene]:
                                gene_position[gene] = start
                        else:
                            gene_position[gene] = start
    print "Number of genes in GTF file = ", len(gene_index), len(gene_position)

    #get coordinates for each scaffolds
    scaffold_index = {}
    scaffold_position = {}
    with open(args.assigned_scaffolds, "r") as infile:
        for line in infile:
            line = line.rstrip()
            scaffold = line.split(",")[0]
            if scaffold == "Scaffold":
                line = line.split(",")
                header = line[1:]
                print header
            else:
                line = line.split(",")
                if line[2].startswith("Chr"):
                    chromosome = line[2]
                    start = float(line[3])
                    scaffold_index[scaffold] = chromosome
                    scaffold_position[scaffold] = start
    print "Number of scaffolds assigned to annotated chromosomes =", len(scaffold_index), len(scaffold_position)
    
    #extract read counts for genes on scaffolds assigned to chromosomes
    count = 0
    with open(args.read_counts_outfile, "w") as outfile_reads:
        with open(args.read_counts, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith("MSTRG"):
                    gene = line.split("\t")[0]
                    scaffold = gene_index[gene]
                    if scaffold in scaffold_index:
                        count += 1
                        outfile_reads.write(line)
                        outfile_reads.write("\n")
                else:
                    outfile_reads.write(line)
                    outfile_reads.write("\n")
    print "Number of genes with count data on scaffolds assigned to annotated chromosomes =", count

    #extract positional information for genes on scaffolds assigned to chromosomes
    count = 0
    with open(args.gene_position_outfile, "w") as outfile_start:
        outfile_start.write("Geneid,Scaffold,Chromosome,Startwithinscaffold,Startwithingenome")
        outfile_start.write("\n")
        with open(args.read_counts, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith("MSTRG"):
                    gene = line.split("\t")[0]
                    scaffold = gene_index[gene]
                    if scaffold in scaffold_index:
                        count += 1
                        chromosome = scaffold_index[scaffold]
                        scaffoldstart = scaffold_position[scaffold]
                        genestartwithinscaffold = gene_position[gene]
                        genestartwithingenome = scaffoldstart+genestartwithinscaffold-1
                        outfile_start.write(gene)
                        outfile_start.write(",")
                        outfile_start.write(scaffold)
                        outfile_start.write(",")
                        outfile_start.write(chromosome)
                        outfile_start.write(",")
                        outfile_start.write(str(genestartwithinscaffold))
                        outfile_start.write(",")
                        outfile_start.write(str(genestartwithingenome))
                        outfile_start.write("\n")
    print "Number of genes with count data on scaffolds assigned to annotated chromosomes =", count

if __name__ == '__main__':
    main()

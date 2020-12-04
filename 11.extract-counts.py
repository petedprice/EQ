#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' EXTRACT read counts
Takes a folder of folders containing read count files from HTSEQ-count.
Prints read counts into one file.
'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
from collections import OrderedDict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder of folders containing expression files from HTSEQ-count")
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
def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if is_number(f)]

def list_file(current_dir):
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if not name.endswith(".DS_Store"):
                if name.endswith(".txt"):
                    f = os.path.join(path, name)
                    return f

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    
    #extract gene names with expression
    print "1. Extracting list of all genes expressed ..."
    folders = list_folder(args.infolder)
    print "No. of folders =", len(folders)
    genenames = []
    for folder in folders:
        file = list_file(folder)
        with open(file, "r") as infile:
            for line in infile:
                line = line.rstrip()
                gene = line.split("\t")[0]
                if gene.startswith("MSTRG"):
                    genenames.append(gene)
    genenames = set(genenames)
    print "Total no. of genes with count data across all samples =", len(genenames)

    #check that same genes present in all files (ie mapped against reference gtf)
    #script will not work if this criteria is not met
    for folder in folders:
        file = list_file(folder)
        samplegenenames = []
        with open(file, "r") as infile:
            for line in infile:
                line = line.rstrip()
                gene = line.split("\t")[0]
                if gene.startswith("MSTRG"):
                    samplegenenames.append(gene)
        if len(samplegenenames) != len(genenames):
            print "ERROR - expression not quantified against reference gtf"

    #extract  expression
    geneexpressiondict = defaultdict(list)
    print "2. Extracting expression ..."
    header = []
    for folder in folders:
        file = list_file(folder)
        sample = int(file.split("/")[-2])
        header.append(sample)
        #get expression values for each sample
        with open(file, "r") as infile:
            sampleexpressiondict = {}
            for line in infile:
                if line.startswith("MSTRG"):
                    line = line.rstrip()
                    gene = line.split("\t")[0]
                    counts = float(line.split("\t")[1])
                    sampleexpressiondict[gene] = counts
        print "No. of genes with count data in", sample, "=", len(sampleexpressiondict)
        #attach expression to each gene
        for gene in genenames:
            counts = float(sampleexpressiondict[gene])
            geneexpressiondict[gene].append(counts)
    print "Total no. of genes for which read count info is extracted =", len(geneexpressiondict)

    print "3. Printing all expresson ..."
    with open(args.outfile, "w") as outfile:
        #print header
        outfile.write("Geneid")
        for sample in header:
            outfile.write("\t")
            outfile.write(str(sample))
        outfile.write("\n")
        #print expression
        count = 0
        for gene in geneexpressiondict:
            count += 1
            outfile.write(gene)
            for g in geneexpressiondict[gene]:
                outfile.write("\t")
                outfile.write(str(g))
            outfile.write("\n")
        print "Number of genes printed =", count
        
if __name__ == '__main__':
    main()

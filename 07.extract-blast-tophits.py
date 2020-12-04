#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' EXTRACTS tophits for blast search
Takes a blast output file format (outfmt "6 qseqid sseqid pident length mismatch gapopen qstart 
qend sstart send evalue bitscore sseq") and identifies the top blast hit for each query.
Top blast hit = minimum 30 pidentity, greatest blast score and greatest pidentity.
If a query has two hits with identical blast score and pidentity once is chosen randomly 
as the tophit. Tophit identity is not used for identifying ncrna. 
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
parser.add_argument("infile", type=str,
                    help="A blast result file. Tabular output")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def read_blastout(source):
    try:
        with open(source, "r") as infile:
            count =0
            blasthits = {}
            blastlist = []
            currentquery = None
            for line in infile:
                count +=1
                sys.stdout.write('%d\r' % (count))
                sys.stdout.flush()
                qseqid = line.rstrip().split()[0]
                if currentquery is None:
                    currentquery = qseqid
                    blastlist.append(line)
                else:
                    if currentquery == qseqid:
                        blastlist.append(line)
                    else:
                        blasthits = read_blastout_bitscore(blastlist, blasthits)
                        blastlist = []
                        currentquery = qseqid
                        blastlist.append(line)
            #final loop for final scaffold
            blasthits = read_blastout_bitscore(blastlist, blasthits)

            print "Number of lines =", count
        return blasthits
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % source
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

def read_blastout_bitscore(blastlist, blasthits):
    if len(blastlist) == 0:
        line = line.rstrip().split()
        qseqid = line[0]
        sseqid = line[1]
        pidentity = float(line[2])
        length = int(line[3])
        sstart = float(line[8])
        send = float(line[9])
        evalue = float(line[10])
        bitscore = float(line[11])
        blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
    else:
        countdel = 0
        check = []
        for line in blastlist:
            line = line.rstrip().split()
            qseqid = line[0]
            sseqid = line[1]
            pidentity = float(line[2])
            length = int(line[3])
            sstart = float(line[8])
            send = float(line[9])
            evalue = float(line[10])
            bitscore = float(line[11])
            # Check for min 30% identity
            if pidentity > 30:
                if qseqid in blasthits:
                    if blasthits[qseqid][1] < bitscore:
                        blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
                    elif blasthits[qseqid][1] == bitscore:
                        if blasthits[qseqid][2] < pidentity:
                            blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
                else:
                    blasthits[qseqid] = (sseqid, bitscore, pidentity, sstart, send)
    return blasthits

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    blastoutput = read_blastout(args.infile)
    print "No of genes with top blasthits on scaffolds =", len(blastoutput)

    outfilename = args.infile+".tophits"
    print "Printing results to ..", outfilename
    with open(outfilename, "w") as outfile:
        for blast in blastoutput:
            outfile.write(blast)
            outfile.write(",")
            outfile.write(blastoutput[blast][0])
            outfile.write(",")
            outfile.write(str(blastoutput[blast][1]))
            outfile.write(",")
            outfile.write(str(blastoutput[blast][2]))
            outfile.write(",")
            outfile.write(str(blastoutput[blast][3]))
            outfile.write(",")
            outfile.write(str(blastoutput[blast][4]))
            outfile.write("\n")

if __name__ == '__main__':
    main()

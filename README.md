# Quantification of Gene Expression and Splicing Variation 
**Preperation pipeline**

1. Filter FASTQ files

2. Map reads

**Expression Analysis Pipeline**

1. GTF generation with StringTie

2. Filter transcriptome

3. Quantify gene expression

4. Obtain reciprocal orthologs

5. Analyse expression

**Splicing Variation Pipeline**

# Preperation Pipeline

## 1. Filter fastqs

Pipeline to quality trim FASTQ files. Assumes illumina naming eg WTCHG_243504_002_1 WTCHG_243504_002_2

* **[FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** - A quality control tool for high throughput sequence data

* **python 01.trimmomatic.py**
This script takes a folder containing paired-end fastq.gz files, finds forward and reverse pairs, and quality trims with Trimmomatic. Must run version 0.35 for correct naming scheme.
The script removes all reads shorter than 95nt for useage within rMATs
  
  [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) - A flexible read trimming tool for Illumina NGS data eg
  >java -jar trimmomatic-0.35.jar PE -qualityscore input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:adaptorfile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

* **python 02.max_length95.py**
This script trims reads to 95nt for input to rMATs

## 2. Map reads

Pipeline to map RNA-seq data to a reference genome, identify transcripts and generate a reference transcriptome.

* **Map reads with [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)** - A fast and sensitive alignment program for mapping next-generation sequencing reads.
  
  * **HISAT2 index**
    >hisat2-build -f reference_genome.fa reference_genome_db

  * **python 03.prep-HISAT2.py**
  This script takes a folder containing paired-end fastq files and prepares .sh scripts to run HISAT2. Makes a new folder for the scripts (called wdpath/scripts), and a new folder for each sample (called wdpath/sample). Assumes naming WTCHG_208_1_forward_paired.fastq

    >hisat2 reference_genome -1 forward.fastq -2 reverse.fastq qualityscore -q -p 12 --no-discordant --no-mixed --no-unal --dta - S samoutput --met-file metfile
  
    note that HISAT2, used with the --dta option, is now the recommended aligner to use for StringTie. Every spliced read alignment (i.e. an alignment across at least one junction) in the input SAM file must contain the tag XS to indicate the genomic strand that produced the RNA from which the read was sequenced. Alignments produced by TopHat and HISAT2 (when ran with --dta option) already include this tag, but if you use a different read mapper you should check that this XS tag is included for spliced alignments.
  
    --dta/--downstream-transcriptome-assembly
    Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computationa and memory usage.
  
* **Coordinate sort sam files with [Samtools](http://www.htslib.org)** - A suite of programs for interacting with high-throughput sequencing data.
  
  * **Samtools Flagstat**
  Simple stats on sam output from HISAT2
    >samtools flagstat file.sam

  * **python 04.prep-coordinate-sort-sam.py**
  This script takes an infolder containing folders of sam files and prepares .sh scripts to run samtools sort. Makes a new folder for the scripts.

    >samtools view -Su file.sam | samtools sort - coord_sorted
  
    StringTie takes as input a binary SAM (BAM) file sorted by reference position. This file contains spliced read alignments and can be produced directly by programs such as TopHat or it can be obtained by converting and sorting the output of HISAT2. We recommend using HISAT2 as it is a fast and accurate alignment program. A text file in SAM format which was produced by HISAT2 must be sorted and converted to BAM format using the samtools program. The file resulted from the below command (alns.sorted.bam) can be used as input for StringTie.
    >samtools view -Su alns.sam | samtools sort - alns.sorted

# Expression Analysis Pipeline

## 1. GTF generation with StringTie
* **Extract and merge gene coordinates for each sample with [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml)** - A fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

  * **python 05.prep-StringTie-nogtf.py**
  This script takes an infolder containing folders of bam files and prepares .sh scripts to run StringTie. Makes a new folder for the scripts (called wdpath/scripts), and a new folder for each sample (called wdpath/sample). Produces a GTF file for each sample. The Gene transfer format (GTF) is a file format used to hold information about gene structure and position.
  
    >stringtie coord_sorted.bam -o StringTie_sample.gtf -p 12 -A StringTie_sample.gene_abund

  * **stringtie merge**
  Transcript merge mode.
    
    This is a special usage mode of StringTie, distinct from the assembly usage mode. In the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. This is the transcriptome.

    >stringtie --merge StringTie_gtfs.list -o StringTie_merged.gtf
    
    >-m <min_len>	minimum input transcript length to include in the merge (default: 50)

## 2. Filter transcriptome

Pipeline to filter GTF file to remove non coding RNA

* **Filter transcriptome to remove non coding RNA**
  
   * **Extract transcript sequences with [Bedtools](http://bedtools.readthedocs.io/en/latest/)**
  
      * Index reference genome
      >samtools faidx reference_genome.fa

      * Extract exon sequences
      
      Extracts sequence for each exon separately
      >bedtools getfasta -fi reference_genome.fa -bed StringTie_merged.gtf -fo StringTie_merged.exons.fa
      
      * Concatenate exon sequences to get transcript sequences
      
        **python 06.merge-exon-sequences.py**
        This script takes fasta files and removes replicate sequences. Prints merged fasta file into new file. (!Note: Assumes that if a gene is present in more than one fasta file the genes sequence is the same in all those files)
        
   * **[BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi) to identify and remove non coding RNA**
      
      Download fasta file of non coding RNA from closely related species from [ENSEMBL](http://www.ensembl.org/info/data/ftp/index.html) eg Oryzias_latipes.MEDAKA1.ncrna.fa

      * BLAST index
      >makeblastdb -in Oryzias_latipes.MEDAKA1.ncrna.fa -input_type fasta -dbtype nucl -title MEDAKA1.ncrna_db -out  MEDAKA1.ncrna_db

      * BLAST transcript sequences to fasta file(s) of known ncRNA
      >blastn -evalue 10e-10 -db MEDAKA1.ncrna_db -query StringTie_merged.transcript.fa -out StringTie_merged.transcripts_MEDAKAncrna.blastout -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"
       
        **python 07.extract-blast-tophits.py**
        This script takes a blast output file format (outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq") and identifies the top blast hit for each query. Top blast hit = minimum 30 pidentity, greatest blast score and greatest pidentity. If a query has two hits with identical blast score and pidentity once is chosen randomly as the tophit. Tophit identity is not used for identifying ncrna.  
      
      * Remove ncRNA from reference GTF file
      
        **python 08.filter-GTF-ncrna.py**
        Takes a folder of blast tophit files and filters the GTF to remove these transcripts. In addition, it removes all other transcripts of a given gene from the GTF file if any one transcript maps to ncRNA.

## 3. Quantify gene expression

Pipeline to quantify and filter gene expression

* **Name sort sam files with [Samtools](http://www.htslib.org)** - A suite of programs for interacting with high-throughput sequencing data.
  
  * **python 09.prep-name-sort-sam.py**
  This script takes an infolder containing folders of sam files and prepares .sh scripts to run samtools sort. Makes a new folder for the scripts.

* **Extract read counts with [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)** - Python package that provides infrastructure to process data from high-throughput sequencing assays
  
  For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the samtools sort function of samtools to sort it. Use this option, with name or pos for <order> to indicate how the input data has been sorted. The default is name.

  If name is indicated, htseq-count expects all the alignments for the reads of a given read pair to appear in adjacent records in the input data. For pos, this is not expected; rather, read alignments whose mate alignment have not yet been seen are kept in a buffer in memory until the mate is found. While, strictly speaking, the latter will also work with unsorted data, sorting ensures that most alignment mates appear close to each other in the data and hence the buffer is much less likely to overflow.
  
   * **python 10.prep-htseq-count.py**
  This script takes an infolder containing folders of sam files and prepares .sh scripts to run htseq-count. Makes a new folder for the scripts (called wdpath/scripts), and a new folder for each sample (called wdpath/sample).

    >htseq-count -f sam -r name -s no <alignment_file> <gff_file> 

* **Extract read counts separately for each tissue** - Recommended to extract read counts and normalise expression for each tissue separately. Then the male:female expression ratio can be calculated separately for each tissue and compared across tissues.
    >Advantages: Tissues often have very different expression profiles and this can mean normalise is not very effective.
   
    >Disadvantages: Cannot now compare male or female expression across tissues only the male:female expression ratio

  * **python 11.extract-counts.py**
  This script takes a folder of folders containing read count files from HTSEQ-count. Prints read counts into one file.

  * **OPTIONAL: python 12.extract-counts-annotated.py**
  This script takes a file of read counts across samples extracted from HTseq-count. Extracts read counts for genes on scaffolds
  assigned to chromosomes and prints read counts into one file. Extracts positional information for genes on scaffolds
  assigned to chromosomes and prints positional information into another file.

* **Filter transcriptome to remove lowly expressed genes**
    Filter transcriptome separately for each tissue
  
 * **python 13.extract-gene-length.py**
  This script takes a file of read counts across samples extracted from HTseq-count. Extracts gene lengths and prints 
  into one file.

 * **python 14.convert-counts-to-rpkm.R**
  Converts read count data to RPKM values with edgeR. Prints RPKM data to a new file.

 * **python 15.filter-expression-2rpkm-halformoreofsamples.py**
    This script takes a file of genes and their rpkm values for each sample as output from edgeR. Filters expression. Gene must be expressed > 2FPKM in half or more of the individuals. The script creates a list of genes that have passed the filtering threshold and outputs a file containing genes that have passed the filtering threshold and their rpkm values. Takes a file of read counts across samples extracted from HTseq-count. (Header starts with "Geneid", and each gene name starts with "MSTRG"). Writes file with read counts for all genes in the list that have passed the 2rpkm filtering threshold.

## 4. Obtain reciprocal orthologs

Pipeline to extract gene sequences for expressed genes and obtain reciprocal orthologs. Reduces redundancy in the gene set. 

* **Extract sequences for expressed genes**

  * **python 16.get-longest-isoforms-StringTie.py**
  This script processes an fasta file and picks the longest isoform for each gene. Outputs a new fasta file with ending in _longest.fasta. Assumes naming follows StringTie format eg transcript = MSTRG.22287.1 gene = MSTRG.22287

  * **python 17.extract-expressed-gene-sequences.py**
  Takes a fasta file with the longest transcript sequence for each gene and extracts genes in the read count file. Assumes naming follows StringTie format eg transcript = MSTRG.22287.1 gene = MSTRG.22287. Assumes one transcript sequences per gene.

  * **OPTIONAL: python 18.merge-fasta-files.py**
  Takes fasta files and removes replicate sequences. Prints merged fasta file into new file

* **Identify reciprocal orthologs**

  Download coding sequences of closely related species from Ensembl. For this section naming of all fasta files must be as follows: Speciesname.fa eg Oryziaslaptipes.fa. Therefore, user might need to change the name of the gene sequence fasta file to follow format.

  * **python 19.get-longest-isoforms-Ensembl.py**
  This script processes an Ensembl fasta file and picks the longest isoform for each gene. Outputs a new fasta file ending in _longest.fasta. This script will only work if fasta file is in Ensembl format eg
  >ENSPFOT00000011287.1 cds:known scaffold:PoeFor_5.1.2:KI519728.1:726488:734097:-1 gene:ENSPFOG00000011302.1 gene_biotype:protein_coding transcript_biotype:protein_coding description:Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:A0A087XZR8]

  * **python 20.run-blastall.py**
  This script automatizes a complete reciprocal besthit blast analysis.
  A list of input genomes is blasted against each other and output is written in provided format. Name of outfile is db name.query name.bla. Run with "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq".

    Input genomes are fasta files containing one isoform per gene. Naming of fasta files must be as follows: Speciesname_longest.fasta eg Oryziaslaptipes_longest.fasta. [BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi) must be compiled and added to the path of the user.
    
    >python 1.run-blastall.py -e 10e-10 -b blastn -p 4 -f "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" species1_longest.fasta species2_longest.fasta

  * **python 21.top-blasthit.py**
  This script takes the outputs blast files of a reciprocal blast. 
  First, identifies tophit for each blast. Minimum 30 pidentity, then picks the tophit with greatest bitscore. If bitscores are identical then the blast hit with greatest pidentity is picked as the tophit. If multiple sequenced have the same bitscore and pidentity, the ortholog is discarded.Second, finds 1-1 reciprocal orthologs. The blastout output format supported is:
  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"

  * **python 22.ortho-cluster.py**
  This script takes the pickle output of reciprocal orthologs and identifies clusters of reciprocal orthologs.

  * **python 23.get-expression-reciprocal-orthologs.py**
  This script takes a file containing read counts and extracts expression for reciprocal orthologs. Prints read counts into one file.

## 5. Analyse expression

Pipeline to normalise and analyse gene expression

* **Normalise read counts with [EdgeR](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) in [R](https://www.r-project.org)** - a Bioconductor package for differential expression analysis of digital gene expression data.
    Normalise separately for each tissue. Only read count data can be used as input for EdgeR TMM. Do not use RPKM values.

  * **24.edgeR-normalisation.R**
  Normalises read count data with [TMM](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25), produces graphs to check normalisation was successful, and prints RPKM data to a new file.

# Splicing variation pipeline

%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%

SegSeq RELEASE 1.0.1

Jan 28 2009

http://www.broad.mit.edu/cancer/pub/solexa_copy_numbers/


----------------------
  TABLE OF CONTENTS
----------------------
A. Introduction
B. Usage
C. Input file formats
D. ChangeLog

==========================
    A. INTRODUCTION
==========================

This Matlab code detects copy-number alterations from short sequence
reads obtained from a test (tumor) sample and its matched normal.  The
algorithm begins with a local changepoint analysis to find candidate
copy-number breakpoints, based on a log-ratio difference statistic.  An
iterative merging procedure that removes the least significant
breakpoint, until a stopping p-value is reached.

For details, please refer to and cite the following publication:
Chiang DY, Getz G et al (2009)  Nature Methods 6: 99-103


==========================
       B. USAGE
==========================

B1. SYNTAX
-------------

->  SegSeq -i sample_info_file -s sample_name [Optional params]

In this syntax, sequence read positions are first mapped to the alignable
portion of the hg18 human genome build (calculated with 36 bp reads and
the ARACHNE aligner).  These alignable positions are saved in a Matlab
file (named [sample_name]_aligned_reads.mat) before segmentation is called
with default parameters (W=400,a=1000,b=10).

->  SegSeq -t tumor_aligned_reads.mat -s sample_name [Optional params]

This syntax can be used to assess the effects of different input parameters
(W,a,b) on segmentation results.  Data from the tumor_aligned_reads.mat file
is loaded for the segmentation procedure only


B2. EXAMPLES
--------------

SegSeq -i H2347_BL2347_solexa_sampleinfo.txt -s H2347

SegSeq -t H2347_aligned_reads.mat -s H2347 -W 400 -a 100 -b 10


B3. OPTIONAL PARAMETERS
-------------------------
[SEGMENTATION PARAMETERS]
-W  Size of local windows (i.e., # of consecutive normal reads)
-a  Number of false positive candidate breakpoints for initialization
-b  Number of false positive segments for termination

The default parameter settings are: W = 400, a = 1000, b = 10

[SEQUENCE READ INFORMATION]
-n  Matlab file with positions of aligned sequence reads in normal replicates
-t  Matlab file with positions of aligned sequence reads in tumor replicates

[COPY-NUMBER RATIO CALCUATIONS]
-c  Chromosome information file (see section C3 below)
-e  Coordinates of non-overlapping genome windows
-r  Flag to calculate copy-number ratios in the genome windows
    specified by -e
-f  Flag to set the median of calculated copy-number ratios to 1
-m  Minimum number of counts in the normal sample in each window
-g  Coordinates of excluded regions, i.e. germline copy-number
variants
-v  Flag to remove excluded regions



==========================
  C. INPUT FILE FORMATS
==========================

C1. SAMPLE INFORMATION FILE
---------------------------------
<Example> HCC1143_BL1143_5550_solexa_sampleinfo.txt

This tab-delimited file denotes the tumor/normal status for each file
of aligned read positions (below).  It has the following format:
[Column 1] Filename with read start positions for a single sequencing
lane
[Column 2] Sample name
[Column 3] Tumor/Normal status.  Please note that "Tumor" and "Normal"
are case-sensitive

Other optional columns may provide additional annotation.
 

C2. ALIGNED READ POSITIONS
---------------------------------
<Example> arachne/5827.8.qltout.txt

This tab-delimited or space-delimited file contains the start positions
of UNIQUELY aligning reads.  Each line has the format:
[Chr] [Pos] [Strand]

These alignments were obtained from short read processing by the ARACHNE
software suite.  The original sequence reads in FASTQ file format can
be retrieved from the NCBI Short Read Archive, with project accession
number SRP000246.

The FASTQ files are annotated with {flowcell.lane}, i.e. 4343:2 in the
SRR002793.fastq file.  The identity of tumor and normal cell lines are
found in the sample information files (*sampleinfo.txt)



C3. CHROMOSOME INFORMATION FILE
---------------------------------
<Example> chromInfo_hg18.txt

This tab-delimited file contains the physical and alignable length
of each chromosome.
[Column 1] Chromosome number
[Column 2] Length in physical coordinates
[Column 3] Length in alignable coordinates




==========================
  D. ChangeLog
==========================

version	1.0.1   2009-01-28   Added documentation in README.txt

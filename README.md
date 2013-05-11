
==========================================================
  aLib: a sequencing pipeline for ancient and modern DNA
==========================================================

aLib is:
- an open-source set of modules to process ancient and modern DNA data from raw intensities to usable BAM files
- It is aimed at mid-size sequencing centers processing ancient DNA using Illumina sequencers
- It performs basic processing that should be common to any sequencing pipeline.
- It features:
 - Online monitoring interface for sequencing runs
 - Online form with analysis requirements 
 - Basecalling either with freeIbis or using the default basecaller Bustard
 - Trimming of sequencing adapters, merging of partially overlapping mate-pairs and flagging of chimeric sequences
 - Demultiplexing based on read group assignment quality 
 - Basic filtering
 - Basic Quality control



To install:

Requirements:
- Bamtools (works fine with 2.2.2) 
- libgab (https://github.com/grenaud/libgab)
- freeIbis (https://github.com/grenaud/freeIbis)
- fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- 

Replace the path in the makefiles with yours and take make. Sorry for the lack of elegance, it will be corrected soon.

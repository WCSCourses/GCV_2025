#!/bin/sh

## Performs normalisation of FASTQs for the GCV WCSC

interleave-reads.py -o interleaved.fastq $1 $2  
normalize-by-median.py -pk31 -C5 -o normalised.fastq interleaved.fastq 
split-paired-reads.py -1 normalised_R1.fastq -2 normalised_R2.fastq normalised.fastq
rm interleaved.fastq 
rm normalised.fastq

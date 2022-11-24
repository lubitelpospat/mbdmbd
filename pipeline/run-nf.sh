#!/usr/bin/bash

genome='/storage/liam/B2_UR.fa.filtlength-125'
sampleA='/storage/liam/re-basecall/p48_anti-adar_guppy6.3/' 
sampleB='/storage/liam/re-basecall/p48_anti-adar_guppy6.3/'

## A
nextflow m6anet-pipeline.nf -resume --mode ncrna --genome $genome  --fastqDirectory $sampleA"fastq_pass" --fast5Directory $sampleA"fast5_pass"

## B
#nextflow m6anet-pipeline.nf -resume --mode ncrna --genome $genome  --fastqDirectory $sampleB"fastq_pass" --fast5Directory $sampleB"fast5_pass"

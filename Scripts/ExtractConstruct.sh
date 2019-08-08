#!/bin/bash -l

#####################
# This script extracts mate pairs that maps to construct and host genome. 
#####################

# $1 = bam 
# $2 = construct_name

echo 'extracts reads from bamfile mapped to construct'
echo 'module load samtools/0.1.19'
module load samtools/0.1.19

# Take only reads mapped to construct and then save the firt 10 fields 
echo 'samtools view PPM1D_M47_BWA_sorted.bam.bam -region $2 | cut -f 1-10 > ${2}.sam' 
samtools view $1 -region $2 | awk '$7 != "="' | cut -f 1-10 > ${2}.sam

echo 'Complete!'


#!/bin/bash -l

#####################
# This script extracts mate pairs that maps to construct and host genome. 
#####################


#$ -S /bin/bash 
#$ -N runConstructExtraction 
#$ -j y
#$ -M vanja.borjesson@gu.se
#$ -m bea
#$ -cwd
#$ -pe mpi 8
#$ -q bfxcore.q@node3-bfx.medair.lcl

bam = $1
construct_name = $2

echo 'extracts reads from bamfile mapped to construct'
echo 'module load samtools/0.1.19'
module load samtools/0.1.19

# Take only reads mapped to construct and then save the firt 10 fields 
echo 'samtools view PPM1D_M47_BWA_sorted.bam.bam -region Wip1 > Wip.sam' 
samtools view PPM1D_M47_BWA_sorted.bam.bam -region $2 | cut -f 1-10 > ${2}.sam

echo 'Complete!'


#!/bin/bash -l

#$ -S /bin/bash 
#$ -N SoftClipEx 
#$ -j y
#$ -M vanja.borjesson@gu.se
#$ -m bea
#$ -cwd
#$ -pe mpi 8
#$ -q bfxcore.q@node3-bfx.medair.lcl

module load samtools/0.1.19

#echo 'samtools view $1 |awk '$6 ~ /S/' |awk '($3 ~ /Wip1/ && $15 !~ /Wip/) || ($3 !~ /Wip1/ && $15 ~ /Wip1/)' |awk '$15 ~ /SA/' | samtools view -Sb - > PPM1D_M47_softclipped_Wip1.bam'
samtools view $1 |awk '$6 ~ /S/' |awk '($3 ~ /Wip1/ && $15 !~ /Wip/) || ($3 !~ /Wip1/ && $15 ~ /Wip1/)' |awk '$15 ~ /SA/' > $2

echo 'Complete!'
 



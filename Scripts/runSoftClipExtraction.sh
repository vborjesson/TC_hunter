#!/bin/bash -l


module load samtools/0.1.19



#echo 'samtools view $1 |awk '$6 ~ /S/' |awk '($3 ~ /Wip1/ && $15 !~ /Wip/) || ($3 !~ /Wip1/ && $15 ~ /Wip1/)' |awk '$15 ~ /SA/' | samtools view -Sb - > PPM1D_M47_softclipped_Wip1.bam'
samtools view $1 |awk '$6 ~ /S/' |awk '($3 ~ /Wip1/ && $15 !~ /Wip/) || ($3 !~ /Wip1/ && $15 ~ /Wip1/)' |awk '$15 ~ /SA/' > $2

echo 'Complete!'
 



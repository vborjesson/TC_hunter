#!/bin/bash -l


module load samtools/0.1.19


echo $3
c=$3
echo $c
echo $1
echo $2

sa_name=SA:Z:${3}
echo $sa_name
#echo "'samtools view $1 $c |awk '$6 ~ /S/' |awk '($3 ~ /$c/ && $15 !~ /$c/) || ($3 !~ /$c/ && $15 ~ /$c/)' |awk '$15 ~ /SA/' > $2'"

#samtools view $1 $c |awk '$6 ~ /S/' |awk -v c=$3 '($3 ~ /c/ && $15 !~ /c/) || ($3 !~ /c/ && $15 ~ /c/)' |awk '$16 ~ /SA/ || $15 ~ /SA/ || $14 ~ /SA/' > $2

#samtools view /jumbo/WorkingDir/B19-001/Intermediate/M4_all/work/8a/04f83696c83f363c5b088feaa9a67b/PPM1D_M42_indexed.bam |awk '$6 ~ /S/' |awk '($3 ~ /T_scan_out/ && $15 !~ /T_scan_out/) || ($3 !~ /T_scan_out/ && $15 ~ /T_scan_out/)' |awk '$15 ~ /SA/' > test.sam
#samtools view $1 |awk '$6 ~ /S/' |awk '($3 ~ /$c/ && $15 !~ /$c/) || ($3 !~ /$c/ && $15 ~ /$c/)' |awk '$15 ~ /SA/' > $2

samtools view $1 $3| awk '$6 ~ /S/' | grep SA: |grep -v $sa_name > $2

echo 'Complete!'
 


 
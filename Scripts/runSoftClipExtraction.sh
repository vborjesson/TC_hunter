#!/bin/bash -l


module load samtools/0.1.19


echo $3
c=$3
echo $c
echo $1
echo $2

sa_name=SA:Z:${3}
echo $sa_name

samtools view $1 $3| awk '$6 ~ /S/' | grep SA: |grep -v $sa_name > $2

echo 'Complete!'
 


 
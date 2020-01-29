#!/bin/bash -l

#$ -S /bin/bash
#$ -N seqtk
#$ -j y
#$ -M vanja.borjesson@gu.se
#$ -m bea
#$ -cwd
#$ -pe mpi 16
#$ -q bfxcore.q@node6-bfx.medair.lcl

seqtk sample -s100 $1 225140734 > M41_div2.fq  
seqtk sample -s100 $2 225140734 > M41_div2.fq 

seqtk sample -s100 $1 112570367 > M41_div2.2.fq
seqtk sample -s100 $2 112570367 > M41_div2.2.fq

seqtk sample -s100 $1 56285184 > M41_div2.2.2.fq
seqtk sample -s100 $2 56285184 > M41_div2.2.2.fq




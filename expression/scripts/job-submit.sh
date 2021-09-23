#!/bin/bash
nam=$1
cd $nam
cp ../hisat2.sub ./
ls *_R1.fq.gz > input.fofn
jobs=$(cat input.fofn |wc -l |awk '{print $1-1}')
sed -i 's/job-name=hisat2/job-name='$nam'/g' hisat2.sub
sbatch --array=0-${jobs} hisat2.sub

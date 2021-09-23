#!/bin/bash
set -x
nam=$1
cd $nam
jobs=$(echo $(cat input.fofn |wc -l) -1 |bc)
sbatch --array=0-${jobs} hisat2-self.sub

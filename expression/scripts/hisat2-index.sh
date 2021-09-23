#!/bin/bash
genome=$1
dir="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/f-rnaseq-bowtie2"
ml hisat2
nam=$(basename $genome |cut -f 1 -d ".")
mkdir -p ${dir}/${nam}/index
cd ${dir}/${nam}/index
hisat2-build -p 36 ${genome} ${nam}

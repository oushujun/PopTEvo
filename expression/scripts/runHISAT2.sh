#!/bin/bash
module load hisat2
nam=$1
# 8DAS_root_MN01011_R1.fq.gz
index="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/f-rnaseq-bowtie2/${nam}/index/${nam}"
read1=$2
read2=$(echo $read1 |sed 's/_R1.fq.gz/_R2.fq.gz/g')
out=$(basename ${read1} | sed 's/_R1.fq.gz/_'${nam}'/g')
hisat2 \
   --threads 9 \
   -k 20 \
   --summary-file ${out}_summary.txt \
   -S ${out}.sam \
   -x ${index} \
   -1 $read1 \
   -2 $read2 || {
echo >&2 hisat2 failed for $out
exit 1
}
SAM=${out}.sam
ml samtools
samtools view --threads 9 -b -o ${SAM%.*}.bam ${SAM} || {
echo >&2 bam conversion failed for $out
exit 1
}
samtools sort --threads 9 -o ${SAM%.*}_sorted.bam -T ${SAM%.*}_temp ${SAM%.*}.bam || {
echo >&2 bam sorting failed for $out
exit 1
}
samtools stats --threads 9 ${SAM%.*}_sorted.bam > ${SAM%.*}_sorted.stats || {
echo >&2 bam stats failed for $out
exit 1
}
rm ${SAM%.*}.bam ${SAM}

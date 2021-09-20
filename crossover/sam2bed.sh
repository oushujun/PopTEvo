#!/bin/bash
ml samtools
ml bedtools2
sam=$1
samtools view -b -o ${sam%.*}.bam ${sam}
samtools sort -o ${sam%.*}_sorted.bam ${sam%.*}.bam
bedtools bamtobed -i ${sam%.*}_sorted.bam > ${sam%.*}.bed

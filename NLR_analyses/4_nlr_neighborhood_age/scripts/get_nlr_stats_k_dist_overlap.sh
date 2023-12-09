#!/usr/bin/env bash

NAMLINE=$1
K=$2
BGD_FILE=$3
FGD_FILE=$4
NAME=$5

# BACKGROUND
bedtools closest -d -k ${K} -a <(grep -v scaf ${BGD_FILE} |sort -k1,1 -k2,2n) -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed|sort -k1,1 -k2,2n) |cut -f1-3,47|  sort -k1,1 -k2,2n -k4,4nr| awk '!seen[$1,$2,$3]++'| awk '{sum+=$4}END{print sum/NR}' | sed "s/^/${NAMLINE}\t/" | sed 's/ /\t/g'| sed "s/$/\tBACKGROUND\t${K}/" >> OUT_${NAME}_${K}_OVERLAP_DIST

# FOREGROUND
bedtools closest -d -k ${K} -a <(grep -v scaf ${FGD_FILE}|sort -k1,1 -k2,2n|awk '$4>1') -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed| sort -k1,1 -k2,2n)| cut -f1-3,47|  sort -k1,1 -k2,2n -k4,4nr| awk '!seen[$1,$2,$3]++'| awk '{sum+=$4}END{print sum/NR}' | sed "s/^/${NAMLINE}\t/" | sed 's/ /\t/g'| sed "s/$/\tNLR\t${K}/" >> OUT_${NAME}_${K}_OVERLAP_DIST


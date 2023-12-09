#!/usr/bin/env bash

NAMLINE=$1
K=$2
BGD_FILE=$3
FGD_FILE=$4
NAME=$5

bedtools closest -k $K -a <(grep -v scaf ${BGD_FILE} |sort -k1,1 -k2,2n) -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed | sort -k1,1 -k2,2n) |sed 's/ /_/g' | cut -f10| sort | uniq -c| sed 's/ *//'| tr ' ' '\t'|sed "s/^/${NAMLINE}\t/"|sed 's/$/\tBACKGROUND/'| sed "s/^/OUT_ALL_ENRICH_100kb_${K}_OVERLAPFAM\t/" >> OUT_ALL_${NAME}_${K}_OVERLAPFAM

bedtools closest -k $K -a <(grep -v scaf ${FGD_FILE} |sort -k1,1 -k2,2n) -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed | sort -k1,1 -k2,2n) |sed 's/ /_/g' | cut -f10| sort | uniq -c| sed 's/ *//'| tr ' ' '\t'|sed "s/^/${NAMLINE}\t/"|sed 's/$/\tNLR/'|sed "s/^/OUT_ALL_ENRICH_100kb_${K}_OVERLAPFAM\t/" >> OUT_ALL_${NAME}_${K}_OVERLAPFAM


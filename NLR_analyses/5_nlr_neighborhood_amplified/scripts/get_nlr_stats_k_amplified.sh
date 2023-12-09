#!/usr/bin/env bash

NAMLINE=$1
K=$2
BGD_FILE=$3
FGD_FILE=$4
NAME=$5

# CALC STATS
bedtools closest -k ${K} -a <(grep -v scaf ${BGD_FILE} |sort -k1,1 -k2,2n) -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed|sort -k1,1 -k2,2n) |cut -f29|  sed 's/ /_/g'|sort| uniq -c| sed 's/^ *//'| tr -s ' '| sed "s/^/${NAMLINE}\t/" | sed 's/ /\t/g'| sed 's/$/\tBACKGROUND/' >> OUT_AMP_${NAME}_${K}

# FOREGROUND
bedtools closest -k ${K} -a <(grep -v scaf ${FGD_FILE}|sort -k1,1 -k2,2n|awk '$4>1') -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed | sort -k1,1 -k2,2n)| cut -f29|sed 's/ /_/g'|sort| uniq -c| sed 's/^ *//'|tr -s ' '| sed "s/^/${NAMLINE}\t/" | sed 's/ /\t/g'| sed 's/$/\tNLR/' >> OUT_AMP_${NAME}_${K}

#
BGD_NOTAMP=`awk '$4=="BACKGROUND" && $3!="Tropical_amplification"' OUT_AMP_${NAME}_${K}|grep $NAMLINE|sed 's/ /_/g'|awk '{sum += $2}END{print sum}'`
BGD_AMP=`awk '$4=="BACKGROUND" && $3=="Tropical_amplification"' OUT_AMP_${NAME}_${K}|grep $NAMLINE|sed 's/ /_/g'|awk '{sum += $2}END{print sum}'`
BGD_ANA=`echo "$BGD_AMP / $BGD_NOTAMP"| bc -l`



FGD_NOTAMP=`awk '$4=="NLR" && $3!="Tropical_amplification"' OUT_AMP_${NAME}_${K}|grep $NAMLINE|sed 's/ /_/g'| awk '{sum += $2}END{print sum}'`
FGD_AMP=`awk '$4=="NLR" && $3=="Tropical_amplification"' OUT_AMP_${NAME}_${K}|grep $NAMLINE|sed 's/ /_/g'| awk '{sum += $2}END{print sum}'`
FGD_ANA=`echo "$FGD_AMP / $FGD_NOTAMP"| bc -l`

echo "$NAME $NAMLINE $K $BGD_ANA $FGD_ANA $BGD_NOTAMP $BGD_AMP $FGD_NOTAMP $FGD_AMP"


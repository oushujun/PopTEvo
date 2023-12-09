#!/usr/bin/env bash

NAMLINE=$1
K=$2
BGD_FILE=$3
FGD_FILE=$4
NAME=$5


# BACKGROUND
bedtools closest -k ${K} -a <(grep -v scaf ${BGD_FILE} |sort -k1,1 -k2,2n) -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed|sort -k1,1 -k2,2n) |cut -f27|  sort| uniq -c| sed 's/^ *//'| tr -s ' '| sed "s/^/${NAMLINE}\t/" | sed 's/ /\t/g'| sed 's/$/\tBACKGROUND/' >> OUT_${NAME}_${K}

# FOREGROUND
bedtools closest -k ${K} -a <(grep -v scaf ${FGD_FILE}|sort -k1,1 -k2,2n|awk '$4>1') -b <(grep -v '^s' ../2_ltr/data/${NAMLINE}_ltr.bed | sort -k1,1 -k2,2n)| cut -f27|sort| uniq -c| sed 's/^ *//'|tr -s ' '| sed "s/^/${NAMLINE}\t/" | sed 's/ /\t/g'| sed 's/$/\tNLR/' >> OUT_${NAME}_${K}

BGD_NOTYOUNG=`awk '$4=="BACKGROUND" && $3!="Young"' OUT_${NAME}_${K}| grep $NAMLINE| awk '{sum += $2}END{print sum}'`
BGD_YOUNG=`awk '$4=="BACKGROUND" && $3=="Young"' OUT_${NAME}_${K}|grep $NAMLINE|awk '{sum += $2}END{print sum}'`
BGD_OLD=`awk '$4=="BACKGROUND" && $3=="Old"' OUT_${NAME}_${K}|grep $NAMLINE|awk '{sum += $2}END{print sum}'`
BGD_YO=`echo "$BGD_OLD / $BGD_YOUNG"| bc -l`

FGD_NOTYOUNG=`awk '$4=="NLR" && $3!="Young"' OUT_${NAME}_${K}|grep $NAMLINE|awk '{sum += $2}END{print sum}'`
FGD_YOUNG=`awk '$4=="NLR" && $3=="Young"' OUT_${NAME}_${K}|grep $NAMLINE|awk '{sum += $2}END{print sum}'`
FGD_OLD=`awk '$4=="NLR" && $3=="Old"' OUT_${NAME}_${K}|grep $NAMLINE|awk '{sum += $2}END{print sum}'`
FGD_YO=`echo "$FGD_OLD / $FGD_YOUNG"| bc -l`

echo "$NAME $K $NAMLINE $BGD_YO $FGD_YO $BGD_YOUNG $BGD_OLD $FGD_YOUNG $FGD_OLD"

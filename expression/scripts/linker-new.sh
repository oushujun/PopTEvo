#!/bin/bash
for fq in /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/f-rnaseq/raw_data/RNASeq_plate_1-6/*fq.gz; do
nam=$(basename ${fq} | cut -f 1 -d "_")
ofq=$(basename ${fq} | cut -f 2- -d "_")
mkdir -p ${nam};
ln -s $fq $nam/${ofq};
done
# post-renaming for consistency
mv KI11 Ki11
mv Mo18w Mo18W
mv Ms71 MS71
mv TX303 Tx303
mv TZI-8 Tzi8
# fix the name
mv CML227 CML277

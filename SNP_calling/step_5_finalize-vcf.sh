#!/bin/bash
#mkdir vcf-files
#mv *.vcf ./vcf-files/
#rm *.idx
#cd vcf-files
for f in *.vcf; do grep -v "^#" $f > $f.1; done
grep "^#" mexicana.v1_10_34000001-36000000.vcf > header
cat *.vcf.1 > body
cat header body > ../merged_mexicana-v1.vcf
cd ..
ml vcftools
ml bcftools
cat merged_mexicana-v1.vcf | vcf-sort -t $TMPDIR -p 36 -c > merged_mexicana-v1_sorted.vcf
bcftools stats merged_mexicana-v1_sorted.vcf > merged_mexicana-v1.vchk
plot-vcfstats merged_mexicana-v1.vchk -p plots/
merged=merged_mexicana-v1
maxdepth=$(grep -oh ";DP=.*;" merged_mexicana-v1_sorted.vcf | cut -d ";" -f 2 | cut -d "="  -f 2 | datamash mean 1 sstdev 1 | awk '{print $1+$2*5}')
echo "$maxdepth" > maxdepth.txt
# maxdepth=1110.97
vcftools --vcf ${merged}_sorted.vcf --keep-only-indels --recode --recode-INFO-all --out ${merged}_sorted-indels.vcf
vcftools --vcf ${merged}_sorted.vcf --remove-indels --recode --recode-INFO-all --out ${merged}_sorted-snps.vcf 
ref="/ptmp/LAS/arnstrm/shujun/mexicana/zmexicana.fa"
gatk --java-options "-Xmx80g -XX:+UseParallelGC" VariantFiltration     --reference $ref     --variant ${merged}_sorted-snps.vcf.recode.vcf     --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 45.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > ${maxdepth}"     --filter-name "FAIL"     --output ${merged}_filtered-snps.vcf
vcftools --vcf ${merged}_filtered-snps.vcf --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out merged-snps_pass-only_maxmis-0.5_maf-3_minq-30
vcftools --vcf ${merged}_filtered-snps.vcf --missing-indv

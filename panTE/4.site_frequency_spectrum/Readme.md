
## process synLTRs


## convert synLTR table to VCF and filter

```bash
# convert to genotypes
perl -nle 's/pan_TE\s+/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT /; print $_ and next if /#CHROM/; $pos++; s/\s+\S+_(intact|truncated)/ 1\/1/g; s/\s+\S+_null/ 0\/0/g; s/\s+NA/ .\/./g; my ($id, $geno) = (split /\s+/, $_, 2); next unless $geno =~ /\//; print "1 $pos $id A C 30 PASS AA= GT $geno"' <(cut -f 1-24 synLTR.txt) > synLTR.geno

# make vcf
cat head.vcf synLTR.geno | perl -nle 's/\s+/\t/g unless /^##/; print $_' > synLTR.vcf

# discard loci with missing ≥ 50%
cat synLTR.vcf | vcf-subset -e -c B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,HP301,P39,Il14H | bcftools view -i 'F_MISSING<0.5' - > synLTR_keep.temp.vcf
cat synLTR.vcf | vcf-subset -e -c CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8 | bcftools view -i 'F_MISSING<0.5' - > synLTR_keep.trop.vcf

# convert missing to ref (assuming no LTR insertion)
perl -nle 's/\.\/\./0\/0/g; print $_' synLTR.vcf > synLTR.convert.vcf
cat synLTR.convert.vcf | vcf-subset -e -c B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,HP301,P39,Il14H > synLTR.convert.temp.vcf
cat synLTR.convert.vcf | vcf-subset -e -c CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8 > synLTR.convert.trop.vcf
```

## calculate synLTR SFS

```
## all synLTR
# unfolded using age
cat <(grep '#' synLTR_temp.vcf) <(perl ~/las/git_bin/EDTA/util/output_by_list.pl 3 synLTR.convert.temp.vcf 1 <(perl -nle 'my @info = (split /\t/, $_); print "$info[0]\t$info[32]" if $info[32] < 20' synLTR.txt))| ~/las/git_bin/SoFoS/sofos -n 10 -u -r -a 1.0 -b 1.0 - > synLTR.convert.temp.post20k.sfs
cat <(grep '#' synLTR_trop.vcf) <(perl ~/las/git_bin/EDTA/util/output_by_list.pl 3 synLTR.convert.trop.vcf 1 <(perl -nle 'my @info = (split /\t/, $_); print "$info[0]\t$info[32]" if $info[32] < 20' synLTR.txt))| ~/las/git_bin/SoFoS/sofos -n 10 -u -r -a 1.0 -b 1.0 - > synLTR.convert.trop.post20k.sfs

# folded
for i in synLTR.convert.*.vcf; do ~/las/git_bin/SoFoS/sofos -n 10 -f -r -a 1.0 -b 1.0 $i > ${i%.*}.folded.sfs & done


## filtered synLTR (missing data rate < 50%)
# unfolded using mexicana
cat <(grep '#' synLTR_keep_temp.vcf) <(perl ~/las/git_bin/EDTA/util/output_by_list.pl 3 synLTR_keep_temp.vcf 1 <(perl -nle 'my @info = (split /\t/, $_); print "$info[0]\t$info[41]" if $info[41] eq "yes"' synLTR_keep.txt))| ~/las/git_bin/SoFoS/sofos -n 10 -u -r -a 1.0 -b 1.0 - > synLTR_keep_temp.mexnull.sfs
cat <(grep '#' synLTR_keep_trop.vcf) <(perl ~/las/git_bin/EDTA/util/output_by_list.pl 3 synLTR_keep_trop.vcf 1 <(perl -nle 'my @info = (split /\t/, $_); print "$info[0]\t$info[41]" if $info[41] eq "yes"' synLTR_keep.txt))| ~/las/git_bin/SoFoS/sofos -n 10 -u -r -a 1.0 -b 1.0 - > synLTR_keep_trop.mexnull.sfs

# folded
~/las/git_bin/SoFoS/sofos -n 10 -f -r -a 1.0 -b 1.0 synLTR_keep_temp.vcf > synLTR_keep_temp.folded.sfs
~/las/git_bin/SoFoS/sofos -n 10 -f -r -a 1.0 -b 1.0 synLTR_keep_trop.vcf > synLTR_keep_trop.folded.sfs
```


## filter SNPs

```bash
# remove SNPs from scaffolds; remove heterozygous, remove multi-allelic sites
cat <(perl -nle 's/IL14H/Il14H/gi; s/MS37W/M37W/gi; s/OH43/Oh43/gi; s/OH7B/Oh7B/gi; s/TX303/Tx303/gi; s/TZi8/Tzi8/gi; print $_' mexicana.header.vcf) \
    <(grep -v -P 'scaf|0/1|1/0' merged_mexicana-v1_filtered-snps.pass.vcf | awk '{if($5!~/,/) print $0}') > merged_mexicana-v1_filtered-snps.biallelic.homo.vcf

# discard loci with missing ≥ 50%
cat merged_mexicana-v1_filtered-snps.biallelic.homo.vcf | \
	vcf-subset -e -c B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,HP301,P39,Il14H | \
	bcftools view -i 'F_MISSING<0.5' - > \
	merged_mexicana-v1_filtered-snps.biallelic.homo.miss50.temp.vcf &

cat merged_mexicana-v1_filtered-snps.biallelic.homo.vcf | \
	vcf-subset -e -c CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8 | \
	bcftools view -i 'F_MISSING<0.5' - > \
	merged_mexicana-v1_filtered-snps.biallelic.homo.miss50.trop.vcf &

# no missing data filter, convert missing to ref (assuming no SNP)
perl -nle 's/:\S+//g; s/\.\/\./0\/0/g; print $_' merged_mexicana-v1_filtered-snps.biallelic.homo.vcf | \
        vcf-subset -e -c B73,B97,Ky21,M162W,Ms71,Oh43,Oh7B,HP301,P39,Il14H > \
	merged_mexicana-v1_filtered-snps.biallelic.homo.convert.temp.vcf &

perl -nle 's/:\S+//g; s/\.\/\./0\/0/g; print $_' merged_mexicana-v1_filtered-snps.biallelic.homo.vcf | \
	 vcf-subset -e -c CML52,CML69,CML103,CML228,CML247,CML277,CML322,CML333,Ki3,Ki11,NC350,NC358,Tzi8 | \
	merged_mexicana-v1_filtered-snps.biallelic.homo.convert.trop.vcf &

```


## calculate folded SNP SFS

```
## all SNPs
for i in merged_mexicana-v1_filtered-snps.biallelic.homo.convert.*.vcf; do \
	sed 's/AC=/AA=/' $i | ~/las/git_bin/SoFoS/sofos -n 10 -f -r -a 1.0 -b 1.0 > ${i%.*}.folded.sfs & \
done

## filtered SNPs
for i in merged_mexicana-v1_filtered-snps.biallelic.homo.miss50.*.vcf ; do \
	sed 's/AC=/AA=/' $i | ~/las/git_bin/SoFoS/sofos -n 10 -f -r -a 1.0 -b 1.0 > ${i%.*}.folded.sfs & \
done
```



## Data download

Wallace et al., 2014 [1] published nearly 7000 GBS markers for the RILs and the data was available through panzea.org (hosted on CyVerse):

```bash
/iplant/home/shared/panzea/genotypes/GBS/v27
```

From that source, we obtained the metadata for all the individuals and the SNPs file (VCF format) for constructing new genetic maps. The files are:

```bash
ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf.gz
AllZeaGBSv2.7_publicSamples_metadata20140411.xlsx
```

## Data cleanup

From the excel sheet, all the individuals related to NAM project were filtered using the Excel `Filter` function on `Project` column and the names were separated to different files using  the `Population` column using bash scripting (see below). The individuals marked as F1 in `Pedigree` column, Blank in the `DNASample` column were excluded.  The parents were identified using the `Population` column, containing `282 maize association mapping panel`.  Since there were 25 replicates for B73, all but one was retained (`B73(PI550473):250027745`). After filtering, [`subset.txt`](assets/subset.txt) with just `FullName` and `Population` columns was retained and processed as follows:

```bash
for namcross in $(cut -f 2 subset.txt |sort |uniq ); do
grep -wE "$namcross$" subset.txt >> ${namcross}.subset.txt;
done
mv Population.subset.txt header.txt
# verify if all data has been split successfully
wc -l subset.txt
wc -l *.subset.txt
# the sum of split files matches exactly to the lines of subset.txt
# parents are in different file and needs to be separated
less 282_maize_association_mapping_panel.subset.txt
for parent in $(cut -f 1 -d "(" 282_maize_association_mapping_panel.subset.txt | sort | uniq); do
grep -E "^$parent" 282_maize_association_mapping_panel.subset.txt >> ${parent}.parents.txt;
done
# also to match the names in the excel file to the vcf file, we had to trim the names a little bit
# from B73(PI550473):62P7LAAXX:4:250027872 in excel to B73(PI550473):250027872 in vcf file
for namcross in *.subset.txt; do
  cut -f 1 $namcross | cut -f 1,4- -d ":" > ${namcross%%.*}.newsubset.txt;
done
# and for parents:
for parents in *.parents.txt; do
  cut -f 1 $parents |cut -f 1,4- -d ":"  >> ${f%%.*}.newparents.txt;
done
# merge files
for nam in *.newsubset.txt; do
  parent=$(echo $nam |sed 's/B73x//g' |sed 's/.newsubset.txt/.newparents.txt/g');
  cat B73.newparents.txt $parent $nam >> ${nam}.full.txt;
done
rename .newsubset.txt.full.txt _full.txt *.newsubset.txt.full.txt
# since we no longer need the second column in these files, we will remove them
for nam in *_full.txt; do
  cut -f 1 $nam > $nam.1;
  mv $nam.1 $nam
done
# names.txt file is the individual names grabbed from the vcf file):
grep "^#CHROM" ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf |cut -f 10- > names.txt
# sanity check
for namcross in *.newsubset.txt; do
  vcfnames=$(grep -wF -f $namcross ../name_filtering/names.txt |wc -l);
  excelnames=$(cat $namcross |wc -l);
  echo -e "$namcross\t$vcfnames\t$excelnames";
done
```

Splitting the VCF file for each NAM cross

```
for namcross in *_full.txt; do
bcftools view \
   --threads 12 \
   --output-type z \
   --output-file ${namcross%.*}.vcf.gz \
   --samples-file $namcross ../ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf.gz;
done
```

Filtering the VCF files

Some filtering is necessary to: (A) reduce the missingness of data, (B) make processing files faster by reducning noise. We will use BCFTools for this purpose. The steps are as follows:

```bash
# retain bi-alleic only SNPs and remove any SNPs that are missing in >50% individuals
module load bcftools
for vcf in *.vcf.gz; do
  vcftools \
     --gzvcf $vcf \
     --min-alleles 2 \
     --max-alleles 2 \
     --max-missing 0.5 \
     --recode --recode-INFO-all \
     --out ${vcf%%.*}_biallelic_only_maxmissing_0.5
  vcftools \
     --vcf ${vcf%%.*}_biallelic_only_maxmissing_0.5.recode.vcf \
     --missing-indv \
     --out ${vcf%%.*}_biallelic_only_maxmissing_0.5
done
```
Clean-up names for the individuals as they have parenthesis and hyphens that are problematic in R

```bash
#!/bin/bash
vcf="$1"
grep -v "^#" $vcf > temp.3
grep "^##" $vcf > temp.1
grep "^#CHROM" $vcf |\
   tr "\t" "\n" |\
   sed 's/:/_/g' |\
   sed 's/(/\t/g' |\
   cut -f 1 |tr "\n" "\t" |\
   sed 's/\t$/\n/g' > temp.2
cat temp.1 temp.2 temp.3 >> ${vcf%.*}-renamed.vcf
rm temp.1 temp.2 temp.3
rename _biallelic_only_maxmissing_0.5.recode-renamed.vcf _cleaned.vcf ${vcf%.*}-renamed.vcf
```
Run it as
```bash
for vcf in *.vcf; do
  ./nameCleaner.sh $vcf;
done
```
To process this in Tassel (converting vcf to ABH format), we will need files with parent names (two text files for each cross, each with single line, listing the parent used in the cross)

```bash
for vcf in *_cleaned.vcf; do
  A=$(grep "^#CHROM" $f |cut -f 10);
  B=$(grep "^#CHROM" $f |cut -f 11);
  echo $A > $A.txt;
  echo $B > $B.txt;
done
```
Compress the files

```bash
gzip *_cleaned.vcf
mkdir orig-vcf
mv *.vcf.gz ./orig-vcf/
```

## Lifting over to NAM coordinates

Download the v4 reference, whcih these VCF files are based on (`Zm-B73-REFERENCE-GRAMENE-4.0.fa`), form [MaizeGDB]().

```bash
# prepare v4
wget https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
gunzip Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
sed 's/^>Chr/>/g' Zm-B73-REFERENCE-GRAMENE-4.0.fa > B73-v4.fa
# prepare v5
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz
hisat2-build -p 36 Zm-B73-REFERENCE-NAM-5.0.fa B73-v5
```


```bash
for f in *.vcf; do grep -v "^##" $f |cut -f 1-4 > $f.1; done
for f in *.1; do grep -v "^#" $f | awk '{print $1"\t"$2-50"\t"$2+50"\t"$3"\t.\t."}' > $f.2; done
rename .vcf.1.2 .bed *.vcf.1.2
for f in *.bed; do bedtools getfasta -fi B73-v4.fa -fo ${f%.*}.fasta -bed $f -name; done
for f in *_cleaned.fasta; do hisat2 -p 12 --mp 1,1 --no-softclip -f -x B73-v5 -U ${f}  1> ${f%.*}_B73-v5.sam 2> ${f%.*}_mapping_stats.txt; done
for f in *.sam; do ./sam2bed.sh $f; done
mkdir process-bed
cd process-bed/
for f in ../*_cleaned_B73-v5.bed; do awk '$3-$2==100' $f > $f.1; done
mv ../*.bed.1 ./
for f in *.1; do awk '{print $1"\t"$2+50"\t"$4}' $f > $f.2; done
for f in *.2; do cut -f 3 $f | sort | uniq -u > ${f%.*}.names; done
for f in ../orig-vcf/*.vcf; do grep -Fw -f $(basename ${f%.*})_B73-v5.bed.1.names $f > ${f%.*}_trim.vcf; done
mv ../orig-vcf/*trim.vcf ./
for f in *.bed.1.2; do awk '{print $3"\t"$0}' $f > $f.3; done
for f in *.vcf; do awk '{print $3"\t"$0}' $f > $f.1; done
for f in *.vcf.1; do g=$(echo $f |sed 's/_trim.vcf.1//g'); awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$0;next}{ print $0, a[$1]}' ${g}_B73-v5.bed.1.2.3 $f | cut -f 5- > ${g}_shuf.txt; done
for f in *_shuf.txt; do awk 'BEGIN{OFS=FS="\t"}{print $(NF-2), $(NF-3), $(NF-1), $0}' $f |rev |cut -f 5- |rev > $f.1; done
for f in ../orig-vcf/*.vcf; do grep "^#" $f > $f.1; done
mv ../orig-vcf/*.1 ./
for f in *_cleaned.vcf.1; do g=$(echo $f |sed 's/.vcf.1//g'); cat $f ${g}_shuf.txt.1 >> ${g}_B73v5.vcf; done
gzip *_cleaned_B73v5.vcf
mkdir ../VCF_files.v5
mv *_cleaned_B73v5.vcf.gz ../VCF_files.v5/
```
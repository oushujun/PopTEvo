## get the outgroup genome and check

```bash
wget https://download.maizegdb.org/Zx-PI566673-REFERENCE-YAN-1.0/Zx-PI566673-REFERENCE-YAN-1.0.fa.gz
perl ./bin/assemblathon_stats.pl Zx-PI566673-REFERENCE-YAN-1.0.fa > Zx-PI566673-REFERENCE-YAN-1.0.fa.stats
```

## get intact LTR

```bash
cat ~/las/git_bin/PopTEvo/TE_annotation/data/*pass.list > NAM.intact.LTR.pass.list
```

## get site list

```bash
perl ~/las/git_bin/PopTEvo/panTE/3.rooting/bin/get_null_site_seq.pl > synLTR_keep.list
```

## get null site seq (~500bp centered on the insertion site)

```bash
for i in `ls ../../NAM_canu1.8/verified/*fasta`; do \
	perl ~/las/git_bin/EDTA/util/call_seq_by_list.pl <(grep yes synLTR_keep.list) -C $i > $(echo $i|sed 's/.*\///').null.500bp.fa & \
done

# check
grep -c \> *bp.fa|sed 's/:/ /'|awk '{i+=$2} END {print i}'
grep yes synLTR_keep.list|wc -l

# aggregate null site seqs
cat *500bp.fa > NAM.LTR.null.500bp.fa
```

## get flank seq of insertions (500bp, 250bp each side)

```bash
perl ./bin/get_flanking_seq.pl

for i in `ls ../../NAM_canu1.8/verified/*fasta`; do \
	perl ~/las/git_bin/EDTA/util/call_seq_by_list.pl synLTR_keep.list.list -C $i > $(echo $i|sed 's/.*\///').flank.250bp.fa & \
done

# aggregate and combine
cat *flank.250bp.fa > NAM.LTR.flank.250bp.fa
perl ./bin/combine_flanking.pl NAM.LTR.flank.250bp.fa > NAM.LTR.flank.cmb500bp.fa
#rm *flank.250bp.fa
```

## blast sites sequences against the outgroup

```bash
# prototype commands:
nohup blastn -query NAM.LTR.flank.cmb500bp.fa -subject Zx-PI566673-REFERENCE-YAN-1.0.fa -out NAM.LTR.flank.cmb500bp.fa.teosinte.out -perc_identity 0.8 -qcov_hsp_perc 0.8 -outfmt 6 -num_alignments 1000 -num_threads 36 &

nohup blastn -query NAM.LTR.null.500bp.fa -subject Zx-PI566673-REFERENCE-YAN-1.0.fa -out NAM.LTR.null.500bp.fa.teosinte.out -perc_identity 0.8 -qcov_hsp_perc 0.8 -outfmt 6 -num_alignments 1000 -num_threads 36 &

# consider to split queries to smaller files and run parallel jobs, eg:
# split to 60 files
split -l 31000 --additional-suffix=.NAM.LTR.flank.cmb500bp.fa NAM.LTR.flank.cmb500bp.fa

# run array jobs
ls x*.NAM.LTR.flank.cmb500bp.fa > split.list
sbatch --array=0-59 ./bin/blast.sub
```

## Identify LTR insertions (null in mexicana, inserted in maize)

```bash
# count full-length copies
perl ./bin/filter_null_sites.pl NAM.LTR.flank.cmb500bp.fa.teosinte.out > NAM.LTR.all.flank.cmb500bp.fa.teosinte.out.fl.count &
perl ./bin/filter_null_sites.pl NAM.LTR.null.500bp.fa.teosinte.out > NAM.LTR.all.null.500bp.fa.teosinte.out.fl.count &

# if you run parallel jobs:
for i in x*out; do \
	perl ./bin/filter_null_sites.pl $i > $i.fl.count & \
done

# only keep single-copy full-length insertion sites
awk '{if ($2==1) print $1}' NAM.LTR.all.flank.cmb500bp.fa.teosinte.out.fl.count NAM.LTR.all.null.500bp.fa.teosinte.out.fl.count|perl -nle 's/.*\|//; print $_' > synLTR_all.mexicana.null.list
```


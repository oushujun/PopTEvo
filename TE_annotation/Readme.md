## Files explained

NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa	Pan-genome TE library (examplars).
*fasta.mod.pass.list	Structurally intact LTRs, LTR_retriever format.
*EDTA.intact.gff3.gz	Structurally intact TEs, GFF3 format.
*fasta.out.gz		Whole-genome TE annotation, Homology-based, generated by RepeatMasker with the pan-TE library. RepeatMasker.out format.
*EDTA.TEanno.gff3.gz	Whole-genome TE annotation, structural + homology, nested annotation (overlapping) contained, GFF3 format.
*EDTA.TEanno.split.gff3.gz	Whole-genome TE annotation, structural + homology, each bp is (almost) uniquely annotated, GFF3 format.
*EDTA.TEanno.sum	Summary of whole-genome TE annotation.
NAM.intact.LTR.genedist.gz	Distance to closest genes (both left and right) for each intact LTR
bin/	Contain scripts used in this section.


## Get pan-genome TE curve


get full length TEs from homo-based masking

```bash
for i in *fasta.out; do
  perl ~/las/git_bin/EDTA/util/find_flTE.pl $i |\
  grep -v -P "CL569186.1|AF013103.1|\)n|cent|Cent|telo|knob|TR-1|osed|sela" > $i.flTE &
done
```

get uniq list of flTEs

```bash
for i in *flTE; do
  awk '{print $10}' $i | grep -v -P 'A-rich|G-rich' | sort -u > $i.list &
done
```

bootstrap pan-TE curve for 1000 times

```bash
for k in {1..100}; do
  for j in {1..10}; do
    for i in `ls *list|grep -v -P 'AB10|Ab10'|shuf`; do
      cat $i >> temp.$j.$k;
      sort -u temp.$j.$k | wc -l;
    done |\
    perl transpose3.pl - > result.$j.$k;
    rm temp.$j.$k;
  done &
done
cat result.* > pan_TE_bootstrap1000.summary26.txt
rm result.*
```

## Identify solo LTRs

Find LTR coordinates from the pan-TE library

```bash
perl ./bin/find_LTR.pl -lib NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa > NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa.LTR.info 
```

Find solo LTR and gather data

```bash
for i in *fasta.out; do 
	perl ./bin/solo_finder.pl -i $i -info NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa.LTR.info > $(echo $i|sed 's/.out//').panTE.solo & 
done

for i in *solo; do 
	perl -i -nle 's/^\s+//; my $id="$1\t" if /^(\S+)_[0-9a-z]+/i; print "${id}$_"' $i & 
done

cat *solo > NAM.EDTA1.9.0.MTEC02052020.TE.v1.anno.solo
```


## Get closest genes of each intact LTR

Get bed files for intact LTR

```bash
for i in ./data/*mod.EDTA.TEanno.gff3; do 
	grep struc $i|grep LTR_retrotransposon | \
	perl -nle 'my ($chr, $str, $end, $info) = (split)[0,3,4,8]; \
		my ($id, $iden) = ($1, $2) if $info =~ /Name=(.*);Classification.*ltr_identity=([0-9.]+);/; \
		print "$chr\t$str\t$end\t$id\t$iden"' \
	> $i.intact.LTR.bed & \
done
```


Find closest genes

```bash
for i in *gff3.intact.LTR.bed; do \
	closest-features --dist <(perl -nle 's/^[0-9a-zA-Z]+_chr/chr/; print $_' $i) \
		<(awk '{print $1" "$2" "$3" "$4" "$6}' ../gene_annotation/data/Zm-$(echo $i|sed 's/\..*//')-*/*noTE) | \
	perl -nle 's/\|/ /g; print $_' > $i.genedist & \
done

for i in *Oh7b*gff3.intact.LTR.bed; do \
	closest-features --dist <(perl -nle 's/^[0-9a-zA-Z]+_chr/chr/; print $_' $i) \
		<(awk '{print $1" "$2" "$3" "$4" "$6}' ../gene_annotation/data/Zm-*Oh7*/*noTE) | \
		perl -nle 's/\|/ /g; print $_' > $i.genedist & \
done
```


Gather gene distance info

```bash
for i in *bed.genedist; do \
	perl -snale '$id=~s/\..*//; print "$id\t$_"' -- -id=$i $i; \
done > NAM.intact.LTR.genedist

perl -i -nle 's/NA NA/NA NA NA NA NA NA/g; s/\s+/\t/g; print $_' NAM.intact.LTR.genedist
```



## Files explained
meth*UMR_chr.bed	UMRs of each genome, obtained from Hufford et al. (2021).
NAM26.intact.lLTR.UMR.nostart	All UM-LTRs originated within 5-LTRs.


Get bed files for left LTR (strand oriented)

```bash
for i in ../TE_annotation/data/*pass.list; do \
	perl -nle 'next if /^#/; my ($id, $in, $dir, $iden)=(split)[0,6,8,7]; \
		my ($chr, $str, $end)=($1, $2, $3) if $id=~/(.*):([0-9]+)\.\.([0-9]+)/; \
		my ($lin, $rin) = ($1-1, $2+1) if $in=~/IN:([0-9]+)\.\.([0-9]+)/; \
		print "$chr\t$str\t$lin\t$id\t$dir\t$iden" if $dir eq "+"; \
		print "$chr\t$rin\t$end\t$id\t$dir\t$iden" if $dir eq "-";' $i \
	> $i.intact.lLTR.bed & \
done
```


Convert genome ids

```bash
mv meth_MS71.ref_MS71.UMR_chr.bed meth_Ms71.ref_Ms71.UMR_chr.bed
mv meth_IL14H.ref_IL14H.UMR_chr.bed meth_Il14H.ref_Il14H.UMR_chr.bed
perl -i -nle 's/Oh7B/Oh7b/g; print $_' meth_Oh7b.ref_Oh7b.UMR_chr.bed
for i in meth_*bed; do 
	perl -nle 'my $id = $1 if /meth_(.*).ref/; print "${id}_$_"' $i > $i.mod & \
done
```


Get UM-LTRs, requires BEDtools `intersect`

```bash
for i in *pass.list.intact.lLTR.bed; do \
	bedtools intersect -a $i -b meth_$(echo $i|sed 's/\..*//')*bed.mod -wo > $i.umr & \
done
```


Get UM LTR inside and outside genes

```bash
for i in *pass.list.intact.lLTR.bed.umr; do 
	perl ~/las/git_bin/EDTA/util/output_by_list.pl 4 <(awk '{if ($13>=200) print $0}' $i) \
		1 <(perl -snale 'my ($chr, $str, $end, $dist1, $dist2)=(split)[0,1,2,10,16]; \
			next unless $chr=~/chr/; $chr = "${1}_${chr}" if $genome=~/([0-9a-z_]+)\..*/i; \
			print "$chr:$str..$end\t$dist1, $dist2" if $dist1 eq 0 or $dist2 eq 0' -- -genome=$i \
			../gene_annotation/data/$(echo $i|sed 's/\..*//').*EDTA.TEanno.gff3.intact.LTR.bed.genedist) | \
	sort -suV > $i.inside & \
done

for i in *pass.list.intact.lLTR.bed.umr; do \
	perl ~/las/git_bin/EDTA/util/output_by_list.pl 4 <(awk '{if ($13>=200) print $0}' $i) \
		1 <(perl -snale 'my ($chr, $str, $end, $dist1, $dist2)=(split)[0,1,2,10,16]; \
			next unless $chr=~/chr/; $chr = "${1}_${chr}" if $genome=~/([0-9a-z_]+)\..*/i; \
			print "$chr:$str..$end\t$dist1, $dist2" if $dist1 ne 0 and $dist2 ne 0' -- -genome=$i \
			../gene_annotation/data/$(echo $i|sed 's/\..*//').*EDTA.TEanno.gff3.intact.LTR.bed.genedist) | \
	sort -suV > $i.outside & \
done
```


Get UM LTRs that are not overlapping with start sites of LTRs

```bash
for i in *pass.list.intact.lLTR.bed.umr.inside *pass.list.intact.lLTR.bed.umr.outside; do \
	perl -snale 'my ($genome, $side) = ($1, $2) if $file=~/([0-9a-z_]+)\..*\.(inside|outside)$/i; \
	my ($lltrs, $lltre, $dir, $umrs, $umre) = (split)[1,2,4,7,8]; \
	print "$genome\t$side\t$_" if ($umrs > $lltrs and $dir eq "+") or ($umre < $lltre and $dir eq "-")' -- -file=$i $i ; \
done > NAM26.intact.lLTR.UMR.nostart
```


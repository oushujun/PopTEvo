## Files explained
*.narrowPeak	Accessible chromatin regions (ACRs) of each NAM genomes obtained from Hufford et al. (2021)
NAM26.intact.LTR.ACR.nostart2	AC-LTRs originated within intact LTR-RTs


Covert seq id

```bash
mv OH7B.allreps_Zm-Oh7B-REFERENCE-NAM-1.0_no-dups-mapq30_peaks.narrowPeak.mod OH7B.allreps_Zm-Oh7b-REFERENCE-NAM-1.0_no-dups-mapq30_peaks.narrowPeak.mod

perl -i -nle 's/OH7B/Oh7b/gi; print $_' OH7B.allreps_Zm-Oh7b-REFERENCE-NAM-1.0_no-dups-mapq30_peaks.narrowPeak.mod

for i in *Peak; do \
	perl -nle 'next unless /^[0-9]+/; my $id = $1 if /-(.*)-REF/; print "${id}_chr$_"' $i > $i.mod & \
done
```

Get LTR coordinate info

```bash
for i in ../TE_annotation/data/*mod.EDTA.TEanno.gff3; do 
	grep struc $i|grep LTR_retrotransposon | \
	perl -nle 'my ($chr, $str, $end, $info) = (split)[0,3,4,8]; \
		my ($id, $iden) = ($1, $2) if $info =~ /Name=(.*);Classification.*ltr_identity=([0-9.]+);/; \
		print "$chr\t$str\t$end\t$id\t$iden"' > $i.intact.LTR.bed & \
done
```


Get LTR overlap with ACRs, requires `intersect` from BEDTools

```bash
for i in *EDTA.TEanno.gff3.intact.LTR.bed; do \
	bedtools intersect -a $i -b *-$(echo $i|sed 's/.*\///; s/\..*//')*narrowPeak.mod -wo > $(echo $i|sed 's/.*\///').acr & \
done
```


Get AC LTR inside and outside genes

```bash
for i in *TEanno.gff3.intact.LTR.bed.acr; do \
	perl ~/las/git_bin/EDTA/util/output_by_list.pl 1 <(awk '{print $1":"$2".."$3"\t"$0}' $i) \
		1 <(perl -snale 'my ($chr, $str, $end, $dist1, $dist2)=(split)[0,1,2,10,16]; \
			next unless $chr=~/chr/; $chr = "${1}_${chr}" if $genome=~/([0-9a-z_]+)\..*/i; \
			print "$chr:$str..$end\t$dist1, $dist2" if $dist1 eq 0 or $dist2 eq 0' -- -genome=$i \
			../TE_annotation/data/$(echo $i|sed 's/\..*//').*EDTA.TEanno.gff3.intact.LTR.bed.genedist) | \
	sort -suV > $i.inside & \
done

for i in *TEanno.gff3.intact.LTR.bed.acr; do 
	perl ~/las/git_bin/EDTA/util/output_by_list.pl 1 <(awk '{print $1":"$2".."$3"\t"$0}' $i) \
		1 <(perl -snale 'my ($chr, $str, $end, $dist1, $dist2)=(split)[0,1,2,10,16]; \
			next unless $chr=~/chr/; $chr = "${1}_${chr}" if $genome=~/([0-9a-z_]+)\..*/i; \
			print "$chr:$str..$end\t$dist1, $dist2" if $dist1 ne 0 and $dist2 ne 0' -- -genome=$i \
			../TE_annotation/data/$(echo $i|sed 's/\..*//').*EDTA.TEanno.gff3.intact.LTR.bed.genedist) | \
	sort -suV > $i.outside & \
done
```

Get AC-LTRs that are not ovelapping with start sites of 5-LTRs

```bash
for i in *acr.inside *acr.outside; do \
	perl -snale 'my ($genome, $side) = ($1, $2) if $file=~/([0-9a-z_]+)\..*\.(inside|outside)$/i; \
		my ($ltrs, $ltre, $acrs, $acre) = (split)[2,3,7,8]; \
		print "$genome\t$side\t$_" if $acrs > $ltrs and $acre < $ltre' -- -file=$i $i; \
done > NAM26.intact.LTR.ACR.nostart2
```



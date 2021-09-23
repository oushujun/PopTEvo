## Files explained
list, A list of genomes for pairwise calculations

*pass.list.count.xls, Pairwise syntenic LTR information


Rename chr names with genome ID

```bash
for i in `cat list`; do 
	perl -i -slane ' s/chr/${ID}_chr/; s/${ID}${ID}/$ID/; print $_ ' -- -ID=$(echo $i|perl -nle 's/\..*//; print $_') $i.pass.list & \
done
```

Make pairwise job scripts, after this step, exectute these scripts in HPC. SLOW (1 month) and HUGE (2T immediate files/pair) warning!
```bash
perl ../bin/make_batch_pairwise_qsub.pl list
```

There are some utility commands in `../bin/memo` for running, checking, cancelling, rerunning, and finalizing a large number of jobs. After this step, you will get `*pass.list.count.xls` files.

## More information about `*pass.list.count.xls` files
1. Site_type: Full/Flanking/Part
              Full: syntenic LTR present with support from both sides
              Flanking: empty site (no synLTR) with support from both flanking sides
              Part: full-flanking or flanking-full situations, uncertain but show some sign of synLTR; no Site_info will be given.
2. Synteny_scaff = Synteny_part = Total_match if Site_type is not Part.
3. Synteny_LTR: Chr/Scaffold information of the syntenic LTR (itself or on the other genome)
4. Syntenic_info: [Flanking://Full:]   provides coordinates of the syntenic locus, not bp-accurate
5. Examples:
`Intact_LTR      Identity        Synteny_scaff   Synteny_part    Synteny_LTR     Total_match     Syntenic_info
Oh7b_chr10:54170076..54178592   0.9669  2       2       Oh7b_chr10,NC350_chr10  2       Sync_Chr:NC350_chr10|Site_type:Full|Site_info:[Flanking://Full:NC350_chr10:98810839..98837298]  
NC350_chr9:126054763..126058744 0.9815  1       1       NC350_chr9      2       Sync_Chr:Oh7b_chr9|Site_type:Flanking|Site_info:[Flanking:Oh7b_chr9:165085357..165085897//Full:]        
NC350_chr1:95665945..95673403   0.9896  2       2       NC350_chr1,Oh7b_chr1    2       Sync_Chr:Oh7b_chr1|Site_type:Full|Site_info:[Flanking://Full:Oh7b_chr1:96630070..96637528]      
Oh7b_chr6:168443232..168455705  0.9926  1       2       Oh7b_chr6       1       Sync_Chr:NC350_chr6|Site_type:Part|Site_info:[Flanking://Full:] 
Oh7b_chr4:130587867..130595046  0.9978  2       2       Oh7b_chr4,NC350_chr4    2       Sync_Chr:NC350_chr4|Site_type:Full|Site_info:[Flanking://Full:NC350_chr4:130875145..130882325]  
Oh7b_chr9:149605376..149611545  0.9668  1       1       Oh7b_chr9       2       Sync_Chr:NC350_chr9|Site_type:Flanking|Site_info:[Flanking:NC350_chr9:107949422..107949874//Full:] 
`

Get perfect syntenic pairs

```bash
for i in *pass.list.count.xls; do 
	grep -v Part $i | \
	perl -nle 'next if /Intact_LTR/; my ($id, $syn) = (split)[0,6]; $syn=~s/.*\[//; $syn=~s/Flanking:\/\///; $syn=~s/\/\/Full:]//; $syn=~s/\]//g; $syn=~s/:/\t/; print "$id\t$syn"' > $i.paird & 
done
```

Modify genome id formats

```bash
for i in *IL14H*paird *MS71*paird; do mv $i $(echo $i|sed 's/IL14H/Il14H/; s/MS71/Ms71/'); done
for i in *AB10*; do mv $i $(echo $i|sed 's/B73_AB10/B73AB10/g'); done
for i in *paird; do perl -i -nle 's/IL14H/Il14H/gi; s/MS71/Ms71/gi; s/B73_AB10/B73AB10/gi; print $_' $i & done
```

Extract intact LTR info in extended bed format

```bash
for i in ../../../TE_annotation/data/*mod.EDTA.TEanno.gff3.gz; do 
	zcat $i | grep struc |grep LTR_retrotransposon | \
	perl -nle 'my ($chr, $str, $end, $info) = (split)[0,3,4,8]; my ($id, $iden) = ($1, $2) if $info =~ /Name=(.*);Classification.*ltr_identity=([0-9.]+);/; print "$chr\t$str\t$end\t$id\t$iden"' > \
	$(echo $i|sed 's/.*\///; s/.gz$//').intact.LTR.bed & \
done
cat *gff3.intact.LTR.bed | perl -nle 's/IL14H/Il14H/gi; s/MS71/Ms71/gi; s/B73_AB10/B73AB10/gi; print $_' > NAM.EDTA1.9.6.intact.LTR.bed
```

Find the exact coordinates of intact LTRs to paired subjects. If there is an intact LTR matching the subject (dist==0, overlapping), then use the intact LTR instead. Pairs with coordinates b/t query and subject >10Mb is discarded. Also add labels for types of informative loci to syntenic LTR pairs. `closest-features` is a function of the SAMtools package.

```bash
for i in *pass.list.count.xls.paird; do 
	closest-features --closest --delim '\t' --dist \
	<(perl -nle 'my $tgt = (split)[2]; my ($chr, $s, $e)=($1, $2, $3) if $tgt=~/^(.*):([0-9]+)\.\.([0-9]+)/; next if $s == '250'; next unless defined $s; print "$chr\t$s\t$e\t$_"' $i|sort -V) \
	NAM.EDTA1.9.6.intact.LTR.bed | \
	perl -nle 'my ($query, $type, $sbjt, $chr, $s, $e, $dist)=(split)[3,4,5,6,7,8,11]; my $query_s = $1 if $query=~/:([0-9]+)\.\./; next if abs($query_s - $s) > 10000000; my $match="$chr:$s..$e"; my $lab = "NA"; $dist = "NA" unless defined $dist; if ($type eq "Full"){ if ($dist eq 0){ $sbjt = $match; $lab = "intact"} else {$lab = "truncated"}} elsif ($type eq "Flanking"){$lab = "null"}; print "$query#intact\t$sbjt#$lab\t$type\t$dist"' > $i.matched.label & \
done
```

Separate two directions, use query_subject as file names, convert pairwise into fully bi-directional genome pairs (A-B and B-A). *.sep files are the final outputs of this step.

```bash
for i in *count.xls.paird.matched.label; do 
	perl -sanle '$genome=~s/_[0-9a-z]+_cmb.*//i; print if /^${genome}_chr/' -- -genome=$i $i | sort -uV > $i.sep;
	perl -sanle '$genome=~s/_[0-9a-z]+_cmb.*//i; print unless /^${genome}_chr/' -- -genome=$i $i | sort -uV > \
	$(echo $i|perl -nle 'my ($id1, $id2, $str)=(split /_/, $_, 3); print "${id2}_${id1}_$str.sep"') & \
done



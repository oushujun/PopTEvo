## Files explained

*.1.gff3.gz, Gene annotations in GFF3 format, obtained from Hufford et al. (2021)

*.1.cds.fa.gz, CDS of each annotated gene in FASTA format, obtained from Hufford et al. (2021)

*.1.cds.fa.TElist, A list of TE-related genes

*.1.gff3.gene.bed.noTE, A list of non-TE related genes in BED format

*_class2.gtf.gz, Transcriptom assembly using Class2 in GTF format, obtained from Hufford et al. (2021)

*_Cufflinks.gtf.gz, Transcriptom assembly using Cufflinks in GTF format, obtained from Hufford et al. (2021)

*_Strawberry.gtf.gz, Transcriptom assembly using Strawberry in GTF format, obtained from Hufford et al. (2021)

*_Stringtie.gtf.gz, Transcriptom assembly using Strawberry in GTF format, obtained from Hufford et al. (2021)

NAM_all.transcripts.intersect.90reciprocal.uniq2, NAM TE annotations that overlap ≥ 90% with at least one transcriptom assembly reciprocally.


## Identify TE-related genes and remove

Identify TE-related genes using the panNAM TE library and `RepeatMasker`

```bash
lib=../TE_annotation/data/NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa
threads=36
for i in *.cds.fa; do \
	RepeatMasker -pa $threads -q -no_is -norna -nolow -div 40 -lib $lib -cutoff 225 -gff $i
	perl ~/las/git_bin/EDTA/util/cleanup_tandem.pl -f $i.masked -Nscreen 0 -nr 0.8 -trf 0 > $i.masked.cln & \
done
```


Identify TE-related genes using `TEsorter`

```bash
threads=36
for i in *.cds.fa; do \
	TEsorter $i -p $threads
	awk '{if ($1!~/#/)print $1}' $i.masked.cleanup $i.rexdb.cls.tsv|sed 's/_.*//'|sort -u > $i.TElist & \
done
rm *.lib *pep *faa *domtbl *dom.tsv *out.gff *dom.gff3
```


Get bed files for genes, requires `gff2bed` from BEDOPS

```bash
for i in ./data/*1.gff3; do \
	gff2bed < $i|grep gene > $i.gene.bed & \
done
```


Remove TEs in genes

```bash
for i in */*gff3.gene.bed; do 
	perl ~/las/git_bin/EDTA/util/output_by_list.pl 1 $i 1 ./data/$(echo $i|sed 's/\/.*//')*TElist -ex > $i.noTE & \
done
```


## Identify TE-related transcript assemblies

Get data from Cyverse

```bash
for i in `ils /iplant/home/shared/NAM/Transcript_Assemblies|awk '{print $2}'`; do \
	nohup iget -r $i & \
done
```


Format genome names

```bash
for i in Il14H/* Ms71/* Il14H/*/* Ms71/*/*; do \
	mv $i $(echo $i|sed 's/IL14H/Il14H/; s/MS71/Ms71/'); \
done

for i in `ls Il14H/* Ms71/* Il14H/*/* Ms71/*/*|grep -v .gz`; do \
	perl -i -nle 's/IL14H/Il14H/g; s/MS71/Ms71/g; print $_' $i & \
done
```


Get transcripts overlapping with TEs, requires `intersect` from BEDTools

```bash
for i in `cat list`; do \
	bedtools intersect -f 0.9 -r -wa -wb -a \
		<(grep -v -P '\s+repeat_region\s+|target_site_duplication|long_terminal_repeat' ../../TE_annotation/$i.*fasta.mod.EDTA.TEanno.gff3 | \
			grep 'chr'|perl -nle 's/^[a-zA-Z0-9_]+_chr/chr/; my $count=(split); print $_ if $count==9') \
		-b <(cat $i/base-assemblies/*.gtf | grep -v -P '#|\s+exon\s+' | awk -F"\t" '!seen[$1, $3, $4, $5]++' | \
			perl -nle 'my $count=(split); print $_ if $count>=9') \
	> $i/$i.base-assemblies.90reciprocal & \
done
```


Combine and get unique entries

```bash
for i in */*base-assemblies.90reciprocal; do \
	perl -slane '$id=~s/\/.*//; print "$id\t$_"' -- -id=$i $i|cut -f1-10,17; \
done | sort -k1,10 -u > NAM_all.transcripts.intersect.90reciprocal.uniq2 &
```


Summarize counts

```bash
grep structural `ls */*base-assemblies.90reciprocal|grep -v AB10`|awk '{print $11}'|sort|uniq -c
grep homology `ls */*base-assemblies.90reciprocal|grep -v AB10`|awk '{print $11}'|sort|uniq -c
grep -v AB10 NAM_all.transcripts.intersect.90reciprocal.uniq2 | grep -c structural
grep -v AB10 NAM_all.transcripts.intersect.90reciprocal.uniq2 | grep -c homology
```



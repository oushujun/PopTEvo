# generate sum files
for i in `cat genome.list`; do perl ~/bin/EDTA/util/count_base.pl $i -s > $i.stats & done

for i in *fasta.mod.EDTA.intact.gff3; do cat $i | perl ~/bin/EDTA/util/gff2bed.pl - structural |   perl -nle 'my ($chr, $s, $e, $anno, $dir, $supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 $chr $s $e NA $dir $anno $supfam"' |     perl ~/bin/EDTA/util/buildSummary.pl -maxDiv 40 -stats ${i%.*.*.*.*}.stats - > ${i%.*.*}.intact.sum 2>/dev/null & done

for i in *mod.EDTA.TEanno.gff3; do cat $i | grep -v structural | perl ~/bin/EDTA/util/gff2bed.pl - homology |   perl -nle 'my ($chr, $s, $e, $anno, $dir, $supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 $chr $s $e NA $dir $anno $supfam"' |     perl ~/bin/EDTA/util/buildSummary.pl -maxDiv 40 -stats ${i%.*.*.*.*}.stats - > ${i%.*.*}.homo.sum 2>/dev/null & done

for i in *mod.EDTA.TEanno.gff3; do cat <(cat $i | grep -v structural | perl ~/bin/EDTA/util/gff2bed.pl - homology | perl -nle 'my ($chr, $s, $e, $anno, $dir, $supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 $chr $s $e NA $dir $anno $supfam"')     <(cat $i | grep structural | perl ~/bin/EDTA/util/gff2bed.pl - structural | perl -nle 'my ($chr, $s, $e, $anno, $dir, $supfam)=(split)[0,1,2,3,8,12]; print "10000 0.001 0.001 0.001 $chr $s $e NA $dir $anno $supfam"') | sort -k5,5 -k6,6n | perl ~/bin/EDTA/util/buildSummary.pl -maxDiv 40 -stats ${i%.*.*.*.*}.stats - > ${i%.*.*}.sum 2>/dev/null & done

# calculate total TE content
for i in `awk '{print $1}' genome.list`; do echo -n "$i "; echo "scale = 4; $(grep 'total interspersed' $i.mod.EDTA.TEanno.sum|awk '{print $4}')/$(grep -P '^All' $i.stats|awk '{print $2}')" | bc; done > Arabidopsis.panEDTA.TE.sum
# awk '{sum+=$2} END {print sum/NR}' Arabidopsis.panEDTA.TE.sum

# aggregate TE sum info
for i in *mod.EDTA.TEanno.sum; do cat <(echo $i|perl -nle 's/\..*//; print "$_\t${_}_cp\t${_}_bp\t${_}_pcnt"')     <(head -32 $i|grep -v -P "\-\-|=|total|masked" | perl -0777 -ne 's/\s+unknown/\nLTR_unknown/; print $_' | grep %) |       perl ~/bin/PopTEvo/TE_annotation/bin/transpose3.pl -; done > Arabidopsis.panEDTA.TE.anno.sum
head -1 Arabidopsis.panEDTA.TE.anno.sum | perl -nle 's/^\S+/Genome/; print $_' > head
cat head <(grep bp Arabidopsis.panEDTA.TE.anno.sum) > Arabidopsis.panEDTA.TE.anno.bp.txt
cat head <(grep pcnt Arabidopsis.panEDTA.TE.anno.sum) > Arabidopsis.panEDTA.TE.anno.pcnt.txt

# sum intact size for LTR TIR and Helitron in each genome
for j in *fasta.mod.EDTA.intact.gff3; do echo -n "$j "; for i in LTR_ /DT Heli; do cat $j | grep ID | grep $i | awk '{print $1"\t"$4"\t"$5"\t"$3}' | perl ~/bin/EDTA/util/combine_overlap.pl - |     perl ~/bin/EDTA/util/count_mask.pl -; done; done | perl -ne 'chomp; print "\n" if /intact/; print "$_\t"' | perl -nle 'next if /^$/; s/\s+$//; s/.*\///; s/\.\S+//; print $_' > Arabidopsis.panEDTA.TE.intact.sum
cat <(echo "Genome LTR TIR Helitron") Arabidopsis.panEDTA.TE.intact.sum > Arabidopsis.panEDTA.TE.intact.sum.txt



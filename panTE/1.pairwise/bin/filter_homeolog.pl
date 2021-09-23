
#############################################################
####Note: this script is developed to identify syntenic intact LTR-RTs in homeologous chromosomes in poliploid genomes
#Two files will be generated: 
#	*.count.xls: syntenic status of each intact LTR that have insertion sites shared among all homeologous chrs, divergence of LTRs can reveal genome evolutionary events (time of the speciation, allo or auto-polyploidy)
#	*.filtered: filtered blast hits that contain syntenic LTR information. For reanalysis, run with -blast *.filtered can speed up the process.


#############################################################
#### Steps ####
#(0) retain uniq intact LTR-RTs that are found in chromosomes
#sort -k7,7 -u genome.fa.pass.list|grep chr > genome.fa.chr.uniq.pass.list

#(1) extract 200 bp sequences surrounding each intact LTR-RT (100 flanking + 100 LTR)
#perl -nle '$_ = (split)[0]; my ($chr, $s, $e)=($1, $2, $3) if /^(.*)\:([0-9]+)\.\.([0-9]+)$/; my ($left_from, $left_to)=($s-100, $s+100-1); my ($right_from, $right_to)=($e-100+1, $e+100); print "$_\t$chr:$left_from..$left_to\n$_\t$chr:$right_from..$right_to"' genome.fa.chr.uniq.pass.list > genome.fa.pass.flanking.list

#(2) extract sequences
#perl LTR_retriever/bin/call_seq_by_list.pl genome.fa.pass.flanking.list -C genome.fa > genome.fa.pass.flanking.list.fa &

#(3) use intact LTR-RT surrounding seq to blast against the genome
#makeblastdb -in genome.fa -out genome.fa -dbtype nucl
#blastn -query genome.fa.pass.flanking.list.fa -db genome.fa -num_threads 20 -max_target_seqs 20 -outfmt=6 > genome.fa.pass.flanking.list.fa.out &

#(4) run this script:
#perl filter_homeolog.pl -pair chr.pairing.list -blast flanking.list.fa.out -pass LTR_retriever.pass.list &

#############################################################
#Author: Shujun Ou (oushujun@msu.edu)
#Last update: 02/11/2019
##start the script:

#!/usr/bin/perl -w
use strict;

my $usage="
perl filter_homeolog.pl \
	-pair	[file]	Pairing information of scaffolds and monoploid chromosomes. Only put the pairs that you want the program to work on. e.g.
				Chr1A	genomeA	1
				Chr1B	genomeB	1
				Chr2A	genomeA	2
				Chr2B	genomeB	2
	-blast	[file]	Blastn results using LTR flanking seqs as query on the genome (-outfmt=6)
	-pass	[file]	The pass.list of intact LTR-RTs derived from LTR_retriever
	-ds	[int]	Filter out insertion sites with left-right extra distance >= this number [100000 (bp)]. Increase -ds for relaxed search.
	-df	[int]	Filter out empty sites with left-right distance >= this number [10 (bp)]. Increase -df for relaxed search.
	-dup	[float]	Filter out syntenic loci duplicated n times [1.5] per chromosome. Strict: 1; Relax: 2
	-ploidy	[int]	Specify the minimum number of syntenic chromosomes for each intact LTR. Default: Estimate from -pair info
	-h		Display this help information
	\n";

my $LTR; #The pass.list of intact LTRs derived from LTR_retriever
my $blast; #blast results input
my $homeo_pair; #pairing of scaffolds and monoploid chromosomes
my $site_distance = 100000; #For insertions, the maximum extra distance (other than the intact LTR) allowed between left and right syntenic site.
my $flank_dist=10; #For empty sites, the maximum allowed distance between left and right flanking sequences.
my $duplicate=1; #duplication level of syntenic sites. 1 means no site duplication (strict). n means a syntenic locus may duplicate n times in a chromosome.Recommend 1-2 for diploids.

my $num_monoploid = 2; #indicate the number of monoploids in the genome. This number will be updated automatically based on the -pair file
my $ploidy = "NA"; #user defined ploidy of the input genome. Can overwrite $num_monoploidy and $duplicate. Use this for high ploidy cases.
my $max_len = ''; #the length of the query sequence for blast. This number will be updated automatically based on the blast input.
my %homeo; #store homeologous pairing of the genome
my %target; #store intact LTR-RTs that contain homeologous information
my %intact; #store matching cases of each intact LTR-RT
my %LTR; #store identity information of each intact LTR-RT
my %pair; #store pairwise homeologous matches
my %genome; #store scaffold id of each (sub)genome
my $filter = 1; #filter out non-informative blast hits

my $k=0;
foreach (@ARGV){
	$blast = $ARGV[$k+1] if /^-blast$/i;
	$homeo_pair = $ARGV[$k+1] if /^-pair$/i;
	$LTR = $ARGV[$k+1] if /^-pass$/i;
	$site_distance = $ARGV[$k+1] if /^-ds$/i;
	$flank_dist = $ARGV[$k+1] if /^-df$/i;
	$duplicate = $ARGV[$k+1] if /^-dup$/i;
	$ploidy = $ARGV[$k+1] if /^-ploidy$/i;
	$max_len = $ARGV[$k+1] if /^-maxlen$/i;
	die $usage if /^-h|-help|--help$/i;
	$k++;
	}

$filter = 0 if $blast =~ /filtered$/; #no more filtering if input is already filtered

open PAIR, "<$homeo_pair" or die "No pairing file!\n$usage"; #pairing of scaffolds and monoploid chromosomes
#open BLAST, "less $blast|" or die "No BLAST file!\n$usage"; #blast results input
open BLAST, "<$blast" or die "No BLAST file!\n$usage"; #blast results input
open LTR, "<$LTR" or die "No LTR info!\n$usage"; #The pass.list of intact LTRs derived from LTR_retriever
open FILTER, ">$LTR.blast.filtered" if $filter == 1; #output filtered blast out lines to reduce reanalyses work load.
open INT, ">$LTR.count.xls" or die $!; #store synteny counts for each intact LTR-RT on homeologous

while (<PAIR>){
	chomp;
	my ($scaff, $genome, $chr)=(split)[0,1,2];
	next if $chr eq "NA";
	$homeo{$scaff}=$chr;
	$genome{$scaff}=$genome;
	}
close PAIR;
#print "$genome{$_} " foreach keys %homeo;exit;

#find the ploidy level of the genome using the -pair file
my %chr = %homeo;
delete $chr{$_} for grep { $chr{$_} eq "NA"} keys %chr; #delete hash with values = NA (no chr info)
my $num_chr = keys %{{ map { $_ => 1 } keys %chr}};
my $num_homeo = keys %{{ map { $_ => 1 } values %chr}};
$num_monoploid = int($num_chr/$num_homeo) if $num_homeo > 0;
$num_monoploid = $ploidy if $ploidy ne "NA";
#print "$num_monoploid\n";exit; #test

while (<LTR>){
	my ($loc, $iden)=(split)[0,7];
	next if $loc =~ /^#/;
	$LTR{$loc}=$iden;
	}
close LTR;
#my $co=keys %LTR;print "$co\n";exit; #test

##process blast results
while (<BLAST>){
	chomp;
	s/^\s+//;
	s/\s+/\t/g;
	my @info = (split);
#	next unless @info == 12;
	#my $count=@info; print "$count\n@info ";exit; #test
#	my ($id, $scaff, $iden, $len, $start, $end, $subject_s, $subject_e, $eval)=(split)[0,1,2,3,6,7,8,9,10];
	my ($id, $scaff, $iden, $len, $start, $end, $subject_s, $subject_e, $eval)=@info[0,1,2,3,6,7,8,9,10];
	my ($id_s, $id_e, $int_s, $int_e)=($1, $2, $3, $4) if $id =~ /^.*:([0-9]+)\.\.([0-9]+)\|.*:([0-9]+)\.\.([0-9]+)$/;
	$max_len = $id_e - $id_s + 1 if $max_len eq ''; #find the maximum length of the input query sequence
	#print "($id, $scaff, $iden, $len, $start, $end, $subject_s, $subject_e, $eval)\n";
	next if $eval>1e-5 or $len < $max_len*0.4; #filter out short alignments and allows full and half alignments. For half it could be LTR or flanking.
	#print "$_\n"; next; #test
	my ($id_scaff, $intact)=($1, $2) if $id=~/^(.*):.*\|(.*)$/;
	my ($intact_s, $intact_e)=($1, $2) if $intact=~/^.*:([0-9]+)\.\.([0-9]+)$/;
	my $strand="+";
	(($subject_s, $subject_e)=($subject_e, $subject_s) and $strand="-") if $subject_s > $subject_e; #align to negative direction
	#print "$id_scaff\n" if $id_scaff=~/Chr12/; #testing
	next if abs($intact_s - $intact_e) > 100000;
	next if abs($subject_s - $subject_e) > 100000;
	my $hit_id="$scaff:$subject_s..$subject_e";

	#basic filtering
	next if $id_scaff eq $scaff; #filter out same-chromosome alignment
	#filter out hits targeting the same (sub)genome
	next unless (defined $genome{$id_scaff} and defined $genome{$scaff});
	next if $genome{$scaff} eq $genome{$id_scaff}; #filter out within (sub)genome alignments
	if (defined $homeo{$id_scaff} and defined $homeo{$scaff}){ #filter out scaffolds not in the -pair file
		next if $homeo{$scaff} ne $homeo{$id_scaff}; #filter out non-homeolog alignments
		}

#	print "$id_scaff ne $scaff\t$homeo{$scaff} eq $homeo{$id_scaff}\t$genome{$scaff} ne $genome{$id_scaff}\n"; #testing
#	print "$_\n"; next;#testing


#design: Assume the sequence is half flanking half LTR. "Full" means the full-length left or right flanking-LTR sequence, solo-LTR is considered as Full.
#ltr_side:  left--------------------right
#	 Flanking-LTR------------LTR-Flanking
#hit_side: left-right		  left-right
#	  Half-Half		  Half-Half
#	    Full		    Full
#print "$id\t$_\n"; next;#test
#print "$id_info, $int_info, $id_s, $id_e, $int_s, $int_e\t$_\n"; next; #test

	#find which side $id_scaff was derived from and if a hit contains flanking sequences
	my ($ltr_side, $hit_side)=('','');
	$ltr_side="left" if $id_s + $max_len/2 == $int_s; #the QUERY is from the left LTR
	$ltr_side="right" if $int_e + $max_len/2 == $id_e; #the QUERY is from the right LTR

	#For partial hits try to keep flanking-contained ones
	if ($len < $max_len*0.75){ #alignment length less then 75% full length is counted as partial alignment
		$hit_side="left" if $start < $max_len*0.25; #the HIT match to the left of the QUERY
		$hit_side="right" if $end > $max_len*0.75; #the HIT match to the right of the QUERY
		next unless $hit_side eq $ltr_side; #exclude hits that has no flanking seq matching to make sure the hit is on the flanking seq but not LTR

		#print out filtered hits to reduce reanalyses work load
		print FILTER "$_\n" if $filter == 1;

		if ($strand eq "-"){ #convert the direction for HITs
			if ($ltr_side eq "left"){
				$ltr_side="right";
				} else {
				$ltr_side="left";
				}
			}
		push @{$target{$intact}{$scaff}{"Flanking"}{$ltr_side}}, $hit_id;
		#print "$intact\t$scaff\tFlanking\t$ltr_side\t$hit_id\n"; #test
		} else {
		#print "full\t$_\n"; #test
		#print out filtered hits to reduce reanalyses work load
		print FILTER "$_\n" if $filter == 1;

		if ($strand eq "-"){ #convert the direction for HITs
			if ($ltr_side eq "left"){
				$ltr_side="right";
				} else {
				$ltr_side="left";
				}
			}
	#Store "Full" hits that contain both flanking and LTR-terminal sequences.
		push @{$target{$intact}{$scaff}{"Full"}{$ltr_side}}, $hit_id;
		#print "$intact\t$scaff\tFull\t$ltr_side\t$hit_id\n"; #test
		}
	}
close BLAST;
close FILTER if $filter == 1;

##find out pairing information of syntanic candidates
foreach my $intact (keys %target){ #each intact LTR-RT
	my $intact_len = abs($2-$1) + 1 if $intact=~/.*:([0-9]+)\.\.([0-9]+)/;
	foreach my $scaff (keys %{$target{$intact}}){ #each matching homeologous chromosome
		foreach my $hit_id_left (@{$target{$intact}{$scaff}{"Flanking"}{"left"}}){ #each partial matching hits to the left flanking seq 
			my ($left_start, $left_end)=($1, $2) if $hit_id_left=~/.*:([0-9]+)\.\.([0-9]+)/;
			#consider flanking-left + flanking-right
			foreach my $hit_id_right (@{$target{$intact}{$scaff}{"Flanking"}{"right"}}){ #try to find the matching-right flanking seq
				my ($right_start, $right_end)=($1, $2) if $hit_id_right=~/.*:([0-9]+)\.\.([0-9]+)/;
				my $distance = $right_start - $left_end + 1; #distance b/t the left and right site
				push @{$target{$intact}{$scaff}{"Flanking"}{"paired"}}, "$scaff:$left_start..$right_end" if $distance <= $flank_dist-5 and $distance >= -$flank_dist-5; #record the paired region which is syntenic to the $id in $scaff. The insertion site is syntenic but no LTR insertion found. Flanking synteny is determined if the flanking sites are in the left-right order with the correct direction. Small overlaps/gaps (+-10 bp) b/t flanking seqs are tolerated.

#test lines
#print "$intact\tFlanking\t$hit_id_left\t$hit_id_right\t$distance\n" if $distance > 50 or $distance < -50; #discarded
#print "$intact\tFlanking\t$hit_id_left\t$hit_id_right\t$distance\n" if ($distance > 10 and $distance <=15) or ($distance >= -15 and $distance < -10); #discarded
#print "$intact\tFlanking\t$hit_id_left\t$hit_id_right\t$distance\n" if $distance <= $flank_dist-5 and $distance >= -$flank_dist-5; #keep, -5 is TSD overlap
#print "$intact\tFlanking\t$hit_id_left\t$hit_id_right\t$distance\n" if $distance <= 50 and $distance > -50; #keep
				}
			#consider flanking-left + full-right ##experimental
			foreach my $hit_id_right (@{$target{$intact}{$scaff}{"Full"}{"right"}}){ #try to find the matching-right full seq
				my ($right_start, $right_end)=($1, $2) if $hit_id_right=~/.*:([0-9]+)\.\.([0-9]+)/;
				my $distance = $right_start - $left_end + 1 + 0.5*$max_len; #distance b/t the left and right site
				push @{$target{$intact}{$scaff}{"Part"}{"paired"}}, "$scaff:$left_start..$right_end" if $distance-0.5*$intact_len <= $site_distance and $distance >= 200; #200 is for the case of solo-LTR, minimum length 200bp
				} 

			}
		foreach my $hit_id_left (@{$target{$intact}{$scaff}{"Full"}{"left"}}){ #each full matching hits to the left full seq
			my ($left_start, $left_end)=($1, $2) if $hit_id_left=~/.*:([0-9]+)\.\.([0-9]+)/;
			my $start=$left_start+$max_len/2; #start of the matched element
			#consider full-left + full-right
			foreach my $hit_id_right (@{$target{$intact}{$scaff}{"Full"}{"right"}}){ #try to find the matching-right full seq
				my ($right_start, $right_end)=($1, $2) if $hit_id_right=~/.*:([0-9]+)\.\.([0-9]+)/;
				my $end=$right_start+$max_len/2; #end of the matched element
				my $distance = $right_start - $left_end + 1 + $max_len; #distance b/t the left and right site
#				my $intact_len = abs($2-$1) + 1 if $intact=~/.*:([0-9]+)\.\.([0-9]+)/;
				push @{$target{$intact}{$scaff}{"Full"}{"paired"}}, "$scaff:$start..$end" if $distance-$intact_len <= $site_distance and $distance >= 200; #keep hits that are pari-end syntenic to $id in $scaff. The syntenic site has the full LTR insertion. Allow some nested TE insertions (default up to 100Kb)  on top of the syntenic LTR. Inclusive for solo-LTR cases (len >= 200bp).

#test lines
#print "$intact\tFull\t$hit_id_left\t$hit_id_right\t$intact_len\t$distance\n" if $distance-$intact_len > $site_distance or $distance < 100; #discarded
#print "$intact\tFull\t$hit_id_left\t$hit_id_right\t$intact_len\t$distance\n" if $distance-$intact_len <= $site_distance and $distance >= 100; #keep
				}

			#consider full-left + flanking-right ##experimental
			foreach my $hit_id_right (@{$target{$intact}{$scaff}{"Flanking"}{"right"}}){ #try to find the matching-right flanking seq
				my ($right_start, $right_end)=($1, $2) if $hit_id_right=~/.*:([0-9]+)\.\.([0-9]+)/;
				my $distance = $right_start - $left_end + 1 + 0.5*$max_len; #distance b/t the left and right site
				push @{$target{$intact}{$scaff}{"Part"}{"paired"}}, "$scaff:$left_start..$right_end" if $distance-0.5*$intact_len <= $site_distance and $distance >= 200; #200 is for the case of solo-LTR, minimum length 200bp

#test lines
#print "$intact\tFull\t$hit_id_left\t$hit_id_right\t$intact_len\t$distance\n" if $distance-0.5*$intact_len <= $site_distance and $distance >= 200; #keep
				}

			}
		}
	}

#goal:
#make sure flanking sites paired correctly, full sites paired correctly

#on the same chr, if full exist, then discard flanking
#multiple flanking/full sites on a chr count as one

#some relaxations
#left full + right flank = insertion
#right flank + left full = insertion (min distance<$site_distance)

#keep targets that has 4 flanking sites
#output number of LTR insertion for each target
#how to deal with repeated counting of syntenic targets?


#count syntenic sites and print out qualified results for each intact LTR-RT
#print INT "Intact_LTR\tIdentity\tSynteny_scaff\tSynteny_ends\n";
print INT "Intact_LTR\tIdentity\tSynteny_scaff\tSynteny_part\tSynteny_LTR\tTotal_match\tSyntenic_info\n";
foreach my $intact (keys %target){
	my $id_scaff=$1 if $intact=~/^(.*):.*$/;
	my %count;
	my $full_scaff = "$id_scaff,";
	my $full_count = 1; #1 is the intact LTR-RT itself
	my $info='';
	my $total_match = 1; #1 is the intact LTR-RT itself
	foreach my $scaff (keys %{$target{$intact}}){ #each matching homeologous chromosome
		$count{$scaff}="Flanking" if exists $target{$intact}{$scaff}{"Flanking"}{"paired"}; #as long as "paired" info is present, it counts
		$count{$scaff}="Part" if exists $target{$intact}{$scaff}{"Part"}{"paired"}; #Part sites will overwright flanking sites if exist, and count as syntenic part #experimental
		($count{$scaff}="Full", $full_scaff.="$scaff," and $full_count++) if exists $target{$intact}{$scaff}{"Full"}{"paired"}; #full sites will overwright flanking sites if exist
		last unless defined $count{$scaff};
		@{$target{$intact}{$scaff}{"Flanking"}{"paired"}} = () unless exists $target{$intact}{$scaff}{"Flanking"}{"paired"};
		@{$target{$intact}{$scaff}{"Full"}{"paired"}} = () unless exists $target{$intact}{$scaff}{"Full"}{"paired"};
		$total_match += @{$target{$intact}{$scaff}{"Flanking"}{"paired"}} + @{$target{$intact}{$scaff}{"Full"}{"paired"}}; #the number of matches for each intact LTR entry
	#	$info.="$intact\t$scaff\t$count{$scaff}: @{$target{$intact}{$scaff}{'Flanking'}{'paired'}}\t||||\t@{$target{$intact}{$scaff}{'Full'}{'paired'}}\n"; #experimental
		$info.="Sync_Chr:$scaff|Site_type:$count{$scaff}|Site_info:[Flanking:@{$target{$intact}{$scaff}{'Flanking'}{'paired'}}//Full:@{$target{$intact}{$scaff}{'Full'}{'paired'}}]\t";
		}

	next if keys %count < $num_monoploid-1; #filter out intact LTR-RTs that have only partial syntenic info (insertion sites not shared).
	next if $total_match > $duplicate*$num_monoploid; #filter out repetitive syntenic sites
	$full_scaff=~s/,$//;

	my $full_part_count=1;
	foreach my $site (values %count){
		$full_part_count++ if $site eq "Part" or $site eq "Full";
		}
#	print "$intact\t$LTR{$intact}\t$full_count\t$count\t$full_scaff\n"; #test
#	print INT "$intact\t$LTR{$intact}\t$full_count\t$full_scaff\n";
#	print INT "$intact\t$LTR{$intact}\t$full_count\t$full_part_count\t$full_scaff\n";
#	print INT "$intact\t$LTR{$intact}\t$full_count\t$full_part_count\t$full_scaff\t$info\n";
	print INT "$intact\t$LTR{$intact}\t$full_count\t$full_part_count\t$full_scaff\t$total_match\t$info\n";
	}


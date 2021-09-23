#!/usr/bin/perl -w
use strict;
use File::Basename;

my $usage = "TBD";

#specify paths
my $script_path = dirname(__FILE__);
my $call_seq = "$script_path/call_seq_by_list.pl"; #program to extract sequence
my $filter = "$script_path/filter_homeolog.pl"; #program to filter out non-informative hits and obtain syntenic LTRs
my $timer = "$script_path/seq_divergence.pl"; #program to estimate sequence divergence time
my $blast = "$script_path/BLAST_worker3.1.pl"; #path to the directory that contains BLAST_worker.pl

#parameters for Genome Timer
my $estimate = 0; #determine whether to de novo estimate LTR divergence. Defalut 0 = no. 1 = estimate.
my $flank_len = "500"; #200bp; consider 500/1000 bp for highly repetitive genomes (eg, maize); use 200 for higher sensitivity
my $min_align = 100; #minimum alignment length for LTRs
my $miu = 3e-8; #mutation rate (per bp per year)
my $threads = 4;

#parameters for filter_homeolog.pl
my $site_distance = 100000;
my $flank_dist = 10;
my $duplicate = 1.5;
my $ploidy = "NA";

#get custom parameters
my $k = 0;
my $genome = '';
my $rand="rand".int(rand(1000000));
foreach (@ARGV){

	#get parameters for Genome Timer
	$estimate = $ARGV[$k+1] if /^-estimate$/i;
	$flank_len = $ARGV[$k+1] if /^-len$/i;
	$blast = $ARGV[$k+1] if /^-blast$/i;
	$threads = $ARGV[$k+1] if /^-t|--threads$/i;
	die $usage if /^-h|-help|--help$/i;

	#get parameters for filter_homeolog.pl
	$site_distance = $ARGV[$k+1] if /^-ds$/i;
	$flank_dist = $ARGV[$k+1] if /^-df$/i;
	$duplicate = $ARGV[$k+1] if /^-dup$/i;
	$ploidy = $ARGV[$k+1] if /^-ploidy$/i;

	#get filename, make combined genome, pass list, and extract pairing info
	if (/^-f$/i){
		my $file = $ARGV[$k+1];
		my $filename = fileparse($file);
		my $id = (split /[._]/, $filename)[0];
	#	my $id = (split /\./, $filename)[0];
	#	my $id = (split /_/, $filename)[0];
		$genome .= $id."_";
		`cat $file >> $rand.cmb.fa`;
		`cat $file.pass.list >> $rand.pass.list`;
		`grep \\> $file | perl -nle 's/>//; my \$chr = "NA"; \$chr = \$1 if /(chr[0-9]+)/i; print "\$_\\t$id\\t\$chr"' >> $rand.pairing`;
		}
	$k++;
	}
$blast = '' unless defined $blast;

die $usage if $genome eq '';

#prepare input files for filter_homeolog.pl
$genome = "$genome"."cmb.fa"; #genome fasta file
my $pairing = "$genome.pairing"; #chromosome pairing info
my $LTR_list = "$genome.pass.list"; #pass.list
`mv $rand.cmb.fa $genome`;
`mv $rand.pass.list $LTR_list`;
`mv $rand.pairing $pairing`;

#estimate divergence time of each intact LTR
if ($estimate == 1){
	#call intact LTR sequence
	`grep -v -P \"Int|#\" $LTR_list | awk '{print \$1\"\\t\"\$1}' > $LTR_list.temp`;
	`perl $call_seq $LTR_list.temp -C $genome > $LTR_list.fa; rm $LTR_list.temp`;

	open Seq, "<$LTR_list.fa" or die $!;
	$/="\n>";
	my %seq;
	while (<Seq>){
		s/>//g;
		my ($id, $seq) = (split /\n/, $_, 2);
		$seq =~ s/\s+//g;
		my ($seq_id, $solo_id) = ($1, $2) if $id =~ /(.*)\|(.*)/;
		push @{$seq{$solo_id}}, "$seq_id\@\@\@$seq";
		}
	$/="\n";
	close Seq;

	open Age, ">$LTR_list.age" or die $!;
#	print Age "qseqid\tsseqid\tqlen\tslen\ttot_len\traw_d\tK2P_d\tJC69_d\traw_T\tK2P_T\tJC69_T\n";
	foreach my $key (sort {$a cmp $b} keys %seq){
		my $seq = @{$seq{$key}}[0];
		next unless defined $seq;
		my $info = `perl $timer -seq $seq -min_aln $min_align -miu $miu -blastplus $blast`;
		print Age "$info";
		}
	close Age;
	$LTR_list = "$LTR_list.age";
	}

#get 500bp flanking seq of each intact LTR
`perl -nle 'next if /^#/; \$_ = (split)[0]; my (\$chr, \$s, \$e)=(\$1, \$2, \$3) if /^(.*)\\:([0-9]+)\\.\\.([0-9]+)\$/; my (\$left_from, \$left_to)=(\$s-$flank_len/2, \$s+$flank_len/2-1); my (\$right_from, \$right_to)=(\$e-$flank_len/2+1, \$e+$flank_len/2); print "\$_\\t\$chr:\$left_from..\$left_to\\n\$_\\t\$chr:\$right_from..\$right_to"' $LTR_list > $genome.pass.blast.${flank_len}bp.list`;
`perl $call_seq $genome.pass.blast.${flank_len}bp.list -C $genome > $genome.pass.blast.${flank_len}bp.fa`;

#get syntenic info with filter_homeolog.pl
#`perl $blast -query $genome.pass.blast.${flank_len}bp.fa -subject $genome -t $threads | perl $filter -pair $pairing -pass $LTR_list -blast - -ds $site_distance -df $flank_dist -dup $duplicate -ploidy $ploidy`;
`perl $blast -query $genome.pass.blast.${flank_len}bp.fa -subject $genome -t $threads > $genome.pass.blast.${flank_len}bp.fa.out`;
#`perl $filter -pair $pairing -pass $LTR_list -blast $genome.pass.blast.${flank_len}bp.fa.out -ds $site_distance -df $flank_dist -dup $duplicate -ploidy $ploidy`;
#`rm $genome.pass.blast.${flank_len}bp.fa.out` if -s "$genome.pass.blast.filtered";



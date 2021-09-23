#!/usr/bin/perl -w
use strict;

##Last Modified: 02/20/2019 Shujun Ou (shujun.ou.1@gmail.com)

my $usage = "#######################

Estimate sequence divergence with JC69 and K2P models.

Usage:	perl seq_divergence.pl -seq [file|string] -seq [file|string]

Inputs:	User can provide either fasta files or sequence strings to the script with -seq.
	For strings, as default (-f 0), the format must follow \$seq_ID@@@\$seq, where \$seq_ID and \$seq are the name and the sequence, respectively.
	For files, only the first sequence will be read. User must also specify (-f 1) to activate the file reading mode.
	No mixture between string and file allows for this version.
	User can perform self-alignment by providing only one sequence, or two sequence alignment by providing two sequences (i.e., -seq A -seq B).
	The one-sequence mode is designed for estimating divergence of terminal repeat in LTR retrotransposons. Only the first hit is considered.

Options:
	-f		[0|1] Indicate input type of -seq. Default, 0 = string. 1 = FASTA file.
	-min_aln	[int] Minimum alignment length. Default 100 (bp).
	-min_olp	[0-1] Minimum alignment coverage of the first sequence. For two sequence alignments only. Default 0.95.
	-miu		[float] Mutation rate. Default 3e-8 (per bp per year).
	-blastplus	[path] Path to the directory that contains blastn. Default read from \$ENV.

Output:	qseqid\tsseqid\tqlen\tslen\ttot_len\traw_d\tK2P_d\tJC69_d\traw_T\tK2P_T\tJC69_T\tn_transition\tn_transversion\n

Credit: Shujun Ou (shujun.ou.1\@gmail.com) 02/20/2019

#######################
";

my $blastplus = '';
my $min_align = 100; #minimal alignment length (bp)
my $min_overlap = 0.95; #minimal alignment portion between two sequences
my $max_off = 20; #maximum terminal offshift for LTR-LTR alignments
my $miu = 3e-8; #mutation rate; per bp per year
my $seq1 = "NA"; #the first sequence
my $seq2 = "NA"; #the second sequence; if only $seq1 is specified, then perform self-alignment
my $mode = "self"; #alignment mode, will automatically adjust based on user inputs
my $file = 0; #default 0, treat -seq1 -seq2 as sequence strings. 1, treat them as fasta files containing 1 sequence each, only read the first sequence

my $k = 0;
foreach (@ARGV){
	$min_align = $ARGV[$k+1] if /^-min_aln$/i;
	$min_overlap = $ARGV[$k+1] if /^-min_olp$/i;
	$miu = $ARGV[$k+1] if /^-miu$/i;
	$seq1 = $ARGV[$k+1] if $seq1 eq "NA" and /^-seq$/i;
	$seq2 = $ARGV[$k+1] if $seq1 ne "NA" and /^-seq$/i;
	$blastplus = $ARGV[$k+1] if /^-blastplus$/i;
	$file = $ARGV[$k+1] if /^-f$/i;
	$k++;
	}
$blastplus = '' unless defined $blastplus;

die $usage if $seq1 eq "NA";

#reformat seq1 and seq2
if ($seq1 eq $seq2){
	$mode = "self";
	} else {
	$mode = "two";
	}

#print "$file\t$mode\n$seq1\n$seq2\n";exit;
my @Blast=();

if ($file == 0){
	$seq1 = ">$1\n$2" if $seq1 =~ /^(.*)@@@(.*)$/;
	$seq2 = ">$1\n$2" if $seq2 =~ /^(.*)@@@(.*)$/;
	my $exec = "${blastplus}blastn -subject <(echo -e \"$seq1\") -query <(echo -e \"$seq2\") -dust no -outfmt \"6 qseqid sseqid sstart send qlen slen length nident btop\" -parse_deflines";
	my $try=0;
	while ($try<10){ #it's possible that sequence wrote in memory is rewritten by other programs and caused blast error, this step will try 10 times to guarantee the blast is run correctly
		@Blast=qx(bash -c '$exec' 2> /dev/null) if defined $exec;
		last if $? == 0;
		$try++;
		}
	} elsif ($file == 1) {
	@Blast = `${blastplus}blastn -subject $seq1 -query $seq2 -dust no -outfmt \"6 qseqid sseqid sstart send slen qstart qend qlen length nident btop\" -parse_deflines`;
	}
#print "$mode\n@Blast";

exit if @Blast < 1 and $mode eq "self";
my $count = 0; #count the number of filtered blast hits
my $info = ''; #store the outcome
foreach (@Blast){
	chomp;
	my ($qseqid, $sseqid, $sstart, $send, $slen, $qstart, $qend, $qlen, $length, $nident, $btop) = (split); #btop=Blast trace-back operations, contains alignment info
	next if $length < $min_align; #control alignment length
	next unless $sstart < 100 or $slen - $send < 100; #make sure the hit located within 100bp from the start/end of input sequence

	if ($mode eq "two"){ #filters for two sequence alignments
		next if $length/$qlen < $min_overlap; #control alignment coverage of the query
		}
	elsif ($mode eq "self"){
		next if $nident eq $qlen; #control self alignment
		next unless ($sstart < $max_off and $qlen - $qend < $max_off) or ($qstart < $max_off and $qlen - $send < $max_off); #requires hit within 20bp of the input start/end
		}
#print "$mode\n$_\n";
#1	838	13774	12216	1305813774	847	797
	$count++;

	$btop =~ s/\d+//g; #remove all matches
	my $len_snp = length $btop;
	next unless $len_snp % 2 == 0; #expect string is even length

	#count transitions and transversions
	my $n_transition = 0; #A<->G; C<->T
	my $n_transversion = 0; #A<->C; A<->T; G<->C; G<->T
	my $n_indel = 0; #[AGCT] <-> -
	while ($btop =~ s/([ATCG-][ATCG-])//i){
		my $snp = $1;
		$n_transition++ if $snp =~ /(AG)|(GA)|(CT)|(TC)/i;
		$n_transversion++ if $snp =~ /(AC)|(CA)|(AT)|(TA)|(GC)|(CG)|(GT)|(TG)/i;
		$n_indel++ if $snp =~ /\-/;
		}

	#estimate evolutionary distance
	my $tot_len = $n_transition + $n_transversion + $nident; #SNP only, indel not counted
	my $raw_d = ($n_transition+$n_transversion) / $tot_len; #percent SNP
	my $JC69_d = 1; #the Jukes-Cantor model K= -3/4*ln(1-4*d/3) adjusts for non-coding sequences, d=$raw_d
	if ($raw_d < 0.66){ #highly diverged sequence could not be adjusted by the JC69 model
		$JC69_d = -3/4*log(1-4*$raw_d/3); #log=ln
		} else {
		$JC69_d = $raw_d;
		}
	my $P = $n_transition / $tot_len; #fraction of transition
	my $Q = $n_transversion / $tot_len; #fraction of transversion
	my $K2P_d = -1/2*log((1-2*$P-$Q)*sqrt(1-2*$Q)); #The Kimura 2-parameter model controls difference b/t transition and transversion rates

	#estimate divergence time T = d/2u
	my $raw_T = sprintf ("%.0f", $raw_d/(2*$miu));
	my $JC69_T = sprintf ("%.0f", $JC69_d/(2*$miu));
	my $K2P_T = sprintf ("%.0f", $K2P_d/(2*$miu));
	$raw_d = sprintf ("%.4f", $raw_d);
	$JC69_d = sprintf ("%.4f", $JC69_d);
	$K2P_d = sprintf ("%.4f", $K2P_d);

#	print "$qseqid\t$sseqid\t$qlen\t$slen\t$tot_len\t$raw_d\t$JC69_d\t$K2P_d\t$raw_T\t$JC69_T\t$K2P_T\n";
#	$info = "$qseqid\t$sseqid\t$qlen\t$slen\t$tot_len\t$raw_d\t$JC69_d\t$K2P_d\t$raw_T\t$JC69_T\t$K2P_T\n" if $info eq '';
#	$info = "$qseqid\t$sseqid\t$qlen\t$slen\t$tot_len\t$raw_d\t$K2P_d\t$JC69_d\t$raw_T\t$K2P_T\t$JC69_T\n" if $info eq '';
	$info .= "$qseqid\t$sseqid\t$qlen\t$slen\t$tot_len\t$raw_d\t$K2P_d\t$JC69_d\t$raw_T\t$K2P_T\t$JC69_T\t$n_transition\t$n_transversion\n";
	}

print "$info";


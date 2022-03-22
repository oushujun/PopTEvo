#!/usr/bin/env perl
use strict;
use warnings;

# find alternative contigs based on minimap2 results
# Shujun Ou (shujun.ou.1@gmail.com)
# 03/18/2022
# minimap2 -a $i.chr $i.scf -t 1 > $i.sam

my $genome = $ARGV[0];
my $info = "$genome.list";
my $match = $ARGV[1]; #sam file from minimap2

my $minq = 60; # if mapping quality lower than this cutoff, it's not a reliable alignment
my $cutoff = 0.95; # if a contig match to the genome over this cutoff, it's an alternative contig
my $maxcount = 1; # if a contig match more than this, it's not an alternative contig

# read seq len info
print STDERR "Read sequence length info...\n";
my %len;
open Info, "<$info" or die;
while (<Info>){
	my ($id, $len) = (split)[0,1];
	$len{$id} = $len;
	}
close Info;

# process sam file
print STDERR "Read SAM file $match and filter alignments...\n";
my %count;
my %match;
open Match, "<$match" or die;
while (<Match>){
	next if /^[#@\[]/;
	next unless (split) >= 9;
	my ($id, $sbj, $sstr, $mapq, $seq) = (split)[0,2,3,4,9];
	next unless $mapq =~ /^[0-9]+$/;
	next if $mapq < $minq;
	my $len = length $seq;
	next unless defined $len{$id} and $len{$id} > 0;
	my $cov = $len / $len{$id};
	next unless $cov >= $cutoff;
	$count{$id}++;
	$match{$id} = [$sbj, $sstr, $mapq, $len, $cov];
#	print "$id\t$sbj\t$len\t$len{$id}\t$cov\n" if $cov >= $cutoff;
	}
close Match;

# process match info
print STDERR "Process matches...\n";
foreach my $id (sort{$a cmp $b} (keys %count)){
	my ($sbj, $sstr, $mapq, $len, $cov) = @{$match{$id}};
	print "$id\t$sbj\t$sstr\t$mapq\t$len\t$len{$id}\t$cov\n" if $count{$id} <= $maxcount;
	}


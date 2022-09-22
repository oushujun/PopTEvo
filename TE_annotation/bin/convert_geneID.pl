#!/usr/bin/env perl
use strict;
use warnings;

# Usage: Read in pangene table, use the first gene or each line as the pangene to index other genes in the same line
# Replace genes on the provided list with pangenes
# If no corresponding pangenes, then use the original gene (unchanged)
# perl convert_geneID.pl pangene.txt genelist.txt > pan_genelist.txt
# Shujun Ou (shujun.ou.1@gmail.com)
# 09/21/2022

my $pangene = $ARGV[0];
my $genedist = $ARGV[1];

open PAN, "<$pangene" or die $!;
my %pangene;
while (<PAN>){
	my ($pan, $others) = (split /\s+/, $_, 2);
	$pangene{$pan} = $pan;
	while ($others =~ s/(\S+)\s+//){
		$pangene{$1} = $pan;
		}
}

open List, "<$genedist" or die $!;
while (<List>){
	s/\s+Zm0/\tgene:Zm0/g;
	while (s/gene:(\S+)/GENE/){
		my $gene = $1;
		my $pan = $gene;
		$pan = $pangene{$gene} if exists $pangene{$gene};
		s/GENE/$pan/;
		}
	print $_;
}
		


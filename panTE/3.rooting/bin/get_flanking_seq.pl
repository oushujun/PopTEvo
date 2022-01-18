#!/usr/bin/env perl
use strict;
use warnings;

my $sites = 'synLTR_keep.list'; #'synLTR_left.list'; #synLTR_keep.list';

# get up- and down-strean 250 bp sequences flanking the given coordinate and combine together
# Shujun Ou (shujun.ou.1@gmail.com)
# 10/13/2021


#Pan_TE_ID_631875        M162W_chr6:184560317..184566978 no      1
#Pan_TE_ID_676258        NA      no      0
#Pan_TE_ID_744722        Ms71_chr4:197807257..197814169  no


open Sites, "<$sites" or die $!;
open List, ">$sites.list" or die $!;

while (<Sites>){
	my ($id, $loc, $null) = (split);
	next if $null eq 'yes';
	next if $loc eq 'NA';
	my ($chr, $str, $end) = ($1, $2, $3) if $loc =~ /(.*):([0-9]+)\.\.([0-9]+)/;
	next unless defined $str;
	my ($up_str, $up_end, $down_str, $down_end) = ($str - 251, $str - 1, $end + 6, $end + 256);
	print List "${id}_up\t$chr:$up_str..$up_end\n";
	print List "${id}_down\t$chr:$down_str..$down_end\n";
	}
close Sites;



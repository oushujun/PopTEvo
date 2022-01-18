#!/usr/bin/env perl
use strict;
use warnings;

my $intactLTR = 'NAM.intact.LTR.pass.list';
my $synLTR = 'synLTR.txt'; #'synLTR_left.txt'; #synLTR_keep.txt';


open SYNLTR, "<$synLTR" or die $!;
print "#Pan_TE\tSite\tNull_site\tIntact#\n";
while (<SYNLTR>){
	next unless /^Pan_TE_ID/;
	my $panID = (split)[0];
	my ($null, $site) = ('no', 'NA');
	($null = 'yes' and $site = $1) if /\s+(\S+)_null\s+/;
	$site = $1 if /\s+(\S+)_intact\s+/ and $site eq 'NA';
	$site =~ s/_merge//gi;
	my $intact_count = 0;
	$intact_count++ while s/intact//;
	print "$panID\t$site\t$null\t$intact_count\n";
	}

# data sample
#Pan_TE_ID_31 B73_chr1:1974551..1984177_intact Tzi8_chr1:2124659..2125153_null Ky21_chr1:2132309..2132803_null M162W_chr1:2384638..2385113_null Ms71_chr1:2185357..2
#Pan_TE_ID_149 B73_chr1:8740888..8748722_intact Tzi8_chr1:9983576..9984083_null Ky21_chr1:9375038..9375541_null M162W_chr1:10331528..10332022_null Ms71_chr1:9349831
#Pan_TE_ID_253 B73_chr1:14757164..14766387_intact Tzi8_chr1:15450881..15451366_null Ky21_chr1:15725636..15726121_null M162W_chr1:16293689..16294174_null Ms71_chr1:1
#Pan_TE_ID_309 B73_chr1:17926441..17933747_intact NA Ky21_chr1:22275283..22275792_null M162W_chr1:19304855..19305365_null Ms71_chr1:18134061..18134585_null NA NA NA
#Pan_TE_ID_480 B73_chr1:26997564..27007157_truncated Tzi8_chr1:28101417..28101908_null Ky21_chr1:30919380..30919871_null M162W_chr1:27766475..27766966_null Ms71_chr
#

#!/usr/bin/env perl
use strict;
use warnings;

my $miniden = 95;
my $mincov = 0.8;

my %query;
while(<>){
	my ($query, $subject, $iden, $len) = (split);
	next if $iden < $miniden;
	my $loc = 'NA';
	my $qlen = 500; # default length is 500bp unless the query has infomation
	$qlen = $2 - $1 + 1 if $query =~ /:([0-9]+)\.\.([0-9]+)\|/;
	next if $len < $qlen*$mincov;
	$query{$query}++;
#	print $_;
	}
foreach (keys %query){
	print "$_\t$query{$_}\n";
	}

#B73_chr1:126761..127255|Pan_TE_ID_1047773       9       88.550  262     29      1       86      346     110970798       110971059       2.27e-84        316
#B73_chr1:167541..168034|Pan_TE_ID_1068130       1       100.000 494     0       0       1       494     75772   76265   0.0     913
#B73_chr1:167541..168034|Pan_TE_ID_1068130       1       100.000 494     0       0       1       494     82248   82741   0.0     913
#B73_chr1:167541..168034|Pan_TE_ID_1068130       1       96.341  246     9       0       195     440     77321   77566   4.71e-111       405
#B73_chr1:189169..189663|Pan_TE_ID_1048070       1       99.420  345     2       0       1       345     96834   97178   9.25e-178       627
#B73_chr1:189169..189663|Pan_TE_ID_1048070       1       98.473  131     1       1       365     495     97102   97231   3.05e-58        230
#B73_chr1:189169..189663|Pan_TE_ID_1048070       8       78.986  138     14      7       5       127     45028120        45027983        3.25e-13        80.5
#B73_chr1:195231..195725|Pan_TE_ID_202055        1       99.798  496     0       1       1       495     102801  103296  0.0     909
#B73_chr1:195231..195725|Pan_TE_ID_202055        1       89.247  93      8       2       61      152     20300542        20300451        8.90e-24        115
#B73_chr1:195231..195725|Pan_TE_ID_202055        3       92.202  218     13      1       215     428     164524758       164524541       4.92e-81        305
#B73_chr1:195231..195725|Pan_TE_ID_202055        3       78.481  395     40      29      1       395     122942851       122943200       2.37e-54        217
#B73_chr1:195231..195725|Pan_TE_ID_202055        3       81.017  295     21      15      112     395     122938587       122938857       6.64e-50        202
#B73_chr1:195231..195725|Pan_TE_ID_202055        3       93.846  130     8       0       1       130     122939427       122939556       3.09e-48        196
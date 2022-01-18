#!/usr/bin/env perl
use strict;
use warnings;
# remove 100% overlap if the element has nested insertions; keep itself if no nesting.
# Shujun Ou (shujun.ou.1@gmail.com) 12/22/2020

my %info;
while (<>){
	chomp;
	my ($genome, $chr, $str, $end, $id, $iden, $chr_n, $str_n, $end_n, $id_n, $clas_n) = (split);
	my ($ele_id, $nest_id) = ("$chr.$str.$end.$id", "$chr_n.$str_n.$end_n.$id_n");
	push @{$info{$ele_id}}, $_;
	}

foreach my $ele_id (keys %info){
	my $count = @{$info{$ele_id}};
	foreach (@{$info{$ele_id}}){
		my ($genome, $chr, $str, $end, $id, $iden, $chr_n, $str_n, $end_n, $id_n, $clas_n) = (split);
		my ($ele_id, $nest_id) = ("$chr.$str.$end.$id", "$chr_n.$str_n.$end_n.$id_n");
		($chr_n, $str_n, $end_n, $id_n, $clas_n) = ('NA', 'NA', 'NA', 'NA', 'NA') if ($ele_id eq $nest_id);
		print "$genome\t$chr\t$str\t$end\t$id\t$iden\t$chr_n\t$str_n\t$end_n\t$id_n\t$clas_n\t$ele_id\n" unless $count > 1 and $ele_id eq $nest_id;
		}
	}
		
	



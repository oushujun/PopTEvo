#!/usr/bin/env perl
use strict;
use warnings;

# Combine up and down-stream flanking sequences into one
# Shujun Ou (shujun.ou.1@gmail.com)
# 10/13/2021

open Seq, "<$ARGV[0]" or die $!;

my %seq;
$/ = "\n>";
while (<Seq>){
	s/>//g;
	my ($id, $seq) = (split /\s+/, $_);
	$seq =~ s/\s+//g;
	$id =~ s/.*\|//;
	my $side = $2 if $id =~ s/(.*)_(up|down)/$1/;
	if (!exists $seq{$id}){
		$seq{$id} = $seq;
		} 
	elsif ($side eq "up"){
		$seq{$id} = $seq.$seq{$id};
		} 
	elsif ($side eq "down") {
		$seq{$id} = $seq{$id}.$seq;
		}
#	print ">$id\t$side\t$seq\n";
	}
close Seq;

foreach (keys %seq){
	print ">$_\n$seq{$_}\n";
	}


#!/usr/bin/perl -w
use strict;
use File::Basename;
# usage: make_batch_pairwise_qsub.pl list

my $jobtype="syntenic";
my $len = 500; #length of flanking-LTR regions for blasting
my $path = "/home/oushujun/las/maize/divergence/pairwise_syntenic"; 
my $mem="20GB"; #nubmer of memory for each job
my $threads=36; #number of CPUs for each job

#job options
my $syntenic = 1; #1 will make artificial tetraploids and find syntenic LTRs

#obtain the exact path for the program location
my $script_path = dirname(__FILE__);

my @file;
while (<>){
	chomp;
	next if /^$/;
	next if /^#/;
	my $file = (split)[0];
	push @file, $file;
	}

my $i=0;
foreach my $file1 (@file){
	for (my $j=$i+1; $j<@file; $j++){
		my $file2=$file[$j];
#		my ($path1, $name1);
#		my ($path2, $name2);
#		($path1, $name1)=($1, $2) if $file1=~/(.*)\/(\w+)_.*/;
#		($path2, $name2)=($1, $2) if $file2=~/(.*)\/(\w+)_.*/;
#                my ($path1, $name1)=($1, $2) if $file1=~/(.*)\/(.*)/;
 #               my ($path2, $name2)=($1, $2) if $file2=~/(.*)\/(.*)/;
		my ($name1, $name2) = ($file1, $file2);
                $name1 = (split /\./, $name1)[0];
                $name2 = (split /\./, $name2)[0];
		next if -e "${name1}_${name2}_cmb.fa.pass.list.count.xls";
		open Qsub, ">worker.$jobtype.$name1.$name2.qsub";
#	open Qsub, ">worker.$jobtype.$path.$file.qsub";
		print Qsub "#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 21-0:00:00
#SBATCH --ntasks-per-node $threads
#SBATCH --mem=$mem

# active your conda environment containing all dependencies
conda activate EDTA

cd $path

echo 'Start:'
date

syntenic=$syntenic

if [ \"\$syntenic\" -eq 1 ]; then

perl $script_path/Genome_Timer3.1.pl -f $file1 -f $file2 -dup 1 -ploidy 2 -len $len -t $threads
rm ${name1}_${name2}_cmb.fa ${name1}_${name2}_cmb.fa.nhr ${name1}_${name2}_cmb.fa.nin ${name1}_${name2}_cmb.fa.nsq

perl $script_path/filter_homeolog.pl -dup 1 -ploidy 2 -pair ${name1}_${name2}_cmb.fa.pairing -pass ${name1}_${name2}_cmb.fa.pass.list -blast ${name1}_${name2}_cmb.fa.pass.list.blast.filtered

fi

echo 'End:'
date

";
		close Qsub;
		}
	$i++;
	}

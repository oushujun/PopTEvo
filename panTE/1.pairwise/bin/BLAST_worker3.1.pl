#!/usr/bin/perl -w
use strict;
use threads;
use Thread::Queue;
use threads::shared;

#Shujun Ou (09/11/2019, shujun.ou.1@gmail.com)
#Update: 10/05/2019

my $usage = "\nAcclerate BLASTn with multithreads.
	Finished queries will be reported in file.fa.fin.list and won't be rerun when resume from interruptions.
	perl BLAST_worker.pl -query file.fa -subject genome.fa -threads 36 > blast.out
	\n";

my $query = '';
my $subject = '';
my $blast = ''; #path to the directory that contains blastn
my $threads = 4;

my $k = 0;
foreach (@ARGV){
	$query = $ARGV[$k+1] if /^-query$/i;
	$subject = $ARGV[$k+1] if /^-subject$/i;
	$blast = $ARGV[$k+1] if /^-blast$/i;
	$threads = $ARGV[$k+1] if /^-t$|^-threads$|^--threads$/i;
	die $usage if /^-h$|^-help$|^--help$/i;
	$k++;
	}

## check files
die "The query file is empty or not existant!\n\n$usage" unless -s $query;
die "The subject file is empty or not existant!\n\n$usage" unless -s $subject;

## make blast database
`${blast}makeblastdb -in $subject -out $subject -dbtype nucl 2> /dev/null` unless -s "$subject.nsq" or -s "$subject.00.nsq";

## read completed list of query sequences
my %skip;
if (-s "$query.fin.list"){
	open Skip, "<$query.fin.list";
	while (<Skip>){
		chomp;
		s/>//g;
		$skip{$_} = $_;
		}
	close Skip;
	}

## Store sequence information
## multi-threading using queue, put query into queue for parallel computation
my $queue = Thread::Queue -> new();
open Query, "<$query" or die $usage;
$/ = "\n>";
while (<Query>){
	chomp;
	s/>//g;
	my ($name, $seq) = (split /\n/, $_, 2);
	next if exists $skip{$name};
	$seq =~ s/\s+//g;
	$seq = uc $seq;
	$queue -> enqueue([$name, $seq]);
	}
$/ = "\n";
$queue -> end();
close Query;

## output finished querys to a file
#open Complete, ">>$query.fin.list" or die $!;

## initiate a number of worker threads and run
foreach (1..$threads){
	threads -> create(\&BLAST);
	}
foreach (threads -> list()){
	$_ -> join();
	}

#close Complete;


# subrotine for blast workers
sub BLAST(){
	while (defined($_ = $queue->dequeue())){
		my ($name, $seq) = (@{$_}[0], @{$_}[1]);
		my $query = ">$name\n$seq";
		my $query_len = length $query;
		my $report = "true";
		my $exec = "timeout 600s ${blast}blastn -db $subject -query <(echo -e \"$query\") -outfmt 6 -word_size 7 -evalue 1e-5 -dust no -max_target_seqs 100000000";
		my $blast = '';
		for (my $try = 0; $try < 100; $try++){ #try 100 times to guarantee the blast is run correctly
			$blast = qx(bash -c '$exec' 2> /dev/null) if defined $query;
			last if $? == 0;
			}
		print "$blast\n";
#		print Complete "$name\n";
		}
	}

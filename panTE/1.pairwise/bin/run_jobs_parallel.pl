#!/usr/bin/perl -w
use strict;
use threads;
use Thread::Queue;
use File::Basename;

my $usage = "
~ ~ ~ Run a list of single-CPU jobs with multi-threading ~ ~ ~

Author: Shujun Ou (shujun.ou.1\@gmail.com)
Date: 11/11/2018

Please modify the job lines in worker{} as indicated in annotated lines
\$target will be the item listed in -list

Usage: perl run_jobs_parallel.pl -list [file] -threads [int]
Options:	-list	[file]	Specify a list of job targets
		-threads [int]	Number of jobs run in parallel (1 CPU each)
";

my $script_path = dirname(__FILE__);
my $target_list=''; #a list of targets for multithreading jobs
my $threads=2;

my $i=0;
foreach (@ARGV){
	$target_list = $ARGV[$i+1] if /^-list$/i;
	$threads = $ARGV[$i+1] if /^-threads$|^-t$/i;
	$i++;
	}

##multi-threading using queue, create a worker module for parallel computation
my $process_q=Thread::Queue->new();
sub worker {
	while (my $target = $process_q -> dequeue()){
		chomp $target;
		$target =~ s/\s+//g;
		#execute the job in a single thread
		`awk '{print \$1}' $target|sort -u > $target.fin.list`;
#		`sh $target`;

		#decide what to do if job fails
		if ($? ne 0){
			print "ERROR: Run job on $target failed!!!\n";
			} else {
		#decide what to do if job succeed
			print "\n".localtime() ." CPU".threads -> self() -> tid().": running on $target succeed\n";
			}
		}
	}


#insert seq names into the worker queue
open List, "<$target_list" or die $usage;
$process_q -> enqueue (<List>);
$process_q -> end(); #stop adding items to the queue
close List;

#work and finish
for (1..$threads){
	threads -> create ( \&worker );
	}
foreach my $thr (threads -> list()){
	$thr -> join();
	}



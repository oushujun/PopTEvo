#combine_count_tables_from_sample_list.pl by sna
use strict; use warnings;

die "usage: perl /home/springer/sna/scripts/combine_count_tables_from_sample_list.pl <samples.txt> <column with nicknames> \n" unless @ARGV == 2;

my $cut = $ARGV[1] - 1;

#instead of manual hash generation, create list of files from samples.txt files
open (my $samples, $ARGV[0]) or die $!;

my ($samp,$cat,$cat_u);
my $count = 0;
my %f; #hash of files
my %u_f; #hash of unique files
my @split;
my @samples = ();

while (my $line = <$samples>){
  $count ++;
  @split = split/\t/, $line;
  $samp = $split[$cut];
  chomp $samp;
  $cat = $samp.''.".counts.txt";
  $cat_u = $samp.''.".unique.counts.txt";
  $f{$count} = "$cat";
  $u_f{$count} = "$cat_u";
  push @samples, $samp;
}

#multi output table

my $counter = 1;
my $family;
my $add;
my @line;
my $header = "family";
my %family_hash;
my %sum_hash;
my %prop_hash; #proportion of reads mapping uniquely
my ($total,$pro);

while ($counter <= $count){
  open (my $file_n,"<",$f{$counter}) or die "Can't open $f{$counter}: $!\n";
  while (my $line = <$file_n>) {
    chomp $line;
    @line = split/\t/, $line;
    $family = $line[0];
    $add = "$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]";
    if ($family =~ m/family/){
      $header = $header."\t".$add;
    }
    else {
      chomp $line[1];
      chomp $line[3];
      $total = $line[1] + $line[3];
      if ($total > 0){
	$pro = $line[1] / $total;
      }
      else {
	$pro = 0;
      }
      $family_hash{$family}{$counter} = $add;
      $sum_hash{$family}{$counter}= $total;
      $prop_hash{$family}{$counter} = $pro;
    }
  }
  close $file_n;
  $counter++;
}

my ($u_file, $m_file,$s_file,$p_file) = ("element_combined_counts.txt", "multi_combined_counts.txt","family_sum_combined_counts.txt","family_prop_unique.txt");
my ($uout, $mout,$sout,$pout);
open($mout,">",$m_file) or die "Couldn't open $m_file: $!";
open($sout,">",$s_file) or die "Couldn't open $s_file: $!";
open($pout,">",$p_file) or die "Couldn't open $p_file: $!";


print $mout "$header\n";
my $sample_line = join("\t",@samples);
print $sout "family\t$sample_line\n";
print $pout "family\t$sample_line\n";

foreach my $fa (sort keys %family_hash){
  print $mout "$fa\t";
  print $sout "$fa\t";
  print $pout "$fa\t";
  for (my $k = 1; $k < $count; $k++){
    if (exists $family_hash{$fa}{$k}){
      print $mout "$family_hash{$fa}{$k}\t";
      print $sout "$sum_hash{$fa}{$k}\t";
      print $pout "$prop_hash{$fa}{$k}\t";
    }
    else {
      print $mout "0\t0\t0\t0\t0\t";
      print $sout "0\t";
      print $pout "0\t";
    }
  }
  
  #last column needs new line character
  if (exists $family_hash{$fa}{$count}){
    print $mout "$family_hash{$fa}{$count}\n";
    print $sout "$sum_hash{$fa}{$count}\n";
    print $pout "$prop_hash{$fa}{$count}\n";
  }
  else {
    print $mout "0\t0\t0\t0\t0\n";
    print $sout "0\n";
    print $pout "0\n";
  }
}



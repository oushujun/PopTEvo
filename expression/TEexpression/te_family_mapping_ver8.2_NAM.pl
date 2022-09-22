#te_family_mapping_ver8.2_NAM.pl by sna
use strict; use warnings;

#This file will read through the output SAM file from HTseq mapping and create an output table of reads mapping to TE families under 4 different conditions. 

#See README_dev_atlas_TEexp.txt for workflow

#EDITED SAM FILE MUST BE SORTED BY READ NAME!!! 

die "usage: <te_family_mapping_ver8_NAM.pl> <HTseq.samout> <sample_id> \n" unless @ARGV == 2;

#create output files
my ($stat_file, $multifam_table, $summary_line,$unique_table) = ("$ARGV[1].stats.txt", "$ARGV[1].counts.txt","te_mapping_summary.txt","$ARGV[1].unique.counts.txt");
my ($statout, $multiout,$summaryline,$uniqueout);

open(my $file, $ARGV[0]) or die $!;

#Create sam output files for unique and all accepted reads
my ($unique_sam, $accepted_sam) = ("$ARGV[1].u.sam","$ARGV[1].accepted.sam");
my ($usamout,$asamout);
open ($usamout,">",$unique_sam) or die "Couldn't open $usamout: $!";
open ($asamout,">",$accepted_sam) or die "Couldn't open $asamout: $!";

my @line;
my $read;
my $famid;
my $hit;
my $teid;
my $geneid;
my %family_hash;
$family_hash{Gene}{unique_1te} = 0; #create this category outside of hash
my %ind_hash;
my $ambig_parts;
my @features;
my ($listnum, $matchnum, $has_te);
my $realhits; 
my $unique_nd = 0; #counter for reads that are unique-mapping but are not defined as gene or te
my $unique_1_gene = 0; #counter for unique reads that hit one gene
my $unique_1_te = 0; #counter for unique reads that hit one te
my $unique_mf_1te = 0; #counter for unique te/gene reads
my $multi_nd = 0; #reads with no feature define and NH>1
my $multi_1fam = 0; #accepted multi-hit reads mapping to only 1 family 
my $multi_g_1fam = 0; #multi-hit reads mapping to a gene and only 1 te family
my ($thisline,$real);
my $lastread = "empty";
my $lastreadlines = "empty";
my ($nh, $ambig_reads, $read_fam, $read_fam_permissive, $unique_tehits, $gene_hit) = ("empty","empty","empty","empty","empty","empty"); #internal categories from former %reads
my ($genename, $family, $unique_tename); #variables from former loop 2
my %multi_unsorttest;
my $ever_nh1 = "FALSE"; # added 9-23-19. This is to prioritize the unique mapping end of a paired-end read when assigning features. 

while(my $line = <$file>){
  @line = split/\t/, $line;
  $read = $line[0];
  $hit = pop @line;

  # Edit 4-2-2021
  # Find the NH tag o determine number of hits
  my @hitsplit;
  foreach my $segs (@line){
    if ($segs =~ m/^NH:/){
      @hitsplit = split/:/,$segs;
      $realhits = $hitsplit[2];
    }
  }

  $thisline = join("\t",@line); #merge fields together that can be converted back to a bam
  chomp $read;
  
  #reset variables for read
  my @matching = (); #set to empty array
  $has_te = "no"; #set to no
  
  if ($lastread !~ m/empty/ and $read !~ m/^$lastread$/){ #if the read is different from the last read, decide what to do with the last read
    #loop 2
    if (exists $multi_unsorttest{$read}){
      die "File must be sorted by read name\n$!\nRepeat entry: $read\n";
    }
    if ($nh == 1 or $ever_nh1 =~ m/TRUE/){ # "or $ever_nh1 =~ m/TRUE/" added 9-23-19. Prior to this change, this would only trigger if the first half of a paired-end read was multi-mapping and the second half was unique mapping. This change improves consistency 
      if ($ambig_reads !~ m/empty/){ #ambiguous reads
	if ($read_fam =~ m/Gene/){ 
	  #Paired-end reads sometimes have one end that is unique to a gene when the other end is ambiguous. In these cases, count read to gene if the ambiguous parts contain the gene id
	  if ($gene_hit !~ m/empty/ and $ambig_reads =~ m/$gene_hit/){ 
	    $unique_1_gene++;
	    print $usamout "$lastreadlines\n";
	    print $asamout "$lastreadlines\n";
	    $family_hash{Gene}{unique_1te}++;
	    $genename = $gene_hit;
	    if (exists $ind_hash{$genename}{unique_1te_teid}){
	      $ind_hash{$genename}{unique_1te_teid}++;
	    }
	    else{
	      $ind_hash{$genename}{unique_1te_teid} = 1;
	    }
	  }
	  elsif ($read_fam_permissive !~ m/empty/ and $read_fam_permissive !~ m/Multifam/){
	    $unique_mf_1te++;
	    print $usamout "$lastreadlines\n";
	    print $asamout "$lastreadlines\n";
	    $family = $read_fam_permissive;
	    if (exists $family_hash{$family}{unique_permissive}){
	      $family_hash{$family}{unique_permissive}++;
	    }
	    else {
	      $family_hash{$family}{unique_permissive} = 1;
	    }
	    if ($unique_tehits !~ m/empty/){
	      $unique_tename =  $unique_tehits;
	      if (exists $ind_hash{$unique_tename}{unique_permissive_teid}){
		$ind_hash{$unique_tename}{unique_permissive_teid}++;
	      }
	      else{
		$ind_hash{$unique_tename}{unique_permissive_teid} = 1;
	      }
	    }
	  }
	  else{ #no read_fam_permissive call, or this is Multifam
	    $unique_nd++;
	  }
	}
	elsif ($read_fam =~ m/Multifam/){ 
	  if ($unique_tehits !~ m/empty/ and $ambig_reads =~ m/$unique_tehits/){ #With paired-end reads, one end could be defined while the other is ambiguous. When this is the case, count towards the defined end
	    $unique_1_te++;
	    print $usamout "$lastreadlines\n";
	    print $asamout "$lastreadlines\n";
	    $unique_tename =  $unique_tehits;
	    if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	      $ind_hash{$unique_tename}{unique_1te_teid}++;
	    }
	    else{
	      $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	    }
	    $family =  $unique_tename;
	    if (exists $family_hash{$family}{unique_1te}){
	      $family_hash{$family}{unique_1te}++; 
	    }
	    else {
	      $family_hash{$family}{unique_1te} = 1;
	    }
	  }
	  else {
	    $unique_nd++; #most Multifam will go here
	  }
	}
	else{ #Non-gene ambiguous reads
	  $family = $read_fam;
	  $unique_1_te++; #1te here means one TE family, not individual element. Should not go to te/gene column of summary file. 2/20/17 changed unique_mf_1te to unique_1_te
	  #this is where a read can be assigned to the unique_1fam column without being assigned to a single TE element
	  print $usamout "$lastreadlines\n";
	  print $asamout "$lastreadlines\n";
	  if (exists $family_hash{$family}{unique_1te}){
	    $family_hash{$family}{unique_1te}++; 
	    if ($unique_tehits !~ m/empty/){ #if one end of a PE is confident in an element, count towards that TE
	      $unique_tename =  $unique_tehits;
	      if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
		$ind_hash{$unique_tename}{unique_1te_teid}++;
	      }
	      else{
		$ind_hash{$unique_tename}{unique_1te_teid} = 1;
	      }
	    }
	  }
	  else{
	    $family_hash{$family}{unique_1te} = 1;
	    if ($unique_tehits !~ m/empty/){
	      $unique_tename =  $unique_tehits;
	      if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
		$ind_hash{$unique_tename}{unique_1te_teid}++;
	      }
	      else{
		$ind_hash{$unique_tename}{unique_1te_teid} = 1;
	      }
	    }
	  }
	}
      }
      else{ #features = 1, non-ambiguous reads
	if ($read_fam !~ m/empty/){
	  if ($read_fam =~ m/Gene/){
	    $unique_1_gene++;
	    print $usamout "$lastreadlines\n";
	    print $asamout "$lastreadlines\n";
	    $family_hash{Gene}{unique_1te}++; #allows for total counts to be used to calculate RPM values
	    if ($gene_hit !~ m/empty/){
	      $genename = $gene_hit;
	      if (exists $ind_hash{$genename}{unique_1te_teid}){
		$ind_hash{$genename}{unique_1te_teid}++;
	      }
	      else{
		$ind_hash{$genename}{unique_1te_teid} = 1;
	      }
	    }
	  }
	  else{
	    unless($unique_tehits =~ m/empty/){ ## 9/26/19
	    $family = $read_fam;
	    $unique_1_te++;
	    print $usamout "$lastreadlines\n";
	    print $asamout "$lastreadlines\n";

	    # This likely also catches reads where first part is multi-mapping and second part is unique mapping to a single element. This results in reads being assigned to "Multifam". Fix this then figure out how to catch the cases where the first part of the read is unique and the second is multi-mapping. 

	    $unique_tename =  $unique_tehits;
	    if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	      $ind_hash{$unique_tename}{unique_1te_teid}++;
	    }
	    else{
	      $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	    }
	    $family = $unique_tename;
	    if (exists $family_hash{$family}{unique_1te}){
	      $family_hash{$family}{unique_1te}++; 
	    }
	    else {
	      $family_hash{$family}{unique_1te} = 1;
	    }

	    
	    # if (exists $family_hash{$family}{unique_1te}){
	    #   $family_hash{$family}{unique_1te}++;
	    #   if ($unique_tehits !~ m/empty/){
	    # 	$unique_tename =  $unique_tehits;
	    # 	if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	    # 	  $ind_hash{$unique_tename}{unique_1te_teid}++;
	    # 	}
	    # 	else{
	    # 	  $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	    # 	}
	    #   }
	    # }
	    # else{
	    #   $family_hash{$family}{unique_1te} = 1;
	    #   if ($unique_tehits !~ m/empty/){
	    # 	$unique_tename =  $unique_tehits;
	    # 	if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	    # 	  $ind_hash{$unique_tename}{unique_1te_teid}++;
	    # 	}
	    # 	else{
	    # 	  $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	    # 	}
	    #   }
	    # }
	  }
	  } ## 9/26/19
	}
	else{
	  $unique_nd++;
	}
      }
    }
    #################################### multi reads
    else{
      if ($read_fam !~ m/empty/){
	if ($read_fam =~ m/Gene/){
	  if ($read_fam_permissive =~ m/empty/){
	    $multi_nd++;
	  }
	  elsif ($read_fam_permissive =~ m/Multifam/){
	    $multi_nd++;
	  }
	  else{
	    $family = $read_fam_permissive;
	    if (exists $family_hash{$family}{fam_ct_gene}){
	      $family_hash{$family}{fam_ct_gene}++;
	    }
	    else{
	      $family_hash{$family}{fam_ct_gene} = 1;
	    }
	    $multi_g_1fam++;
	    print $asamout "$lastreadlines\n";
	  }
	}
	elsif ($read_fam =~ m/Multifam/){
	  $multi_nd++;
	}
	else{
	  $family = $read_fam;
	  if(exists $family_hash{$family}{fam_ct_1}){
	    $family_hash{$family}{fam_ct_1}++;
	  }
	  else{
	    $family_hash{$family}{fam_ct_1} = 1;
	  }
	  $multi_1fam++;
	  print $asamout "$lastreadlines\n";
	}
      }
      else{
	$multi_nd++;
      }
    }
    #Create test for unsorted file
    if ($nh > 15){ #want to check for unsorted without saving a ton of reads
      $multi_unsorttest{$lastread} = $nh;
    }
    #end loop with resetting variables for family
    $famid = "error"; #set to error to find problems 15 Feb 17
    $teid = "error"; #set to error to find problems 19 Feb 17
    $geneid = "error"; #set to error to find problems 19 Feb 17
    ($nh, $ambig_reads, $read_fam, $read_fam_permissive, $unique_tehits, $gene_hit) = ("empty","empty","empty","empty","empty","empty");
    $lastreadlines = "empty";
    $ever_nh1 = "FALSE"; # added 9-23-19
  }
####################################################################################################################
  #loop 1, execute regardless of whether the read matches the last one or not
  if ($lastreadlines =~ m/empty/){
    $lastreadlines = $thisline;
  }
  else {
    $lastreadlines = $lastreadlines."\n".$thisline;
  }
  $nh = $realhits;
  if($nh == 1){ # added 9-23-19
    $ever_nh1 = "TRUE"; # added 9-23-19
  } # added 9-23-19
  chomp $hit;
  if ($hit =~ m/no_feature/){
  }
  elsif ($hit =~ m/too_low_aQual/){
    print "error: too_low_aQual\n";
  }
  elsif ($hit =~ m/not_aligned/){
#    print "error: not aligned\n";
  }
  elsif ($hit =~ m/alignment_not_unique/){
#    print "error: alignment_not_unique\n";
  }
  elsif ($hit =~ m/ambiguous/){
    #if hit is ambiguous, go through each thing hit and decide if it is a TE or gene. 
    $ambig_parts = substr $hit,17,-1;
      @features = split/\+/,$ambig_parts;
    $ambig_reads = $hit;
    if ($hit =~ m/gene/){ #in gff3, genes listed as gene:Zm.....
      $read_fam = "Gene"; #if this read hits a gene, it will count in the gene/te column
      foreach my $it1 (@features){
	if ($it1 !~ m/gene/){
	  $has_te = "yes"; 
	  $famid = $it1;
	  if ($famid =~ m/^XF/){      #26 May 21
	    $famid = substr $famid, 5;
	  }
	  if ($famid =~ m/LTR/){ # Added 26 May 21
	    $famid = substr $famid, 0, -4;
	  }
	  if ($famid =~ m/INT/){
	    $famid = substr $famid, 0, -4;
	  }
	  $teid = $famid; # Only really famid
	}
      }
      #If multiple TEs are hit, see if they all match. Specifically, see if all TEs match the last TE listed.
      my $hitstoTEs = 0;
      foreach my $it2 (@features) {
	if (($has_te =~ m/yes/) and ($it2 =~ m/$famid/ or $it2 =~ m/gene/)){ #edited 15 Feb 17 to add $has_te check
	  push(@matching, $it2);
	  if ($it2 !~ m/gene/){
	    $hitstoTEs++;
	  }
	}
      }
      $listnum = @features;
      $matchnum = @matching;
      if ($listnum == $matchnum and $has_te =~ m/yes/){  #edited 15 Feb 17 to add $has_te check
	unless($read_fam_permissive =~ m/empty/){
	  if ($read_fam_permissive !~ m/$famid/){
	    $read_fam_permissive = "Multifam";
	  }
	}
	else{
	  $read_fam_permissive = $famid;
	  if ($realhits == 1 and $hitstoTEs == 1){
	    $unique_tehits = $famid; #only count fams
	  }
	}
      }
      else{
	$read_fam_permissive = "Multifam";
      }
    }
    else {
	#ambiguous hits that have no reads hitting genes will always hit >1 TE, so do not add to unique_tehits table
      $famid = $features[0];
      if ($famid =~ m/^XF/){      #26 May 21
	$famid = substr $famid, 5;
      }
      if ($famid =~ m/LTR/){ # Added 26 May 21
	$famid = substr $famid, 0, -4;
      }
      if ($famid =~ m/INT/){
	$famid = substr $famid, 0, -4;
      }
      foreach my $it3 (@features) {
	if ($it3 =~ m/$famid/){
	  push(@matching, $it3);
	}
      }
      $listnum = @features;
      $matchnum = @matching;
      if ($listnum == $matchnum){
	unless ($read_fam =~ m/empty/){
	  if ($read_fam !~ m/$famid/){
	    unless ($read_fam =~ m/Gene/){
		$read_fam = "Multifam";
	      }
	  } 
	}
	else{
	  if ($famid =~ m/^XF/){      #26 May 21
	    $famid = substr $famid, 5;
	  }
	  $read_fam = $famid;
	}
      }
      else{
	unless ($read_fam =~ m/empty/){
	  unless ($read_fam =~ m/Gene/){
	      $read_fam = "Multifam";
	    }
	}
	  else{
	    $read_fam = "Multifam";
	  }
      }
    }
  }
    #non-ambiguous hits
  else {
    if ($hit =~ m/gene/){
      $read_fam = "Gene";
      $geneid = $hit;
      if ($realhits == 1){
	$gene_hit = $geneid; 
      }
    }
    else{
      $famid =  $hit;
      if ($famid =~ m/^XF/){      #26 May 21
	$famid = substr $famid, 5;
      }
      if ($famid =~ m/LTR/){
	$famid = substr $famid, 0, -4;
      }
      if ($famid =~ m/INT/){
	$famid = substr $famid, 0, -4;
      }
      $teid =  $famid;
      unless ($read_fam =~ m/empty/){
	if ($read_fam !~ m/$famid/){
	  unless ($read_fam =~ m/Gene/){
	    $read_fam = "Multifam";
	  }
	} 
      }
      else{
	$read_fam = $famid;
      }
      unless ( $read_fam_permissive =~ m/empty/){
	if ($read_fam_permissive !~ m/$famid/){
	  $read_fam_permissive = "Multifam";
	}
      }
      else{
	$read_fam_permissive = $famid;
      }
      if ($realhits == 1){
	$unique_tehits = $teid;
      }
    }
  }

  
  $lastread = $read;
}
  
close $file;

########################################
#Assign last read
########################################
    #loop 2
if ($nh == 1){
  if ($ambig_reads !~ m/empty/){ #ambiguous reads
    if ($read_fam =~ m/Gene/){ 
      #Paired-end reads sometimes have one end that is unique to a gene when the other end is ambiguous. In these cases, count read to gene if the ambiguous parts contain the gene id
      if ($gene_hit !~ m/empty/ and $ambig_reads =~ m/$gene_hit/){ 
	$unique_1_gene++;
	print $usamout "$lastreadlines\n";
	print $asamout "$lastreadlines\n";
	$family_hash{Gene}{unique_1te}++;
	$genename = $gene_hit;
	if (exists $ind_hash{$genename}{unique_1te_teid}){
	  $ind_hash{$genename}{unique_1te_teid}++;
	}
	else{
	  $ind_hash{$genename}{unique_1te_teid} = 1;
	}
      }
      elsif ($read_fam_permissive !~ m/empty/ and $read_fam_permissive !~ m/Multifam/){
	$unique_mf_1te++;
	print $usamout "$lastreadlines\n";
	print $asamout "$lastreadlines\n";
	$family = $read_fam_permissive;
	if (exists $family_hash{$family}{unique_permissive}){
	  $family_hash{$family}{unique_permissive}++;
	}
	else {
	  $family_hash{$family}{unique_permissive} = 1;
	}
	if ($unique_tehits !~ m/empty/){
	  $unique_tename =  $unique_tehits;
	  if (exists $ind_hash{$unique_tename}{unique_permissive_teid}){
	    $ind_hash{$unique_tename}{unique_permissive_teid}++;
	  }
	  else{
	    $ind_hash{$unique_tename}{unique_permissive_teid} = 1;
	  }
	}
      }
      else{ #no read_fam_permissive call, or this is Multifam
	$unique_nd++;
      }
    }
    elsif ($read_fam =~ m/Multifam/){ 
      if ($unique_tehits !~ m/empty/ and $ambig_reads =~ m/$unique_tehits/){ #With paired-end reads, one end could be defined while the other is ambiguous. When this is the case, count towards the defined end
	$unique_1_te++;
	print $usamout "$lastreadlines\n";
	print $asamout "$lastreadlines\n";
	$unique_tename =  $unique_tehits;
	if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	  $ind_hash{$unique_tename}{unique_1te_teid}++;
	}
	else{
	  $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	}
	$family =  $unique_tename;
	if (exists $family_hash{$family}{unique_1te}){
	  $family_hash{$family}{unique_1te}++; 
	}
	else {
	  $family_hash{$family}{unique_1te} = 1;
	}
      }
      else {
	$unique_nd++; #most Multifam will go here
      }
    }
    else{ #Non-gene ambiguous reads
      $family = $read_fam;
      $unique_1_te++; #1te here means one TE family, not individual element. Should not go to te/gene column of summary file. 2/20/17 changed unique_mf_1te to unique_1_te
      #this is where a read can be assigned to the unique_1fam column without being assigned to a single TE element
      print $usamout "$lastreadlines\n";
      print $asamout "$lastreadlines\n";
      if (exists $family_hash{$family}{unique_1te}){
	$family_hash{$family}{unique_1te}++; 
	if ($unique_tehits !~ m/empty/){ #if one end of a PE is confident in an element, count towards that TE
	  $unique_tename =  $unique_tehits;
	  if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	    $ind_hash{$unique_tename}{unique_1te_teid}++;
	  }
	  else{
	    $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	  }
	}
      }
      else{
	$family_hash{$family}{unique_1te} = 1;
	if ($unique_tehits !~ m/empty/){
	  $unique_tename =  $unique_tehits;
	  if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	    $ind_hash{$unique_tename}{unique_1te_teid}++;
	  }
	  else{
	    $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	  }
	}
      }
    }
  }
  else{ #features = 1, non-ambiguous reads
    if ($read_fam !~ m/empty/){
      if ($read_fam =~ m/Gene/){
	$unique_1_gene++;
	print $usamout "$lastreadlines\n";
	print $asamout "$lastreadlines\n";
	$family_hash{Gene}{unique_1te}++; #allows for total counts to be used to calculate RPM values
	if ($gene_hit !~ m/empty/){
	  $genename = $gene_hit;
	  if (exists $ind_hash{$genename}{unique_1te_teid}){
	    $ind_hash{$genename}{unique_1te_teid}++;
	  }
	  else{
	    $ind_hash{$genename}{unique_1te_teid} = 1;
	  }
	}
      }
      else{
	$family = $read_fam;
	$unique_1_te++;
	print $usamout "$lastreadlines\n";
	print $asamout "$lastreadlines\n";
	if (exists $family_hash{$family}{unique_1te}){
	  $family_hash{$family}{unique_1te}++;
	  if ($unique_tehits !~ m/empty/){
	    $unique_tename =  $unique_tehits;
	    if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	      $ind_hash{$unique_tename}{unique_1te_teid}++;
	    }
	    else{
	      $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	    }
	  }
	}
	else{
	  $family_hash{$family}{unique_1te} = 1;
	  if ($unique_tehits !~ m/empty/){
	    $unique_tename =  $unique_tehits;
	    if (exists $ind_hash{$unique_tename}{unique_1te_teid}){
	      $ind_hash{$unique_tename}{unique_1te_teid}++;
	    }
	    else{
	      $ind_hash{$unique_tename}{unique_1te_teid} = 1;
	    }
	  }
	}
      }
    }
    else{
      $unique_nd++;
    }
  }
}
#################################### multi reads
else{
  if ($read_fam !~ m/empty/){
    if ($read_fam =~ m/Gene/){
      if ($read_fam_permissive =~ m/empty/){
	$multi_nd++;
      }
      elsif ($read_fam_permissive =~ m/Multifam/){
	$multi_nd++;
      }
      else{
	$family = $read_fam_permissive;
	if (exists $family_hash{$family}{fam_ct_gene}){
	  $family_hash{$family}{fam_ct_gene}++;
	}
	else{
	  $family_hash{$family}{fam_ct_gene} = 1;
	}
	$multi_g_1fam++;
	print $asamout "$lastreadlines\n";
      }
    }
    elsif ($read_fam =~ m/Multifam/){
      $multi_nd++;
    }
    else{
      $family = $read_fam;
      if(exists $family_hash{$family}{fam_ct_1}){
	$family_hash{$family}{fam_ct_1}++;
      }
      else{
	$family_hash{$family}{fam_ct_1} = 1;
      }
      $multi_1fam++;
      print $asamout "$lastreadlines\n";
    }
  }
  else{
    $multi_nd++;
  }
}


########################################
#Create Output Files
########################################
#test to make sure there are actually reads
if ($unique_1_gene < 10){
#  die "Error Reading File\n$!\n";
}

#open stat out and multi-hit out files
open($multiout,">",$multifam_table) or die "Couldn't open $multifam_table: $!";
open($summaryline,">>",$summary_line) or die "Couldn't open $summary_line: $!";
#open ($uniqueout,">",$unique_table) or die "Couldn't open $unique_table: $!";  #no unique hits table required 6/4/2021

#add to summary file following columns: sample unique_1gene unique_1te unique_speculative_te unique_not_determined multi_1fam multi_speculative_1fam multi_not_determined total_reads
my $total_read_num = $unique_1_gene + $unique_1_te + $unique_mf_1te + $unique_nd + $multi_1fam + $multi_g_1fam + $multi_nd;
print $summaryline "$ARGV[1]\t$unique_1_gene\t$unique_1_te\t$unique_mf_1te\t$unique_nd\t$multi_1fam\t$multi_g_1fam\t$multi_nd\t$total_read_num\n";

#print multi-hit family table

print $multiout "family\t$ARGV[1]_u_te.fam\t$ARGV[1]_u_te.g\t$ARGV[1]_m_te.fam\t$ARGV[1]_m_te.g\t$ARGV[1]_total\n";

my $perfamtot;

foreach my $f (sort keys %family_hash){
  unless (exists $family_hash{$f}{unique_1te}){
    $family_hash{$f}{unique_1te} = 0;
  }
  unless (exists $family_hash{$f}{unique_permissive}){
    $family_hash{$f}{unique_permissive} = 0;
  }
  unless (exists $family_hash{$f}{fam_ct_1}){
    $family_hash{$f}{fam_ct_1} = 0;
  }
  unless (exists $family_hash{$f}{fam_ct_gene}){
    $family_hash{$f}{fam_ct_gene} = 0;
  }
  $perfamtot = $family_hash{$f}{unique_1te} + $family_hash{$f}{unique_permissive} + $family_hash{$f}{fam_ct_1} + $family_hash{$f}{fam_ct_gene};
  print $multiout "$f\t$family_hash{$f}{unique_1te}\t$family_hash{$f}{unique_permissive}\t$family_hash{$f}{fam_ct_1}\t$family_hash{$f}{fam_ct_gene}\t$perfamtot\n";
}


#print out unique hit table #no unique hits table required 6/4/2021

#print $uniqueout "element\t$ARGV[1]_unique\t$ARGV[1]_unique_te.g\n";

# foreach my $individuals (sort keys %ind_hash){
#   unless (exists $ind_hash{$individuals}{unique_1te_teid}){
#     $ind_hash{$individuals}{unique_1te_teid} = 0;
#   }
#   unless (exists $ind_hash{$individuals}{unique_permissive_teid}){
#     $ind_hash{$individuals}{unique_permissive_teid} = 0;
#   }
#   print $uniqueout "$individuals\t$ind_hash{$individuals}{unique_1te_teid}\t$ind_hash{$individuals}{unique_permissive_teid}\n";
# }

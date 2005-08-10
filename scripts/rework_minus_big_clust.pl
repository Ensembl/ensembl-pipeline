#!/usr/local/ensembl/bin/perl -w

use strict;

my $kill_list = shift @ARGV;
my $input_dir = shift @ARGV;

my %dodgy_gene_id;

open(KILLLIST, $kill_list) or die;

while (<KILLLIST>){
  chomp;
  $dodgy_gene_id{$_}++
}

close(KILLLIST);

####

my @output_files;

opendir(DIR, $input_dir) or die "Cant open [$input_dir].";

foreach my $item (readdir DIR){

  next if ($item eq '.' || $item eq '..');

  #print $item . "\n";

  if (-d ($input_dir . '/' . $item)){

    opendir(RECDIR, "$input_dir/$item");# or die "Cant open recursed dir [$input_dir/$item]."
	
      foreach my $recurse_item (readdir RECDIR) {
	next if ($recurse_item eq '.' || $recurse_item eq '..');
	#print "  $recurse_item\n";

	if (-d "$input_dir/$item/$recurse_item"){
	  die "Found a level of recursion in the " . 
	    "output directories that wasnt expected."
	  }

	push @output_files, "$input_dir/$item/$recurse_item";
      }

    closedir(RECDIR)
  } else {
    push @output_files, "$input_dir/$item";
  }
}


closedir(DIR);

my %tally;

while (my $output_file = shift @output_files){

  open(INFILE, $output_file);# or die "Cant open file [$output_file]."

    while (<INFILE>){

      next if /No homologous matches/;

      my ($dn,
	  $ds,
	  $n,
	  $s,
	  $lnl,
	  $threshold_on_ds,
	  $query_gene_id,
	  $query_transcript_id,
	  $query_cigar_line,
	  $query_start,
	  $query_end,
	  $query_cov,
	  $query_identity,
	  $query_similarity,
	  $match_gene_id,
	  $match_transcript_id,
	  $match_cigar_line,
	  $match_start,
	  $match_end,
	  $match_cov,
	  $match_identity,
	  $match_similarity) = split /\t/, $_;

#      next 
#       if ($n == 0 and $s == 0); # As opposed to 0.000 - zero (0) denotes a dodgy computation.

      next
	if ($dodgy_gene_id{$query_gene_id} || 
	    $dodgy_gene_id{$match_gene_id});

      next
	if ($dn > 1) or ($ds > 1) or ($dn < 0) or ($ds < 0);

     next
	if (($query_identity < 60)or($match_identity < 60));

      next
	if (($query_cov < 80) or ($match_cov < 80));

      next 
	if ($tally{$query_gene_id . $match_gene_id} or 
	    $tally{$match_gene_id . $query_gene_id});

      $tally{$query_gene_id . $match_gene_id}++;
      $tally{$match_gene_id . $query_gene_id}++;

      print $_
    }
}

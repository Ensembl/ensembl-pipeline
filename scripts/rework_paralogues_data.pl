#!/usr/bin/env perl


# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use warnings ;
use strict;

my $input_dir = shift @ARGV;

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
	if ($dn > 2) or ($ds > 2);

#      print $_ . "\n"
#	if (($query_identity < 60)or($match_identity < 60));

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

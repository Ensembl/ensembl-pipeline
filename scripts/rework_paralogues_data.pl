#! /usr/local/ensembl/bin/perl -w

use strict;

my %tally;

while (<>){

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

#  next 
#    if ($n == 0 and $s == 0); # As opposed to 0.000 - zero (0) denotes a dodgy computation.

  next
    if ($dn > 2) or ($ds > 2);

  print $_ . "\n"
    if (($query_identity < 60)or($match_identity < 60));

#  next
#    if (($query_cov < 80) or ($match_cov < 80));

  next 
    if ($tally{$query_gene_id . $match_gene_id} or 
	$tally{$match_gene_id . $query_gene_id});

  $tally{$query_gene_id . $match_gene_id}++;
  $tally{$match_gene_id . $query_gene_id}++;

 # print $_
}

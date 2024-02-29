#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

my $MAX_CLUSTER_SIZE = 50;

my %monster_hash;

while (<>){
#0.0000	0.0000	0	0	0	0.6	ENSG00000188551	ENSP00000342676	235M	1	235	100.00	100.00	100.00	ENSG00000175756	ENSP00000319778	235M	1	235	100.00	100.00	100.00
  my ($dn,
      $ds,
      $n,
      $s,
      $lnl,
      $threshold_on_ds,
      $query_gene_id,
      $query_translation_id,
      $query_cigar_line,
      $query_cigar_start,
      $query_cigar_end,
      $query_perc_cov,
      $query_perc_id,
      $query_perc_pos,
      $match_gene_id,
      $match_translation_id,
      $match_cigar_line,
      $match_cigar_start,
      $match_cigar_end,
      $match_perc_cov,
      $match_perc_id,
      $match_perc_pos) = split /\t/, $_;

  push @{$monster_hash{$query_gene_id}},  [$match_gene_id, 
      $dn,
      $ds,
      $n,
      $s,
      $lnl,
      $threshold_on_ds,
      $query_gene_id,
      $query_translation_id,
      $query_cigar_line,
      $query_cigar_start,
      $query_cigar_end,
      $query_perc_cov,
      $query_perc_id,
      $query_perc_pos,
      $match_translation_id,
      $match_cigar_line,
      $match_cigar_start,
      $match_cigar_end,
      $match_perc_cov,
      $match_perc_id,
      $match_perc_pos];

  push @{$monster_hash{$match_gene_id}}, [$query_gene_id,
      $dn,
      $ds,
      $n,
      $s,
      $lnl,
      $threshold_on_ds,
      $query_gene_id,
      $query_translation_id,
      $query_cigar_line,
      $query_cigar_start,
      $query_cigar_end,
      $query_perc_cov,
      $query_perc_id,
      $query_perc_pos,
      $match_gene_id,
      $match_translation_id,
      $match_cigar_line,
      $match_cigar_start,
      $match_cigar_end,
      $match_perc_cov,
      $match_perc_id,
      $match_perc_pos];
}

my %seen_ids;
my @clusters;

foreach my $id (keys %monster_hash) {
  next if $seen_ids{$id};
  $seen_ids{$id}++;

  my @cluster;
  foreach my $match (@{$monster_hash{$id}}){
    my @array = ($id, @$match);
    push @cluster, \@array;
    $seen_ids{$match->[0]}++
  }
  push @clusters, \@cluster;
}

my %already_printed;

for (my $i = 0; $i < scalar @clusters; $i++) {
  next unless ((scalar @{$clusters[$i]} + 1) > $MAX_CLUSTER_SIZE);


  foreach my $match (@{$clusters[$i]}) {

    if (! $already_printed{$match->[0]}){
      print $match->[0] . "\n";
      $already_printed{$match->[0]}++
    }
    if (! $already_printed{$match->[1]}) {
      print $match->[1] . "\n";
      $already_printed{$match->[1]}++
    }
  }
}

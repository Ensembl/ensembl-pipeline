#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

# script to collect the output lines
# from the compare_isoforms runs.
# It collects two types of output:
# 'TRANPAIR' and 'GENEPAIR' lines
# each one containing summary information
# and transcript and gene levels, respectively

my $tran_file            = 'dp_transcript_pairs.out';
my $gene_file            = 'dp_gene_pairs.out';
my $blastz_trans_file    = 'blastz_transcript_pairs.out';
my $exact_trans_file     = 'exact_transcript_pairs.out';
my $semiexact_trans_file = 'semiexact_transcript_pairs.out';
my $blastz_gene_file     = 'blastz_gene_pairs.out';
my $CDS_pair_file        = 'CDS_pairs_comparison_file';

unless ( $ARGV[0] ){
  print STDERR "Script to collect the results from the run of compare_isoforms.pl\n";
  print STDERR "it prints the results into two files transcripts.out and genes.out\n";
  print STDERR "Usage: $0 /dir/with/stderr/results/\n";
  exit(0);
}

my $dir = $ARGV[0];
#print STDERR "dir = $dir\n";
open (LS, "ls $dir |") or die ("Can't open pipe from ls : $!");

open(TRAN, ">$tran_file" )                  or die ("cannot open $tran_file");
open(GENE, ">$gene_file" )                  or die ("cannot open $gene_file");
open(TRANS, ">$blastz_trans_file" )         or die ("cannot open $blastz_trans_file");
open(EXACT, ">$exact_trans_file" )          or die ("cannot open $exact_trans_file");     
open(SEMIEXACT, ">$semiexact_trans_file" )  or die ("cannot open $semiexact_trans_file");
open(BLASTZ_GENE, ">$blastz_gene_file" )    or die ("cannot open $blastz_gene_file" );  
open(CDS, ">$CDS_pair_file" )    or die ("cannot open $CDS_pair_file" );  

FILES:
while(<LS>){
  chomp;
  my $file_name = $_;
  
  #print STDERR "checking $file_name\n";
  open(IN, "<$dir/$file_name") or die "Can't open file [$file_name]\n";
  
 THISFILE:
  while (<IN>){
    chomp; 
    my $line = $_;
    if ( $line =~ /TRANPAIR/i ){
      print TRAN "$line\n";
    }
    elsif ( $line =~ /GENEPAIR/ && !($line =~ /BLASTZ_GENEPAIR/) ){
      print GENE "$line\n";
    }
    elsif ( $line =~ /TRANS_PAIR/ ){
      print TRANS "$line\n";
    }
    elsif ( $line =~ /BLASTZ_GENEPAIR/ ){
	print BLASTZ_GENE "$line\n";
    }
    elsif ( $line =~ /TRANS_SEMI_MATCH/ ){
        print SEMIEXACT "$line\n";
    }
    elsif ( $line =~ /TRANS_EXACT_MATCH/ ){
        print EXACT "$line\n";
    }
    elsif ( $line =~ /CDS_PAIR/ ){
      print CDS "$line\n";
    }
  }
  close(IN);
}
close(CDS);
close(EXACT);
close(SEMIEXACT);
close(BLASTZ_GENE);
close(TRANS);
close(GENE);
close(TRAN);
close(LS);


#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute
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

my $file            = 'pseudogene_label.out';

unless ( defined $ARGV[0] ){
  print STDERR "Script to collect the results from the run of label_pseudogenes.pl\n";
  print STDERR "it prints the results into a file pseudogene_label.out\n";
  print STDERR "Usage: $0 /dir/with/results/\n";
  exit(0);
}

my $dir = $ARGV[0];
#print STDERR "dir = $dir\n";
open (LS, "ls $dir |") or die ("Can't open pipe from ls : $!");

open(TRAN, ">$file" )                  or die ("cannot open $file");

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
    if ( $line =~ /RESULT/ ){
      print TRAN "$line\n";
    }
  }
  close(IN);
}
close(TRAN);
close(LS);


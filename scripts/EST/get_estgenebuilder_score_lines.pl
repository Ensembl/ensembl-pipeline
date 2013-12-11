#!/usr/bin/env perl


# Copyright [1999-2013] Genome Research Ltd. and the EMBL-European Bioinformatics Institute
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

open LS, 'ls -R |' or die "Can't open pipe from ls : $!";

my $dir;
my %time_dist;

my $trans_file   = "transcript_scores";
my $cluster_file = "gene_scores";
my $site_file    = "site_scores";

open( TRANS, ">$trans_file" ) or die ("cannot open file $trans_file");
open( CLUST, ">$cluster_file" ) or die ("cannot open file $cluster_file");
open( SITE, ">$site_file" ) or die ("cannot open file $site_file");

FILES:
while(<LS>){
  chomp;
  my $file_name = $_;
  
  ## filter jobs error files look like 1.161000001-162000000.err
  if ( $file_name =~/\/(\S+):$/ ){
    $dir = $1;
  }
  next unless ( $file_name  =~/(\S+\.\d+-\d+)\.err$/ );
  my $name = $1;
  open(IN, "<$dir/$file_name") or die "Can't open file [$file_name]\n";
  #print STDERR "checking $file_name\n";
    
 THISFILE:
  while (<IN>){
    chomp;
    my $line = $_;
    
    if ( $line=~/TRAN/ ){
      print TRANS "$line\n";
    }
    if ( $line=~/GENE/ ){
      print CLUST "$line\n";
    }
    if ( $line=~/SITE/ ){
      print SITE "$line\n";
    }
  }
  close(IN)
}

close(SITE);
close(CLUST);
close(TRANS);

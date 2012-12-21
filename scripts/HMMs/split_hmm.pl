#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/HMMs/split_hmm.pl,v $
# $Revision: 1.2 $

use warnings ;
use strict;


$/ = '//'; 


SEQFETCH: 
while(<>){
  my @lines = split /\n/;
  unless ( $lines[0] ){
    shift @lines;
  }
  
  my $name;
  foreach my $line (@lines) {
    if ( $line=~/NAME/ ){
      my @entry = split '\s+', $line;
      $name = $entry[1];
      print STDERR "NAME: $name\n";
      last;
    }
  }

  if ( $name ){
    open( HMM, ">$name.hmm") or die("could not open $name.hmm"); 
    foreach my $line (@lines){
      print HMM $line."\n";
    }
    close HMM;
  }
}

$/ = "\n";

#!/usr/local/ensembl/bin/perl -w

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

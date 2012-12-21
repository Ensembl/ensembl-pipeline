#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/DataConversion/tetraodon/annotation/genoscope2simpleGFF.pl,v $
# $Revision: 1.2 $
 
use warnings ;
use strict;
 
use Getopt::Long;
 
# source: 
# EP1_H, EP1_S        (ecore)

my (@features);

while(<>) {
  /^\#/ and next;

  my @l = split /\t/;

  $l[0] =~ s/^chr//;

  if ($l[1] =~ /^EP1/) {    
    my ($fid) = $l[8] =~ /$l[2]\s+(\S+)/;      

    $l[8] = "display_id \"$fid\";";

    print join("\t", @l), "\n";

  }
} 

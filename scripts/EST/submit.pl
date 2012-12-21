#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/EST/submit.pl,v $
# $Revision: 1.3 $

use warnings ;
use strict;

while( <> ) {
  my $input = $_;
  chomp ($input);
  system("$input");
#print "$input\n";

    system("sleep 1");
}



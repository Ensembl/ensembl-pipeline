#!/usr/local/bin/perl -w

use strict;

while( <> ) {
  my $input = $_;
  chomp ($input);
  system("$input");
#print "$input\n";

    system("sleep 1");
}



#!/usr/local/bin/perl -w

use strict;

while( <> ) {
  my $input = $_;
  chomp ($input);
  system("$input");
  print "$input\n";

  print "sleeping for 5 seconds\n";
  system("sleep 15 ");
}

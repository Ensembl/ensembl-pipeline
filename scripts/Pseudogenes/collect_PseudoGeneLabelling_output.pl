#!/usr/local/ensembl/bin/perl -w

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


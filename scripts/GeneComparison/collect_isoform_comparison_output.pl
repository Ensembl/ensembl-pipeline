#!/usr/local/ensembl/bin/perl -w

use strict;

# script to collect the output lines
# from the compare_isoforms runs.
# It collects two types of output:
# 'TRANPAIR' and 'GENEPAIR' lines
# each one containing summary information
# and transcript and gene levels, respectively

my $tran_file = 'transcripts.out';
my $gene_file = 'genes.out';

unless ( $ARGV[0] ){
  print STDERR "Script to collect the results from the run of compare_isoforms.pl\n";
  print STDERR "it prints the results into two files transcripts.out and genes.out\n";
  print STDERR "Usage: $0 /dir/with/stderr/results/\n";
  exit(0);
}

my $dir = $ARGV[0];
#print STDERR "dir = $dir\n";
open (LS, "ls $dir |") or die ("Can't open pipe from ls : $!");

open(TRAN, ">$tran_file" ) or die ("cannot open $tran_file");
open(GENE, ">$gene_file" ) or die ("cannot open $gene_file");

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
    if ($line =~ /GENEPAIR/i){
      print GENE "$line\n";
    }
  }
  close(IN);
}

close(GENE);
close(TRAN);
close(LS);

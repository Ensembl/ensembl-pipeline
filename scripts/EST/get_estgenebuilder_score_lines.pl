#!/usr/local/ensembl/bin/perl -w
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

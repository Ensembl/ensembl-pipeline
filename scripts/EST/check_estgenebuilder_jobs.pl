#!/usr/local/ensembl/bin/perl -w
use strict;

### usage: /dir/where the /chr/dirs/sit/> ~/ensembl-scripts/check_estgenebuilder_jobs.pl

use Getopt::Long;

my $exit;
my $memory;
my $exception;

&GetOptions( 
	    'exit'     => \$exit,
	    'memory'   => \$memory,
	    'exception'=> \$exception,
	   );

unless ( $exit || $memory || $exception ){
  print STDERR "$0 -exit -exception -memory\n";
  exit(0);
}

open LS, 'ls -R |' or die "Can't open pipe from ls : $!";

my $dir;

FILES:
while(<LS>){
  chomp;
  my $file_name = $_;
  
  ## filter jobs error files look like 1.161000001-162000000.err
  if ( $file_name =~/\/(\S+):$/ ){
    $dir = $1;
  }
  next unless ( $file_name  =~/(\S+\.\d+-\d+)\./ );
  

  #unless ( $exception || $memory ){
  # next unless ( $file_name  =~/(\S+\.\d+-\d+)\.out$/ );
  #}
  #unless ( $exit ){
  # next unless ( $file_name  =~/(\S+\.\d+-\d+)\.err$/ );
  #}
  my $name = $1;

  open(IN, "<$dir/$file_name") or die "Can't open file [$file_name]\n";
  
  #print STDERR "checking $file_name\n";
  
 THISFILE:
  while (<IN>){
    chomp;
    my $line = $_;
    
    if ( $exit && $line =~ /exit/i ){
      #print "$file_name: $line\n";
      $file_name =~/(\S+\.\d+-\d+)\.out/;
      my $input_id = $1;
      print "$input_id Exited\n";
    }
    if ( $exception && $line =~ /exception/i ){
      $file_name =~/(\S+\.\d+-\d+)\./; 
      my $input_id = $1; 
      print "$input_id Exception\n"; 
    }
    if ( $memory && $line =~ /out of memory/i ){
      $file_name =~/(\S+\.\d+-\d+)\./;  
      my $input_id = $1;  
      print "$input_id Out of memory\n";  
    }
  }
  close(IN);
#  if ($exit >0 ){
#    print STDERR "$file_name: Exited\n";
#  }
  
#  if ( $exception != 0 && $wrote_gene != 0){
#      print STDERR "$file_name : found $exception exceptions and $wrote_gene written genes\n";
#  }
#  elsif ( $exception !=0 && $wrote_gene == 0 ){
#      print STDERR ">>> $file_name <<< : found $exception exceptions\n";
#  }
  #elsif ( $exception == 0 && $wrote_gene == 0){
  #    print STDERR "$file_name : no genes found\n";
  #}
#  else{
#      next;
#  }
}

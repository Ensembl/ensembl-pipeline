#!/usr/local/ensembl/bin/perl -w
use strict;
use Getopt::Long;

my $job_file;

&GetOptions( 
	    'job_file:s'     => \$job_file,
	   );

unless ( $job_file ){
  print STDERR "$0 -job_file file_with_bsubs < input_ids\n";
  exit(0);
}


while(<>){
  chomp;
  my $input_id = $_;
  $input_id    =~ /(\S+\.\d+-\d+)/;
  #print $input_id."\n";
  
  #open (JOB,"<$job_file") or die("cannot open file $job_file");
  #while(<JOB>){
  #  chomp;
  #  my $line = $_;
  #  if ( $line =~/\s+$input_id/ ){
  #    print $line."\n";
  #  }
  #}
  #close (JOB);
  my $target = " $1";
    
  #my $command = "grep \"$target\" $job_file";
  #print $command."\n";;
  my $line = `grep "$target" $job_file`;
  print $line;

}

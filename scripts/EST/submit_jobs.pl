#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Env; 

###############
#### Usage ####
###############
#
# submit_jobs.pl [ -jobs -sleep -slowdown -limit ] bsub_commands_file.jobs
#

####################
#### Parameters ####
####################

# how many jobs we submit at once
my $jobs     = 100;

# how many seconds we sleep between every <$jobs> jobs
my $sleep    = 150;

# how many jobs we allow to be in the queue before we panic and go to sleep again
my $slowdown = 200; 

# how many jobs we allow in the queue before we STOP the whole thing to prevent the farm to melt
my $panic    = 10000;

# the file with the bsub commands
my $file;

&GetOptions( 
            'jobs:i'     => \$jobs,
            'sleep:i'    => \$sleep,
            'slowdown:i' => \$slowdown,
            'panic:i'    => \$panic,
	    'file:s'     => \$file,
	   );

unless ( $file ){
  print "Usage: submit_jobs.pl  [ -jobs N1 -sleep N2 -slowdown N3 -panic N4 ] -file bsub_commands_file\n";
  print "       where sensible values are N1..N4 integers, N4>N3>N1 and 600<N2<1200\n";
  exit(0);
}

# open the job file
open(FILE,"<$file");

# loop over the job file
my $counter = 1;
while( <FILE> ) {
  
  # this is supposed to be reading from the bsub_jobs file
  my $input = $_;
  chomp ($input);
  system("$input");
  print "$counter: $input\n";
  
  if ($counter%$jobs == 0){
    
    # every job takes T minutes in one processor
    # with N processors, if we submit say 400 jobs, that will take 400*T/N minutes
    # we then sleep for S = 60*400*T/N seconds every 400 jobs
    # IDEALLY number_jobs = number_processors

    # for ecsnodes N = 32        --> for T=10 and N=32  => S= 125min * 60 = 7500 (2 hours and 5 min)
    # for the whole farm N = 382 --> for T=10 and N=382 => S= 10.5min* 60 = 630 

    # T is CPU time, not real time, we could give to sleep a slightly higher number

    my $minutes = $sleep/60;
    print "sleeping for $sleep seconds ( $minutes minute(s) )\n";
    system("sleep $sleep");
    
    ## wake up and check howmany jobs are left to run
    #system("bjobs -w -u eae -q acari | wc -l > number_file");
    my $result = &read_jobs;

    
    # bjobs' output contains a header with column names, hence...
    print "jobs still running (or pending) in the farm: $result\n";
    
    if ( $result > $panic ){
      print "We have too many jobs\n";
      print "counter: $counter\n";
      print "last job submitted: $input\n";
      die("LSF may go nuts!");
    }
    
    while ( $result > $slowdown ){
      my $slow    = $sleep;
      my $minutes = $slow/60;
      print "sleeping for another $slow seconds ( $minutes minute(s) )\n";
      system("sleep $slow");
      $result = &read_jobs;
      print "jobs still running (or pending) in the farm: $result\n";
    }
  } 
  $counter++; 
}

close( FILE );


############################################################


sub read_jobs{
  #system("busers | grep eae > number_file");
  system("bjobs -w | wc -l > number_file");
  open(IN,"<number_file");
  my $result;
  while(<IN>){
    chomp;
    #my @entry = split;
    #$result = $entry[5];
    $result = $_;
  }
  close(IN);
  return $result;
}


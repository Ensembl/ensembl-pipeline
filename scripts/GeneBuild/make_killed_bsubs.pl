#!/usr/local/ensembl/bin/perl -w


#this is a script which takes the file produced by RuleManager_Genebuild.pl of killed jobs 
#and creates bsubs lines to run them by hand

use strict;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (GB_DBNAME
							     GB_DBHOST
							     GB_DBUSER
							     GB_DBPASS
							    );
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw (GB_KILLED_RUNNER
							   GB_QUEUE
							   GB_KILLED_OUTPUT_DIR
							   GB_KILLED_INPUT_IDS
							   GB_KILLED_BSUB_FILE
							   GB_SPLIT_JOB_SIZE
							  );





if((!$GB_KILLED_INPUT_IDS) || (!$GB_SPLIT_JOB_SIZE) || (!$GB_KILLED_OUTPUT_DIR) || (!$GB_KILLED_RUNNER)){
  print STDERR "Need these options filled in a Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts\n";
  print STDERR "GB_KILLED_INPUT_IDS location of killed file produced by RuleManager\n";
  print STDERR "GB_KILLED_OUTPUT_DIR location of directory to write the LSF output\n";
  print STDERR "GB_KILLED_RUNNER location of runner script \n";
  print STDERR "GB_SPLIT_JOB_SIZE size of jobs to be split to\n";
}
open(FH, $GB_KILLED_INPUT_IDS) or die("couldn't open ".$GB_KILLED_INPUT_IDS." $!");
my $queue = $GB_QUEUE;
if(!$queue){
  $queue = 'acari';
  print STDERR "GB_QUEUE not defined assuming queue is acari\n";
}
my $fh;
if($GB_KILLED_BSUB_FILE){
  open(KILL, ">".$GB_KILLED_BSUB_FILE) or die "couldn't open ".$GB_KILLED_BSUB_FILE." $!";
  $fh = \*KILL;
}else{
  print STDERR "you haven't specified a file to write the bsub lines to writing to stdout\n";
  $fh = \*STDERR;
}

while(<FH>){
  chomp;
  my ($input_id, $logic_name, $runnable) = split;
  #print STDERR "have ".$input_id." analysis ".$logic_name." runnable ".$runnable."\n";
  my ($outfile, $errfile) = &make_filenames($GB_KILLED_OUTPUT_DIR, $input_id, $logic_name); 
  
  my $runner = $GB_KILLED_RUNNER;
  my $command = "bsub -q $queue -R alpha -C0 -o $outfile -e $errfile -E \"$runner -check -runnable  Bio::EnsEMBL::Pipeline::RunnableDB::$runnable -analysis $logic_name\"";
  $command .= "  $runner ";
  $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable -analysis $logic_name -split_size $GB_SPLIT_JOB_SIZE";
  $command .= " -input_id $input_id ";
  $command .= " -write"; 
  
  print $fh $command."\n";
  
}




sub make_filenames {
  my ($dir, $input_id, $logic_name) = @_;
  
  my $num = int(rand(10));

  my $full_dir = $dir . "/$num/";
  if( ! -e $full_dir ) {
    system( "mkdir $full_dir" );
  }

  my $stub = $input_id.".";
  $stub .= $logic_name.".";
  $stub .= time().".".int(rand(1000));

  my $stdout = $full_dir."/".$stub.".out";
  my $stderr = $full_dir."/".$stub.".err";
 
  return($stdout, $stderr);
}

#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer;
use Getopt::Long;

$| = 1;

print STDERR "starttime: ",time,"\n";

my $host   = 'ecs1b';
my $port   = undef;
my $dbname = 'abel_compara_humanNCBI28_mouse';
my $dbuser = 'ensadmin';
my $pass   = 'ensembl';
my $alnprog = 'bl2seq';
my $min_score = 0.01;

&GetOptions('host:s' => \$host,
	    'port:i' => \$port,
	    'dbname:s' => \$dbname,
	    'dbuser:s' => \$dbuser,
	    'pass:s' => \$pass,
	    'alnprog:s' => \$alnprog,
	    'min_score:f' => \$min_score);

my $input_file = shift;

if (! defined $input_file) {
  die "Must call run_compara with an input file!";
}


my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor (-host => $host,
						      -user => $dbuser,
						      -pass => $pass,
						      -dbname => $dbname );


open INPUT, "$input_file" || die "could not open $input_file, $!\n";

my $log_file = $input_file.".log";
my %already_ran;

if (-e $log_file) {
  open LOG, "$log_file";
  while (defined (my $line =<LOG>)) {
    chomp $line;
    $already_ran{$line} = 1;
  }
  close LOG;
}

open LOG, ">> $log_file";

my $exit_jobs_file = "exit_jobs_file";

while (defined (my $line = <INPUT>)) {

  if (-e $exit_jobs_file) {
    print STDERR "Job properly ended before all pairwise comparisons were completed\n";
    print STDERR "endtime: ",time,"\n";
    exit 0;
  }

  if ($line =~ /^\S+:\S+::\S+:\S+$/) {
    chomp $line;
  } else {
    die "line input in  $input_file should be dbname:contig_id::dbname:contig_id\n";
  }

  next if (defined $already_ran{$line});

  my $rundb = new Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer (-alnprog => $alnprog,
								     -dbobj => $db,
								     -input_id => $line,
								     -min_score => $min_score);
  # run the runnabledb
  
  print STDERR "$line\n";
  
  $rundb->fetch_input();
  $rundb->run();
#  if ($rundb->write_output()) {
  $rundb->write_output();
    $already_ran{$line} = 1;
    print LOG "$line\n";
#  }
}

close LOG;
close INPUT;

print STDERR "endtime: ",time,"\n";

#!/usr/local/bin/perl

=head1 NAME

  make_bsubs.pl

=head1 SYNOPSIS
 
  make_bsubs.pl

=head1 DESCRIPTION

=head1 OPTIONS
  
  Options are to be set in GeneBuild config files
  The important ones for this script are:

     GeneBuild::Databases::GB_DBNAME   
     GeneBuild::Databases::GB_DBHOST   
     GeneBuild::Databases::GB_DBUSER   
     GeneBuild::Databases::GB_DBPASS   

     GeneBuild::Scripts::GB_RUNNER
     GeneBuild::Scripts::GB_QUEUE
     GeneBuild::Scripts::GB_OUTPUT_DIR
     GeneBuild::Scripts::GB_LENGTH_RUNNABLES
     GeneBuild::Scripts::GB_PM_OUTPUT
     GeneBuild::Scripts::GB_SIZE
							   
     GeneCombinerConf::RUNNER
     GeneCombinerConf::GENECOMBINER_RUNNABLES
     GeneCombinerConf::SLICE_SIZE
     GeneCombinerConf::QUEUE
     GeneCombinerConf::OUTPUT_DIR

=cut

use strict;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (GB_DBNAME
							     GB_DBHOST
							     GB_DBUSER
							     GB_DBPASS
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw (GB_RUNNER
							   GB_QUEUE
							   GB_OUTPUT_DIR
							   GB_LENGTH_RUNNABLES
							   GB_PM_OUTPUT
							   GB_SIZE
							  );

use Bio::EnsEMBL::Pipeline::GeneCombinerConf qw (RUNNER
						 GENECOMBINER_RUNNABLES
						 SLICE_SIZE
						 QUEUE
						 OUTPUT_DIR
						);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

if($GB_DBUSER eq 'ensadmin' && $GB_DBPASS eq ''){
  print "You cannot have dbuser set to ensadmin with no dbpass set!\nPlease correct the entries in GeneBuild config files\n";
  exit(1);
}

foreach my $arg($GB_RUNNER, $GB_DBNAME, $GB_DBHOST, $GB_DBUSER, $GB_QUEUE, $GB_OUTPUT_DIR){
  if ($arg eq '' ){
    print "You need to set various parameters in GeneBuild config files\n" .  
      "Here are your current values for required settings: \n" .
      "runner      => $GB_RUNNER\n" .
      "dbname      => $GB_DBNAME\n" .
      "dbhost      => $GB_DBHOST\n" .
      "dbuser      => $GB_DBUSER\n" .
      "dbpass      => $GB_DBPASS\n" .
      "queue       => $GB_QUEUE\n" .		  
      "output_dir      => $GB_OUTPUT_DIR\n" ;
    
    exit(1);
  }
}

my %chrhash;

&get_chrlengths;

### read runnable and analysis logic names ###

foreach my $length_runnable_list (@{$GB_LENGTH_RUNNABLES}) {
  my $analysis = $length_runnable_list->{analysis};
  my $runnable = $length_runnable_list->{runnable};
  make_lbsubs($runnable, $analysis) unless $runnable eq '';
}

foreach my $genecombiner_runnable ( @{$GENECOMBINER_RUNNABLES} ){
  my $analysis =  $genecombiner_runnable->{analysis};
  my $runnable = $genecombiner_runnable->{runnable};
  make_genecombiner_bsubs( $runnable, $analysis ) unless $runnable eq '';
}



### SUBROUTINES ###

sub get_chrlengths{

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $GB_DBHOST,
					      -user   => $GB_DBUSER,
					      -pass   => $GB_DBPASS,
					      -dbname => $GB_DBNAME,
					     );

  my $q = "SELECT c.name, max(a.chr_end) 
           FROM   chromosome c, assembly a
           WHERE  c.chromosome_id = a.chromosome_id
           GROUP BY c.name";

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  
}

sub make_lbsubs {
  my ($runnable,$analysis_logic_name) = @_;
  
  my $runner      = $GB_RUNNER;
  my $dbname      = $GB_DBNAME;
  my $dbhost      = $GB_DBHOST;
  my $dbuser      = $GB_DBUSER;
  my $dbpass      = $GB_DBPASS;
  my $queue       = $GB_QUEUE;
  my $size        = $GB_SIZE;
  my $dir         = $GB_OUTPUT_DIR . "/$runnable";
  print STDERR "outputdir = ".$dir."\n";

  # check them!
  foreach my $arg($size, $GB_OUTPUT_DIR){
    if ($arg eq '' ){
      print "You need to set various parameters in GeneBuild config files\n" .  
	"Here are your current values for required settings: \n" .
	"GB_SIZE   => $size\n" .
	"GB_OUTPUT_DIR => $GB_OUTPUT_DIR\n";
      
      exit(1);
    }
  }
  

  # ought to do some checking here - does it exist, are there files in it, etc
  system("mkdir $dir") unless opendir(DIR, $dir);
  closedir(DIR);

  my $outf  = "$GB_OUTPUT_DIR/$runnable.jobs.dat";
  open(OUTF, ">$outf") or die "Can't open $outf\n";
  
  foreach my $chr(keys %chrhash) {
    my $length = $chrhash{$chr};

    my $chrdir = $dir . "/$chr";
    
    # needs checks
    system("mkdir $chrdir") unless opendir(DIR, $chrdir);
    closedir(DIR);
    
    my $count = 1;
    
    while ($count < $length) {
      my $start = $count;
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      my $outfile  = $chrdir . "/$input_id.out";
      my $errfile  = $chrdir . "/$input_id.err";
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable  Bio::EnsEMBL::Pipeline::RunnableDB::$runnable -analysis $analysis_logic_name\"";
      $command .= "  $runner ";
      $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable -analysis $analysis_logic_name";
      $command .= " -input_id $input_id ";
      $command .= " -write";      
      print OUTF "$command\n";
      
      $count = $count + $size;
    }
  }
  close OUTF;
}
	     
sub make_genecombiner_bsubs {
  my ($runnable,$analysis) = @_;
  
  my $runner      = $RUNNER;
  my $queue       = $QUEUE;
  my $size        = $SLICE_SIZE;
  my $dir         = $OUTPUT_DIR . "/$runnable";
  print STDERR "outputdir = ".$dir."\n";
  
  # check them!
  foreach my $arg ($size, $OUTPUT_DIR){
    if ($arg eq '' ){
      print "You need to set various parameters in GeneBuild config files\n" .  
	"Here are your current values for required settings: \n" .
	"GB_SIZE       => $SLICE_SIZE\n" .
	"GB_OUTPUT_DIR => $OUTPUT_DIR\n";
      
      exit(1);
    }
  }
  

  # ought to do some checking here - does it exist, are there files in it, etc
  system("mkdir $dir") unless opendir(DIR, $dir);
  closedir(DIR);
  
  my $outf  = "$OUTPUT_DIR/$runnable.jobs.dat";
  open(OUTF, ">$outf") or die "Can't open $outf\n";
  
  foreach my $chr(keys %chrhash) {
    my $length = $chrhash{$chr};
    
    my $chrdir = $dir . "/$chr";
    
    # needs checks
    system("mkdir $chrdir") unless opendir(DIR, $chrdir);
    closedir(DIR);
    
    my $count = 1;
    
    while ($count < $length) {
      my $start = $count;
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      my $outfile  = $chrdir . "/$input_id.out";
      my $errfile  = $chrdir . "/$input_id.err";
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable  Bio::EnsEMBL::Pipeline::RunnableDB::$runnable -analysis $analysis\"";
      $command .= "  $runner ";
      $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable  -analysis $analysis";
      $command .= " -input_id $input_id ";
      $command .= " -write";      
      print OUTF "$command\n";
      
      $count = $count + $size;
    }
  }
  close OUTF;
}

#!/usr/local/bin/perl
use strict;

use Bio::EnsEMBL::Pipeline::GeneConf qw (GB_RUNNER
					 GB_DBNAME
					 GB_DBHOST
					 GB_DBUSER
					 GB_DBPASS
					 GB_QUEUE
					 GB_OUTPUT_DIR
					 GB_LENGTH_RUNNABLES
					 GB_TARGETTED_RUNNABLES
					 GB_PM_OUTPUT
					 GB_SIZE
					);

use Bio::EnsEMBL::Pipeline::GeneCombinerConf qw (RUNNER
						 RUNNABLE
						 SLICE_SIZE
						 QUEUE
						 OUTPUT_DIR
						);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

if($GB_DBUSER eq 'ensadmin' && $GB_DBPASS eq ''){
  print "You cannot have dbuser set to ensadmin with no dbpass set!\nPlease correct the entries in GeneConf.pm\n";
  exit(1);
}

foreach my $arg($GB_RUNNER, $GB_DBNAME, $GB_DBHOST, $GB_DBUSER, $GB_QUEUE, $GB_OUTPUT_DIR){
  if ($arg eq '' ){
    print "You need to set various parameters in GeneConf.pl\n" .  
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

foreach my $other_runnable_list (@{$GB_TARGETTED_RUNNABLES}) {
  my $analysis = $other_runnable_list->{analysis};
  my $runnable = $other_runnable_list->{runnable};
  make_tbsubs($runnable, $analysis) unless $runnable eq '';
}

make_genecombiner_bsubs( $RUNNABLE ) if ($RUNNABLE);



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

sub make_tbsubs {
  my ($runnable,$analysis_logic_name) = @_;
  
  my $runner      = $GB_RUNNER;
  my $dbname      = $GB_DBNAME;
  my $dbhost      = $GB_DBHOST;
  my $dbuser      = $GB_DBUSER;
  my $dbpass      = $GB_DBPASS;
  my $queue       = $GB_QUEUE;
  my $output_dir  = $GB_OUTPUT_DIR . "/$runnable";
  
  my $pm_out      = $GB_PM_OUTPUT;

  $pm_out     .= "pm_best.out";
  # check them!
  foreach my $arg($pm_out){
    if ($arg eq '' ){
      print "You need to set various parameters in GeneConf.pl\n" .  
	"Here are your current values for required settings: \n" .
	"pm_output  => $GB_PM_OUTPUT\n" ;

      exit(1);
    }
  }

  # parse pmatch results
  my %pm_ids;
  
  # set up jobfile & output dirs
  system("mkdir $output_dir") unless opendir(DIR, $output_dir);
  closedir(DIR);
  my $outf  = "$GB_OUTPUT_DIR/$runnable.jobs.dat";

  open(OUTF, ">$outf") or die "Can't open outfile $outf $!\n";

  # generate bsubs, one per protein
  open(PM, "<$pm_out") or die "Can't open pmoutfile $pm_out\n";
  my $tracker = 0;
  my $resdir = $output_dir . "/jobs0";
  
  while(<PM>){
    # do we need to start another output directory? Limit the files in each so we can parse them easily
    if($tracker%100 == 0){
      $resdir = $output_dir . "/jobs" . $tracker/100;
      print STDERR "$resdir\n";
      system("mkdir $resdir") unless opendir(DIR, $resdir);
      closedir(DIR);   
    }
    $tracker++;
    # exact format of input_id varies depending on which runnable we're using ...
    chomp;
    my $input_id;
    if($runnable eq "TargettedGeneWise"){
      $input_id = $_;
    }
    else {
      #  ctg15907:34857075,34858997:NP_005336.2:1,641
      my @cols = split /:/;
      $input_id =  $cols[0] . ":" . $cols[1] . ":" . $cols[2] . ":";
    }

    my $outfile  = $resdir . "/$input_id.out";
    my $errfile  = $resdir . "/$input_id.err";
    my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable -analysis $analysis_logic_name\"";

    $command .= "  $runner ";
    $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable -analysis $analysis_logic_name";
    $command .= " -input_id $input_id -write";      
    print OUTF "$command\n";
    
  }

  close PM;
  close OUTF;
  
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
      print "You need to set various parameters in GeneConf.pl\n" .  
	"Here are your current values for required settings: \n" .
	  "size   => $size\n" .
	    "tmpdir => $GB_OUTPUT_DIR\n";
      
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
  my ($runnable) = @_;
  
  my $runner      = $RUNNER;
  my $queue       = $QUEUE;
  my $size        = $SLICE_SIZE;
  my $dir         = $OUTPUT_DIR . "/$runnable";
  

  # check them!
  foreach my $arg ($size, $OUTPUT_DIR){
    if ($arg eq '' ){
      print "You need to set various parameters in GeneConf.pl\n" .  
	"Here are your current values for required settings: \n" .
	  "size   => $SLICE_SIZE\n" .
	    "tmpdir => $OUTPUT_DIR\n";
      
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
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable  Bio::EnsEMBL::Pipeline::RunnableDB::$runnable \"";
      $command .= "  $runner ";
      $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable ";
      $command .= " -input_id $input_id ";
      $command .= " -write";      
      print OUTF "$command\n";
      
      $count = $count + $size;
    }
  }
  close OUTF;
}

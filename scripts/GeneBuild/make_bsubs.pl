#!/usr/local/bin/perl
use strict;

use Bio::EnsEMBL::Pipeline::GeneConf qw (GB_RUNNER
					 GB_DBNAME
					 GB_DBHOST
					 GB_DBUSER
					 GB_DBPASS
					 GB_GOLDEN_PATH
					 GB_QUEUE
					 GB_TMPDIR
					 GB_LENGTH_RUNNABLES
					 GB_TARGETTED_RUNNABLES
					 GB_PM_OUTPUT
					 GB_SIZE
					);

use Bio::EnsEMBL::DBSQL::DBAdaptor;

if($GB_DBUSER eq 'ensadmin' && $GB_DBPASS eq ''){
  print "You cannot have dbuser set to ensadmin with no dbpass set!\nPlease correct the entries in GeneConf.pm\n";
  exit(1);
}

foreach my $arg($GB_RUNNER, $GB_DBNAME, $GB_DBHOST, $GB_DBUSER, $GB_QUEUE, $GB_TMPDIR){
    if ($arg eq '' ){
      print "You need to set various parameters in GB_conf.pl\n" .  
	"Here are your current values for required settings: \n" .
	"runner      => $GB_RUNNER\n" .
	"dbname      => $GB_DBNAME\n" .
	"dbhost      => $GB_DBHOST\n" .
	"dbuser      => $GB_DBUSER\n" .
	"dbpass      => $GB_DBPASS\n" .
	"queue       => $GB_QUEUE\n" .
	"tmpdir      => $GB_TMPDIR\n" .
	"golden_path => $GB_GOLDEN_PATH ( empty string will use UCSC )\n" ;

      exit(1);
    }
  }

my %chrhash;

&get_chrlengths;

foreach my $lr(@{$GB_LENGTH_RUNNABLES}) {
  make_lbsubs($lr) unless $lr eq '';
}

foreach my $tr(@{$GB_TARGETTED_RUNNABLES}) {
  make_tbsubs($tr) unless $tr eq '';
}


### SUBROUTINES ###

sub get_chrlengths{

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $GB_DBHOST,
					      -user   => $GB_DBUSER,
					      -pass   => $GB_DBPASS,
					      -dbname => $GB_DBNAME,
					     );

  my $q = "SELECT chr_name,max(chr_end) FROM static_golden_path GROUP BY chr_name";

  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  
}

sub make_tbsubs {
  my ($runnable) = @_;
  
  my $runner      = $GB_RUNNER;
  my $dbname      = $GB_DBNAME;
  my $dbhost      = $GB_DBHOST;
  my $dbuser      = $GB_DBUSER;
  my $dbpass      = $GB_DBPASS;
  my $queue       = $GB_QUEUE;
  my $golden_path = $GB_GOLDEN_PATH;
  my $dir         = $GB_TMPDIR . "/$runnable";
  my $pm_out      = $GB_PM_OUTPUT;

  $pm_out     .= "pm_best.out";
  $golden_path = 'UCSC' unless (defined $golden_path && $golden_path ne '');
  # check them!
  foreach my $arg($pm_out){
    if ($arg eq '' ){
      print "You need to set various parameters in GB_conf.pl\n" .  
	"Here are your current values for required settings: \n" .
	"pm_output  => $GB_PM_OUTPUT\n" ;

      exit(1);
    }
  }

  # parse pmatch results
  my %pm_ids;
  
  # set up jobfile & output dirs
  system("mkdir $dir") unless opendir(DIR, $dir);
  closedir(DIR);
  my $outf  = "$runnable.jobs.dat";

  open(OUTF, ">$outf") or die "Can't open outfile $outf\n";

  # generate bsubs, one per protein
  open(PM, "<$pm_out") or die "Can't open pmoutfile $pm_out\n";
  my $tracker = 0;
  my $resdir = $dir . "/jobs0";
   
  while(<PM>){
    # do we need to start another output directory? Limit the files in each so we can parse them easily
    if($tracker%100 == 0){
      $resdir = $dir . "/jobs" . $tracker/100;
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
    my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable\"";

    $command .= "  $runner ";
    $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable ";
    $command .= " -input_id $input_id -write";      
    print OUTF "$command\n";
    
  }

  close PM;
  close OUTF;
  
}


sub make_lbsubs {
  my ($runnable) = @_;
  
  my $runner      = $GB_RUNNER;
  my $dbname      = $GB_DBNAME;
  my $dbhost      = $GB_DBHOST;
  my $dbuser      = $GB_DBUSER;
  my $dbpass      = $GB_DBPASS;
  my $golden_path = $GB_GOLDEN_PATH;
  my $queue       = $GB_QUEUE;
  my $size        = $GB_SIZE;
  my $dir         = $GB_TMPDIR . "/$runnable";
  
  $golden_path = 'UCSC' unless (defined $golden_path && $golden_path ne '');

  # check them!
  foreach my $arg($size, $GB_TMPDIR){
    if ($arg eq '' ){
      print "You need to set various parameters in GeneConf.pl\n" .  
	"Here are your current values for required settings: \n" .
	"size   => $size\n" .
	"tmpdir => $GB_TMPDIR\n";

      exit(1);
    }
  }


  # ought to do some checking here - does it exist, are there files in it, etc
  system("mkdir $dir") unless opendir(DIR, $dir);
  closedir(DIR);

  my $outf  = "$runnable.jobs.dat";
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
	     

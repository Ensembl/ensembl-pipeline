#!/usr/local/bin/perl
use strict;

require "Bio/EnsEMBL/Pipeline/GB_conf.pl";

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my %conf    = %::scripts_conf;
my %db_conf = %::db_conf;

if($db_conf{'dbuser'} eq 'ensadmin' && $db_conf{'dbpass'} eq ''){
  print "You cannot have dbuser set to ensadmin with no dbpass set!\nPlease correct the entries in GB_conf.pl\n";
  exit(1);
}

foreach my $arg($conf{'runner'}, $db_conf{'dbname'}, $db_conf{'dbhost'}, $db_conf{'dbuser'}, $conf{'queue'},$conf{'tmpdir'}){
    if ($arg eq '' ){
      print "You need to set various parameters in GB_conf.pl\n" .  
	"Here are your current values for required settings: \n" .
	"runner     => $conf{'runner'}\n" .
	"dbname     => $db_conf{'dbname'}\n" .
	"dbhost     => $db_conf{'dbhost'}\n" .
	"dbuser     => $db_conf{'dbuser'}\n" .
	"dbpass     => $db_conf{'dbpass'}\n" .
	"queue      => $conf{'queue'}\n" .
	"tmpdir     => $conf{'tmpdir'}\n" .
	"golden_path => $db_conf{'golden_path'} ( empty string will use UCSC )\n" ;

      exit(1);
    }
  }

my %chrhash;

&get_chrlengths;

foreach my $lr(@{$conf{'length_runnables'}}) {
  make_lbsubs($lr) unless $lr eq '';
}

foreach my $tr(@{$conf{'targetted_runnables'}}) {
  make_tbsubs($tr) unless $tr eq '';
}


### SUBROUTINES ###

sub get_chrlengths{

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $db_conf{'dbhost'},
					      -user   => $db_conf{'dbuser'},
					      -pass   => $db_conf{'dbpass'},
					      -dbname => $db_conf{'dbname'},
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
  
  my $runner      = $conf{'runner'};
  my $dbname      = $db_conf{'dbname'};
  my $dbhost      = $db_conf{'dbhost'};
  my $dbuser      = $db_conf{'dbuser'};
  my $dbpass      = $db_conf{'dbpass'};
  my $queue       = $conf{'queue'};
  my $golden_path = $db_conf{'golden_path'};
  my $dir         = $conf{'tmpdir'} . "/$runnable";
  my $pm_out      = $conf{'pm_output'};
  my $cdnas       = $conf{'cdna_pairs'};

  $pm_out .= "pm_best.out";
  $golden_path = 'UCSC' unless (defined $golden_path && $golden_path ne '');
  # check them!
  foreach my $arg($pm_out, $cdnas){
    if ($arg eq '' ){
      print "You need to set various parameters in GB_conf.pl\n" .  
	"Here are your current values for required settings: \n" .
	"pm_output  => $conf{'pm_output'}\n" .
	"cdna_pairs => $cdnas\n";

      exit(1);
    }
  }

  # parse pmatch and cdna results
  my %pm_ids;
  my %cdnas;
  
  open(CDNA, "<$cdnas") or die "Can't open cdnafile $cdnas\n";

  while(<CDNA>){
    #  NP_002923.1 : L26953
    if(/(\S+)\s+:\s+(\S+)/){
      if(!defined($cdnas{$1})){ $cdnas{$1} = $2; }
      else { die "already seen $1 with $cdnas{$1}\n"  }
    }
  }
  close CDNA;

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

      if($cdnas{$cols[2]}){
	$input_id .= $cdnas{$cols[2]};      
      }
    }

    my $outfile  = $resdir . "/$input_id.out";
    my $errfile  = $resdir . "/$input_id.err";
    my $command = "bsub -q $queue -o $outfile -e $errfile -E \"$runner -check -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable\"";
    $command .= "  $runner ";
    $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable ";
    $command .= " -dbuser $dbuser -pass $dbpass -dbname $dbname -host $dbhost ";
    $command .= " -input_id $input_id -parameters golden_path=$golden_path -write";      
    print OUTF "$command\n";
    
  }

  close PM;
  close OUTF;
  
}


sub make_lbsubs {
  my ($runnable) = @_;
  
  my $runner      = $conf{'runner'};
  my $dbname      = $db_conf{'dbname'};
  my $dbhost      = $db_conf{'dbhost'};
  my $dbuser      = $db_conf{'dbuser'};
  my $dbpass      = $db_conf{'dbpass'};
  my $golden_path = $db_conf{'golden_path'};
  my $queue       = $conf{'queue'};
  my $size        = $conf{'size'};
  my $dir         = $conf{'tmpdir'} . "/$runnable";
  
  $golden_path = 'UCSC' unless (defined $golden_path && $golden_path ne '');

  # check them!
  foreach my $arg($size, $conf{'tmpdir'}){
    if ($arg eq '' ){
      print "You need to set various parameters in GB_conf.pl\n" .  
	"Here are your current values for required settings: \n" .
	"size   => $size\n" .
	"tmpdir => $conf{'tmpdir'}\n";

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
      my $command = "bsub -q $queue -o $outfile -e $errfile -E \"$runner -check -runnable  Bio::EnsEMBL::Pipeline::RunnableDB::$runnable \"";
      $command .= "  $runner ";
      $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable ";
      $command .= " -dbuser $dbuser -pass $dbpass -dbname $dbname -host $dbhost ";
      $command .= " -input_id $input_id -parameters golden_path=$golden_path -write";      
      print OUTF "$command\n";
      
      $count = $count + $size;
    }
  }
  close OUTF;
}
	     

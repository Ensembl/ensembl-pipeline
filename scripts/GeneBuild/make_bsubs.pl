#!/usr/local/bin/perl -w
use strict;

BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
  require "GB_conf.pl";
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my %gbc = %::GB_conf;

my %chrhash;

&get_chrlengths;

foreach my $lr(@{$gbc{'length_runnables'}}) {
  make_lbsubs($lr);
}

foreach my $tr(@{$gbc{'targetted_runnables'}}) {
  make_tbsubs($tr);
}


### SUBROUTINES ###

sub get_chrlengths{
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $gbc{'dbhost'},
					      -user   => $gbc{'dbuser'},
					      -dbname => $gbc{'dbname'},
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
  
  my $runner = $gbc{'runner'};
  my $dbname = $gbc{'dbname'};
  my $dbhost = $gbc{'dbhost'};
  my $dbuser = $gbc{'dbuser'};
  my $queue  = $gbc{'queue'};
  my $dir    = $gbc{'tmpdir'} . "/$runnable";
  my $pm_out = $gbc{'pm_output'};
  $pm_out .= "pm_best.out";
  my $cdnas  = $gbc{'cdna_pairs'};

  # parse pmatch and cdna results
  my %pm_ids;
  my %cdnas;
  
  open(CDNA, "<$cdnas") or die "Can't open $cdnas\n";

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

  open(OUTF, ">$outf") or die "Can't open $outf\n";

  # generate bsubs, one per protein
  open(PM, "<$pm_out") or die "Can't open $pm_out\n";
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
    my $command = "bsub -q $queue -o $outfile -e $errfile -E \"$runner -check \"";
    $command .= "  $runner ";
    $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable ";
    $command .= " -dbuser $dbuser -dbname $dbname -host $dbhost ";
    $command .= " -input_id $input_id -write";      
    print OUTF "$command\n";
    
  }

  close PM;
  close OUTF;
  
}


sub make_lbsubs {
  my ($runnable) = @_;
  
  my $runner = $gbc{'runner'};
  my $dbname = $gbc{'dbname'};
  my $dbhost = $gbc{'dbhost'};
  my $dbuser = $gbc{'dbuser'};
  my $queue  = $gbc{'queue'};
  my $size   = $gbc{'size'};
  my $dir    = $gbc{'tmpdir'} . "/$runnable";
  
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
      my $command = "bsub -q $queue -o $outfile -e $errfile -E \"$runner -check \"";
      $command .= "  $runner ";
      $command .= " -runnable Bio::EnsEMBL::Pipeline::RunnableDB::$runnable ";
      $command .= " -dbuser $dbuser -dbname $dbname -host $dbhost ";
      $command .= " -input_id $input_id -write";      
      print OUTF "$command\n";
      
      $count = $count + $size;
    }
  }
  close OUTF;
}
	       
		  

#!/usr/local/bin/perl -w

=head1 NAME

  make_bsubs.pl

=head1 SYNOPSIS
 
  make_bsubs.pl
  Makes bsub entries for exonerate_ests.pl and filter_and_e2g.pl
  bsubs can be submitted using submit.pl - they\'re not automatically 
  done from here as it\'s better to submit a few and check they come 
  back OK before sending a genome worth.

  Makes sure all the various scratch subdirectories needed are in place, 
  and makes them if necessary.

=head1 DESCRIPTION


=head1 OPTIONS

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
require "Bio/EnsEMBL/Pipeline/EST_conf.pl";

my %chrhash;

# declare these here so we can refer to them later
my $ex_resultsdir       = "exonerate_est/results";
my $ex_bsubdir          = "exonerate_est/bsub";
my $filterdir           = "filter_and_e2g";
my $est_genebuilder_dir = "est_genebuilder";

# get chr info from the database where you have contig, clone and static_golden_path
&get_chrlengths();

# make output directories
&make_directories();

# create jobs file for Exonerate
&make_exonerate_bsubs();

# create jobs file for FilterESTs_and_E2G
&make_filter_bsubs();

# create jobs file for EST_GeneBuilder
&make_EST_GeneBuilder_bsubs();

=head2 make_directories

  Title   : make_directories
  Usage   : make_directories
  Function: makes sure needed output directories exist, and if not, makes them
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_directories {
  my $scratchdir =  $::scripts_conf{'tmpdir'};

  my @resdirs = split /\//, $ex_resultsdir;

  # exonerate_ests
  my $exoneratedir = $scratchdir . "/" . $resdirs[0] . "/";
  makedir($exoneratedir);

  # exonerate output directories
  my $exdir = $exoneratedir . $resdirs[1] . "/";  
  my $exerr = $exdir . "stderr/";
  my $exout = $exdir . "stdout/";
  makedir($exdir);
  makedir($exerr);
  makedir($exout);
  
  # bsub output directories
  my $bsubdir = $scratchdir . "/" . $ex_bsubdir . "/";
  my $bsuberr = $bsubdir . "stderr/";
  my $bsubout = $bsubdir . "stdout/";
  makedir($bsubdir);
  makedir($bsuberr);
  makedir($bsubout);

  # filter_and_e2g
  my $filter = $scratchdir . "/" . $filterdir . "/";
  makedir($filter);
  
  foreach my $chr(keys %chrhash){
    my $chrdir = $filter . $chr . "/";
    makedir($chrdir);
  }

  # EST_GeneBuilder output directories
  my $estbuilder_path = $scratchdir . "/" . $est_genebuilder_dir . "/";
  makedir($estbuilder_path);
  
  foreach my $chr(keys %chrhash){
    my $chrdir = $estbuilder_path . $chr . "/";
    makedir($chrdir);
  }

}

=head2 get_chrlengths

  Title   : get_chrlengths
  Usage   : get_chrlengths
  Function: gets lengths of all chromosomes from reference database
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub get_chrlengths{
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $::db_conf{'refdbhost'},
					      -user   => $::db_conf{'refdbuser'},
					      -dbname => $::db_conf{'refdbname'},
					     );
  my $q = "SELECT chr_name,max(chr_end) FROM static_golden_path GROUP BY chr_name";
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  
}

=head2 make_exonerate_bsubs

  Title   : make_exonerate_bsubs
  Usage   : make_exonerate_bsubs
  Function: makes bsubs to run exonerate_est.pl
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_exonerate_bsubs {
  my $jobfile = $::scripts_conf{'exonerate_bsubsfile'};
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  my $queue         = $::scripts_conf{'queue'};
  my $scriptdir     = $::scripts_conf{'scriptdir'};
  my $check         = $scriptdir . "/check_node.pl";
  my $exonerate_est = $scriptdir . "/exonerate_est.pl";
  my $bsuberr       = $::scripts_conf{'tmpdir'} . "/" . $ex_bsubdir . "/stderr/";
  my $bsubout       = $::scripts_conf{'tmpdir'} . "/" . $ex_bsubdir . "/stdout/";
  
  my $estfile = $::scripts_conf{'estfile'}; # may be a full path
  my @path = split /\//, $estfile;
  $estfile = $path[$#path];
  $estfile .= "_chunk_";

  my $numchunks = $::scripts_conf{'estchunknumber'};

  for(my $i = 0; $i < $numchunks; $i++){
    my $num = $i;
    while (length($num) < 7){
      $num = "0" . $num;
    }
    
    my $chunk     = $estfile . $num;
    my $outfile   = $bsubout . $chunk;
    my $errfile   = $bsuberr . $chunk;
    
    my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$check $chunk\" $exonerate_est -chunkname $chunk";
    print OUT "$command\n";
  }
  
  close (OUT) or die (" Error closing $jobfile: $!");
}

=head2 make_filter_bsubs

  Title   : make_filter_bsubs
  Usage   : make_filter_bsubs
  Function: makes bsubs to run filter_and_e2g.pl
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_filter_bsubs {
  my $jobfile    = $::scripts_conf{'filter_bsubsfile'};
  my $scratchdir = $::scripts_conf{'tmpdir'};
  my $queue      = $::scripts_conf{'queue'};
  my $scriptdir  = $::scripts_conf{'scriptdir'};
  my $filter_e2g = $scriptdir . "/filter_and_e2g.pl";

  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  my $filter = $scratchdir . "/" . $filterdir . "/";
  my $size   = $::scripts_conf{'filter_chunksize'};
  my $runner = $::scripts_conf{'runner'};

  foreach my $chr(keys %chrhash) {
    my $length = $chrhash{$chr};
    
    my $chrdir = $filter . "/$chr";
    my $count = 1;

    while($count < $length) {
      my $start = $count;
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      my $outfile  = $chrdir . "/$input_id.out";
      my $errfile  = $chrdir . "/$input_id.err";
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G\" $filter_e2g -input_id $input_id";
      print OUT "$command\n";
      
      $count = $count + $size;
    }
  }

  close (OUT) or die (" Error closing $jobfile: $!");
}



=head2 make_EST_GeneBuilder_bsubs

  Title   : make_EST_GeneBuilder_bsubs
  Function: makes bsubs to run EST_GeneBuilder

=cut

sub make_EST_GeneBuilder_bsubs{
  my $jobfile    = $::scripts_conf{'EST_GeneBuilder_bsubsfile'};
  my $scratchdir = $::scripts_conf{'tmpdir'};
  my $queue      = $::scripts_conf{'queue'};
  my $scriptdir  = $::scripts_conf{'scriptdir'};

  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  # where out and err files go
  my $filter = $scratchdir . "/" . $est_genebuilder_dir . "/";

  # genomic size for each job
  my $size   = $::scripts_conf{'est_genebuilder_chunksize'};
  
  my $runner = $::scripts_conf{'runner'};

  foreach my $chr(keys %chrhash) {
    my $length = $chrhash{$chr};
    
    my $chrdir = $filter . "/$chr";
    my $count = 1;

    while($count < $length) {
      my $start = $count;
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      my $outfile  = $chrdir . "/$input_id.out";
      my $errfile  = $chrdir . "/$input_id.err";

      # if you don't want it to write to the database, eliminate the -write option
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder\" $runner -runnable Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder -input_id $input_id -write";
      print OUT "$command\n";
      
      $count = $count + $size;
    }
  }

  close (OUT) or die (" Error closing $jobfile: $!");
}

=head2 makedir

  Title   : makedir
  Usage   : makedir
  Function: checks to see if the given directory exists, and makes it if it doesn\'t
  Returns : nothing, but will kill script if mkdir fails
  Args    : $dir - directory to be created

=cut

sub makedir{
  my ($dir) = @_;
  if(opendir(DIR, $dir)){ closedir(DIR); }
  else{ system("mkdir $dir") == 0 or die "error creating $dir\n"; }
}

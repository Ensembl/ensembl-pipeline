#!/usr/local/bin/perl -w

BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  use lib $script_dir;
  require "EST_conf.pl";
}

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

my %conf =  %::EST_conf; # from EST_conf.pl
my %chrhash;

# declare these here so we can refer to them later
my $ex_resultsdir = "exonerate_est/results";
my $ex_bsubdir    = "exonerate_est/bsub";
my $filterdir     = "filter_and_e2g";

&get_chrlengths();
&make_directories();
#&make_exonerate_bsubs();
&make_filter_bsubs();


=head2 make_directories

  Title   : make_directories
  Usage   : make_directories
  Function: makes sure needed output directories exist, and if not, makes them
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_directories {
  my $scratchdir =  $conf{'tmpdir'};

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

}

=head2 get_chrlengths

  Title   : get_chrlengths
  Usage   : get_chrlengths
  Function: gets lengths of all chromosomes from reference database
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub get_chrlengths{
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $conf{'refdbhost'},
					      -user   => $conf{'refdbuser'},
					      -dbname => $conf{'refdbname'},
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
  my $jobfile = $conf{'exonerate_bsubsfile'};
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  my $queue   = $conf{'queue'};
  my $check   = "check_node.pl";
  my $bsuberr = $conf{'tmpdir'} . "/" . $ex_bsubdir . "/stderr/";
  my $bsubout = $conf{'tmpdir'} . "/" . $ex_bsubdir . "/stdout/";

  my $estfile = $conf{'estfile'}; # may be a full path
  my @path = split /\//, $estfile;
  $estfile = $path[$#path];
  $estfile .= "_chunk_";

  my $numchunks = $conf{'estchunknumber'};

  for(my $i = 0; $i < $numchunks; $i++){
    my $num = $i;
    while (length($num) < 7){
      $num = "0" . $num;
    }
    
    my $chunk     = $estfile . $num;
    my $outfile   = $bsubout . $chunk;
    my $errfile   = $bsuberr . $chunk;
    my $chunkfile = $conf{'estfiledir'} . "/" . $chunk;
    
    my $command = "bsub -q $queue -o $outfile -e $errfile -E \"$check $chunk\" exonerate_est.pl -chunkname $chunkfile";
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
  my $jobfile = $conf{'filter_bsubsfile'};
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  my $filter = $scratchdir . "/" . $filterdir . "/";
  my $size = $conf{'filter_chunksize'};
  my $runner = $conf{'runner'};

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
      my $command = "bsub -q $queue -o $outfile -e $errfile -E \"$runner -check \" filter_and_e2g.pl -input_id $input_id"
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

#!/usr/local/bin/perl -w

=head1 NAME

  make_firstef_bsubs.pl

=head1 SYNOPSIS
 
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
use Bio::EnsEMBL::Pipeline::Config::FirstEF qw (
						FEF_TMPDIR
						FEF_REFDBHOST
						FEF_REFDBUSER
						FEF_REFDBNAME
						FEF_BSUB_FILE
						FEF_QUEUE
						FEF_CHUNKSIZE
						FEF_RUN_SCRIPT
					       );

my %chrhash;

# get chr info from the database
&get_chrlengths();

# make output directories
&make_output_directories();

# create jobs file for EST_GeneBuilder
&make_firstef_bsubs();


=head2 make_directories

  Title   : make_directories
  Usage   : make_directories
  Function: makes sure needed output directories exist, and if not, makes them
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_output_directories {

  makedir($FEF_TMPDIR . "/firstef/");
  
  foreach my $chr(keys %chrhash){
    makedir($FEF_TMPDIR . "/firstef/" . $chr . "/");
  }
  
}


=head2 get_chrlengths

  Function: gets lengths of all chromosomes from reference database

=cut

sub get_chrlengths{
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $FEF_REFDBHOST,
					      -user   => $FEF_REFDBUSER,
					      -dbname => $FEF_REFDBNAME,
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


=head2 make_firstef_bsubs

  Function: makes bsubs to run firstef

=cut

sub make_firstef_bsubs{

  open (OUT, ">$FEF_BSUB_FILE") or die ("Can't open $FEF_BSUB_FILE for writing bsub commands: $!");

  # Name of the script that runs the analysis for each slice.
  my $runner   = $FEF_RUN_SCRIPT;

  foreach my $chr (keys %chrhash) {
    my $length = $chrhash{$chr};
    
    my $chrdir = $FEF_TMPDIR . "/firstef/" . $chr;
    my $count = 1;

    while($count < $length) {
      my $start = $count;
      my $end   = $count + $FEF_CHUNKSIZE -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      my $outfile  = $chrdir . "/$input_id.out";
      my $errfile  = $chrdir . "/$input_id.err";
     

      my $command = "bsub -q $FEF_QUEUE -C0 -o $outfile -e $errfile $FEF_RUN_SCRIPT -input_id $input_id -write";
      print OUT "$command\n";
      
      $count = $count + $FEF_CHUNKSIZE;
    }
  }

  close (OUT);
}



=head2 makedir

  Title   : makedir
  Usage   : makedir
  Function: Creates a directory if it doesnt already exist.
  Returns : none
  Args    : directory path

=cut

sub makedir{
  my ($dir) = @_;
  unless(-e $dir && -d $dir){
    system("mkdir $dir") == 0 or die "Error creating $dir\n" 
  }
}

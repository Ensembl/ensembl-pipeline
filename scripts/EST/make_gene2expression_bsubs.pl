#!/usr/local/bin/perl -w

=head1 NAME

  make_gene2expression_bsubs.pl

=head1 SYNOPSIS
 
  make_gene2expression_bsubs.pl
  Makes bsub entries for MapGenesToExpression
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
use Bio::EnsEMBL::Pipeline::Config::cDNAs_ESTs::GenesToExpression qw (
								      EST_TMPDIR
								      EST_REFDBHOST
								      EST_REFDBUSER
								      EST_REFDBNAME
                                                                      EST_REFDBPORT
								      EST_QUEUE
								      EST_TMPDIR
								      EST_EXPRESSION_RUNNER
								      EST_EXPRESSION_CHUNKSIZE
								      EST_EXPRESSION_BSUBS
								      EST_EXPRESSION_RUNNABLE
								      EST_EXPRESSION_ANALYSIS    
								     );

my %chrhash;

# declare these here so we can refer to them later
my $gene2expression_dir   = "gene2ests";

# get chr info from the database where you have contig, clone and static_golden_path
&get_chrlengths();

# make output directories
&make_directories();

# create jobs file MapGeneToExpression 
&make_MapGeneToExpression_bsubs();

=head2 make_directories

  Title   : make_directories
  Usage   : make_directories
  Function: makes sure needed output directories exist, and if not, makes them
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_directories {
  my $scratchdir =  $EST_TMPDIR;

  # MapGeneToExpression directories
  my $gene2expression_path = $scratchdir . "/" . $gene2expression_dir . "/";
  makedir($gene2expression_path);
  
  foreach my $chr (keys %chrhash){
    my $chrdir = $gene2expression_path . $chr . "/";
    makedir($chrdir);
  }
}

############################################################

=head2 get_chrlengths

  Function: gets lengths of all chromosomes from reference database

=cut

sub get_chrlengths{
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $EST_REFDBHOST,
					      -user   => $EST_REFDBUSER,
					      -dbname => $EST_REFDBNAME,
                                              -port   => $EST_REFDBPORT,
					     );
#  my $q = "SELECT c.name, max(a.chr_end) 
#           FROM   chromosome c, assembly a
#           WHERE  c.chromosome_id = a.chromosome_id
#           GROUP BY c.name";
 
  my $q = "select distinct sr.name,sr.length 
           from seq_region_attrib sra, 
                attrib_type at,
                seq_region sr 
           where sra.seq_region_id = sr.seq_region_id and 
                 sra.attrib_type_id = at.attrib_type_id and 
                 at.code = \'toplevel\'";


  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  
}

=head2 make_MapGeneToExpression_bsubs

    Function: makes bsubs to run MapGeneToExpression

=cut

sub make_MapGeneToExpression_bsubs{
  my $jobfile    = $EST_EXPRESSION_BSUBS;
  my $scratchdir = $EST_TMPDIR;
  my $queue      = $EST_QUEUE;
  my $script     = $EST_EXPRESSION_RUNNER;
  
  my $runnable   = $EST_EXPRESSION_RUNNABLE;
  my $analysis   = $EST_EXPRESSION_ANALYSIS;
  
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");
  
  # where out and err files go
  my $filter = $scratchdir . "/" .  $gene2expression_dir . "/";
  
  # genomic size for each job
  my $size   = $EST_EXPRESSION_CHUNKSIZE;

  foreach my $chr(keys %chrhash) {
    my $length = $chrhash{$chr};
    
    my $chrdir = $filter . "$chr";
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
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$script -check -runnable $runnable -analysis $analysis\" $script -runnable $runnable -analysis $analysis -input_id $input_id -write";
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

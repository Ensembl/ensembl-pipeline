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
use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_TMPDIR
					EST_REFDBHOST
					EST_REFDBUSER
					EST_REFDBNAME
					EST_EXONERATE_BSUBS
					EST_FILTER_BSUBS
					EST_QUEUE
					EST_SCRIPTDIR
					EST_CHUNKNUMBER
					EST_FILTERCHUNKSIZE
					EST_GENEBUILDER_BSUBS
					EST_GENEBUILDER_CHUNKSIZE
					EST_FILE
					EST_RUNNER
                                        EST_GENE_RUNNER
					MAP_GENES_TO_ESTS
					EST_EXPRESSION_RUNNER
					EST_EXPRESSION_CHUNKSIZE
					EST_EXPRESSION_BSUBS
				       );



my %chrhash;

# declare these here so we can refer to them later
my $ex_resultsdir       = "exonerate_est/results";
my $ex_bsubdir          = "exonerate_est/bsub";
my $filterdir           = "filter_and_e2g";
my $est_genebuilder_dir = "est_genebuilder";
my $gene2expression_dir = "gene2ests";

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

if ( $MAP_GENES_TO_ESTS ){
  &make_MapGeneToExpression_bsubs();
}

=head2 make_directories

  Title   : make_directories
  Usage   : make_directories
  Function: makes sure needed output directories exist, and if not, makes them
  Returns : none - uses globals
  Args    : none - uses globals

=cut

sub make_directories {
  my $scratchdir =  $EST_TMPDIR ;

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

  # MapGeneToExpression directories
  if ( $MAP_GENES_TO_ESTS ){
    my $gene2expression_path = $scratchdir . "/" . $gene2expression_dir . "/";
    makedir($gene2expression_path);
    
    foreach my $chr (keys %chrhash){
      my $chrdir = $gene2expression_path . $chr . "/";
      makedir($chrdir);
    }
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
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $EST_REFDBHOST,
					      -user   => $EST_REFDBUSER,
					      -dbname => $EST_REFDBNAME,
					     );
  my $q = "SELECT chr_name,max(chr_end) FROM static_golden_path GROUP BY chr_name";
  
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
    print STDERR "$chr -> $length\n";
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
  my $jobfile = $EST_EXONERATE_BSUBS;
  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  my $queue         = $EST_QUEUE;
  my $scriptdir     = $EST_SCRIPTDIR;
  my $check         = $scriptdir . "/check_node.pl";
  my $exonerate_est = $scriptdir . "/exonerate_est.pl";
  my $bsuberr       = $EST_TMPDIR . "/" . $ex_bsubdir . "/stderr/";
  my $bsubout       = $EST_TMPDIR . "/" . $ex_bsubdir . "/stdout/";
  
  my $estfile = $EST_FILE; # may be a full path
  my @path = split /\//, $estfile;
  $estfile = $path[$#path];
  $estfile .= "_chunk_";

  my $numchunks = $EST_CHUNKNUMBER;

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
  my $jobfile    = $EST_FILTER_BSUBS;
  my $scratchdir = $EST_TMPDIR;
  my $queue      = $EST_QUEUE;
  my $scriptdir  = $EST_SCRIPTDIR;
  my $filter_e2g = $scriptdir . "/filter_and_e2g.pl";

  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  my $filter = $scratchdir . "/" . $filterdir . "/";
  my $size   = $EST_FILTERCHUNKSIZE;
  my $runner = $EST_RUNNER;

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
  my $jobfile    = $EST_GENEBUILDER_BSUBS;
  my $scratchdir = $EST_TMPDIR;
  my $queue      = $EST_QUEUE;
  my $scriptdir  = $EST_SCRIPTDIR;

  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  # where out and err files go
  my $filter = $scratchdir . "/" . $est_genebuilder_dir. "/";

  # genomic size for each job
  my $size   = $EST_GENEBUILDER_CHUNKSIZE;
  
  my $runner = $EST_GENE_RUNNER;

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
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder\" $runner -runnable Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder -input_id $input_id -write";
      print OUT "$command\n";
      
      $count = $count + $size;
    }
  }

  close (OUT) or die (" Error closing $jobfile: $!");
}


=head2 make_MapGeneToExpression_bsubs

    Function: makes bsubs to run MapGeneToExpression

=cut

sub make_MapGeneToExpression_bsubs{
  my $jobfile    = $EST_EXPRESSION_BSUBS;
  my $scratchdir = $EST_TMPDIR;
  my $queue      = $EST_QUEUE;
  my $scriptdir  = $EST_SCRIPTDIR;

  open (OUT, ">$jobfile") or die ("Can't open $jobfile for writing: $!");

  # where out and err files go
  my $filter = $scratchdir . "/" .  $gene2expression_dir . "/";

  # genomic size for each job
  my $size   = $EST_EXPRESSION_CHUNKSIZE;
  my $runner = $EST_EXPRESSION_RUNNER;

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
      my $command = "bsub -q $queue -C0 -o $outfile -e $errfile -E \"$runner -check -runnable Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression\" $runner -runnable Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression -input_id $input_id -write";
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

# Cared for by Steve Searle  <searle@sanger.ac.uk>
#
# Copyright Steve Searle 
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Hsperate

=head1 SYNOPSIS

  # To run a hsperate job from scratch do the following.


  my $hsperate =  Bio::EnsEMBL::Pipeline::Runnable::Hsperate->new 
    ('-queries'   => $queries,
     '-program'   => 'hsperate',
     '-database'  => 'swall',
     '-threshold' => 150,
     '-filter'    => $filter,
     '-options'   => '-l 5');

  $hsperate->run();

  @featurepairs = $hsperate->output();

  foreach my $fp (@featurepairs) {
      print $fp->gffstring . "\n";
  }

  # Additionally if you have hsperate runs lying around that need parsing
  # you can do

  open(HSPERATE,"<hsperate.out");

  my @featurepairs = Bio::EnsEMBL::Pipeline::Runnable::Hsperate->parse_results(\*HSPERATE);

  close(HSPERATE);

=head1 DESCRIPTION

Hsperate takes a Bio::Seq (or Bio::PrimarySeq) object and runs hsperate with, the
output is parsed by BPLite and stored as Bio::EnsEMBL::FeaturePairs. 

The output features can be filtered by probability using the -threshold option.
Other options can be passed to the hsperate program using the -options method

=head1 CONTACT

Describe contact details here

=head1 APPENDIX


=cut

package Bio::EnsEMBL::Pipeline::Runnable::Hsperate;


use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::RootI;

use FileHandle;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter;
use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::Runnable::HsperateParser;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Hsperate->new 
    (-query    => $seq,
     -program  => 'hsperate',
     -database => 'swall',
     -threshold => 150,
     -filter    => $filter,
     -options   => '-l 5');

    Function:   Initialises Hsperate object
    Returns :   a Hsperate Object
    Args    :   A Bio::Seq object
                The hsperate executable (-HSPERATE) and database (-DB).

=cut

sub new {
  my ($class,@args) = @_;

  my $self = $class->SUPER::new(@args);    

  $self->{'_queries'}   = undef;     # ref to ARRAY of Bio::Seq objects
  $self->{'_program'}   = undef;     # location of Hsperate
  $self->{'_database'}  = undef;     # name of database
  $self->{'_threshold'} = undef;     # Threshold for hit filterting
  $self->{'_options'}   = undef;     # arguments for hsperate
  $self->{'_filter'}    = 0;         # Do we filter features?
  $self->{'_fplist'}    = [];        # an array of feature pairs (the output)

  $self->{'_workdir'}   = undef;     # location of temp directory
  $self->{'_filename'}  = undef;     # file to store Bio::Seq object
  $self->{'_results'}   = undef;     # file to store results of analysis
  $self->{'_prune'}     = 1;         # 
  $self->{'_coverage'}  = 10;

  # Now parse the input options and store them in the object

  my ($queries, $program, $database, $threshold, $threshold_type, 
      $filter,$coverage,$prune,$options) = 
            $self->_rearrange([qw(QUERIES 
                                  PROGRAM 
                                  DATABASE 
                                  THRESHOLD
                                  THRESHOLD_TYPE
                                  FILTER 
                                  COVERAGE
                                  PRUNE
                                  OPTIONS)], 
                              @args);
  if ($queries) {
    $self->queries($queries);
  } else {
    $self->throw("No query sequence array input.");
  }

  $self->program($self->find_executable($program));
  if ($database) {
    $self->database($database);
  } else {
    $self->throw("No database input");
  }
    
  if ($options) {
    $self->options($options);
  } else {
    $self->options(' ');  
  }
    
  if (defined($threshold)) {
    $self->threshold($threshold);
  }

  if (defined($threshold_type)) {
    $self->threshold_type($threshold_type);
  }

  if (defined($filter)) {
    $self->filter($filter);
  }
    
  if (defined($prune)) {
    $self->prune($prune);
  }
  if (defined($coverage)) {
    $self->coverage($coverage);
  }
  return $self; # success - we hope!
}


###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs hsperate and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
  my ($self, $dir) = @_;

  my $seqs = $self->queries || $self->throw("Query seqs required for Hsperate\n");

  $self->workdir('/tmp') unless ($self->workdir($dir));
  $self->checkdir();

  #write sequences to file
  $self->write_queries(); 
  $self->run_analysis();

  #parse output and create features
  $self->parse_results;
  $self->deletefiles();
}

sub write_queries {
  my ($self) = @_;
 
  my $queries = $self->queries;

  foreach my $query (@$queries) {
    $query->isa("Bio::PrimarySeqI")
      || $self->throw("Input isn't a Bio::PrimarySeqI");
  }
 
  my $filename = $self->filename;
  my $queryOutput =
    Bio::SeqIO->new(-file => ">$filename", '-format' => 'Fasta')
    or $self->throw("Can't create new Bio::SeqIO from $filename '$' : $!");
 
  foreach my $qseq (@$queries) {
    $self->throw("problem writing query seqeunce to $filename\n")
      unless $queryOutput->write_seq($qseq);
  }
  $queryOutput->close();
}     

sub databases {
  my ($self,@dbs) = @_;

  if (!defined($self->{_databases})) {
     $self->{_databases} = [];
  }
  if (defined(@dbs)) {
     push(@{$self->{_databases}},@dbs);
  }
  return @{$self->{_databases}};
}

=head2 run_analysis

    Title   :   run_analysis
    Usage   :   $obj->run_analysis
    Function:   Runs the hsperate query
    Returns :   nothing
    Args    :   none

=cut

sub run_analysis {
  my ($self) = @_;

  # This routine expands the database name into $db-1 etc for
  # split databases

  my @databases = $self->fetch_databases;
  my @hits;

  $self->databases(@databases);

  foreach my $database (@databases) {
    my $db = $database;
    $db =~ s/.*\///;
    my $command = $self->program ;
    $command .= ' '.$self->options;
    $command .= ' '.$self->filename;
    $command .= ' '.$database;
    $command .= ' | sort -k 1,1 -k 9,9 -k 10,10n '; 
    $command .= ' > '.$self->results . ".$db";

    print "Command " . $command . "\n";
    $self->throw("Failed during hsperate run $!\n") unless (system ($command) == 0) ;
  }
}


=head2 fetch_databases

    Title   :   fetch_databases
    Usage   :   $obj->fetch_databases
    Function:   For split databases this
                method checks whether the database in $self->database
                exists.  If it doesn\'t it checks whether the
                database has been split into $database-1, $database-2 etc.
    Returns :   Array of database names
    Args    :   none

=cut

sub fetch_databases {
  my ($self) = @_;
    
  my @databases;
    
  my $fulldbname;

  # If we have passed a full path name don't append the $BLASTDB
  # environment variable.

  if ($self->database =~ /\//) {
    $fulldbname = $self->database;
  } else {
    $fulldbname = $ENV{BLASTDB} . "/" .$self->database;
  }

  # If the expanded database name exists put this in
  # the database array.
  #
  # If it doesn't exist then see if $database-1,$database-2 exist
  # and put them in the database array

  if (-e $fulldbname) {
    push(@databases,$self->database);
  } else {
    my $count = 1;

    while (-e $fulldbname . "-$count") {
      push(@databases,$fulldbname . "-$count");
      $count++;
    }
  }

  if (scalar(@databases) == 0) {
    $self->throw("No databases exist for " . $self->database);
  }

  return @databases;
}



sub get_parsers {
  my ($self)  = @_;

  my @parsers;

  foreach my $db ($self->databases) {
    $db =~ s/.*\///;

    my $fh = new FileHandle;
    $fh->open("<" . $self->results . ".$db");
    
    my $parser = new Bio::EnsEMBL::Pipeline::Runnable::HsperateParser ('-fh' => $fh);
    
    push(@parsers,$parser);
  } 

  return @parsers;
}


=head2 parse_results

    Title   :   parse_results
    Usage   :   $obj->parse_results
    Function:   Parses the hsperate results and stores the output in
                an array of feature pairs.  The output will be filtered
                by the score.
                If we input a filehandle to this method the results 
                will be read from this file rather than the filename 
                stored in the object.  This is handy for processing
                hsperate runs that have been run previously.
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   Optional filehandle. 
                 

=cut

sub parse_results {
  my ($self,$fh) = @_;

  my @parsers;

  if (defined($fh)) {
    my $parser = new Bio::EnsEMBL::Pipeline::HsperateParser(-fh => $fh);
    push(@parsers,$parser);
  } else {
    @parsers = $self->get_parsers;
  }

  foreach my $parser (@parsers) {
    print STDERR "New parser\n";

    my $fh = $parser->fh;

    HSP: while (<$fh>) {

      my ($queryid,$type,$source,$qstart,$qend,$score,$qstrand,$qphase,$hitid,
          $hstart,$hend,$hstrand,$hphase,$length) = split /\t/;
 
# NOTE  percent not implemented yet
      my $percent = 0;
        
      my $name = $hitid;          

      # print STDERR "Name " . $name . "\n";

      my ($ug) = $name =~ m{/ug=(.*?)\ };
  
      if ($name =~ /\|UG\|(\S+)/) {
        # scp - unigene ID 'patch'
        # there must be a better way of doing this...
        if (length $ug > 0) { # just in case "/ug=" not in header
          $name = $ug;
        } else {
          $name = $1;
        }
      } elsif ($name =~ /\S+\|(\S+)\|\S+/) {
        $name = $1;
      } elsif ($name =~ /^(\S+) (\S+)/) {
        $name = $1 || $2;
      }

      if ($self->threshold_type eq "PID") {
        next HSP if ($length < $self->threshold);
      } elsif ($self->threshold_type eq "SCORE") {
        if ($score < $self->threshold) {
          # print "Below SCORE threshold \n";
          next HSP 
        }
      }
      # $self->_addHSPasFeaturePair($hsp,$name);

      my $analysis = new Bio::EnsEMBL::Analysis(-db              => $self->database,
                                            -db_version      => 1,
                                            -program         => $source,
                                            -program_version => 1,
                                            -gff_source      => $source,
                                            -gff_feature     => 'similarity',
                                            -logic_name      => 'hsperate');
      my $fp = $self->_makeFeaturePair($queryid, $qstart,$qend,$qstrand,
                                   $hstart,$hend,$hstrand,$score,
                                   $percent,0,$name,$analysis);
      $self->growfplist($fp);                                         
    }
  }

  my @allfeatures = $self->output;
  if ($self->threshold_type eq "PID") {
    @allfeatures = sort {$b->percent_id <=> $a->percent_id} @allfeatures;
  } else {
    @allfeatures = sort {$a->score <=> $b->score} @allfeatures;
  }
  if ($self->filter) {
    my $search = new Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter(
                                                -prune    => $self->prune,
                                                -coverage => $self->coverage
                                                                    );

    my @pruned = $search->run(@allfeatures);

    print STDERR "dbg", scalar(@allfeatures), " ", scalar(@pruned), "\n";
    $self->output(@pruned);
  }

  print "Output size = " . $self->output;
  return $self->output;
}

sub prune {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_prune} = $arg;
  }
  return $self->{_prune};
}

sub coverage {
  my($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_coverage} = $arg;
  }
  return $self->{_coverage};
}

    
=head2 _makeFeaturePair

    Title   :   _makeFeaturePair
    Usage   :   $obj->_makeFeaturePair
    Function:   Internal function that makes feature pairs
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   
                 

=cut

sub _makeFeaturePair {
  my ($self,$qid,$qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$score,
      $pid,$evalue,$name,$analysis)  = @_;

  my $source = $self->program;             
  $source =~ s/\/.*\/(.*)/$1/;

  my $feature1 = new Bio::EnsEMBL::SeqFeature(-seqname     => $qid,
                                              -start       => $qstart,
                                              -end         => $qend,
                                              -strand      => $qstrand * $hstrand,
                                              -source_tag  => $source,
                                              -primary_tag => 'similarity',
                                              -analysis    => $analysis,
                                              -score       => $score);
        
  $feature1->percent_id($pid);
  $feature1->p_value($evalue);

  my $feature2 = new Bio::EnsEMBL::SeqFeature(-seqname     => $name,
                                              -start       => $hstart,
                                              -end         => $hend,
                                              -strand      => $hstrand * $qstrand,
                                              -source_tag  => $source,
                                              -primary_tag => 'similarity',
                                              -analysis    => $analysis,
                                              -score       => $score);
    
    
  my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
                                         -feature2 => $feature2);
   
  $feature2->percent_id($pid);
  $feature2->p_value($evalue);
 
  return $fp;
}

##############
# input/output methods
#############

=head2 output

    Title   :   output
    Usage   :   obj->output()
    Function:   Returns an array of feature pairs
    Returns :   Returns an array of feature pairs
    Args    :   none

=cut

sub output {
  my ($self, @arg) = @_;

  if (@arg) {
    @{$self->{'_fplist'}} = @arg;
  }

  if (!defined($self->{'_fplist'})) {
    $self->{'_fplist'} = [];
  }
  return @{$self->{'_fplist'}};
}


#################
# get/set methods 
#################

sub queries {
  my ($self, $value) = @_;
 
  if ($value) {
    $self->throw ("queries should be a ref to an ARRAY") unless (ref($value) eq 'ARRAY');
    $self->{'_queries'} = $value;

    $self->filename("hsperate_".$$. ".seq");
    $self->results($self->filename.".hsperate.out"); 
  }

 
  #NB ref to an array of Bio::Seq
  return $self->{'_queries'};
}     


=head2 program

    Title   :   program
    Usage   :   $obj->program('/usr/local/pubseq/bin/hsperate');
    Function:   Get/set method for the location of the hsperate flavour
    Returns :   string
    Args    :   string

=cut

sub program {
  my ($self, $location) = @_;
  
  if ($location) {
    $self->throw("executable not found at $location: $!\n")  
      unless (-e $location && -x $location);
    $self->{'_program'} = $location ;
  }
  return $self->{'_program'};
}

=head2 database

    Title   :   database
    Usage   :   $obj->database('dbEST');
    Function:   Get/set method for the location of database
    Args    :   none

=cut

sub database {
  my ($self, $db) = @_;

  if (defined($db)) {
    $self->{'_database'} = $db ;
  }
  return $self->{'_database'};
}

=head2 options

    Title   :   options
    Usage   :   $obj->options(' -I ');
    Function:   Get/set method for hsperate arguments
    Args    :   File path (optional)

=cut

sub options {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{'_options'} = $args ;
  }
  return $self->{'_options'};
}

sub filter {
  my ($self,$args) = @_;

  if (defined($args)) {
    if ($args != 0 && $args != 1) {
      $self->throw("Filter option must be 0 or 1");
    }
    $self->{'_filter'} = $args;
  }
  return $self->{'_filter'};
}

sub get_threshold_types {
  my ($self) = @_;

  return ("PID","SCORE");
}

sub threshold_type {
  my ($self,$type) = @_;

  my @types = $self->get_threshold_types;
  
  if (defined($type)) {
    my $found = 0;
    foreach my $allowed_type ($self->get_threshold_types) {
      if ($type eq $allowed_type) {
        $found = 1;
      }
    }
    if ($found == 0) {

      $self->throw("Type [$type] is not an allowed type.  Allowed types are [@types]");
    } else {
      $self->{_threshold_type} = $type;
    }
  }
  return $self->{_threshold_type} || $types[0];
}


sub get_pars {
  my ($self) = @_;

  if (!defined($self->{_hits})) {
     $self->{_hits} = [];
  }
  return @{$self->{_hits}};
}

1;

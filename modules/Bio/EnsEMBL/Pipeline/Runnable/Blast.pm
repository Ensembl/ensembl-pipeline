#!/usr/local/bin/perl -w

#
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Blast

=head1 SYNOPSIS

  # To run a blast job from scratch do the following.

  my $query = new Bio::Seq(-file   => 'somefile.fa',
			   -format => 'fasta');

  my $blast =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new (-query    => $query,
							     -program  => 'blastp',
							     -database => 'swir',
							     -threshold => 1e-6,
							     -options   => 'V=1000000');

  $blast->run();

  @featurepairs = $blast->output();

  foreach my $fp (@featurepairs) {
      print $fp->gffstring . "\n";
  }

  # Additionally if you have blast runs lying around that need parsing
  # you can do

  open(BLAST,"<blast.out");

  my @featurepairs = Bio::EnsEMBL::Pipeline::Runnable::Blast->parse_results(\*BLAST);

  close(BLAST);

=head1 DESCRIPTION

Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs blast with, the
output is parsed by BPLite and stored as Bio::EnsEMBL::FeaturePairs. 

The output features can be filtered by probability using the -threshold option.
Other options can be passed to the blast program using the -options method

=head1 CONTACT

Describe contact details here

=head1 APPENDIX


=cut

package Bio::EnsEMBL::Pipeline::Runnable::Blast;


use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis;
use Bio::PrimarySeq; 
use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::Object;
use BPlite;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::Object);

=head2 new

    Title   :   new
    Usage   :   my obj =  Bio::EnsEMBL::Pipeline::Runnable::Blast->new (-query    => $seq,
									-program  => 'blastp',
								        -database => 'swir',
								        -threshold => 1e-6,
								        -options   => 'V=1000000');
    Function:   Initialises Blast object
    Returns :   a Blast Object
    Args    :   A Bio::Seq object
                The blast executable (-BLAST) and database (-DB).

=cut

sub _initialize {
    my ($self,@args) = @_;

    my $make = $self->SUPER::_initialize(@_);    

    $self->{_query}     = undef;     # location of Bio::Seq object
    $self->{_program}   = undef;     # location of Blast
    $self->{_database}  = undef;     # name of database
    $self->{_threshold} = undef;     # Threshold for hit filterting
    $self->{_options}   = undef;     # arguments for blast

    $self->{_fplist}    = [];        # an array of feature pairs (the output)

    $self->{_workdir}   = undef;     # location of temp directory
    $self->{_filename}  = undef;     # file to store Bio::Seq object
    $self->{_results}   = undef;     # file to store results of analysis

      
    # Now parse the input options and store them in the object

    my( $query, $program, $database, $threshold, $options) = 
	$self->_rearrange(['QUERY', 'PROGRAM', 'DATABASE','THRESHOLD','OPTIONS'], @args);
   
    if ($query) {
      $self->clone($query);
    } else {
      $self->throw("No query sequence input.");
    }

    if ($program =~ m!/!) { #path to blast is provided 
      $self->program($program);
    } elsif ($program =~ /blastn|blastx|blastp|tblastn|tblastx/) {
      $self->program($self->locate_executable($program));  
    } else {
      $self->throw("Path to blast executable required [$program]\n");
    }
    
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

    return $self; # success - we hope!
  }


###########
# Analysis methods
##########

=head2 run

    Title   :  run
    Usage   :   $obj->run()
    Function:   Runs blast and BPLite and creates array of feature pairs
    Returns :   none
    Args    :   none

=cut

sub run {
    my ($self, $dir) = @_;

    my $seq = $self->clone || $self->throw("Query seq required for Blast\n");

    $self->workdir('/tmp') unless ($self->workdir($dir));
    $self->checkdir();

    #write sequence to file
    $self->writefile(); 
    $self->run_analysis();

    #parse output and create features
    $self->parse_results();
    $self->deletefiles();
  }

=head2 run_analysis

    Title   :   run_analysis
    Usage   :   $obj->run_analysis
    Function:   Runs the blast query
    Returns :   nothing
    Args    :   none

=cut

sub run_analysis {
    my ($self) = @_;
    
    print STDERR ("Running blast and BPlite:\n " . $self->program  . ' ' .
                  		                   $self->database . ' ' . 
		                                   $self->filename . ' ' .
		                                   $self->options  . ' > ' .
		                                   $self->results."\n");

    # This routine expands the database name into $db-1 etc for
    # split databases

    my @databases = $self->fetch_databases;

    foreach my $database (@databases) {
	$self->throw("Failed during blast run $!\n")
	    
	    unless (system ($self->program . ' ' . 
			    $database      . ' ' .
			    $self->filename. ' ' .
			    $self->options . ' > ' . 
			    $self->results) == 0) ;
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

=head2 parse_results

    Title   :   parse_results
    Usage   :   $obj->parse_results
    Function:   Parses the blast results and stores the output in
                an array of feature pairs.  The output will be filtered
                by the probability threshold.
                If we input a filehandle to this method the results 
                will be read from this file rather than the filename 
                stored in the object.  This is handy for processing
                blast runs that have been run previously.
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   Optional filehandle. 
                 

=cut


sub parse_results {
  my ($self,$fh) = @_;

  # If we have input a filehandle use that. Otherwise use
  # the results file stored in the object
  my $filehandle;

  if (defined($fh)) {
      $filehandle = $fh;
  } elsif (ref ($self->results) !~ /GLOB/) {
      open (BLAST, "<".$self->results)
	  or $self->throw ("Couldn't open file ".$self->results.": $!\n");
      $filehandle = \*BLAST;
  } else {
      $filehandle = $self->results;
  }    

  unless (<$filehandle>) {
    print "No hit found with blast \n";
    return;
  }

  # BPlite does most of the work

  my $parser = new BPlite ($filehandle);

  # Loop over each blast hit
  while(my $sbjct = $parser->nextSbjct)  {

      my $name = $sbjct->name ;

      $name =~ s/^>(\S+).*/$1/;

      if ($name =~ /\|UG\|(\S+)/) {
         $name = $1;
      } elsif ($name =~ /\S+\|\S+\|(\S+)/) {
         $name = $1;
      }
      
    HSP: while (my $hsp = $sbjct->nextHSP) {
	next HSP if ($hsp->P > $self->threshold);
	
	# Each HSP is a gapped alignment.
	# This method splits the gapped alignment into
	# ungapped pieces

	$self->split_HSP($hsp,$name);

    }
  } 
  return $self->output;
}


=head2 split_HSP

    Title   :   split_HSP
    Usage   :   $obj->split_HSP
    Function:   Takes a gapped blast HSP 
                and turns it into an array of ungapped feature pairs.
    Returns :   Nothing
    Args    :   BPlite::HSP,name string
                 

=cut

sub split_HSP {
    my ($self,$hsp,$name) = @_;

    # First of all some jiggery pokery to find out what sort of alignment
    # we have - (dna-dna, dna-pep, pep-dna etc).
    # We also work out which strand each sequence is on.
    # 
    # Both these variables are needed to work out the increment
    # values for the query and the hit sequences.  As we split
    # up the gapped alignment we need to increment the coordinates.
    #
    # We have to take care with dna-pep alignments to increment the dna
    # sequences by 3 as there are 3 bases in a codon.
    #
    # For dna-dna alignments this will be 1 or -1 depending on
    # the strand.  
    #
    # For pep-dna alignments and vice-versa the increments will be +-3 
    # for the dna sequence and +- 1 for the peptide sequence. 

    my ($qtype,  $htype)   = $self->_findTypes     ($hsp);
    my ($qstrand,$hstrand) = $self->_findStrands   ($hsp);
    my ($qinc,   $hinc)    = $self->_findIncrements($hsp,$qstrand,$hstrand,$qtype,$htype);

    print STDERR "Alignment q : " . $hsp->queryBegin . "\t" . $hsp->queryEnd . "\t" . $hsp->queryAlignment . "\n";
    print STDERR "Alignment s : " . $hsp->sbjctBegin . "\t" . $hsp->sbjctEnd . "\t" . $hsp->sbjctAlignment . "\n";

    print STDERR "types (increments) $qtype ($qinc) : $htype ($hinc)\n";

    # We split the alignment strings into arrays of one char each.  
    # We then loop over this array and when we come to a gap
    # in either the query sequence or the hit sequence we make 
    # a new feature pair. 
    #
    # Before making a new feature pair we check whether we have something to make
    # a feature pair out of - i.e. we aren't just on the 2nd or greater base in a gapped
    # region.
    #
    # As we loop over the array we need to keep track of both the query coordinates
    # and the hit coordinates.  We track the last start of a feature pair and
    # also the current end of a feature pair.  The feature pair start is reset every time 
    # we hit a gapped position and the feature pair end is reset every time we hit 
    # an aligned region.

    my @gap;

    my @qchars = split(//,$hsp->queryAlignment);  # split alignment into array of char
    my @hchars = split(//,$hsp->sbjctAlignment);  # ditto for hit sequence
    
    my $qstart = $hsp->queryBegin;                # Start off the feature pair start
    my $hstart = $hsp->sbjctBegin;                # ditto
    
    my $qend   = $hsp->queryBegin;                # Set the feature pair end also
    my $hend   = $hsp->sbjctBegin;                # ditto
    
    my $count = 0;                                # counter for the bases in the alignment
    my $found = 0;                                # flag saying whether we have a feature pair
    
    my $source = $self->program;             
    $source =~ s/\/.*\/(.*)/$1/;

    my $analysis = new Bio::EnsEMBL::Analysis(-db              => $self->database,
					      -db_version      => 1,
					      -program         => $source,
					      -program_version => 1,
					      -gff_source      => $source,
					      -gff_feature     => 'similarity');
    
    # Here goes...

    while ($count <= $#qchars) {
	# We have hit an ungapped region.  Increase the query and hit counters.
	# and flag that we have a feature pair.

	if ($qchars[$count] ne '-' &&
	    $hchars[$count] ne '-') {

	    $qend += $qinc;
	    $hend += $hinc;
	    
	    $found = 1;
	} else {

	    # We have hit a gapped region.  If the feature pair flag is set ($found)
	    # then make a feature pair, store it and reset the start and end variables.

	    if ($found == 1) {
		my $fp = $self->_convert2FeaturePair($qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$qinc,$hinc,$hsp,$name,$analysis);
		
		$self->growfplist($fp);                             	    
	    }
	
	    # We're in a gapped region.  We need to increment the sequence that
	    # doesn't have the gap in it to keep the coordinates correct.
	    # We also need to reset the current end coordinates.

	    if ($qchars[$count] ne '-') {
		$qstart = $qend   + $qinc;
	    } else {
		$qstart = $qend;
	    }
	    if ($hchars[$count] ne '-') {
		$hstart = $hend   + $hinc;
	    } else {
		$hstart = $hend;
	    }
	    
	    $qend = $qstart;
	    $hend = $hstart;

	    $found = 0;
	}
	$count++;
    }			     

    # Remember the last feature
    if ($found == 1) {
	my $fp = $self->_convert2FeaturePair($qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$qinc,$hinc,$hsp,$name,$analysis);
	$self->growfplist($fp);                             	    
    }

}

=head2 _convert2FeaturePair

    Title   :   convert2FeaturePair
    Usage   :   $obj->convert2FeaturePair
    Function:   Internal function taking a set of coords and
                converts them into a feature pair
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   int,int,int,int,int,int,int,int,BPlite::HSP
                 

=cut

sub _convert2FeaturePair {
    my ($self,$qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$qinc,$hinc,$hsp,$name,$analysis) = @_;

    # The actual end of the alignment is the previous character.
    
    my $tmpqend = $qend; $tmpqend -= $qinc;
    my $tmphend = $hend; $tmphend -= $hinc;
    
    my $tmpqstart = $qstart;
    my $tmphstart = $hstart;
    
    # This is for dna-pep alignments.  The actual end base
    # will be +- 2 bases further on.
    if (abs($qinc) > 1) {
	$tmpqend += $qstrand * 2;
    }
    if (abs($hinc) > 1) {
	$tmphend += $hstrand * 2;
    }
    
    # Make sure start is always < end
    if ($tmpqstart > $tmpqend) {
	my $tmp    = $tmpqstart;
	$tmpqstart = $tmpqend;
	$tmpqend   = $tmp;
    }
    if ($tmphstart > $tmphend) {
	my $tmp    = $tmphstart;
	$tmphstart = $tmphend;
	$tmphend   = $tmp;
    }
    
    print "Creating feature pair " . $tmpqstart . "\t" . $tmpqend . "\t" . $qstrand . "\t" . $tmphstart . "\t" . $tmphend . "\t" . $hstrand . "\n";

    my $fp = $self->_makeFeaturePair($tmpqstart,$tmpqend,$qstrand,$tmphstart,$tmphend,$hstrand,$hsp->score,
				     $hsp->percent,$hsp->P,$name,$analysis);

    return $fp;
}

=head2 _makeFeaturePair

    Title   :   _makeFeaturePair
    Usage   :   $obj->_makeFeaturePair
    Function:   Internal function that makes feature pairs
    Returns :   Bio::EnsEMBL::FeaturePair
    Args    :   int,int,int,int,int,int,int,int,BPlite::HSP
                 

=cut

sub _makeFeaturePair {
    my ($self,$qstart,$qend,$qstrand,$hstart,$hend,$hstrand,$score,$pid,$evalue,$name,$analysis)  = @_;

    my $source = $self->program;             
    $source =~ s/\/.*\/(.*)/$1/;

    my $feature1 = new Bio::EnsEMBL::SeqFeature(-seqname     => $self->clone->id,
						-start       => $qstart,
						-end         => $qend,
						-strand      => $qstrand * $hstrand,
						-source_tag  => $source,
						-primary_tag => 'similarity',
						-analysis    => $analysis,
						-score       => $score);
	
    $feature1->percent_id($pid);
    $feature1->p_value($evalue);

    my $feature2 = new Bio::EnsEMBL::SeqFeature(-seqname => $name,
						-start   => $hstart,
						-end     => $hend,
						-strand  => $hstrand * $qstrand,
						-source_tag  => $source,
						-primary_tag => 'similarity',
						-analysis => $analysis,
						-score    => $score);
    
    
    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
					   -feature2 => $feature2);
   
    $feature2->percent_id($pid);
    $feature2->p_value($evalue);
 
    return $fp;
}

sub _findIncrements {
    my ($self,$hsp,$qstrand,$hstrand,$qtype,$htype) = @_;

    my $qinc   = 1 * $qstrand;
    my $hinc   = 1 * $hstrand;
    
    if ($qtype eq 'dna' && $htype eq 'pep') {
	$qinc = 3 * $qinc;
    } 
    if ($qtype eq 'pep' && $htype eq 'dna') {
	$hinc = 3 * $hinc;
    }

    return ($qinc,$hinc);
}

sub _findStrands {
    my ($self,$hsp) = @_;

    my $qstrand;
    my $hstrand;

    if ($hsp->queryBegin < $hsp->queryEnd) {
	$qstrand = 1;
    } else {
	$qstrand = -1;
    }
    if ($hsp->sbjctBegin < $hsp->sbjctEnd) {
	$hstrand = 1;
    } else {
	$hstrand = -1;
    }

    return ($qstrand,$hstrand);
}

sub _findTypes {
    my ($self,$hsp) = @_;

    my $type1;
    my $type2;

    my $len1 = abs($hsp->queryEnd - $hsp->queryBegin) + 1;
    my $len2 = abs($hsp->sbjctEnd - $hsp->sbjctBegin) + 1;


    if ($len1/$len2 > 2) {
	$type1 = 'dna';
	$type2 = 'pep';
    } elsif ($len2/$len1 > 2) {
	$type1 = 'pep';
	$type2 = 'dna';
    } else {
	$type1 = 'dna';
	$type2 = 'dna';
    }

    return ($type1,$type2);
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
    my ($self) = @_;

    if (!defined($self->{'_fplist'})) {
       $self->{'_fplist'} = [];
    }
    return @{$self->{'_fplist'}};
  }


#################
# get/set methods 
#################

sub clone {
    my ($self, $seq) = @_;
    if ($seq) {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
	$self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
      }

      $self->{_query} = $seq ;

      $self->filename($self->clone->id.".$$.seq");
      $self->results($self->filename.".blast.out");
      
    }
    return $self->{_query};
}

=head2 program

    Title   :   program
    Usage   :   $obj->program('/usr/local/pubseq/bin/blastn');
    Function:   Get/set method for the location of the blast flavour
    Returns :   string
    Args    :   string

=cut

sub program {
  my ($self, $location) = @_;
  
  if ($location) {
    $self->throw("executable not found at $location: $!\n") 	unless (-e $location && -x $location);
    $self->{_program} = $location ;
  }
  return $self->{_program};
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
      $self->{_database} = $db ;
    }
    return $self->{_database};
  }

=head2 options

    Title   :   options
    Usage   :   $obj->options(' -I ');
    Function:   Get/set method for blast arguments
    Args    :   File path (optional)

=cut

sub options {
  my ($self, $args) = @_;
  
  if (defined($args)) {
    $self->{_options} = $args ;
  }
  return $self->{_options};
}

1;


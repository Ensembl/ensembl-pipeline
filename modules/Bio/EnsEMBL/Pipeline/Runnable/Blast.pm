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

  #create and fill Bio::Seq object
  my $dnafile   = '/nfs/disk65/mq2/temp/bA151E14.seq';
  my $seqstream = Bio::SeqIO->new(-file => $clonefile, -fmt => 'Fasta');
  my $seq       = $seqstream->next_seq();

  #create Bio::EnsEMBL::Pipeline::Runnable::Blast object
  my $blast = Bio::EnsEMBL::Pipeline::Runnable::Blast->new (-query     => $seq,
							    -program   => 'blastn',
							    -database  => 'dbest',
							    -threshold => 100);
  $blast->run();
  @featurepairs = $blast->output();

=head1 DESCRIPTION

Blast takes a Bio::Seq (or Bio::PrimarySeq) object and runs blast with, the
output is parsed by BPLite and stored as Bio::EnsEMBL::FeaturePairs. 


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

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

use Data::Dumper;

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
    $self->{_database}  = undef;     # name and location of database
    $self->{_threshold} = undef;     # Threshold for hit filterting
    $self->{_options}   = undef;     # arguments for blast

    $self->{_fplist}    = [];        # an array of feature pairs (the output)

    $self->{_workdir}   = undef;     # location of temp directory
    $self->{_filename}  = undef;     # file to store Bio::Seq object
    $self->{_results}   = undef;     # file to store results of analysis

      
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

    #set directory if provided (MC this should be in a config file somewhere)
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

sub fetch_databases {
    my ($self) = @_;
    
    my @databases;
    
    my $fulldbname;

    if ($self->database =~ /\//) {
	$fulldbname = $self->database;
    } else {
	$fulldbname = $ENV{BLASTDB} . "/" .$self->database;
    }

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
    Function:   Runs the blast query
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   optional filename

=cut


sub parse_results {
  my ($self,$fh) = @_;

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

  my $parser = new BPlite ($filehandle);

  $parser->query;
  $parser->database;

  my $threshold = $self->threshold;

  while(my $sbjct = $parser->nextSbjct)  {
    HSP: while (my $hsp = $sbjct->nextHSP) {
	next HSP if ($hsp->P > $self->threshold);
	
	$self->split_HSP($sbjct,$hsp);

    }
  } 
  return $self->output;
}


sub split_HSP {
    my ($self,$sbjct,$hsp) = @_;

    my $name = $sbjct->name ;
    $name =~ s/^>(\S+).*/$1/;

    my $type1;
    my $type2;

    my $len1 = abs($hsp->queryEnd - $hsp->queryBegin) + 1;
    my $len2 = abs($hsp->sbjctEnd - $hsp->sbjctBegin) + 1;


    # Find the strands first

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

    print STDERR "Alignment q : " . $hsp->queryBegin . "\t" . $hsp->queryEnd . "\t" . $hsp->queryAlignment . "\n";
    print STDERR "Alignment s : " . $hsp->sbjctBegin . "\t" . $hsp->sbjctEnd . "\t" . $hsp->sbjctAlignment . "\n";

    my @features;
    
    my $qinc   = 1 * $qstrand;
    my $sinc   = 1 * $hstrand;
    
    if ($type1 eq 'dna' && $type2 eq 'pep') {
	$qinc = 3 * $qinc;
    } 
    if ($type1 eq 'pep' && $type2 eq 'dna') {
	$sinc = 3 * $sinc;
    }

    print STDERR "types $type1 ($qinc) : $type2 ($sinc)\n";
    my @gap;
    
    my @qchars = split(//,$hsp->queryAlignment);
    my @schars = split(//,$hsp->sbjctAlignment);
    
    my $qstart = $hsp->queryBegin;
    my $sstart = $hsp->sbjctBegin;
    
    my $qend   = $hsp->queryBegin;
    my $send   = $hsp->sbjctBegin;
    
    my $count = 0;
    my $found = 0;
    
    my $source = $self->program;
    $source =~ s/\/.*\/(.*)/$1/;

    my $analysis = new Bio::EnsEMBL::Analysis(-db              => $self->database,
					      -db_version      => 1,
					      -program         => $source,
					      -program_version => $1,
					      -gff_source      => $source,
					      -gff_feature     => 'similarity');
    
    while ($count <= $#qchars) {
	if ($qchars[$count] ne '-' &&
	    $schars[$count] ne '-') {

	    $qend += $qinc;
	    $send += $sinc;
	    
	    $found = 1;
	} else {
	    if ($found == 1) {

		my $tmpqend = $qend; $tmpqend -= $qinc;
		my $tmpsend = $send; $tmpsend -= $sinc;

		my $tmpqstart = $qstart;
		my $tmpsstart = $sstart;

		if (abs($qinc) > 1) {
		    $tmpqend += $qstrand * 2;
		}
		if (abs($sinc) > 1) {
		    $tmpsend += $hstrand * 2;
		}

		if ($tmpqstart > $tmpqend) {
		    my $tmp    = $tmpqstart;
		    $tmpqstart = $tmpqend;
		    $tmpqend   = $tmp;
		}
		if ($tmpsstart > $tmpsend) {
		    my $tmp    = $tmpsstart;
		    $tmpsstart = $tmpsend;
		    $tmpsend   = $tmp;
		}

		print "Creating feature pair " . $tmpqstart . "\t" . $tmpqend . "\t" . $qstrand . "\t" . $tmpsstart . "\t" . $tmpsend . "\t" . $hstrand . "\n";

		my $feature1 = new Bio::EnsEMBL::SeqFeature(-seqname     => $self->clone->id,
							    -start       => $tmpqstart,
							    -end         => $tmpqend,
							    -strand      => $qstrand * $hstrand,
							    -source_tag  => $source,
							    -primary_tag => 'similarity',
							    -analysis    => $analysis,
							    -score       => $hsp->score);
		
		my $feature2 = new Bio::EnsEMBL::SeqFeature(-seqname => $name,
							    -start   => $tmpsstart,
							    -end     => $tmpsend,
							    -strand  => $hstrand * $qstrand,
							    -source_tag  => $source,
							    -primary_tag => 'similarity',
							    -analysis => $analysis,
							    -score    => $hsp->score);
		
		
		my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
						       -feature2 => $feature2);
		
		print $fp->source_tag . "\t" . $fp->primary_tag . "\t" . $fp->score . "\n";
		push(@features,$fp);
		$self->growfplist($fp);                             
		
	    }
	    if ($qchars[$count] ne '-') {
		$qstart = $qend   + $qinc;
	    } else {
		$qstart = $qend;
	    }
	    if ($schars[$count] ne '-') {
		$sstart = $send   + $sinc;
	    } else {
		$sstart = $send;
	    }
	    
	    $qend = $qstart;
	    $send = $sstart;

	    $found = 0;
	}
	$count++;
    }			     

    # Remember the last feature
    if ($found == 1) {
	my $tmpqend = $qend; $tmpqend -= $qinc;
	my $tmpsend = $send; $tmpsend -= $sinc;
	
	my $tmpqstart = $qstart;
	my $tmpsstart = $sstart;
	
	if (abs($qinc) > 1) {
	    $tmpqend += $qstrand * 2;
	}
	if (abs($sinc) > 1) {
	    $tmpsend += $hstrand * 2;
	}
	
	if ($tmpqstart > $tmpqend) {
	    my $tmp = $tmpqstart;
	    $tmpqstart = $tmpqend;
	    $tmpqend   = $tmp;
	}
	if ($tmpsstart > $tmpsend) {
	    my $tmp = $tmpsstart;
	    $tmpsstart = $tmpsend;
	    $tmpsend   = $tmp;
	}
	
	print "Creating feature pair " . $tmpqstart . "\t" . $tmpqend . "\t" . $qstrand . "\t" . $tmpsstart . "\t" . $tmpsend . "\t" . $hstrand . "\n";
	
	my $feature1 = new Bio::EnsEMBL::SeqFeature(-seqname     => $self->clone->id,
						    -start       => $tmpqstart,
						    -end         => $tmpqend,
						    -strand      => $qstrand,
						    -source_tag  => $source,
						    -primary_tag => 'similarity',
						    -analysis    => $analysis,
						    -score       => $hsp->score,
						    );
	
	my $feature2 = new Bio::EnsEMBL::SeqFeature(-seqname => $name,
						    -start   => $tmpsstart,
						    -end     => $tmpsend,
						    -strand  => $hstrand,
						    -source_tag  => $source,
						    -primary_tag => 'similarity',
						    -analysis => $analysis,
						    -score    => $hsp->score);
	
	my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $feature1,
					       -feature2 => $feature2);
	
	push(@features,$fp);
	$self->growfplist($fp);                             
	
    }

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


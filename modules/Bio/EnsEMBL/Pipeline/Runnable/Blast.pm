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
      $self->query($query);
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

    my $seq = $self->query || $self->throw("Query seq required for Blast\n");

    #set directory if provided (MC this should be in a config file somewhere)
    $self->workdir('/work2/michele') unless ($self->workdir($dir));
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
    
    $self->throw("Failed during blast run $!\n")

      unless (system ($self->program.' '.$self->database.' '.$self->filename
		      .' '.$self->options .' > '.$self->results) == 0) ;
  }

=head2 parse_results

    Title   :   parse_results
    Usage   :   $obj->parse_results
    Function:   Runs the blast query
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   optional filename

=cut


sub parse_results {
  my ($self) = @_;

  my $filehandle;

  if (ref ($self->results) !~ /GLOB/) {
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
      
      my (%feat1, %feat2);

      if ($self->query->id) {
	$feat1{name}     = $self->query->id;
      } else {
	$self->results =~ m!/.+/(.+)|(.+)!; #extract filename
	$feat1{name} = ($1) ?  $1 :  $2;
      }
            
      $feat1 {score}       = $hsp->score;
      $feat2 {score}       = $hsp->score;
      $feat1 {percent}     = $hsp->percent;
      $feat2 {percent}     = $hsp->percent;
      $feat1 {p}           = $hsp->P;
      $feat2 {p}           = $hsp->P;
            
      if ($hsp->queryBegin < $hsp->queryEnd) {
	$feat1 {start}   = $hsp->queryBegin;
	$feat1 {end}     = $hsp->queryEnd;
	$feat1 {strand}  = 1;
      } else {
	$feat1 {start}   = $hsp->queryEnd;
	$feat1 {end}     = $hsp->queryBegin;
	$feat1 {strand}  = -1;
      }
      
      $sbjct->name =~ /[\||\s|:](\w+)[\||\s|:]/; #extract subjectname
      $feat2 {name}    = $1;
      
      #MC 29/11/00 The new viersion of BPlite deals with strand 
      # so we won't have to do this.
      if ($hsp->sbjctBegin < $hsp->sbjctEnd) {
	$feat2 {start}   = $hsp->sbjctBegin;
	$feat2 {end}     = $hsp->sbjctEnd;
	$feat2 {strand}  = 1;
      } else {
	$feat2 {start}   = $hsp->sbjctEnd;
	$feat2 {end}     = $hsp->sbjctBegin;
	$feat2 {strand}  = -1;
      }

      # We force the feat1 strand always to be 1.  This is not necessarily correct
      my $tmpstrand = $feat1{strand};
      
      $feat1{strand} = 1;
      $feat2{strand} = $tmpstrand * $feat2{strand};
      
      if ($self->database) {
	$feat2 {db} = $self->database;
      } else {
	$feat2 {db} = 'unknown';
      }
                        
      if ($self->program) {
	$self->program =~ m!/.+/(.+)|(.+)!; #extract executable name
	if ($1)  {
	  $feat2 {program} = $1; 
	} elsif ($2)  {
	  $feat2 {program} = $2; 
	}
      } else {
	$feat2 {program} = 'unknown';
      }

      $feat2 {p_version}   = '1';
      $feat2 {db_version}  = '1';
      $feat1 {primary}     = 'similarity';
      $feat1 {source}      = $feat2{program};
      $feat2 {primary}     = 'similarity';
      $feat2 {source}      = $feat2{program};
      
      #if alignments contain gaps, the feature needs to be split
      $feat1 {alignment} = $hsp->queryAlignment;
      $feat2 {alignment} = $hsp->sbjctAlignment;
            
      if ($feat1 {alignment} =~ /-/ or $feat2 {alignment} =~ /-/) {
	$self->split_gapped_feature(\%feat1, \%feat2); 
      } else {                    
	$self->createfeaturepair(\%feat1, \%feat2); 
      }
    }
  } 
}


#This function creates mini-features from a gapped feature. 
#The gaps are discarded and the mini alignments have the attributes of the parent feature. 
sub split_gapped_feature {
  my ($self, $feat1, $feat2) = @_;
    
  my $type1;
  my $type2;

  my $len1 = $feat1->{end} - $feat1->{start} + 1;
  my $len2 = $feat2->{end} - $feat2->{start} + 1;

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
  
  #There's a strange bug in wublastn where a gap is inserted at the end of an alignment
  #The alignment is trimmed of any terminal gaps so that the alignment ends on a valid feature
  while ($feat1->{'alignment'} =~ /-$/ || $feat2->{'alignment'} =~ /-$/) {
    $feat1->{'alignment'} =~ s/.$//;
    $feat2->{'alignment'} =~ s/.$//;
  }

  #The terminal gap bug is so odd, it's probably worth checking for an initial gap
  $self->throw("Alignment starts with a gap! F1\n"
	       .$feat1->{'alignment'}."\nF2\n".$feat2->{'alignment'}."\n")
    if ($feat1->{'alignment'} =~ /^-/ || $feat2->{'alignment'} =~ /^-/);
    
  my (@masked_f1, @masked_f2);

  #replace bases and gaps with positions and mask number
  if ($type1 = 'pep' && $type2 == 'dna') {
    @masked_f1 = $self->mask_alignment($feat1->{'start'}, $feat1->{'strand'}, $feat1->{'alignment'},1);
    @masked_f2 = $self->mask_alignment($feat2->{'start'}, $feat2->{'strand'}, $feat2->{'alignment'},3);
  } elsif ($type1 = 'dna' && $type2 == 'pep') {
    @masked_f1 = $self->mask_alignment($feat1->{'start'}, $feat1->{'strand'}, $feat1->{'alignment'},3);
    @masked_f2 = $self->mask_alignment($feat2->{'start'}, $feat2->{'strand'}, $feat2->{'alignment'},1);
  }
  $self->throw("Can't split feature where alignment lengths don't match: F1 ("
	       .scalar(@masked_f1).") F2 (".scalar(@masked_f2).")\n")
    if (scalar(@masked_f1) != scalar(@masked_f2)); 
    
  my $building_feature;
  my $mask_len = scalar(@masked_f1);
  my ($f1_start, $f2_start);
    for (my $index =0; $index < $mask_len; $index++)
    {
        
        if ($masked_f1[$index] == -1 || $masked_f2[$index] == -1 || $index == $mask_len -1)
        {
            #One of the alignments contains an insertion.
            if ($building_feature)
            {            
                #feature ended at previous position unless alignment end
                my $f1_end     = ($index == $mask_len -1) 
                                        ? $masked_f1[$index] : $masked_f1[$index-1];
                my $f2_end     = ($index == $mask_len -1) 
                                        ? $masked_f2[$index] :$masked_f2[$index-1];
                
                if ($feat1->{'strand'} == 1)
                {
                    $feat1->{'start'}   = $f1_start;
                    $feat1->{'end'}     = $f1_end;
                }
                else
                {
                    $feat1->{'start'}   = $f1_end;
                    $feat1->{'end'}     = $f1_start;
                }
                if ($feat2->{'strand'} == 1)
                {
                    $feat2->{'start'}   = $f2_start;
                    $feat2->{'end'}     = $f2_end;
                }
                else
                {
                    $feat2->{'start'}   = $f2_end;
                    $feat2->{'end'}     = $f2_start;
                }
                
                #print STDERR "Subfeat1: ".$feat1->{'start'}." - ".$feat1->{'end'}
                #             ."\tSubfeat2: ".$feat2->{'start'}." - ".$feat2->{'end'}."\n";
                
                my $f1_len = $feat1->{'end'} - $feat1->{'start'} +1;
                my $f2_len = $feat2->{'end'} - $feat2->{'start'} +1; 
                
#                $self->throw("FeaturePair lengths don't match! ".
#                        "F1 $f1_start - $f1_end ($f1_len) F2 $f2_start - $f2_end ($f2_len)\n") 
#                        if ( $f1_len != $f2_len );
                 
                        
                $self->createfeaturepair($feat1, $feat2);
                $building_feature = 0;
            }
            
        }
        else
        {
            #Alignment of two bases found
            if (!$building_feature)
            {
                $f1_start = $masked_f1[$index];
                $f2_start = $masked_f2[$index];
                $building_feature =1;
            }
        }
    }
  }

#Fills gapped alignment with base position number or -1 for insertions. 
sub mask_alignment {
    my ($self, $start, $strand, $alignment,$inc) =@_;
    my @masked_array;
    my @array = split (//,$alignment);
    $_ = $alignment;
    my $valid_bases = tr/A-Za-z//;
    
    #print STDERR "Start: $start Strand: $strand Len ".scalar(@array)." Valids $valid_bases\n";
    
    my $base_count = ($strand == 1) ? $start : $start + ($valid_bases -1);
    
    foreach my $base (@array)
    {
        if ($base ne '-')
        {
            push (@masked_array, $base_count);
	      $base_count = ($strand == 1) ? $base_count + $inc  : $base_count - $inc; 

        }
        else
        {
            push (@masked_array, -1);
        }   
    }
    
    #print STDERR "@masked_array\n";
    
    return @masked_array;
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

sub query {
    my ($self, $seq) = @_;
    if ($seq) {
      unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq")) {
	$self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
      }

      $self->{_query} = $seq ;

      $self->filename($self->query->id.".$$.seq");
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

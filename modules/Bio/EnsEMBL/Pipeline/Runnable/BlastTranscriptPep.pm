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

Bio::EnsEMBL::Pipeline::Runnable::BlastTranscriptPep

=head1 SYNOPSIS

  my $blasttran = Bio::EnsEMBL::Pipeline::Runnable::BlastTranscriptPep->new ( 
    -genomic    => $genomic,
    -transcript => $tran,
    -database   => 'swir',
    -threshold  => 1e-6,
    -options    => 'B=1000')
                                                    );
  $blastgenscan->run();
  $blastgenscan->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on the 
peptide of the given transcript against a protein database.
The resulting blast hits are stored as DnaPepAlignFeatures.

=head1 CONTACT

B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BlastTranscriptPep;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::DnaPepAlignFeature;
#use Bio::EnsEMBL::Transcript;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-genomic     => $genomic
                           -transcript  => $tran,
                           -program     => 'blastp',
			   -database    => 'swir',
			   -threshold   => 1e-6);
			   
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::Runnable::BlastTranscriptPep object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::BlastTranscriptPep object
    Args    :   -genomic : Bio::Seq containing genomic dna
                -transcript     : Bio::EnsEMBL::Transcript
                -parameters : ref. to a hash of parameters

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_featurepairs'}= [];
  
  $self->{'_transcript'}     = undef;
  $self->{'_genomic'}     = undef;
  $self->{'_program'}     = undef;
  $self->{'_database'}    = undef;
  $self->{'_threshold'}   = undef;
  $self->{'_threshold_type'}   = undef;
  $self->{'_options'}     = undef;
    
  # Read the input parameters and set them in the object

  my ( $genomic,$transcript,$program,$database,$threshold,$threshold_type,$options) = 
    $self->_rearrange ([qw(GENOMIC TRANSCRIPT PROGRAM DATABASE THRESHOLD THRESHOLD_TYPE
			   OPTIONS)], @args);
  
  if (($genomic) && $genomic->isa("Bio::PrimarySeqI")) {
    $self->genomic($genomic);
  } else {
    $self->throw("No genomic sequence input");
  }
  
  if (($transcript) && $transcript->isa("Bio::EnsEMBL::Transcript")) {
      $self->transcript($transcript);
  } elsif ($transcript) {
      $self->throw("[$transcript] is not a Bio::EnsEMBL::Transcript");
  } else {
    $self->throw("No transcript input");
  }
  
  if ($program) {
    $self->program($self->find_executable($program));
  } else {
    $self->throw ("No program input");
  }
  
  if ($database) {
    $self->database($database);
  } else {
    $self->throw("No database ");
  }

  if ($threshold) {
    $self->threshold($threshold);
  } else {
    $self->threshold(0);
  }

  if ($threshold_type) {
    $self->threshold_type($threshold_type);
  } else {
    $self->threshold_type('PVALUE');
  }
  
  if ($options) {
    $self->options($options);
  } 
  
  return $self;
}


=head2 run

    Args       : none
    Description: runs the blast program and creates features
    Returntype : none
    Exceptions : none
    Caller     : general

=cut

sub run {
  my ($self) = @_;
  
  my $transcript = $self->transcript;

  if (!$transcript) {
    $self->throw("No transcript input");
  }
  
  if (not $transcript->translation) {
    print "Transcript has no translation; skipping\n";
    return;
  }

  my $peptide = Bio::PrimarySeq->new(-id         => 'Transcipt_translation',
				     -seq        => $transcript->translate->seq,
				     -moltype    => 'protein' );


  if ($peptide->length < 3) {
    print "Peptide too short (min length is 3); skipping transcript\n";
    return;
  }
  
  my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Blast
      (-query     => $peptide,
       -program   => $self->program,
       -database  => $self->database,
       -threshold => $self->threshold,
       -threshold_type => $self->threshold_type,
       -ungapped => 0,
       -options   => $self->options,
       -filter    => 1);

  $runnable->run();
  $self->align_hits_to_contig($runnable->output);

}


=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::BlastTranscriptPep->output()
    Returns :   An array of Bio::EnsEMBL::FeaturePair objects
    Args    :   none

=cut

sub output {
    my ($self) = @_;

    return $self->featurepairs();  
}


# This function creates a hash which is used to map between the exon genomic position
# and a position within the genscan predicted peptide. The hash is then matched
# against blast peptide hits to return a set of featurepairs of exons and blast
# peptides

sub align_hits_to_contig {
  my ( $self, @features )  = @_;

  for my $feature ( @features ) {
    my %exon_hash = ();
    # for each ungapped piece in it
    foreach my $ugFeature ( $feature->ungapped_features() ) {
      my $cdna_total = 1;
      #convert peptide coords to genomic coords
      my @split = $self->transcript->pep2genomic($ugFeature->start(),
					      $ugFeature->end());
      foreach my $gcoord ( @split ) {
	if($gcoord->isa('Bio::EnsEMBL::Mapper::Gap')) {
	  $cdna_total += $gcoord->end - $gcoord->start + 1;
	  next;
	}

	my $gstart = $gcoord->start;
	my $gend   = $gcoord->end;
	my $gstrand = $gcoord->strand;
	my $cdna_start = $cdna_total;
	my $cdna_end = $cdna_start + $gend - $gstart;
	$cdna_total += $gend - $gstart + 1;

	#determine which exon this genomic coordinate overlaps
	my $exon;
	foreach my $e (@{$self->transcript->get_all_Exons}) {
	  if($gstart >= $e->start && $gend <= $e->end) {
	    $exon = $e;
	    last;
	  }
	}

	# first, eat away non complete codons from start
	while(( $cdna_start - 1 ) % 3 != 0 ) {
	  $cdna_start++;
	  if( $gstrand == 1 ) {
	    $gstart++;
	  } else {
	    $gend--;
	  }
	}

	# and from end
	while( $cdna_end  % 3 != 0 ) {
	  $cdna_end--;
	  if( $gstrand == 1 ) {
	    $gend--;
	  } else {
	    $gstart++;
	  }
	}

	if( $cdna_end <= $cdna_start ) {
	  next;
	}

	my $pep_start = ($cdna_start+2)/3;
	my $pep_end = ( $cdna_end / 3 );

	my $fp = Bio::EnsEMBL::FeaturePair->new;

	$fp->start($gstart);
	$fp->end($gend);
	$fp->strand($gstrand);
	$fp->score($feature->score);
	$fp->p_value($feature->p_value);
	$fp->percent_id($feature->percent_id);
	$fp->analysis($feature->analysis);
	$fp->seqname($feature->seqname);

	$fp->hseqname($feature->hseqname);
	$fp->hstart($pep_start + $ugFeature->hstart() - 1);
	$fp->hend($pep_end + $ugFeature->hstart() - 1);
	$fp->slice($feature->slice);

	push( @{$exon_hash{$exon}}, $fp );
      }
    }

    # Take the pieces for each exon and make gapped feature
    foreach my $ex ( keys %exon_hash ) {
      my $dna_align_feature = Bio::EnsEMBL::DnaPepAlignFeature->new
	(-features => $exon_hash{$ex});
      $self->featurepairs($dna_align_feature);
    }
  }
}



# This function takes a blast transcript feature hit and a set of matching exons and
# creates a set of featurepairs aligned to genomic coordinates. It will split
# features if they cross exon boundaries
sub create_peptide_featurepairs {
    my ($self, $fp, @aligned_exons) = @_;
    #create featurepairs
    my @output_features;

    foreach my $ex_align (@aligned_exons)
    {
        my ($ex_start, $ex_end, $pep_start, $pep_end, $start_phase, $end_phase);

      #This splits features across multiple exons and records phases

      
      if ($ex_align->{'pep_start'}  < $fp->start) {
	#feature starts inside current exon
	if ($ex_align->{strand} == 1) {
	  $ex_start   = $ex_align->{'gen_start'}
	  + (($fp->start - $ex_align->{'pep_start'})*3);
	} else {
	  $ex_end     = $ex_align->{'gen_end'}
	  - (($fp->start - $ex_align->{'pep_start'})*3);
	}

	$start_phase = 0;
	$pep_start   = $fp->hstart;
      } else {
	#feature starts in a previous exon or absolute start of current exon

	if ($ex_align->{strand} == 1) {
	  $ex_start   = $ex_align->{'gen_start'};
	} else {
	  $ex_end     = $ex_align->{'gen_end'};
	}
	$start_phase = 0;
	$pep_start  = $fp->hstart + ($ex_align->{'pep_start'} - $fp->start);
      }
        
      if ($$ex_align{'pep_end'}    > $fp->end) {
	#feature ends in current exon
	if ($ex_align->{strand} == 1) {
	  $ex_end     = $ex_align->{'gen_start'}
	  + (($fp->end -  $ex_align->{'pep_start'})*3)+2;
	} else {
	  $ex_start   = $ex_align->{'gen_end'}
	  - (($fp->end -  $ex_align->{'pep_start'})*3)-2;
	}
	$end_phase  = 0;
	$pep_end    = $fp->hend;
      } else {
	#feature ends in a later exon or absolute end of current exon
	if ($ex_align->{strand} == 1) {
	  $ex_end     = $ex_align->{'gen_end'};
	} else {
	  $ex_start   = $ex_align->{'gen_start'};
	}
	$end_phase  = $ex_align->{'end_phase'};
	$pep_end    = $fp->hstart + ($ex_align->{'pep_end'} - $fp->start);
	
      }
        
      my $start_frac = $ex_align->{'phase'} + 1;
      my $end_frac   = (( 3 - $$ex_align{'end_phase'})%3) + 1;
      my $dna_feat = Bio::EnsEMBL::SeqFeature->new (-seqname    =>  $ex_align->{'name'},
						    -start      =>  $ex_start, 
						    -end        =>  $ex_end,
						    -strand     =>  $ex_align->{'strand'},
						    -score      =>  $fp->score,
						    -p_value    =>  $fp->p_value,
						    -percent_id =>  $fp->percent_id,
						    -analysis   =>  $fp->analysis,
						    -primary_tag=>  $fp->primary_tag,
						    -source_tag =>  $fp->source_tag, 
						    -phase      =>  $start_phase,  
						    -end_phase  =>  $end_phase );
      
	my $pep_feat = Bio::EnsEMBL::SeqFeature->new (-seqname    =>  $fp->hseqname,
							-start      =>  $pep_start,
							-end        =>  $pep_end,
							-strand     =>  $ex_align->{'strand'},
							-start_frac =>  $start_frac,
							-end_frac   =>  $end_frac,
							-score      =>  $fp->score,
							-p_value    =>  $fp->p_value,
							-percent_id =>  $fp->percent_id,
							-analysis   =>  $fp->analysis,
							-primary_tag=>  $fp->primary_tag,
							-source_tag =>  $fp->source_tag );
      
      my $featurepair = Bio::EnsEMBL::FeaturePair->new (-feature1   => $dna_feat,
							-feature2   => $pep_feat );

	$featurepair->slice($self->genomic);
	$featurepair->{'_exon_align'} = $ex_align;
	push( @output_features, $featurepair );
      }
    return @output_features;
}



# Get/set functions for the data follow here

=head2 genomic

    Title   :   genomic
    Usage   :   $self->genomic
    Function:   Get/set for the genomic dna sequence
    Returns :   Bio::PrimarySeqI
    Args    :   Bio::PrimarySeqI

=cut

sub genomic {
    my($self,$seq) = @_;
    
    if ($seq) {
      if (!($seq->isa("Bio::PrimarySeqI"))) {
        $self->throw("[$seq] is not a Bio::PrimarySeqI");
      }

      $self->{'_genomic'} = $seq;
    }

    return $self->{'_genomic'};
}

=head2 transcript

    Title   :   transcript
    Usage   :   $self->transcript($pep)
    Function:   Get/set for the transcript
    Returns :   Bio::EnsEMBL::Transcript
    Args    :   Bio::EnsEMBL::Transcript

=cut

sub transcript {
    my($self,$seq) = @_;
    
    if ($seq) {
      if (!($seq->isa("Bio::EnsEMBL::Transcript"))) {
	$self->throw("[$seq] is not a Bio::EnsEMBL::Transcript");
      }
      $self->{'_transcript'} = $seq;
    }

    return $self->{'_transcript'};
}

=head2 program

    Title   :   program
    Usage   :   $self->program('blastp');
    Function:   Get/set for the flavour of blast to run
    Returns :   String
    Args    :   String

=cut

sub program {
    my($self,$arg) = @_;
    
    if ($arg) {
      $self->{'_program'} = $arg;
    }

    return $self->{'_program'};
}

=head2 database

    Title   :   database
    Usage   :   $self->database('swir');
    Function:   Get/set for the database to search against
    Returns :   String
    Args    :   String

=cut

sub database {
    my($self,$arg) = @_;
    
    if ($arg) {
      $self->{'_database'} = $arg;
    }

    return $self->{'_database'};
}


=head2 featurepairs

    Title   :   featurepairs
    Usage   :   @out = $obj->featurepairs
    Function:   Returns the output featurepairs after they\'ve
                been converted back to genomic coords
                Also adds feature pairs to the array
    Returns :   @Bio::EnsEMBL::FeaturePair
    Args    :   Bio::EnsEMBL::FeaturePair

=cut


sub featurepairs {
    my ($self, $fp) = @_;
    if ($fp)
    {
        $self->throw("Input isn't a Bio::EnsEMBL::FeaturePair") 
                unless $fp->isa("Bio::EnsEMBL::FeaturePair");
        push (@{$self->{'_featurepairs'}}, $fp);
    }

    return @{$self->{'_featurepairs'}};
}

1;

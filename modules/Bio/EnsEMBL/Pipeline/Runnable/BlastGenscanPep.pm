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

Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep

=head1 SYNOPSIS

  my $blastgenscan = Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep->new ( 
                                                    -genomic    => $genomic,
                                                    -peptide    => $pep,
                                                    -database   => 'swir',
                                                    -threshold  => 1e-6,
                                                    -options    => 'B=1000'
                                                    );
  $blastgenscan->run();
  $blastgenscan->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides
constructed from assembling genscan predicted features to peptide
sequence. The resulting blast hits are stored as DnaPepAlignFeatures.

=head1 CONTACT

B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Pipeline::BioperlDBConf qw (
					      BIOPERLDB
					     );

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

=head2 new

    Title   :   new
    Usage   :   $self->new(-genomic     => $genomic
                           -peptide     => $peptide
                           -program     => 'blastp',
			   -database    => 'swir',
			   -threshold   => 1e-6);
			   
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep object
    Args    :   -genomic : Bio::Seq containing genomic dna
                -pep     : Bio::EnsEMBL::Transcript
                -parameters : ref. to a hash of parameters

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'_featurepairs'}= [];
  
  $self->{'_peptide'}     = undef;
  $self->{'_genomic'}     = undef;
  $self->{'_program'}     = undef;
  $self->{'_database'}    = undef;
  $self->{'_threshold'}   = undef;
  $self->{'_threshold_type'}   = undef;
  $self->{'_options'}     = undef;
    
  
  # Read the input parameters and set them in the object

  my ( $genomic,$peptide,$program,$database,$threshold,$threshold_type,$options) = 
    $self->_rearrange ([qw(GENOMIC PEPTIDE PROGRAM DATABASE THRESHOLD THRESHOLD_TYPE
			   OPTIONS)], @args);
  
  if (defined($genomic) && $genomic->isa("Bio::PrimarySeqI")) {
    $self->genomic($genomic);
  } else {
    $self->throw("No genomic sequence input");
  }
  
  if (defined($peptide) && $peptide->isa("Bio::EnsEMBL::PredictionTranscript")) {
      $self->peptide($peptide);
  } elsif (defined($peptide)) {
      $self->throw("[$peptide] is not a Bio::EnsEMBL::PredictionTranscript");
  } else {
    $self->throw("No peptide input");
  }
  
  if (defined($program)) {
    $self->program($self->find_executable($program));
  } else {
    $self->throw ("No program input");
  }
  
  if (defined($database)) {
    $self->database($database);
  } else {
    $self->throw("No database defined");
  }

  if (defined($threshold)) {
    $self->threshold($threshold);
  } else {
    $self->threshold(0);
  }

  if (defined($threshold_type)) {
    $self->threshold_type($threshold_type);
  } else {
    $self->threshold_type('PVALUE');
  }
  
  if (defined($options)) {
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
  
  my $transcript = $self->peptide;

  if (!defined($transcript)) {
    $self->throw("No peptide input");
  }
  
  print STDERR "Creating BioPrimarySeq for peptide ".$transcript->translate()."\n";

  my $peptide = Bio::PrimarySeq->new(-id         => 'Genscan_prediction',
				     -seq        => $transcript->translate->seq,
				     -moltype    => 'protein' );

  #print STDERR "Peptide length: ", $peptide->length, "\n";
  if ($peptide->length < 3) {
    print "Peptide too short (min length is 3); skipping transcript\n";
    return;
  }
  
  #print ::LOG "New BlastGenscanPep Runnable\n";
  #print ::LOG $transcript->_dump();
  #print ::LOG "\n";

  my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Blast  (-query     => $peptide,
							       -program   => $self->program,
							       -database  => $self->database,
							       -threshold => $self->threshold,
							       -threshold_type => $self->threshold_type,
							       -ungapped => 0,
							       -options   => $self->options,
                                                               -filter    => 1);
  $runnable->run();

  print "Output " . $runnable->output. "\n"; 

  $self->align_hits_to_contig($runnable->output);
  print "Output " . $runnable->output. "\n"; 
}


=head2 output

    Title   :   output
    Usage   :   $self->output();
    Function:   Runs Bio::EnsEMBL::Pipeline::Runnable::BlastGenscanPep->output()
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
      my @split = $self->peptide->pep2genomic($ugFeature->start(),
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
	foreach my $e (@{$self->peptide->get_all_Exons}) {
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
	$fp->contig($feature->contig);

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



# This function takes a blast peptide feature hit and a set of matching exons and
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

	$featurepair->attach_seq($self->genomic);
	$featurepair->{'_exon_align'} = $ex_align;
	push( @output_features, $featurepair );
      }
    return @output_features;
}

sub check_features {
  my ($self,$pep,@f) = @_;

  #print STDERR "Peptide is " . $pep . "\n";

  my %seqhash;

  foreach my $f (@f) {
    if (!defined($seqhash{$f->hseqname})) {
      my $seq = $self->get_Sequence($f->hseqname);
      $seqhash{$f->hseqname} = $seq;
    }
    my $seq = $seqhash{$f->hseqname};
    my $fdna = $f->seq;
    my $fpep = $fdna->translate('*','X',$f->phase);
    my $hpep = substr($seq->seq,$f->hstart-1,($f->hend - $f->hstart + 1));
    my $pos  = index($pep,$fpep->seq);

    #print "\tFeature " . $f->start . "\t" . $f->end . "\t" . $f->strand . "\t" . $f->phase . "\t" . $f->hstart . "\t" . $f->hend . " " . $pos . "\n" ;
     # print $fpep->seq . "\n" . $hpep . "\n";

  }
}

sub get_Sequence {
    my ($self,$id) = @_;


    next ID unless defined($id);

    #print(STDERR "Sequence id :  is [$id]\n");
    my $seq;
    if ($BIOPERLDB) {
      
      my $bpDBAdaptor = $self->bpDBAdaptor;
      my $seqfetcher = $bpDBAdaptor->fetch_BioSeqDatabase_by_name($self->database);
      eval {
	$seq = $seqfetcher->get_Seq_by_acc($id);
      };
      if ($@) {
	print STDERR ("Couldn't find sequence for [$id] in BlastGenscanPep");
	return;
      }
      
    }
    else {
      open(IN,"pfetch -q $id |") || $self->throw("Error fetching sequence for id [$id]");
	
      my $seqstr;
	
      while (<IN>) {
	chomp;
	$seqstr .= $_;
      }
    
    

      if (!defined($seqstr) || $seqstr eq "no match") {
        print STDERR ("Couldn't find sequence for [$id] in BlastGenscanPep");
	return;
      }
      
      $seq = new Bio::Seq(-id  => $id,
			     -seq => $seqstr);
    }

    #print (STDERR "Found sequence for $id [" . $seq->length() . "]\n");

    return $seq;
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
    
    if (defined($seq)) {
      if (!($seq->isa("Bio::PrimarySeqI"))) {
	$self->throw("[$seq] is not a Bio::PrimarySeqI");
      }

      $self->{'_genomic'} = $seq;
    }

    return $self->{'_genomic'};
}

=head2 peptide

    Title   :   peptide
    Usage   :   $self->peptide($pep)
    Function:   Get/set for the peptide transcript
    Returns :   Bio::EnsEMBL::Transcript
    Args    :   Bio::EnsEMBL::Transcript

=cut

sub peptide {
    my($self,$seq) = @_;
    
    if (defined($seq)) {
      if (!($seq->isa("Bio::EnsEMBL::PredictionTranscript"))) {
	$self->throw("[$seq] is not a Bio::EnsEMBL::PredicitionTranscript");
      }
      $self->{'_peptide'} = $seq;
    }

    return $self->{'_peptide'};
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
    
    if (defined($arg)) {
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
    
    if (defined($arg)) {
      $self->{'_database'} = $arg;
    }

    return $self->{'_database'};
}

=head2 threshold

    Title   :   threshold
    Usage   :   $self->threshold(100)
    Function:   Get/set for the score threshold to filter with
    Returns :   int
    Args    :   int

=cut


=head2 options

    Title   :   options
    Usage   :   $self->options('V=1000000');
    Function:   Get/set for the options to pass to the blast command
    Returns :   string
    Args    :   string

=cut


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

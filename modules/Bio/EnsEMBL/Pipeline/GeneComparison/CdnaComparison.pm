=head1 NAME - Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison.pm

=head1 DESCRIPTION

Takes a list of sequences and maps it to entries in a Blast nucleotide
database. Written with the aim of creating a mapping between Unigene
cluster IDs and Ensembl transcript IDs, this module may also prove
useful for other applications. Sequence comparison is performed using
Blastn.

=head1 SYNOPSIS

  Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison->new(
                  -CDNA_DATABASE       => blast_dbname,
                  -TRANSCRIPT_ARRAYREF => \@transcript_seqs);
  $cdna_comp->run_transcript_mapping;
  my %name_map = $cdna_comp->get_transcript_mapping;
  foreach my $transcript_id (keys %name_map) {
    print "$transcript_id\t", $name_map{$transcript_id}, "\n";
  }

=head1 CONTACT

Ensembl: ensembl-dev@ebi.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison;

use strict;

use vars qw(@ISA);
use IO::File;
use POSIX qw(tmpnam);
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Runnable::Blast;

@ISA = qw(Bio::EnsEMBL::Root);

=head2 new

    Title   :   new
    Usage   :   my $cdna_comp =
                Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison->new(
		  -CDNA_DATABASE       => blast_dbname,
                  -TRANSCRIPT_ARRAYREF => \@transcript_seqs);
    Function:   Initialises CdnaComparison object
    Returns :   A CdnaComparison object
    Args    :   Reference to array of PrimarySeq objects
                representing transcripts (-TRANSCRIPTS);
		name of Blast cDNA database (-CDNA_DATABASE);
		minimum percent identity of a FeaturePair for it
		to be considered (-PERCENT_ID);
		minimum threshold for coverage of the transcript
		by FeaturePairs (each of at least PERCENT_ID
		percent identity), given as total percentage of
		the transcript involved in such FeaturePairs
		(-TRANSCRIPT_PERCENT_COVERAGE, defaults to 97)

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($transcript_arrayref,
      $cdna_database,
      $transcript_percent_coverage,
      $percent_id)
    = $self->_rearrange([qw(TRANSCRIPT_ARRAYREF
                            CDNA_DATABASE
			    TRANSCRIPT_PERCENT_COVERAGE
			    PERCENT_ID)],
		        @args);

  $self->transcript_arrayref($transcript_arrayref);
  $self->cdna_database($cdna_database);
  $transcript_percent_coverage = 97
    unless defined $transcript_percent_coverage;
  $self->transcript_percent_coverage($transcript_percent_coverage);
  $percent_id = 97
    unless defined $percent_id;
  $self->percent_id($percent_id);

  return $self; # success - we hope!
}

=head2 transcript_arrayref

    Title   :   transcript_arrayref
    Usage   :   $cdna_comp->transcript_arrayref( \@transcript_seqs);
    Function:   get/set for reference to array of transcript
                PrimarySeqs
    Returns :   reference to array of PrimarySeq
    Args    :   optional reference to array of PrimarySeq

=cut

sub transcript_arrayref {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    $self->{_cdna_comparison_transcript_arrayref} = $value;
  }
  return $self->{_cdna_comparison_transcript_arrayref};
}

=head2 cdna_database

    Title   :   cdna_database
    Usage   :   $cdna_comp->cdna_database('unigene.seq');
    Function:   get/set for cDNA blast database file name
    Returns :   string
    Args    :   string

=cut

sub cdna_database {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    $self->{_cdna_comparison_cdna_database} = $value;
  }
  return $self->{_cdna_comparison_cdna_database};
}

=head2 get_transcript_mapping

    Title   :   get_transcript_mapping
    Usage   :   my $results = $cdna_comp->get_transcript_mapping;
    Function:   returns results of mapping (to be called after
                run_transcript_mapping)
    Returns :   hash (keys: transcript IDs; values: cDNA IDs)
    Args    :   none

=cut

sub get_transcript_mapping {
  my $self = shift;
  return %{$self->{_cdna_comparison_transcript_mapping_hashref}};
}

=head2 transcript_percent_coverage

    Title   :   transcript_percent_coverage
    Usage   :   $cdna_comp->transcript_percent_coverage(97.5);
    Function:   get/set for minimum acceptable percentage of a
                transcript involved in FeaturePairs with a cDNA
		(considering only those FeaturePairs with
		percent_id >= percent_id)
    Returns :   floating point value in the range [1..100]
    Args    :   optional floating point value in the range
                [1..100]

=cut

sub transcript_percent_coverage {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    $self->throw("expected a number in the range 1 to 100")
      if ($value < 1 or $value > 100);
    $self->{_cdna_comparison_transcript_percent_coverage} = $value;
  }
  return $self->{_cdna_comparison_transcript_percent_coverage};
}

=head2 percent_id

    Title   :   percent_id
    Usage   :   $cdna_comp->percent_id(97.5);
    Function:   get/set for minimum acceptable percentage of a
                transcript involved in FeaturePairs with a cDNA
    Returns :   floating point value in the range [1..100]
    Args    :   optional floating point value in the range
                [1..100]

=cut

sub percent_id {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    $self->throw("expected a number in the range 1 to 100")
      if ($value < 1 or $value > 100);
    $self->{_cdna_comparison_percent_id} = $value;
  }
  return $self->{_cdna_comparison_percent_id};
}

=head2 run_transcript_mapping

    Title   :   run_transcript_mapping
    Usage   :   my %mapping_by_transcripts =
                  $cdna_comp->run_transcript_mapping;
    Function:   Generate a mapping of transcript IDs to cDNA IDs.
                In this mapping, each transcript is represented
		either zero or one times, but each cDNA is
		represented zero, one, or more times. For each
		transcript, we choose the cDNA whose same-strand
		FeaturePairs cover the hit by
		transcript_percent_coverage or more. Where there
		is more than one such sequence, we chose the one
		of lowest p-value. Use get_transcript_mapping to
		access the results.
    Returns :   none
    Args    :   none

=cut

sub run_transcript_mapping {
  my ( $self ) = shift;

  my %mapping;
  foreach my $transcript_seq (@{$self->transcript_arrayref}) {
    # run sequence comparison
    my $blast_obj
      = Bio::EnsEMBL::Pipeline::Runnable::Blast->new(
                                  -query  => $transcript_seq,
				  -database => $self->cdna_database,
				  -program => 'wublastn',
				  -prune => '1');
    $blast_obj->run;
    my @features = $blast_obj->output;
    
    # store mapping for transcript if best coverage is good enough and
    # p-value is lower than other such sequences
    print STDERR "XXX run_transcript_mapping: ", scalar(@features),
      " features\n";
    if (@features) {	# Blast may have found none
      my %coverage_by_hits =
        %{$self->_get_coverage($transcript_seq->length, @features)};
      my %p_value_by_hits = %{$self->_get_min_p(@features)};
      my $mapped_name;
      my $min_p_for_hit = 9999999;	# arbitrary huge value
      my $coverage = 0;
      foreach my $hit (keys %coverage_by_hits) {
        if ($coverage_by_hits{$hit} >= $self->transcript_percent_coverage
          and $p_value_by_hits{$hit} < $min_p_for_hit) {
  	  $coverage = $coverage_by_hits{$hit};
          $mapped_name = $hit;
          $min_p_for_hit = $p_value_by_hits{$hit};
        }
      }
      if ($coverage >= $self->transcript_percent_coverage) {
        $mapping{$transcript_seq->display_id} = $mapped_name;
      }
    }
  }

  #store results
  $self->{_cdna_comparison_transcript_mapping_hashref} = \%mapping;
}

=head2 _get_min_p

    Title   :   _get_min_p
    Usage   :
    Function:   Records the minimum p-value for each cDNA,
                considering only same-strand FeaturePairs.
    Returns :   reference to a hash (keys: hit sequence name,
                values: floating point p-value)
    Args    :   array of FeaturePair objects

=cut

sub _get_min_p {
  my ( $self, @featurepairs ) = @_;
  $self->throw('interface fault') if (@_ < 2);

  my %hit_hash;	# keys: hit names
                # values: minimum p-values

  foreach my $featurepair (@featurepairs) {
    my $hseqname = $featurepair->hseqname;
    if (! $hit_hash{$hseqname}) {
      $hit_hash{$hseqname} = 999999;	# arbitrary huge value
    }
    $hit_hash{$hseqname} = $featurepair->p_value
         if $featurepair->strand == $featurepair->hstrand
        and $featurepair->p_value < $hit_hash{$hseqname};
  }

  return \%hit_hash;
}

=head2 _get_coverage

    Title   :   _get_coverage
    Usage   :
    Function:   Calculates percentage of a transcript covered by
                FeaturePairs involving one or more hit sequences,
		considering only same-strand FeaturePairs with
		percent IDs of at least percent_id. Each hit
		sequence may be involved in one or more
		FeaturePairs. Must be supplied with precisely all
		the FeaturePairs for a single transcript, across
		the cDNA database in use.  Results will not be
		distorted if the same area of transcript sequence
		is involved in many FeaturePairs, since we work
		with simple presence/absence of coverage for each
		base in the transcript. But if the same area of a
		cDNA is involved in many FeaturePairs, this could
		lead to high "coverage". (Not sure if this
		behaviour is correct or not.)
    Returns :   reference to a hash (keys: hit sequence name,
                values: floating point percentage coverage of
		transcript)
    Args    :   length of transcript (integer), array of
                FeaturePair objects for that transcript

=cut

sub _get_coverage {
  my ( $self, $cdna_length, @original_featurepairs ) = @_;
  $self->throw('interface fault') if (@_ < 3);

  # ignore FeaturePairs of too-low percent ID
  my @featurepairs;
  foreach my $featurepair (@original_featurepairs) {
    if ($featurepair->percent_id >= $self->percent_id) {
      print STDERR "XXX _get_coverage: using feature of percent_id ",
        $featurepair->percent_id, "\n";
      push @featurepairs, $featurepair;
    } else {
      print STDERR "XXX _get_coverage: ignoring feature of percent_id ",
        $featurepair->percent_id, "\n";
    }
  }

  my %hit_hash;	# keys: hit names
                # values: refs to coverage vectors (in which undef
		# means 'not covered')

  # create coverage vectors
  foreach my $featurepair (@featurepairs) {
    my $hseqname = $featurepair->hseqname;
    if (! $hit_hash{$hseqname}) {
      $hit_hash{$hseqname} = [];
    }
    if ($featurepair->strand == $featurepair->hstrand) {
      for (my $i = $featurepair->start; $i <= $featurepair->end; $i++) {
        $hit_hash{$hseqname}->[$i] = 1;
      }
    }
  }
  
  my %coverage;	# keys: hit names
                # values: percentage of cDNA bases involved in HSPs
		# with the named hit

  # calculate coverage percentages
  foreach my $hseqname (keys %hit_hash) {
    my $bases_covered = 0;
    my $coverage_vector_arr_ref = $hit_hash{$hseqname};
    for (my $i = 1; $i <= $cdna_length; $i++) {
      $bases_covered += 1 if defined $$coverage_vector_arr_ref[$i];
    }
    my $cdna_percent = 100 * $bases_covered / $cdna_length;
    $coverage{$hseqname} = $cdna_percent;
  }

  return \%coverage;
}

1;

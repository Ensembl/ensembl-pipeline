=head1 NAME - Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison.pm

=head1 DESCRIPTION

Takes two lists of sequences and maps one to the other using sequence
similarity at the DNA level. Written with the aim of assigning
Unigene cluster IDs to Ensembl transcript IDs, this module may also
prove useful for other applications. Sequence comparison is
performed using Exonerate with intron modelling switched off.

=head1 SYNOPSIS

  my $cdna_comp =
    Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison->new(
      -TRANSCRIPTS            => \@ensembl_transcript_seqs,
      -CDNAS                  => \@unigene_cluster_representative_seqs,
      -TRANSCRIPT_PERCENT_COVERAGE => 97.5);

  $cdna_comp->run_transcript_mapping;
  my %mapping_by_transcripts = $cdna_comp->transcript_mapping;

=head1 CONTACT

Ensembl: ensembl-dev@ebi.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison;

use strict;

use constant EXECUTABLE => '/acari/work2/gs2/gs2/bin/exonerate-0.3d';
use constant ARGS => '-a 400 -p no -w 14 -t 65 -H 150 -D 5 -m 768';

use vars qw(@ISA);
use IO::File;
use POSIX qw(tmpnam);
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;

@ISA = qw(Bio::EnsEMBL::Root);

=head2 new

    Title   :   new
    Usage   :   my $cdna_comparison
                = Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison->new(
                  -TRANSCRIPT_ARRAYREF => \@ensembl_transcript_seqs,
                  -CDNA_ARRAYREF       => \@unigene_representative_seqs,
                  -TRANSCRIPT_PERCENT_COVERAGE => 97.5 );

    Function:   Initialises CdnaComparison object
    Returns :   A CdnaComparison object
    Args    :   Reference to array of PrimarySeq objects representing
                transcripts (-TRANSCRIPTS);
		reference to array of PrimarySeq objects representing
		cDNAs (-CDNAS);
		threshold for coverage of the transcript by
		featurepairs, given as total percentage of the
		transcript involved in any featurepairs
		(-TRANSCRIPT_PERCENT_COVERAGE, defaults to 95)

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($transcript_arrayref,
      $cdna_arrayref,
      $transcript_percent_coverage)
    = $self->_rearrange([qw(TRANSCRIPT_ARRAYREF
                            CDNA_ARRAYREF
			    TRANSCRIPT_PERCENT_COVERAGE)],
		        @args);

  $self->transcript_arrayref($transcript_arrayref);
  $self->cdna_arrayref($cdna_arrayref);
  $transcript_percent_coverage = 95
    unless defined $transcript_percent_coverage;
  $self->transcript_percent_coverage($transcript_percent_coverage);

  return $self; # success - we hope!
}

=head2 transcript_arrayref

    Title   :   transcript_arrayref
    Usage   :   $cdna_comp->transcript_arrayref(
                  \@ensembl_transcript_seqs);
    Function:   get/set for reference to array of transcript
                PrimarySeqs

=cut

sub transcript_arrayref {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    $self->{_cdna_comparison_transcript_arrayref} = $value;
  }
  return $self->{_cdna_comparison_transcript_arrayref};
}

=head2 cdna_arrayref

    Title   :   cdna_arrayref
    Usage   :   $cdna_comp->cdna_arrayref(\@cdna_seqs);
    Function:   get/set for reference to array of cdna PrimarySeqs

=cut

sub cdna_arrayref {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    $self->{_cdna_comparison_cdna_arrayref} = $value;
  }
  return $self->{_cdna_comparison_cdna_arrayref};
}

=head2 transcript_mapping

    Title   :   transcript_mapping
    Usage   :   my $results = $cdna_comp->transcript_mapping;
    Function:   returns results of run_transcript_mapping
    Returns :   hash (keys: transcript IDs; values: cDNA IDs)

=cut

sub transcript_mapping {
  my $self = shift;
  return %$self->_cdna_comparison_transcript_mapping_hashref;
}

=head2 transcript_percent_coverage

    Title   :   transcript_percent_coverage
    Usage   :   $cdna_comp->transcript_percent_coverage(97.5);
    Function:   get/set for threshold percentage of a transcript
                involved in featurepairs with a cDNA

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

=head2 _get_coverage

    Title   :   _get_coverage
    Usage   :
    Function:   Calculates percentage of a transcript covered by
                featurepairs involving one or more hit sequences,
		where each hit sequence may be involved in one or
		more featurepairs. Intended to give correct
		coverage even with overlapping hits, but must be
		supplied with precisely all the FeaturePairs for
		a single genomic transcript across the cDNA
		database in use.
    Returns :   reference to a hash (keys: hit sequence name,
                values: floating point percentage coverage of
		transcript)
    Args    :   length of transcript (integer), array of
                FeaturePair objects for that transcript

=cut

sub _get_coverage {
  my ( $self, $cdna_length, @featurepairs ) = @_;
  $self->throw('interface fault') if (@_ < 3);

  my %hit_hash;	# keys: hit names
                # values: refs to coverage vectors (in which undef
		# means 'not covered')

  # create coverage vectors
  foreach my $featurepair (@featurepairs) {
    my $hseqname = $featurepair->hseqname;
    if (! $hit_hash{$hseqname}) {
      $hit_hash{$hseqname} = [];
    }
    for (my $i = $featurepair->start; $i <= $featurepair->end; $i++) {
      $hit_hash{$hseqname}->[$i] = 1;
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

=head2 run_transcript_mapping

    Title   :   run_transcript_mapping
    Usage   :   my %mapping_by_transcripts =
                  $cdna_comp->run_transcript_mapping;
    Function:   Generate a mapping of transcript IDs to cDNA IDs.
                In this mapping, each transcript is represented
		either zero or one times, but each cDNA is
		represented zero, one, or more times.
    Returns :   none (use get_transcript_mapping to obtain results)
    Args    :   none

=cut

sub run_transcript_mapping {
  my ( $self ) = shift;

  # obtain +- unique 'genomic' (actually, transcript) file name
  my $genomic_fnam;
  my $genomic_fh;
  do { $genomic_fnam = tmpnam() }
    until $genomic_fh =
      IO::File->new($genomic_fnam, O_RDWR|O_CREAT|O_EXCL);
  close $genomic_fh or die "error closing file $genomic_fnam";
  unlink $genomic_fnam or die "error deleting file $genomic_fnam";

  # write temporary (but uniquely named) file of all cDNAs
  # we only have to do this once
  my $cdna_fnam;
  my $cdna_fh;
  do { $cdna_fnam = tmpnam() }
    until $cdna_fh = IO::File->new($cdna_fnam, O_RDWR|O_CREAT|O_EXCL);
  my $cdna_stream = Bio::SeqIO->new(-fh     => $cdna_fh,
                                    -format => 'Fasta');
  my $cdna_seqs_arrayref = $self->cdna_arrayref;
  foreach my $cdna_seq (@$cdna_seqs_arrayref) {
    $self->throw("error writing to $cdna_fnam")
      unless $cdna_stream->write_seq($cdna_seq);
  }
  close $cdna_fh or die "error closing file $cdna_fnam";

  # generate mapping

  my %mapping;
  foreach my $genomic_seq (@$self->transcript_arrayref) {

    # create the 'genomic' file containing the current transcript
    $genomic_fh = IO::File->new($genomic_fnam, O_RDWR|O_CREAT|O_EXCL);
    $self->throw("temporary file name $genomic_fnam not unique")
      unless $genomic_fh;
    my $genomic_stream = Bio::SeqIO->new(-fh     => $genomic_fh,
                                         -format => 'Fasta');
    $self->throw("error writing to $genomic_fnam")
      unless $genomic_stream->write_seq($genomic_seq);
    close $genomic_fh or die "error closing file $genomic_fh";

    # run sequence comparison
    my $exonerate_obj
      = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
                                           -genomic   => $genomic_fnam,
                                           -est       => $cdna_fnam,
					   -arguments => ARGS
                                          );
    $exonerate_obj->exonerate(EXECUTABLE);
    $exonerate_obj->run();
    my @features = $exonerate_obj->output;
    
    # store mapping for transcript if best coverage is good enough
    my %coverage_by_hits = %$self->_get_coverage($genomic_seq->length,
                                                 @features);
    my $max_coverage = 0.0;
    my $mapped_name;
    foreach my $hit (keys %coverage_by_hits) {
      if ($coverage_by_hits{$hit} > $max_coverage) {
        $mapped_name = $hit;
        $max_coverage = $coverage_by_hits{$hit};
      }
    }
    if ($max_coverage >= $self->transcript_percent_coverage) {
      $mapping{$genomic_seq->name} = $mapped_name;
    }

    unlink $genomic_fnam or die "error deleting file $genomic_fnam";
  }

  unlink $cdna_fnam or die "error deleting file $cdna_fnam";

  #store results
  $self->_cdna_comparison_transcript_mapping_hashref = \%mapping;
}

1;

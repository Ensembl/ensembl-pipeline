=head1 NAME - Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison.pm

=head1 DESCRIPTION

Takes two lists of sequences and maps one to the other using sequence similarity at
the DNA level. Useful for, among other things, mapping Unigene cluster IDs to
Ensembl transcript IDs.

=head1 SYNOPSIS

  my $cdna_comparison
    = Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison->new(
      -TRANSCRIPTS        => \@ensembl_transcript_seqs,
      -CDNAS              => \@unigene_cluster_representative_seqs,
      -TRANSCRIPT_PERC_ID => 99);

  my %mapping_by_transcripts = $cdna_comparison->get_transcript_mapping;
  # my %mapping_by_cdnas = $cdna_comparison->get_cdna_mapping;

  my $unmapped_ensembl_transcripts = $cdna_comparison->get_unmapped_transcripts;
  my $unmapped_ensembl_transcripts = $cdna_comparison->get_unmapped_cdnas;

=head1 CONTACT

Ensembl: ensembl-dev@ebi.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison;

use vars qw(@ISA);
use strict;
use POSIX::tmpnam;
use IO::File;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;

@ISA = qw(Bio::EnsEMBL::Root);

=head2 new

    Title   :   new
    Usage   :   my $cdna_comparison
                = Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison->new(
                 -TRANSCRIPT_ARRAYREF => \@ensembl_transcript_seqs,
                 -CDNA_ARRAYREF       => \@unigene_representative_seqs,
                 -TRANSCRIPT_PERC_COVERAGE => 99 );

    Function:   Initialises CdnaComparison object
    Returns :   A CdnaComparison object
    Args    :   Reference to array of PrimarySeq objects representing
                transcripts (-TRANSCRIPTS);
		reference to array of PrimarySeq objects representing
		cDNAs (-CDNAS);
		threshold for coverage of the transcript by HSPs, given as
		total percentage of the transcript involved in any HSPs
		(-TRANSCRIPT_PERC_COVERAGE).

=cut

=head2 transcript_arrayref

    Title   :   transcript_arrayref
    Usage   :   $cdna_comp->transcript_arrayref(\@ensembl_transcript_seqs);
    Function:   get/set for reference to array of transcript PrimarySeqs

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

=head2 transcript_perc_coverage

    Title   :   transcript_perc_coverage
    Usage   :   $cdna_comp->transcript_perc_coverage(99);
    Function:   get/set for 

=cut

sub transcript_perc_coverage {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    $self->throw("expected a number in the range 1 to 100")
      if ($value < 1 or $value > 100);
    $self->{_cdna_comparison_transcript_perc_coverage} = $value;
  }
  return $self->{_cdna_comparison_transcript_perc_coverage};
}

=head2 mapping_by_cdnas

    Title   :   mapping_by_cdnas
    Usage   :   N/A
    Function:   Not yet implemented

=cut

sub mapping_by_cdnas {
  my $self = shift;
  $self->throw('not yet implemented!');
}

sub get_transcript_mapping {
  my $self = shift;

  # create uniquely-named genomic file
  my $genomic_fnam;
  my $genomic_fh;
  do { $genomic_fnam = POSIX::tmpnam() }
    until $genomic_fh = IO::File->new($genomic_fnam, O_WR|O_CREAT|O_EXCL);
  close $genomic_fh or die "error closing file $genomic_fnam";

  # create and open uniquely named cDNA file
  my $cdna_fnam;
  my $cdna_fh;
  do { $cdna_fnam = POSIX::tmpnam() }
    until $cdna_fh = IO::File->new($cdna_fnam, O_WR|O_CREAT|O_EXCL);

  # write cDNA file
  my $cdna_stream = Bio::SeqIO->new(-fh => $cdna_fh, -format => 'Fasta');
  my $cdna_seqs_arrayref = $self->cdna_arrayref;
  foreach my $cdna_seq (@$cdna_seqs_arrayref) {
    $self->throw "error writing to $cdna_fnam"
      if ! $cdna_stream->write_seq($cdna_seq);
  }
  close $cdna_fh or die "error closing file $cdna_fnam";

  my $genomic_seqs_arrayref = $self->transcript_arrayref;
  foreach my $genomic_seq (@$genomic_seqs_arrayref) {
    $genomic_fh = open ">$genomic_fnam";
    my $genomic_stream = Bio::SeqIO->new(-fh => $genomic_fh, -format => 'Fasta');
    $self->throw "error writing to $genomic_fnam"
      if ! $genomic_stream->write_seq($genomic_seq);
    my $exonerate_obj = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
                                              -genomic => $genomic_fnam,
                                              -est     => $cdna_fnam
                                              );
    $exonerate_obj->run();

    my @features = $exonerate_obj->output;


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

Bio::EnsEMBL::Pipeline::RunnableDB::BlastTranscriptPep

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBAdaptor->new($locator);
my $btpep     = Bio::EnsEMBL::Pipeline::RunnableDB::BlastTranscriptPep->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );

$btpep->fetch_input();
$btpep->run();
$btpep->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Pipeline::Runnable::Blast on peptides
obtained by translating a representative transcript from each gene
in the region. The resulting blast hits are written back as
DnaPepAlignFeature's.

The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _'

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::BlastTranscriptPep;

use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastTranscriptPep;
use Bio::PrimarySeq;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);


=head2 fetch_input

  Args       : none
  Example    : $runnable->fetch_input
  Description: Fetches input data for BlastTranscriptPep and makes runnable
  Returntype : none
  Exceptions : $self->input_id is not defined
  Caller     : run_RunnableDB, Bio::EnsEMBL::Pipeline::Job

=cut

sub fetch_input {
  my($self) = @_;

  $self->throw("No input id") unless defined($self->input_id);
  
  $self->fetch_sequence;
  
  my (@representative_trans);
  my $genes = $self->db->get_GeneAdaptor->fetch_all_by_Slice($self->query);
  # get longest transcript; only sensible to do it for that
  foreach my $gene (@$genes) {
    my ($longest_pep_tran, $longest_pep_tran_len);
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      if ($tran->translation) {
        my $pep_string = $tran->translate->seq;
        
        if (not defined $longest_pep_tran or length($pep_string) > $longest_pep_tran_len) {
          $longest_pep_tran = $tran;
          $longest_pep_tran_len = length($pep_string);
        }
      }
    }
    push @representative_trans, $longest_pep_tran;
  }
    
  my ($thr, $thr_type);
  my %p = $self->parameter_hash;
  
  if (defined $p{-threshold} && defined $p{-threshold_type}) {
    $thr      = $p{-threshold};
    $thr_type = $p{-threshold_type};
  }
  else {
    $thr_type = 'PVALUE';
    $thr      = 0.001;
  }
  
  foreach my $t (@representative_trans) {
    foreach my $db (split ',', ($self->analysis->db_file)) {
      $self->runnable(Bio::EnsEMBL::Pipeline::Runnable::BlastTranscriptPep
                        ->new(
                              -genomic        => $self->query,
                              -transcript     => $t,
                              -database       => $db,
                              -program        => $self->analysis->program_file,
                              -options        => $self->arguments || undef,
                              -threshold      => $thr,
                              -threshold_type => $thr_type
                             ));
    }
  }
  return 1;
}


1;

=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

AlignmentSquizz - 

=head1 SYNOPSIS


=head1 DESCRIPTION

  Simple hacked module for looking at alignment quality, especially
  coverage and identity.  Designed for finding odd alignments.
  Use the script ensembl-pipeline/scripts/alignment_squiz.pl to 
  drive this module.

=head1 METHODS

=cut


package AlignmentSquizz;

sub new {
  my ($class, $align) = @_;

  my $self = bless {}, $class;

  unless (defined $align && $align->isa("Bio::EnsEMBL::Pipeline::Alignment")){
    die "Alignment object is not an Ensembl pipeline " .
      "alignment object, it is a [$align]";
  }

  $self->{_alignment} = $align;

  return $self
}


sub appraise {
  my $self = shift;

  my $align_seqs = $self->{_alignment}->fetch_AlignmentSeqs;

  foreach my $align_seq (@$align_seqs) {

    if ($align_seq->name eq 'genomic_sequence'){
      $self->{_genomic_sequence} = $align_seq;
      next
    } elsif ($align_seq->name eq 'exon_sequence') {
      $self->{_exon_sequence} = $align_seq;
      next
    } elsif ($align_seq->name eq 'translated_exon_sequence') {
      $self->{_translated_exon_sequence} = $align_seq;
      next
    }

    $self->_sequence_comparison($align_seq);

  }

  return 1
}


sub _sequence_comparison {
  my ($self, $align_seq) = @_;

  my $master_seq;

  if ($align_seq->type eq 'nucleotide'){
    $master_seq = $self->{_exon_sequence}
  } elsif ($align_seq->type eq 'protein') {
    $master_seq = $self->{_translated_exon_sequence}
  } else {
    die "Unrecognised sequence type [".$align_seq->type."]"
  }

  my $alignment_length = $master_seq->length;
  my $matched_bases = 0;
  my $unmatched_bases = 0;
  my $identical_bases = 0;
  my $total_evidence_bases = 0;


  for (my $i = 1; $i <= $alignment_length; $i++){

    my $master_base = substr($master_seq->seq, $i, 1);
    my $aligned_base = substr($align_seq->seq, $i, 1);

    if ($aligned_base ne '-') {
      $total_evidence_bases++;
      if ($master_base ne '-') {
	$matched_bases++;
	if ($master_base eq $aligned_base){
	  $identical_bases++
	}
      } else {
	$unmatched_bases++
      }
    }
  }

  my $evidence_length = $matched_bases + $unmatched_bases;
  my $evidence_coverage = sprintf("%3.2f",$matched_bases/$evidence_length);
  my $uncovered_evidence = sprintf("%3.2f",$unmatched_bases/$evidence_length);
  my $identity = sprintf("%3.2f",$identical_bases/$evidence_length);

  push @{$self->{_output}}, [$align_seq->name,
			     $align_seq->type,
			     $evidence_length,
			     $matched_bases,
			     $evidence_coverage,
			     $unmatched_bases,
			     $uncovered_evidence,
			     $identical_bases,
			     $identity];

}


sub print_tabulated_output {
  my $self = shift;

  my @output = @{$self->{_output}};

  while (my $row = shift @output){
    print STDOUT join("\t", @$row) . "\n";
  }

  return 1;
}

sub output_array {
  return $_[0]->{_output}
}



1;


# Simple hacked module for looking at alignment quality, especially
# coverage and identity.  Designed for finding odd alignments.

# Use the script ensembl-pipeline/scripts/alignment_squiz.pl to 
# drive this module.


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

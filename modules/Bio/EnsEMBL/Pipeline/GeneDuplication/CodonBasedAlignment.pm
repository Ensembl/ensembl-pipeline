package Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::Clustalw;

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  return $self;
}

sub filename {
  my $self = shift;

  $self->{filename} = shift if @_;

  return $self->{filename}
}

sub sequences_from_file {
  my $self = shift;

  my $seqio = Bio::SeqIO->new(-file   => $self->filename,
			      -format => 'fasta');

  my @seqs;

  while (my $seq = $seqio->next_seq){
    push @seqs, $seq;
  }

  $self->sequences(\@seqs);

  return 1
}

sub sequences {
  my $self = shift;

  if (@_) {
    $self->{_sequences} = shift;
  }

  return $self->{_sequences}
}


### Alignment Stuff ###

sub run_alignment {
  my ($self) = @_;

  my $clustalw = Bio::Tools::Run::Alignment::Clustalw->new();

  my $seqs = $self->sequences;

  unless (scalar @$seqs > 1) {
    die "It takes more than one sequence to clustalw";
    return 0
  }

  my %nt_seqs;

  foreach my $nt_seq (@$seqs){
    $nt_seqs{$nt_seq->display_id} = $nt_seq; 
  }

  my @aa_seqs;
  foreach my $seq (@$seqs) {
    push (@aa_seqs, $seq->translate(undef, undef, undef, 1));  ### !!! UNIVERSAL TRANSLATION !!! ###
#    push (@aa_seqs, $seq->translate(undef, undef, undef, 2));  ### !!! VERT MITO TRANSLATION !!! ###
  }

  my $alignment = $clustalw->align(\@aa_seqs);

  my @aligned_seqs;

  my $longest_seq = 0;
  foreach my $aligned_seq ($alignment->each_seq){
#print "This seq is " . $aligned_seq->length . " long\n";
    push(@aligned_seqs, $aligned_seq);
    $longest_seq = $aligned_seq->length if $aligned_seq->length > $longest_seq;
#print "Longest seq : " . $longest_seq . "\n";
  }

print "Have " . scalar @aligned_seqs . " aligned sequences.\n";

  foreach my $aligned_seq (@aligned_seqs) {
#print "This seq is " . scalar $aligned_seq->length . " long.\n";
    my $nt_work_seq = $nt_seqs{$aligned_seq->display_id}->seq;
    $nt_work_seq =~ s/(...)/$1:/g;
    my @nt_seq_array = split /:/, $nt_work_seq;
    my @aa_seq_array = split //, $aligned_seq->seq;

    while (scalar @aa_seq_array < $longest_seq){
      push(@aa_seq_array, '-');
      print "Adding gap to make sequences same length.\n";
    }

#    $self->throw("nt and aa sequences are not of corresponding lengths.") 
#      if scalar @nt_seq_array != scalar @aa_seq_array;

    my $aligned_nt_string = '';

    foreach my $aa (@aa_seq_array) {      
      if ($aa eq '.' || $aa eq '-'){
	$aligned_nt_string .= '---';
      } else {
	$aligned_nt_string .= shift @nt_seq_array;
      }
    }

    my $aligned_nt_bioseq = Bio::Seq->new(-display_id => $aligned_seq->display_id,
					  -seq        => $aligned_nt_string);

    $self->_alignment('add', $aligned_nt_bioseq);
  }

  return $self->_alignment;
}

sub _alignment {
  my ($self, $action, $seq) = @_;

  if (defined $action && $action eq 'add') {
    push (@{$self->{_alignment}}, $seq);
  } elsif (defined $action && $action eq 'clear') {
    $self->{_alignment} = [];
  } elsif (defined $action) {
    $self->warn("HomologTool::_alignment - unknown action specified ".
		"[$action].  Must be \'add\' or \'clear\'");
  }
  
  return $self->{_alignment};
}

return 1;

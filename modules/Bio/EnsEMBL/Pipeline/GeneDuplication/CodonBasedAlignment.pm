# Cared for by Dan Andrews <dta@sanger.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment

=head1 SYNOPSIS

use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;

# Create our object, specifying the genetic code needed to 
# translate our sequences.  The commonest codes are: 
#   universal                => 1 
#   vertebrate mitochondrial => 2
# These numbers are the same as those you need to specify to 
# translate any Bio::Seq object.  Hence if you need a truly 
# oddball genetic code, check the Bio::Seq module documentation.
my $cba =
  Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment->new(
    -genetic_code => 1);

# Add some sequence to our alignment module.  Pass a reference to
# an array of Bio::Seq objects.  Everything is going to fall apart
# if you pass amino acid sequences here, some make sure they are
# nucleotide sequences.

$cba->sequences(\@seqs);

# Run the actual alignment process.  The return value is an 
# array of Bio::Seq objects, that when dumped to a multiple 
# fasta file or the like, will display a multiple alignment.

my $align = $cba->run_alignment;

foreach my $seq (@$align){
  my $string = $seq->seq;
  $string =~ s/(.{60})/$1\n/g;

  print STDOUT ">" . $seq->display_id . "\n$string";
}



=head1 DESCRIPTION

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut



package Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;

use strict;
use Bio::EnsEMBL::Root;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($genetic_code
     ) = rearrange([qw(GENETIC_CODE
		      )], @args);

  unless (defined $genetic_code){
    warning("Genetic code not specified.  Defaulting " . 
	    "to the so-called universal code.");
    $genetic_code = 1;
  }

  $self->genetic_code($genetic_code);

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
    push (@aa_seqs, $seq->translate(undef, undef, undef, $self->genetic_code));
  }

  my $alignment = $clustalw->align(\@aa_seqs);

  my @aligned_seqs;

  my $longest_seq = 0;
  foreach my $aligned_seq ($alignment->each_seq){
    push(@aligned_seqs, $aligned_seq);
    $longest_seq = $aligned_seq->length if $aligned_seq->length > $longest_seq;
  }

  foreach my $aligned_seq (@aligned_seqs) {
    my $nt_work_seq = $nt_seqs{$aligned_seq->display_id}->seq;
    $nt_work_seq =~ s/(...)/$1:/g;
    my @nt_seq_array = split /:/, $nt_work_seq;
    my @aa_seq_array = split //, $aligned_seq->seq;

    while (scalar @aa_seq_array < $longest_seq){
      push(@aa_seq_array, '-');
      print "Adding gap to make sequences same length.\n";
    }

#    throw("nt and aa sequences are not of corresponding lengths.") 
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
    warning("CodonBasedAlignment::_alignment - unknown action specified ".
		"[$action].  Must be \'add\' or \'clear\'");
  }
  
  return $self->{_alignment};
}

sub genetic_code {
  my $self = shift;

  if (@_){
    $self->{_genetic_code} = shift;
  }
 
  throw("Genetic code not specified.") 
    unless $self->{_genetic_code};

  return $self->{_genetic_code}
}

return 1;

#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new(
							      -executable => $exe
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using pfetch.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::FileIndex;

use strict;
use Bio::EnsEMBL::Root;
use Bio::DB::RandomAccessI;
use Bio::SeqIO;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($seqfile) = $self->_rearrange(['SEQFILE'], @args);

  if (defined($seqfile)) {
    $self->seqfile($seqfile);
  } else {
    $self->throw("Need a sequence file");
  }
  
  return $self; # success - we hope!
}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_eq_by_acc($accession);
  Function: Does the sequence retrieval
  Returns : Bio::Seq
  Args    : 

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;
  
  if (!defined($acc)) {
    $self->throw("No accession input");
  }  
  
  if (defined($self->{_seqhash}{$acc})) {
     return $self->{_seqhash}{$acc};
  } else {
     $self->throw("Could not fetch sequence for [$acc]\n");
  } 
}

sub seqfile {
  my ($self,$file) = @_;


  if (defined($file)) {
    open (INDEX,"<$file") || $self->throw("Can't open $file");
    
    my $seqio = new Bio::SeqIO(-fh => \*INDEX, -format => 'fasta');
    
    while (my $seq = $seqio->next_seq) {
      $self->{_seqhash}{$seq->id} = $seq;
    }
    
    close(INDEX);
  } else {
    $self->throw("Must supply a seqfile to Bio:EnsEMBL::Pipeline::SeqFetcher::FileIndex");
  }
}
  
sub list_all_ids {
  my ($self) = @_;

  my @ids = keys %{$self->{_seqhash}};

  return \@ids;
}

1;

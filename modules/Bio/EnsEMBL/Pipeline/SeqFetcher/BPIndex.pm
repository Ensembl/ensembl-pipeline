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

  Bio::EnsEMBL::Pipeline::SeqFetcher::BPIndex

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::BPIndex->new(
							       '-index'  => $index,
							       '-format' => 'Fasta',
							      );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

    Object to retrieve sequences as Bio::Seq, from a bioperl index. 
    The index is not made by this module; instead, the absolute path 
    to the bioperl index must be set using $obj->bp_index. The format 
    of the database must be set using bp_format.


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::BPIndex;

use strict;
use Bio::EnsEMBL::Root;
use Bio::DB::RandomAccessI;
use Bio::Seq;
# change these uses as new indices are released in bioperl - these are for 0.6.2
use Bio::Index::Fasta;
use Bio::Index::EMBL;
use Bio::Index::SwissPfam;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  
  my($index, $format) = $self->_rearrange(['INDEX',
					   'FORMAT'], @args);
  
  if (!defined $index) {
    $self->throw("No bioperl indexfile provided to BPIndex\n");
  }
  $self->bp_index($index);
  
  
  if (!defined $format) {
    $self->throw("No bioperl index format provided to BPIndex\n");
  }
  $self->bp_format($format);


  return $self; # success - we hope!
}

=head2 bp_index

    Title   :   bp_index
    Usage   :   $self->bp_index('/usr/local/ensembl/data/bp.inx')
    Function:   Get/set for a bioperl index
    Returns :   path to bioperl index
    Args    :   path to bioperl index

=cut

sub bp_index {
  
  my ($self, $inx) = @_;
  if ($inx)
    {
      $self->{'_inx'} = $inx;
    }
  return $self->{'_inx'};  
  
}

=head2 bp_format

    Title   :   bp_format
    Usage   :   $self->bp_format('Fasta')
    Function:   Get/set for a bioperl index format
    Returns :   String representing format. NOTE - bp_format is used in run to identify the type of Bio::Index module to make - so case is crucial. eg Fasta, EMBL, SwissPfam
    Args    :   String representing format

=cut

sub bp_format {
  
  my ($self, $format) = @_;
  if ($format)
    {
      $self->{'_format'} = $format;
    }
  return $self->{'_format'};  
  
}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_Seq_by_acc($accession);
  Function: Does the sequence retrieval
  Returns : Bio::Seq
  Args    : 

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;
  
  my $inx    = $self->bp_index;
  my $format = $self->bp_format;
  

  if (!defined($acc)) {
    $self->throw("No accession input");
  }  

  if (!defined($inx)) {
    $self->throw("No search index specified; cannot run");
  }  

  if (!defined($format)) {
    $self->throw("No index format specified");
  }  
  
  my $type = 'Bio::Index::' . $format;
  my $index;

  eval {
    $index = $type->new($inx);
  };

  if ($@) {
    my $tmp = $@; # for some reason, warn empties out $@ ...
    $self->warn("Problem opening the index [$inx] - check you have supplied the right format!");
    $self->throw ("[$tmp]!");
  }
  
  # get the sequence
  my $seq;
  eval{
    $seq = $index->fetch($acc); # Returns Bio::Seq object
  };

  $self->throw("Could not fetch sequence for [$acc]") unless defined $seq;
  return $seq; 
}

1;

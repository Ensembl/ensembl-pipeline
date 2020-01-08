=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::SeqFetcher::Getz -

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Getz->new(
							    '-executable' => $exe,
							    '-library'        => $lib,
							   );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using getz.

=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Getz;

use warnings ;
use strict;
use Bio::DB::RandomAccessI;
use Bio::Seq;
use Bio::SeqIO;

use vars qw(@ISA);

@ISA = qw(Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $lib) = $self->_rearrange([
				       'EXECUTABLE',
				       'LIBRARY'], @args);

  if (!defined $exe) {
    $exe = 'getz';
  }
  $self->executable($exe);


  if (defined $lib) {
    $self->library($lib);
  }
  return $self; # success - we hope!
}

=head2 executable

  Title   : executable
  Usage   : $self->executable('/path/to/executable');
  Function: Get/set for the path to the executable being used by the module. If not set, the executable is looked for in $PATH.
  Returns : string
  Args    : string

=cut

sub executable {
  my ($self, $exe) = @_;
  if ($exe)
    {
      $self->{'_exe'} = $exe;
    }
  return $self->{'_exe'};
}

=head2 library

  Title   : library
  Usage   : $self->library('embl');
  Function: Get/set for a library/libraries to search in
  Returns : string
  Args    : string

=cut

sub library {
  my ($self, $lib) = @_;
  if ($lib) {
      $self->{'_lib'} = $lib;
    }
  return $self->{'_lib'};
}


=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_Seq_by_acc($accession);
  Function: Does the sequence retrieval via getz
  Returns : Bio::Seq
  Args    :

=cut

sub get_Seq_by_acc {
  my ($self, $acc) = @_;
  my $libs = $self->library;

  if (!defined($acc)) {
    $self->throw("No id input");
  }

  if (!defined($libs)) {
    $self->throw("No search libs specified");
  }

  my $seqstr;
  my $seq;
  my $getz     = $self->executable;

  open(IN, "$getz  -d -sf fasta '[libs={$libs}-AccNumber:$acc]' |")
    or $self->throw("Error running getz for id [$acc]: $getz");

  my $format = 'fasta';

  my $fh = Bio::SeqIO->new(-fh   => \*IN, "-format"=>$format);

  $seq = $fh->next_seq();
  close IN;

  $self->throw("Could not getz sequence for [$acc]\n") unless defined $seq;
  $seq->display_id($acc);
  $seq->accession_number($acc);

  return $seq;
}
=head2 get_Seq_by_id

  Title   : get_Seq_by_id
  Usage   : $self->get_Seq_by_id($id);
  Function: Does the sequence retrieval via getz
  Returns : Bio::Seq
  Args    :

=cut

sub get_Seq_by_id {
  my ($self, $id) = @_;
  my $libs = $self->library;

  if (!defined($id)) {
    $self->throw("No id input");
  }

  if (!defined($libs)) {
    $self->throw("No search libs specified");
  }

  my $seqstr;
  my $seq;
  my $getz     = $self->executable;

  open(IN, "$getz  -d -sf fasta '[libs={$libs}-ID:$id]' |")
    or $self->throw("Error running getz for id [$id]: $getz");

  my $format = 'fasta';

  my $fh = Bio::SeqIO->new(-fh   => \*IN, "-format"=>$format);

  $seq = $fh->next_seq();
  close IN;

  $self->throw("Could not getz sequence for [$id]\n") unless defined $seq;
  $seq->display_id($id);
  $seq->accession_number($id);

  return $seq;
}

=head2 batch_fetch

  Title   : batch_fetch
  Usage   : $self->batch_fetch(@accs);
  Function: Retrieves batches of sequences
  Returns : listref of Bio::Seq
  Args    :

=cut

# This method is not yet couth.  getz can retrieve multiple sequences at
# once - this method should be reimplemented to do things this way, rather
# than this five-for-the-price-of-five approach.

sub batch_fetch {
  my ($self, @accs) = @_;

  my @seqs;

  while (@accs){

    push @seqs, $self->get_Seq_by_acc(shift @accs);

  }

  return \@seqs
}


1;

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Pipeline::SeqFetcher::Efetch -

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Efetch->new(
							      '-executable' => $exe,
							      '-lib'        => $lib,
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using efetch.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Efetch;

use warnings ;
use strict;
use Bio::DB::RandomAccessI;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $lib) = $self->_rearrange([
				       'EXECUTABLE',
				       'LIBRARY'], @args);

  if (!defined $exe) {
    $exe = 'efetch';
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
  Usage   : $self->library('sw');
  Function: Get/set for a library to search in - eg efetch wil not search swall without the prefix sw:
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
  Function: Does the sequence retrieval via efetch
  Returns : Bio::Seq
  Args    :

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;

  if (!defined($acc)) {
    $self->throw("No accession input");
  }

  my $lib = $self->library;
  if (defined $lib && $lib ne ''){
    $acc = $lib . ":" . $acc;
  }

  my $seqstr;
  my $seq;
  my $efetch = $self->executable;

  open(IN,"$efetch -q $acc |") or $self->throw("Error running efetch for acc [$acc]: $efetch");
  $seqstr = <IN>;
  close IN;

  if(defined $seqstr && $seqstr ne "no match") {
    chomp($seqstr);
    $seq = new Bio::Seq('-seq'               => $seqstr,
			'-accession_number'  => $acc,
			'-display_id'        => $acc);
  }

  $self->throw("Could not efetch sequence for [$acc]\n") unless defined $seq;
  return $seq;

}

1;

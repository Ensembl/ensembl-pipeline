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

  Bio::EnsEMBL::Pipeline::SeqFetcher::Efetch

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Efetch->new(
							      '-executable' => $exe,
							      '-lib'        => $lib,
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using efetch.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Efetch;

use strict;
use Bio::EnsEMBL::Root;
use Bio::DB::RandomAccessI;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root Bio::DB::RandomAccessI);

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

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
package Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

use strict;
use Bio::Root::RootI;
use Bio::DB::RandomAccessI;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my($exe) = $self->_rearrange(['EXECUTABLE'], @args);

  if (!defined $exe) {
    $exe = 'pfetch';
  }
  $self->executable($exe);
  
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


=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_eq_by_acc($accession);
  Function: Does the sequence retrieval via pfetch
  Returns : Bio::Seq
  Args    : 

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;
  
  if (!defined($acc)) {
    $self->throw("No accession input");
  }  
  
  my $seqstr;
  my $seq;
  my $pfetch = $self->executable;
  
  open(IN,"$pfetch -q $acc |") or $self->throw("Error running pfetch for accession [$acc]: $pfetch");
  $seqstr = <IN>;
  close IN;
  
  chomp($seqstr);
  if(defined $seqstr && $seqstr ne "no match") {
    $seq = new Bio::Seq('-seq'               => $seqstr,
			'-accession_number'  => $acc,
			'-display_id'        => $acc);
  }

  $self->throw("Could not pfetch sequence for [$acc]\n") unless defined $seq;

  return $seq;
}

1;

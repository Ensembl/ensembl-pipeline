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

Bio::EnsEMBL::Pipeline::SeqFetcher::xdget

=head1 SYNOPSIS

  my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::xdget->new(
      -executable => $exe
  );
  my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

Object to retrieve sequences using xdget (Wash U).
Database must be formatted with xdformat. Sequence type (protein
or nucleotide) is guessed, based on file extensions of the database
files. Returns sequence as Bio::Seq.
Additional options for xdget can be specifed though no checking
is performed for compatibility.

=head1 CONTACT

B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::SeqFetcher::xdget;

use strict;
use Bio::EnsEMBL::Root;
use Bio::DB::RandomAccessI;
use Bio::SeqIO;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($exe, $options, $db) = $self->_rearrange([qw(
    EXECUTABLE OPTIONS DB
  )], @args);

  $exe ||= 'xdget';
  $self->executable($exe);
  $self->options($options) if defined $options;
  $self->db($db) if defined $db;
  
  return $self;
}

=head2 executable

  Title   : executable
  Usage   : $self->executable('/path/to/executable');
  Function: Get/set for the executable to be used by the module.
  Returns : string
  Args    : string

=cut

sub executable {
  my ($self, $exe) = @_;
  if ($exe) {
      $self->{'_exe'} = $exe;
  }
  return $self->{'_exe'};  
}

=head2 db

  Title   : db
  Usage   : $self->db('/path/to/db');
  Function: Get/set for the db to be used by the module.
  Returns : string
  Args    : string

=cut

sub db {
  my ($self, $db) = @_;
  if ($db) {
      $self->{'_db'} = $db;
      if (glob "$db.xp?") {
	  $self->_moltype('p');
      }
      elsif (glob "$db.xn?") {
	  $self->_moltype('n');
      }
      else {
	  $self->throw("XDF database appears to be missing files");
      }
  }
  return $self->{'_db'};  
}

=head2 options

  Title   : options
  Usage   : $self->options('-r');
            Returns reverse complement of nucleotide sequence
  Function: Get/set for options to xdget
  Returns : string
  Args    : string

=cut

sub options {

  my ($self, $options) = @_;
  if ($options) {
      $self->{'_options'} = $options;
  }
  return $self->{'_options'};  

}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_eq_by_acc($accession);
  Function: retrieves sequence via xdget
  Returns : Bio::Seq
  Args    : Sequence identifier string

=cut

sub get_Seq_by_acc {
  my ($self, $acc) = @_;
  
  $self->throw("No accession input") unless $acc;
  $self->throw("No database defined") unless $self->db;
  
  my $seqstr;
  my $seq;
  my $xdget   = $self->executable;
  my $options = $self->options;
  my $db      = $self->db;

  # maybe should have some checking here to see if -n/-p have
  # already been specified in options
  if ($self->_moltype eq 'n') {
    $options .= " -n";
  }
  else {
    $options .= " -p";
  }

  my $command = "$xdget $options $db $acc";

  local *FH;
  open FH, "$command |" or $self->throw("Error retrieving $acc from $db with $xdget");

  my $seq = Bio::SeqIO->new(
    -format => 'fasta',
    -fh     => \*FH
  )->next_seq;

  close FH or $self->throw("Error retrieving $acc from $db with $xdget");

  return $seq;
}

sub _moltype {
    my ($self, $type) = @_;

    if ($type) {
	$self->{'_type'} = $type;
    }
    return $self->{'_type'};
}


1;

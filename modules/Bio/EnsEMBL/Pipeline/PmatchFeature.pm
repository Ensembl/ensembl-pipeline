=head1 NAME

Bio::EnsEMBL::Analysis::PmatchFeature

=head1 SYNOPSIS

    

=head1 DESCRIPTION


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX


=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::PmatchFeature;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Exon;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($protein_id,$cdna_id,$chr_name,$start,$end,$coverage) = $self->_rearrange([qw(PROTEIN_ID
										 CDNA_ID
										 CHR_NAME
										 START
										 END
										 COVERAGE
										)],@args);
  
  
  $self->protein_id($protein_id);
  $self->chr_name($chr_name);
  $self->start($start);
  $self->end($end);
  $self->cdna_id($cdna_id);
  $self->coverage($coverage);

  return $self;
}

sub  cdna_id {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_cdna_id} = $arg;
  }

  return $self->{_cdna_id};
}

sub coverage {
  my ($self,$arg) = @_;

  if (defined($arg)) {

    if ($arg < 0 || $arg > 100) {
      $self->throw("Coverage must be betwee 0 and 100.  Trying to set it to $arg");
    }
    
    $self->{_coverage} = $arg;
  }

  return $self->{_coverage};
}

sub protein_id {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_protein_id} = $arg;
  }

  return $self->{_protein_id};
}

sub chr_name{
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_chr_name} = $arg;
  }
  return $self->{_chr_name};
}

sub start {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_start} = $arg;
  }
  return $self->{_start};
}

sub end {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_end} = $arg;
  }
  return $self->{_end};
}

1;


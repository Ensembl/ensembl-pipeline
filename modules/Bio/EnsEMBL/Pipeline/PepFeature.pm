=head1 NAME

Bio::EnsEMBL::Analysis::PepFeature

=head1 SYNOPSIS

    
=head1 DESCRIPTION


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX


=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::PepFeature;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::FeaturePair;

@ISA = qw(Bio::EnsEMBL::FeaturePair);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($qgenenum,$qtrans,$qpos,$qres,$qchr,$hgenenum,$htrans,$hpos,$hres,$hchr) = 

    $self->_rearrange([qw(QGENENUM
			  QTRANS
			  QPOS
			  QRES
			  QCHR
			  HGENENUM
			  HTRANS
			  HPOS
			  HRES
			  HCHR
			 )],@args);


  $self->qgenenum($qgenenum);
  $self->qtrans($qtrans);
  $self->qpos($qpos);
  $self->qres($qres);
  $self->qchr($qchr);

  $self->hgenenum($hgenenum);
  $self->htrans($htrans);
  $self->hpos($hpos);
  $self->hres($hres);
  $self->hchr($hchr);
  
  return $self;
}


sub qgenenum {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_qgenenum} = $arg;
  }
  return $self->{_qgenenum};
}

sub hgenenum{
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_hgenenum} = $arg;
  }
  return $self->{_hgenenum};
}

sub qtrans {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_qtrans} = $arg;
  }
  return $self->{_qtrans};
}

sub htrans {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_htrans} = $arg;
  }
  return $self->{_htrans};
}

sub qpos {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_qpos} = $arg;
  }
  return $self->{_qpos};
}

sub hpos {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_hpos} = $arg;
  }
  return $self->{_hpos};
}

sub qres {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_qres} = $arg;
  }
  return $self->{_qres};
}

sub hres {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_hres} = $arg;
  }
  return $self->{_hres};
}


sub qchr {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_qchr} = $arg;
  }
  return $self->{_qchr};
}

sub hchr{
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_hchr} = $arg;
  }
  return $self->{_hchr};
}

sub print {
  my ($self) = @_;

  print $self->gffstring . "\t" . $self->p_value . "\t" . $self->percent_id . "\t" . $self->qtrans . "\t" . $self->qpos . "\t" .$self->qres . "\t" . $self->qgenenum . "\t" . $self->qchr . "\t" . $self->htrans . "\t" . $self->hpos . "\t" .$self->hres . "\t" . $self->hgenenum . "\t" . $self->hchr . "\n";
}

1;


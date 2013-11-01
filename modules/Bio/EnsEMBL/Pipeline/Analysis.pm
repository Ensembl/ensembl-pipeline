=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::Analysis - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Pipeline::Analysis;

use warnings ;
use  vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( deprecate warning throw );

@ISA = qw(Bio::EnsEMBL::Analysis);





sub new {
  my($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);
  
  my ($type) = rearrange([qw(INPUT_ID_TYPE)], @args);

  $self->input_id_type($type);

  return $self;

}


sub input_id_type{
  my ($self, $type) = @_;
  if($type){
    $self->{'type'} = $type;
  }
  
  return $self->{'type'};
}

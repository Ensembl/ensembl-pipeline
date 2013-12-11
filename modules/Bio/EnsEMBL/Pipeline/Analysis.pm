=head1 LICENSE
# Copyright [1999-2013] Genome Research Ltd. and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


# Mar 6, 2006 5:19:41 PM
#
# Created by Mustapha Larbaoui <ml6@sanger.ac.uk>

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor

=head1 SYNOPSIS

my $dbobj = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(
	-host   => $dbhost,
	-dbname => $dbname,
	-user   => $dbuser,
	-pass   => $dbpass,
	-port   => $dbport,
);

=head1 DESCRIPTION

Interface for the connection to the analysis database

=head1 FEEDBACK

=head1 AUTHOR - Mustapha Larbaoui

Mustapha Larbaoui E<lt>ml6@sanger.ac.ukE<gt>

=head1 CONTACT

Post general queries to B<anacode@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;

use warnings ;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

@ISA = qw(Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor);

sub get_available_adaptors {
  my ($self) = @_;

  my $pairs = $self->SUPER::get_available_adaptors();

  $pairs->{'Job'}                = 'Bio::EnsEMBL::Pipeline::DBSQL::Finished::JobAdaptor';
  $pairs->{'StateInfoContainer'} = 'Bio::EnsEMBL::Pipeline::DBSQL::Finished::StateInfoContainer';
  $pairs->{'HitDescription'} = 'Bio::EnsEMBL::Pipeline::DBSQL::Finished::HitDescriptionAdaptor';
  $pairs->{'DnaAlignFeature'} = 'Bio::EnsEMBL::Pipeline::DBSQL::Finished::DnaAlignFeatureAdaptor';
  $pairs->{'ProteinAlignFeature'} = 'Bio::EnsEMBL::Pipeline::DBSQL::Finished::ProteinAlignFeatureAdaptor';

  return $pairs; 
}


1;

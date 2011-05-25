=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

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

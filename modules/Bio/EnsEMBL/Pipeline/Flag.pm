=head1 LICENSE

# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::Flag - 

=head1 SYNOPSIS

  my $flag = Bio::EnsEMBL::Pipeline::Flag->new(
					       '-type'         => 'gene',
					       '-ensembl_id'   => 12345,
					       '-goalAnalysis' => $analysis_object,
					       );

=head1 DESCRIPTION

Object to allow flagging of sequences for further analysis.
The flag consists of the dbID of the object to be flagged, the name of the
table in which it is stored and the analysis that you would like to
pass the object to.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Flag;

use warnings ;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use strict;


@ISA = qw( Bio::EnsEMBL::Root );

=head2 new

  Title   : new
  Usage   : ...Flag->new($analysis);
  Function: Constructor for Flag object
  Returns : Bio::EnsEMBL::Pipeline::Flag
  Args    : A Bio::EnsEMBL::Analysis object.
            adaptor and dbid only used from the adaptor.

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_goal'} = {}; # Container for analysis object
  $self->{'_dbID'} = 0;  # Container for scalar dbID of the flag
  $self->{'_id'} = 0;    # Container for scalar dbID of object to be flagged
  $self->{'_type'} = ""; # Name of the table the object is stored in ie: gene 
  $self->{'_adaptor'} = {}; # Container for FlagAdaptor object

  my ( $goal,$type, $ensembl_id, $adaptor, $dbID ) =
    Bio::EnsEMBL::Utils::Argument->rearrange( [ qw ( GOALANALYSIS
			      TYPE
			      ENSEMBL_ID
			      ADAPTOR
			      DBID
			     ) ], @args );
  $self->throw( "Wrong parameter" ) unless $goal->isa( "Bio::EnsEMBL::Analysis");
  $self->type( $type );
  $self->goalAnalysis( $goal );
  $self->adaptor( $adaptor );
  $self->dbID( $dbID );
  $self->ensembl_id( $ensembl_id  );
  return $self;
}

=head2 goalAnalysis

  Title   : goalAnalysis
  Usage   : $self->goalAnalysis($analysis);
  Function: Get/set method for the goal analysis object of this flag.
  Returns : Bio::EnsEMBL::Analysis
  Args    : Bio::EnsEMBL::Analysis

=cut

sub goalAnalysis {
  my ($self,$arg) = @_;

  ( defined $arg ) &&
    ( $self->{'_goal'} = $arg );
  return $self->{'_goal'};
}

=head2 dbID

  Title   : dbID
  Usage   : $self->dbID($dbID);
  Function: Get/set method for the dbID of this flag.
  Returns : Scalar
  Args    : Scalar

=cut

sub dbID {
  my ( $self, $dbID ) = @_;
  ( defined $dbID ) &&
    ( $self->{'_dbID'} = $dbID );
  return $self->{'_dbID'};
}

=head2 ensembl_id

  Title   : ensemblID
  Usage   : $self->ensembl_id($ensemblID);
  Function: Get/set method for the dbID of the object to be flagged
  Returns : Scalar
  Args    : Scalar

=cut

sub ensembl_id {
  my ( $self, $id) = @_;
  ( defined $id ) &&
    ( $self->{'_id'} = $id );
  return $self->{'_id'};
}

=head2 type

  Title   : type
  Usage   : $self->type('table_name');
  Function: Get/set method for the name of the table that the flagged object sits in
  Returns : Scalar
  Args    : Scalar

=cut

sub type {
  my ( $self, $type) = @_;
  ( defined $type ) &&
    ( $self->{'_type'} = $type );
  $self->{'_type'};
}

=head2 adaptor

  Title   : adaptor
  Usage   : $self->adaptor($adaptor);
  Function: Get/set method for the flag adaptor set is private
  Returns : Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor
  Args    : Bio::EnsEMBL::Pipeline::DBSQL::FlagAdaptor

=cut

sub adaptor {
  my ( $self, $adaptor ) = @_;
  ( defined $adaptor ) &&
    ( $self->{'_adaptor'} = $adaptor );
  return $self->{'_adaptor'};
}

1;




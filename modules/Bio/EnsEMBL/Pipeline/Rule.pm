# Perl module for Bio::EnsEMBL::Pipeline::Rule
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 10.09.2000
# Last modified : 10.09.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Rule 

=head1 SYNOPSIS


=head1 DESCRIPTION
  

=head1 CONTACT

    Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
    Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Rule;
use vars qw(@ISA);
use strict;


@ISA = qw( Bio::Root::RootI );

=head2 Constructor

  Title   : new
  Usage   : ...Rule->new($analysis);
  Function: Constructor for Rule object
  Returns : Bio::EnsEMBL::Pipeline::Rule
  Args    : A Pipeline::Analysis object. Conditions are added later,
            adaptor and dbid only used from the adaptor.

=cut


sub new {
  my $class = shift;
  my $self = bless {},$class;

  my ( $goal, $adaptor, $dbID ) = 
    $self->_rearrange( [ qw ( GOAL
			      ADAPTOR
			      DBID
			     ) ], @_ );
  $self->throw( "Wrong parameter" ) unless
    $goal->isa( "Bio::EnsEMBL::Pipeline::Analysis" );
  $self->dbID( $dbID );
  $self->goalAnalysis( $goal );
  $self->adaptor( $adaptor );
				
  return $self;
}

=head2 add_condition

  Title   : add_condition
  Usage   : $self->add_conditon($logic_name);
  Function: Add method for conditions for the rule.
  Returns : nothing
  Args    : a string describing a testable condition on an object
            For now a logic_name of an analysis.

=cut


sub add_condition {
  my $self = shift;
  my $condition = shift;

  push( @{$self->{_conditions}}, $condition );
}

=head2 list_conditions

  Title   : list_conditions
  Usage   : $self->list_conditions();
  Function: Give a list of all conditions you have to fulfill to make this
            Rule evaluate to true.
  Returns : a list of strings which are probably logic_names
  Args    : -

=cut



sub list_conditions {
  my $self = shift;

  return @{$self->{_conditions}};
}

=head2 goalAnalysis

  Title   : goalAnalysis
  Usage   : $self->goalAnalysis($anal);
  Function: Get/set method for the goal analysis object of this rule.
  Returns : Bio::EnsEMBL::Pipeline::Analysis
  Args    : bio::EnsEMBL::Pipeline::Analysis

=cut

sub goalAnalysis {
  my ($self,$arg) = @_;
  
  ( defined $arg ) &&
    ( $self->{_goal} = $arg );
  $self->{_goal};
}

sub dbID {
  my ( $self, $dbID ) = @_;
  ( defined $dbID ) &&
    ( $self->{_dbID} = $dbID );
  $self->{_dbID};
}

sub adaptor {
  my ( $self, $adaptor ) = @_;
  ( defined $adaptor ) &&
    ( $self->{_adaptor} = $adaptor );
  $self->{_adaptor};
}


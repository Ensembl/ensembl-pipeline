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
use Bio::Root::RootI;
use strict;


@ISA = qw( Bio::Root::RootI );

=head2 Constructor

  Title   : new
  Usage   : ...Rule->new($analysis);
  Function: Constructor for Rule object
  Returns : Bio::EnsEMBL::Pipeline::Rule
  Args    : A Bio::EnsEMBL::Analysis object. Conditions are added later,
            adaptor and dbid only used from the adaptor.

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $goal, $adaptor, $dbID ) =
    $self->_rearrange( [ qw ( GOAL
			      ADAPTOR
			      DBID
			     ) ], @args );
  $self->throw( "Wrong parameter" ) unless
    $goal->isa( "Bio::EnsEMBL::Analysis" );
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

  push( @{$self->{'_conditions'}}, $condition );
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

  my @conditions = @{$self->{'_conditions'}};
  if (! scalar (@conditions) ) {
      $self->throw("No conditions found for this Rule");
  }
  return @conditions;
}

=head2 goalAnalysis

  Title   : goalAnalysis
  Usage   : $self->goalAnalysis($anal);
  Function: Get/set method for the goal analysis object of this rule.
  Returns : Bio::EnsEMBL::Analysis
  Args    : Bio::EnsEMBL::Analysis

=cut

sub goalAnalysis {
  my ($self,$arg) = @_;

  ( defined $arg ) &&
    ( $self->{'_goal'} = $arg );
  $self->{'_goal'};
}


# return 0 if nothing can be done or $goalAnalysis,
# if it should be done.

sub check_for_analysis {
  my $self = shift;
  my @analist = @_;
  my %anaHash;

  # reimplement with proper identity check!
  my $goal = $self->goalAnalysis->dbID;

  # print STDERR "My goal is " . $goal . "\n";

  for my $analysis ( @analist ) {
      # print STDERR " Analysis " . $analysis->logic_name . " " . $analysis->dbID . "\n";
    $anaHash{$analysis->logic_name} = $analysis;
    if( $goal == $analysis->dbID ) {
      # already done
      return 0;
    }
  }

  for my $cond ( $self->list_conditions ) {
    if( ! defined $anaHash{$cond} ) {
      return 0;
    }
  }
  return $self->goalAnalysis;
}





sub dbID {
  my ( $self, $dbID ) = @_;
  ( defined $dbID ) &&
    ( $self->{'_dbID'} = $dbID );
  $self->{'_dbID'};
}

sub adaptor {
  my ( $self, $adaptor ) = @_;
  ( defined $adaptor ) &&
    ( $self->{'_adaptor'} = $adaptor );
  $self->{'_adaptor'};
}

1;

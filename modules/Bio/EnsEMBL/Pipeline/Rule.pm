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

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Rule;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw();

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
  my $self = bless {},$class;

  my ( $goal, $adaptor, $dbID ) =
    rearrange( [ qw ( GOALANALYSIS
                      ADAPTOR
                      DBID
                    ) ], @args ); 
  if ( $goal ) { 
    throw( "Wrong parameter" ) unless
     $goal->isa( "Bio::EnsEMBL::Analysis" ); 
  }
  $self->dbID( $dbID );
  $self->goalAnalysis( $goal );
  $self->adaptor( $adaptor );
  $self->{'_condition_type_cache'} = {};
				
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

  my @conditions;
  if ($self->{'_conditions'}) {
    @conditions = @{$self->{'_conditions'}};
  }
  if (! scalar (@conditions) ) {
      throw("No conditions found for this Rule");
  }
  return \@conditions;
}


sub has_condition_of_input_id_type {
  my $self = shift;
  my $id_type = shift;
  my $ana_adaptor = $self->goalAnalysis->adaptor;

  if (!scalar (@{$self->{'_conditions'}})) {
      throw("No conditions found for this Rule");
  }
  
  return 1 if (exists($self->{'_condition_type_cache'}{$id_type}));

  foreach my $cond (@{$self->{'_conditions'}}) {
      my $cond_anal = $ana_adaptor->fetch_by_logic_name($cond);
      if ($cond_anal->input_id_type eq $id_type) {
          #print " Condition of " . $cond_anal->logic_name . " is of type $id_type\n";
          $self->{'_condition_type_cache'}{$id_type} = 1;
          return 1;
      }
  }
  return 0;
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

=head2 check_for_analysis

 -args: [analysis list], 'input id type', {completed accumulators}, verbose
 -returns: Either bits for status if nothing can be done;
           1 - Failed Input_Id_Type Check.
           2 - Failed Already Complete Check [so is complete].
           4 - Failed Condition Check.
           Or;
           $goalAnalysis if it should be done
=cut


sub check_for_analysis {
  my $self = shift;
  my ($analist, $input_id_type, $completed_accumulator_href, $verbose) = @_;
  my %anaHash;
  my $return = 0;

  # reimplement with proper identity check!
  my $goal_anal    = $self->goalAnalysis;
  my $goal         = $goal_anal->dbID;
  my $goal_id_type = $goal_anal->input_id_type;

  print "\nHave goal type ".$goal_id_type." and input id type ".$input_id_type."\n" if($verbose);

#This id isn't of the right type so doesn't satify goal
  if ($goal_id_type ne 'ACCUMULATOR' &&
      $goal_id_type ne $input_id_type) {
    print "In check_for_analysis failed input_id_type check as goal input_id type ".
      "isn't the same as the input_id type\n" if ($verbose);
    $return += 1;
  }


  print " My goal is " . $goal_anal->logic_name . "\n" if($verbose);
  print " Completed anals on this id:\n" if ($verbose);

  for my $analysis ( @$analist ) {

    print "  Analysis " . $analysis->logic_name . " " . $analysis->dbID . "\n" if($verbose);
    $anaHash{$analysis->logic_name} = $analysis;

    if ($goal == $analysis->dbID) {
      # already done
      print $goal_anal->logic_name." already done\n" if ($verbose);
      $return += 2;
    }
  }

#the completed_accumulator_href contains input_id_type ACCUMULATOR anals that have completed
  print " Checking conditions:\n" if ($verbose);
  for my $cond ( @{$self->{'_conditions'}} ) {
   # remove any trailing whitespace because the matches won't happen if there is whitespace in, eg. your rule_conditions table.
   $cond =~ s/\s+//g; 
    if ( ! exists $anaHash{$cond} && ! exists $completed_accumulator_href->{$cond}) {
      print "   failed condition check for $cond\n" if ($verbose);
      $return += 4;
    }
  }
  if ($return < 4) {
    print "  All conditions satisfied\n" if ($verbose);
  }

  return $return if $return;
  return $goal_anal;
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




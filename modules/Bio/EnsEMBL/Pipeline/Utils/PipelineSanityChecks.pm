=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks -

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;

use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Config::General qw( BIN_DIR LIB_DIR DATA_DIR );
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;


sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = bless({}, $class);

  $self->{'db'} = undef;

  my ($db)= rearrange([qw(DB)], @_);

  $self->db($db) if($db);

  throw("you need to pass at least a DBAdaptor to an PipelineSanityChecks") unless($self->db);

  return $self;
}



sub db{
  my $self = shift;

  if(@_){
    $self->{'db'} = shift;
  }

  return $self->{'db'};
}

sub db_sanity_check{
  my ($self) = @_;

  my ($query, $msg);
  my $warn = 1;
  #check all rules in the rule_goal table have existing analyses
  $query = qq{SELECT COUNT(DISTINCT g.rule_id)
              FROM rule_goal g
              LEFT JOIN analysis a ON g.goal = a.analysis_id
              WHERE a.analysis_id IS NULL};
  $msg = "Some of your goals in the rule_goal table don't seem".
    " to have entries in the analysis table";
  $self->execute_sanity_check($query, $msg);
  #check all rules in the rule_condition table have existing analyses
  $query = qq{SELECT COUNT(DISTINCT c.rule_id)
              FROM rule_conditions c
              LEFT JOIN analysis a ON c.rule_condition = a.logic_name
              WHERE a.logic_name IS NULL};
  $msg = "Some of your conditions in the rule_condition table don't" .
    " seem to have entries in the analysis table";
  $self->execute_sanity_check($query, $msg);
  #check all the analyses have types
  $query = qq{SELECT COUNT(DISTINCT(a.analysis_id))
              FROM analysis a
              LEFT JOIN input_id_type_analysis t
              ON a.analysis_id = t.analysis_id
              WHERE t.analysis_id IS NULL};
  $msg = "Some of your analyses don't have entries in the".
    " input_id_type_analysis table";
  $self->execute_sanity_check($query, $msg, $warn);
  #check that all types which aren't accumulators have entries in
  #input__id_analysis table
  $query = qq{SELECT DISTINCT(t.input_id_type)
              FROM input_id_analysis i
              LEFT JOIN input_id_type_analysis t
              ON i.input_id_type = t.input_id_type
              WHERE t.input_id_type IS NULL
                && t.input_id_type != 'ACCUMULATOR'};
  $msg = "Some of your types don't have entries in the".
    " input_id_type_analysis table";
  $self->execute_sanity_check($query, $msg);

  $query = qq{SELECT count(input_id)
              FROM input_id_analysis
              WHERE input_id_type = ''};
  $msg = "Some of your input_ids don't have a type in the input_id_analysis ".
    "table";
  $self->execute_sanity_check($query, $msg);
}

sub execute_sanity_check{
    my ($self, $query, $msg, $warn) = @_;
    my $db = $self->db;
    my $sth = $db->prepare($query);
    $sth->execute();
    if($warn){
      warning($msg) if $sth->fetchrow();
    }else{
      throw($msg) if $sth->fetchrow();
    }
}


sub accumulator_sanity_check{
  my ($self, $rules, $accumulators) = @_;

  my $sic = $self->db->get_StateInfoContainer;
  my $aa = $self->db->get_AnalysisAdaptor;

 RULE:foreach my $rule(@$rules){
    if($rule->goalAnalysis->input_id_type eq 'ACCUMULATOR'){
      #print STDERR "dealing with rule ".$rule->goalAnalysis->logic_name."\n";
      my @conditions = @{$rule->list_conditions};
      my %input_id_type;
      foreach my $c(@conditions){
        #print STDERR "have condition ".$c."\n";
        my $analysis = $aa->fetch_by_logic_name($c);
        if(!$input_id_type{$analysis->input_id_type}){
          $input_id_type{$analysis->input_id_type} = [];
        }
        push(@{$input_id_type{$analysis->input_id_type}}, $c);
      }
      TYPE:foreach my $type(keys(%input_id_type)){
          #print STDERR "have type ".$type."\n";
        my @ids = @{$sic->list_input_ids_by_type($type)};
          #print STDERR "have ".@ids." ids\n";
        if(!@ids){
          my $logic_names = join(",", @{$input_id_type{$type}});
          print STDERR "can't run with accumulators on as ".
            $rule->goalAnalysis->logic_name." depends on $logic_names with type ".
              $type." which has no entries in the input_id_".
                "analysis table\n";
          $accumulators = 0;
        }else{
          next TYPE;
        }
      }
    }
else{
      next RULE;
    }

  }
    return $accumulators;
}



sub rule_type_sanity{
  my ($self, $rules, $verbose) = @_;

  my $aa = $self->db->get_AnalysisAdaptor;
 RULE:foreach my $rule(@$rules){
    my $type = $rule->goalAnalysis->input_id_type;
    if($type eq 'ACCUMULATOR'){
      next RULE;
    }
  CONDITION:foreach my $name(@{$rule->list_conditions}){
      my $condition = $aa->fetch_by_logic_name($name);
      if(!$condition){
        my $msg = "Can't depend on an analysis which doesn't exist $name";
        throw($msg);
      }
      if($condition->input_id_type eq 'ACCUMULATOR'){
        print STDERR "Skipping ".$name." is an accumulator\n" if($verbose);
        next CONDITION;
      }
      if($condition->input_id_type ne $type){
        my $msg = $rule->goalAnalysis->logic_name."'s type ".$type.
                     " doesn't match condition ".$condition->logic_name.
                     "'s type ".$condition->input_id_type;
        throw($msg);
      }
    }
  }
}


sub config_sanity_check {
  my ($self) = @_;

  my $ok = 1;
  unless ($QUEUE_MANAGER) {
    print "Need to specify QUEUE_MANAGER in Config/BatchQueue.pm\n";
    $ok = 0;
  }
  unless ($LIB_DIR) {
    print "Need to specify LIB_DIR in Config/General.pm\n";
    $ok = 0;
  }
  unless ($DATA_DIR) {
    print "Need to specify DATA_DIR in Config/General.pm\n";
    $ok = 0;
  }
  unless ($BIN_DIR) {
    print "Need to specify BIN_DIR in Config/General.pm\n";
    $ok = 0;
  }

  if(!$ok){
    throw("The pipeline config isn't sane");
  }
}


1;

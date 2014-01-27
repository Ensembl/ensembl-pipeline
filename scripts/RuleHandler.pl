#!/usr/bin/env perl


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


=head1 NAME

RuleHandler.pl - handles insertion, deletion and listing of rules in a database

=head1 SYNOPSIS

RuleHandler.pl -optio

=head1 DESCRIPTION

  This script allows the user to list all analysisprocesses and rules
  in a database, insert rules and delete rules. Checks for circular
  dependencies and fulfilment of conditions are made before rules are
  inserted. Rules cannot be deleted if other rules depend upon them.

=head1 OPTIONS

    -dbhost     host name for database (gets put as host= in locator)

    -dbport     For RDBs, what port to connect to (port= in locator)

    -dbname     For RDBs, what name to connect to (dbname= in locator)

    -dbuser     For RDBs, what username to connect as (dbuser= in locator)

    -dbpass     For RDBs, what password to use (dbpass= in locator)

    -help       Displays script documentation with PERLDOC

    -analyses   Lists all available analysisprocesses

    -rules      Lists all available rules together with conditions and goals

    -delete     Deletes a rule based on value of id. Requires
                additional qualifier, -ruleId.

    -ruleId     ID of rule to be deleted (required by -delete)

    -insert     Inserts a rule that has the given goal and conditions
                (requires additional qualifiers -goal and -condition

     -goal      Goal analysis for the rule to be inserted (goal and
                conditions must be logic_names)

     -condition Conditions required for rule to be inserted. A rule
                can have more than one condition; each must be prefaced by
                -condition, eg -condition condition1 -condition condition2

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::Rule;

my $goal;
my @conditions;
my $insert;
my $delete;
my $ruleId;
my $show_analyses;
my $show_rules;
my $host     = 'localhost';
my $dbname   = 'pipeline';
my $dbuser   = 'ensadmin';
my $pass     = '';
my $port     = '3306';
my $help; 
my $add_rule ; 
my %rule_goals; # hash to link up goals with rules for speed later

GetOptions(
  'host|dbhost|h:s'     => \$host,
  'dbname|db|D:s'     => \$dbname,
  'user|dbuser|u:s'     => \$dbuser,
  'pass|dbpass|p:s'     => \$pass,
  'port|dbport|P:s'     => \$port,
  'goal=s'       => \$goal,
  'condition=s@' => \@conditions,
  'insert'       => \$insert,
  'delete'       => \$delete,
  'analyses'     => \$show_analyses,
  'rules'        => \$show_rules,
  'ruleId=i'     => \$ruleId,
  'help'       => \$help,
  'auto!' =>
    \$add_rule, # automatically create the condtion out of supplied logic_name
);
if (!($insert||$delete||$show_analyses||$show_rules)) {
  $help = 1;
}
if ($help) {
    exec('perldoc', $0);
}
print STDERR $host." ".$dbname." ".$pass ." ".$dbuser." ".$port."\n";

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new
( -host   => $host,
  -dbname => $dbname,
  -user   => $dbuser,
  -pass   => $pass,
  -port   => $port
  );

my $anaAdaptor     = $db->get_AnalysisAdaptor();
my $ruleAdaptor    = $db->get_RuleAdaptor; 


if ( scalar(@conditions) eq 0 && defined $add_rule ) {
  push @conditions, "Submit_" . $goal;
  # check if submission analysis exists in database
  my $sa = $anaAdaptor->fetch_by_logic_name("Submit_$goal");
  if ( !defined $sa ) {
    my $goal_analysis = $anaAdaptor->fetch_by_logic_name($goal);
    my $analysis      = new Bio::EnsEMBL::Pipeline::Analysis;
    $analysis->logic_name("Submit_$goal");
    $analysis->module("Dummy");
    $analysis->input_id_type( $goal_analysis->input_id_type );
    $anaAdaptor->store($analysis);
  }
}

my $basis          = "SubmitContig";

my @existing_rules = $ruleAdaptor->fetch_all();

$rule_goals{$basis} = 1;

foreach my $rule(@existing_rules) {
  my $goal = $rule->goalAnalysis()->logic_name();
  $rule_goals{$goal} = $rule;
}

if ($insert) { &insert_rule(); }
elsif ($delete) { &delete_rule(); }
elsif ($show_analyses) { &show_analyses(); }
elsif ($show_rules) { &show_rules(); }
else { exec('perldoc', $0); }

=head2 insert_rule

  Title   : insert_rule
  Usage   : insert_rule()
  Function: Inserts a new rule into the user specified database, if dependency conditions are satisfied
  Returns : Nothing
  Args    : None, though the command line options -goal and -condition must have been used

=cut

sub insert_rule {

  $goal or die "Cannot insert a rule without a goal logic_name\n";
  @conditions
    or die "Cannot insert a rule without at least one condition logic_name\n";

  # check the goal is valid
  my $analysis = $anaAdaptor->fetch_by_logic_name($goal);
  if ( !$analysis ) {
    die "The goal condition $goal is not found in the database\n";
  }

  # no need to insert a rule for SubmitContig
  if ( $goal eq $basis ) {
    die "There's no need to insert a rule for $basis\n";
  }

  # check the input conditions are unique and valid
  my %checked;
  foreach my $condition (@conditions) {
    die "$condition is not found in the database\n"
      unless scalar( $anaAdaptor->fetch_by_logic_name($condition) );
    $checked{$condition} = 1;
  }

  # check that the goal is not the same as any of the conditions
  die "Sorry, you can't insert a rule "
    . "if the goal is the same as one of the conditions\n"
    if exists $checked{$goal};

  # make and store the rule
  my $rule = Bio::EnsEMBL::Pipeline::Rule->new( -goalanalysis => $analysis );
  foreach my $condition ( keys %checked ) {
    $rule->add_condition($condition);
  }
  if ( check_dependencies($rule) && check_duplications($rule) ) {
    $ruleAdaptor->store($rule);
  }

} ## end sub insert_rule

=head2 delete_rule

  Title   : delete_rule
  Usage   : delete_rule()
  Function: Deletes a rule from the database as long as no other rules depend on it.
  Returns : Nothing
  Args    : None, though the command line option -ruleId must be used

=cut
sub delete_rule {

  $ruleId or die "Cannot remove a rule without the ruleId\n";
  my $rule = $ruleAdaptor->fetch_by_dbID($ruleId);
  if ( ! defined $rule ) { die "There's no rule in the database with id $ruleId\n";  }

  # check dependencies
  my $goal = $rule->goalAnalysis()->logic_name();
  foreach my $old_rule(@existing_rules) {
    foreach my $cond(@{$old_rule->list_conditions()}) {
      die "Sorry, another rule depends on this one\n" unless $cond ne $goal;
    }
  }

  $ruleAdaptor->remove($rule);
}

=head2 show_analyses

  Title   : show_analyses
  Usage   : show_analyses
  Function: Lists details of the analysis processes in the database
  Returns : Nothing
  Args    : None

=cut
sub show_analyses {

  my @analyses = @{ $anaAdaptor->fetch_all() };
  scalar(@analyses) or die "There are no analyses in the database\n";
  foreach my $analysis (@analyses) {
    print "analysisId:\t", $analysis->dbID(), "\n",
      "created:\t",        $analysis->created(),         "\n",
      "logic_name:\t",     $analysis->logic_name(),      "\n",
      "db:\t\t",           $analysis->db(),              "\n",
      "db_version:\t",     $analysis->db_version(),      "\n",
      "db_file:\t",        $analysis->db_file(),         "\n",
      "program:\t",        $analysis->program(),         "\n",
      "program_version:",  $analysis->program_version(), "\n",
      "program_file:\t",   $analysis->program_file(),    "\n",
      "parameters:\t",     $analysis->parameters(),      "\n",
      "module:\t\t",       $analysis->module(),          "\n",
      "module_version:\t", $analysis->module_version(),  "\n",
      "gff_source:\t",     $analysis->gff_source(),      "\n",
      "gff_feature:\t",    $analysis->gff_feature(),     "\n\n",
      ;
  }
}

=head2 show_rules

  Title   : show_rules
  Usage   : show_rules
  Function: Lists details of the rules already in the database
  Returns : Nothing
  Args    : None

=cut

sub show_rules {

  scalar(@existing_rules) or die "There are no rules in the database\n";
  print "\nConditions must be fulfilled "
      . "before the goal analysis can be completed\n\n";
  my @sorted = sort { $a->dbID <=> $b->dbID } @existing_rules;
  foreach my $rule (@sorted) {
    print "Rule ID:", $rule->dbID(), "\n", "    conditions:\t";
    my @conditions = @{ $rule->list_conditions() };
    foreach my $cond (@conditions) {
      print "$cond (", get_analysis_id($cond), ") ";
    }
    print "\n    goal:\t"
        . $rule->goalAnalysis()->logic_name()
        . " (". $rule->goalAnalysis()->dbID, ")\n\n";
  }
}

sub get_analysis_id {
  my ($logic) = shift;

  my $ana = $anaAdaptor->fetch_by_logic_name($logic);

  return $ana->dbID;
}


=head2 check_dependencies

  Title   : check_dependencies
  Usage   : check_dependencies
  Function: Checks a new rule for circular rule dependencies and unfulfilled conditions
  Returns : True if the rule is "clean"
  Args    : A Bio::EnsEMBL::Pipeline::Rule to be evaluated

=cut

sub check_dependencies {

  my ($new_rule) = shift;
  my @new_conditions = @{ $new_rule->list_conditions() };

  # check for unfulfilled conditions - these do not prevent rule insertion
  foreach my $cond (@new_conditions) {
    if ( !$rule_goals{$cond} ) {
      warning( "You should also add a rule that has $cond as its goal,"
             . "\nor this rule will never have its conditions fulfilled.\n" );

      return 1;
    }
  }

  # check for multi stage circles. EEEEK! recursion
  foreach my $old_rule (@existing_rules) {
    check_circles( $new_rule->goalAnalysis()->logic_name(), $old_rule );
  }

  return 1;
}

=head2 check_duplications

  Title   : check_duplications
  Usage   : check_duplications
  Function: Checks to see if a new rule duplicates an exisiting rule
  Returns : True if the rule is "clean"
  Args    : A Bio::EnsEMBL::Pipeline::Rule to be evaluated

=cut
sub check_duplications {
  my ($new_rule) = shift;
  my $goal = $new_rule->goalAnalysis()->logic_name();

RULE: foreach my $old_rule (@existing_rules) {
    my %count;
    my $diff = 0;

    # are the goals different?
    next RULE unless ( $goal eq $old_rule->goalAnalysis()->logic_name() );

    # check the conditions lists
    foreach my $cond ( @{ $new_rule->list_conditions() },
                       @{ $old_rule->list_conditions() } )
    {
      $count{$cond}++;
    }

    foreach my $cond ( keys %count ) {
      next RULE unless $count{$cond} == 2;
    }
    # if we get to here, the two conditions lists are the same
    my $conditions = join ' ', @{ $old_rule->list_conditions() };
    die "Sorry, your new rule duplicates an existing rule:\nRule ID: ",
      $old_rule->dbID, "\n conditions: $conditions\n goal: $goal\n";
  }
  return 1;
} ## end sub check_duplications

=head2 check_circles

  Title   : check_circles
  Usage   : check_circles
  Function: Recursive function called by check_dependencies to look for circular dependencies.
  Returns : True if the rule is "clean"
  Args    : A Bio::EnsEMBL::Pipeline::Rule to be evaluated; A string representing the goal of 
            the original rule

=cut
# recursive routine to check for circular dependencies of one rule on another
sub check_circles {

  my ( $goal, $rule ) = @_;
  my @conditions = @{ $rule->list_conditions() };
  #print STDERR "checking rule ".$rule->dbID." ".$rule->goalAnalysis->logic_name."\n";
  #  die "Sorry, your rule will introduce a circular dependency\n"
  #   if $rule->goalAnalysis()->logic_name() eq $goal;
LOOP: while (@conditions) {
    my $curr_cond = pop @conditions;
    next LOOP if $curr_cond eq $basis;
    # get Rule that has this condition as its goal
    my $new_rule = $rule_goals{$curr_cond};
    if ( !defined($new_rule) ) {
      if ( $curr_cond eq $goal ) {
        die "Sorry, your rule will introduce a circular dependency\n";
      }
      return 1;
    }

    # and recurse
    if ( check_circles( $goal, $new_rule ) ) { return 1; }
  }

  return 1;
} ## end sub check_circles

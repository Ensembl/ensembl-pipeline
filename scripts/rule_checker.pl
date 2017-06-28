#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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


use warnings ;
use strict;
use Bio::EnsEMBL::Pipeline::Utils::Node;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);


my $host = '';
my $port = '3306';
my $pass = '';
my $user = 'ensro';
my $dbname = '';

GetOptions(
  'host|dbhost|h:s'    => \$host,
  'user|dbuser|u:s'    => \$user,
  'dbname|db|D:s'  => \$dbname,
  'pass|dbpass|p:s'    => \$pass,
  'port|dbport|P:n'    => \$port,
);


my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $host,
    -dbname => $dbname,
    -user   => $user,
    -pass   => $pass,
    -port   => $port,
);


check_rules($db);


# Make rules and analyses into graph of nodes
# Find nodes which have no parents (these should be the submits)
# Recursively walk done all the paths from those nodes maintaining a list of the previous nodes
#  At each node on the walk check if we've already seen the node (analysis)
#   If so its a circular walk

sub check_rules {
  my $ra  = $db->get_RuleAdaptor;
  my $aa  = $db->get_AnalysisAdaptor;
  
  my @rules = $ra->fetch_all;
  
  my %nodehash;
  
  RULE: for my $rule (@rules)  {
    my $goal_anal = $rule->goalAnalysis;
  
    if (!exists($nodehash{$goal_anal->logic_name})) {
      my $node = new Node( -anal => $goal_anal );
  
      $nodehash{$goal_anal->logic_name} = $node;
    }
  
    my $node =  $nodehash{$goal_anal->logic_name};
  
    foreach my $cond (@{$rule->list_conditions}) {
      if (exists($nodehash{$cond})) {
        $node->add_parent($nodehash{$cond});
        $nodehash{$cond}->add_child($node);
      } else {
        my $parent_node = new Node( -anal => $aa->fetch_by_logic_name($cond));
  
        $nodehash{$cond} = $parent_node;
  
        $node->add_parent($parent_node);
        $parent_node->add_child($node);
      }
    }
  }
  
  my @tops;
  foreach my $key (keys %nodehash) {
    if (!scalar(@{$nodehash{$key}->parents})) {
      push @tops,$nodehash{$key};
    }
  }
  
  foreach my $top (@tops) {
    my $seenhash = {};
  
    look_for_circles($top, $seenhash);
  }
}

sub look_for_circles {
  my ($node, $seenhash) = @_;

  if (exists($seenhash->{$node})) {
    if ($seenhash->{$node}==2 ) {
      print "ERROR: LOOP in Rules ending at " . $node->anal->logic_name . "\n";
      return;
    }
    print "Loop element " . $node->anal->logic_name . "\n";
  }
  $seenhash->{$node}++;

  foreach my $child (@{$node->children}) {
    look_for_circles($child, $seenhash);
  }
  
  delete($seenhash->{$node});
  return;
}

#!/usr/local/ensembl/bin/perl -w

use strict;
use Node;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Getopt::Long;


my $host = "ecs4";
my $port = 3350;
my $pass = undef;
my $user = "ensro";
my $dbname = "steve_test_pipe";

&GetOptions(
  'host:s'    => \$host,
  'user:s'    => \$user,
  'dbname:s'  => \$dbname,
  'pass:s'    => \$pass,
  'port:n'    => \$port,
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

package FeatureComparison;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Root;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [1]   : FeatureComparison
  Arg [2]   : arrayref, array of features 
  Arg [3]   : arrayref, array of features 
  Arg [4]   : int, binary switch for verbosity
  Function  : create a FeatureComparison module
  Returntype: FeatureComparison
  Exceptions: 
  Example   : my $feature_comp = FeatureComparison->
  (
   -QUERY => $query,
   -TARGET => $target,
   -VERBOSE => $verbose,
  );

=cut


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($query, $target, $verbose) = rearrange( ['QUERY', 'TARGET',
                                               'VERBOSE'], @args);
  $self->query($query);
  $self->target($target);
  $self->verbosity($verbose);
  return $self;
}

=head2 verbosity

  Arg [1]   : FeatureComparison
  Arg [2]   : int
  Function  : toggle as to be verbose or not
  Returntype: int
  Exceptions: 
  Example   : 

=cut

sub verbosity{
  my $self = shift;
  $self->{'verbose'} = shift if(@_);
  return $self->{'verbose'};
}

=head2 query/target

  Arg [1]   : FeatureComparison
  Arg [2]   : arrayref, array of features
  Function  : container for the arrayref of features
  Returntype: arrayref
  Exceptions: throws if not passed an arrayref
  Example   : 

=cut


sub query{
  my ($self, $query) = @_;
  if($query){
    throw("query ust be an array ref not a ".$query) 
      unless(ref($query) eq 'ARRAY'); 
    $self->{'query'} = $query;
  }
  return $self->{'query'};
}

sub target{
  my ($self, $target) = @_;
  if($target){
    throw("target must be an array ref not a ".$target) 
      unless(ref($target) eq 'ARRAY');
    $self->{'target'} = $target;
  }
  return $self->{'target'};
}

sub tables_to_compare{
  my ($self, $target) = @_;
  if($target){
    throw("tables to compare must be an array ref not a ".$target) 
      unless(ref($target) eq 'ARRAY');
    $self->{'tables'} = $target;
  }
  return $self->{'tables'};
}

=head2 fast_sort

  Arg [1]   : FeatureComparison
  Arg [2]   : arrayref array of features
  Function  : sorts the features based on their start coord very quickly
  Returntype: arrayref
  Exceptions: none
  Example   : 

=cut

#this complex map sort map is used as the feature list could potetially
#be very long and this provides a significant speed up from just doing
#sort{$a->start <=> $b->start} as it means the start method is called
#significantly fewer times

sub fast_sort{
  my ($self, $array) = @_;
  my @array = map { $_->[1] } sort { $a->[0] <=> $b->[0] } 
    map { [$_->start, $_] } @$array;
  return \@array;
}

=head2 compare

  Arg [1]   : FeatureComparison
  Function  : Compares the to array of featues noting which features are
  not the same. If the verbose flag is on it will also print out the
  coordinates of the features which are the same
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub compare{
  my ($self) = @_;
  my @query = @{$self->query};
  my @target = @{$self->target};
  my @q_not_matched;
  @query = @{$self->fast_sort(\@query)};
  @target = @{$self->fast_sort(\@target)};
  print "\nThere are ".@query." query features\n";
  print "There are ".@target." target features\n\n";
  my $num_t_feats = scalar(@target);
  my $start_ind = 0;
 QUERY:foreach my $feat(@query){
    if($start_ind > $num_t_feats){
      push(@q_not_matched, $feat);
      next QUERY;
    }
  TARGET:for (my $i = $start_ind; $i< $num_t_feats; $i++){
      my $t_feat = $target[$i];
      next TARGET if(!$t_feat);
      if($t_feat->start == $feat->start &&
         $t_feat->end == $feat->end &&
         $t_feat->strand == $feat->strand){
        if($self->verbosity){
          print "Got match\n";
          print "t_feat = " . $t_feat->start . " " . $t_feat->end . " ".
            $t_feat->strand."\n";
          print "feat   = " . $feat->start . " " . $feat->end . " ".
            $feat->strand."\n";
        }
        $target[$i] = undef;
        next QUERY;
      }elsif($t_feat->start > $feat->start){
        push(@q_not_matched, $feat);
        next QUERY;
      }elsif($t_feat->end < $feat->start){
        $start_ind++;
      }
    }
  }
  print "Unmatched query features\n";
  $self->feature_info(\@q_not_matched);
  print "\nUnmatched target feature\n";
  $self->feature_info(\@target);
}


sub pmatch_compare{
  my ($self) = @_;
  my @query = @{$self->query};
  my @target = @{$self->target};
  my @q_not_matched;
  @query = @{$self->fast_sort(\@query)};
  @target = @{$self->fast_sort(\@target)};
  print "\nThere are ".@query." query features\n";
  print "There are ".@target." target features\n\n";
  my $num_t_feats = scalar(@target);
  my $start_ind = 0;
 QUERY:foreach my $feat(@query){
    if($start_ind > $num_t_feats){
      push(@q_not_matched, $feat);
      next QUERY;
    }
  TARGET:for (my $i = $start_ind; $i< $num_t_feats; $i++){
      my $t_feat = $target[$i];
      next TARGET if(!$t_feat);
      if($t_feat->start == $feat->start &&
         $t_feat->end == $feat->end){
        if($self->verbosity){
          print "Got match\n";
          print "t_feat = " . $t_feat->start . " " . $t_feat->end . "\n";
          print "feat   = " . $feat->start . " " . $feat->end . "\n";
        }
        $target[$i] = undef;
        next QUERY;
      }elsif($t_feat->start > $feat->start){
        push(@q_not_matched, $feat);
        next QUERY;
      }elsif($t_feat->end < $feat->start){
        $start_ind++;
      }
    }
  }
  print "Unmatched query features\n";
  $self->feature_info(\@q_not_matched);
  print "\nUnmatched target feature\n";
  $self->feature_info(\@target);
}

=head2 feature_info

  Arg [1]   : FeatureComparison
  Arg [2]   : arrayref of features
  Function  : prints the start, end and strand of each feature passed in
  Returntype: none
  Exceptions: none
  Example   : 

=cut

sub feature_info{
  my ($self, $array) = @_;
  my $name;
  foreach my $f(@$array){
    if(!$f){
      next;
    }
    $name = $f->analysis->logic_name if(!$name);
    print $name." ".$f->slice->seq_region_name." ".$f->start." ".$f->end.
      " ".$f->strand."\n";
  }
}



sub pipeline_compare{
  my($self, $table_names, $query_db, $ref_db) = @_;
  if(!$table_names){
    $table_names = $self->tables_to_compare;
  }
  $self->db_type_check($query_db);
  $self->db_type_check($ref_db);
  if(@$table_names == 0){
    print "FeatureComparison:pipeline_compare you passed in no tables ".
      "we can't make any comparisons\n";
    return;
  }
  my $query_numbers = {};
  my $ref_numbers = {};
  foreach my $table(@$table_names){
    my $sql = "select logic_name, count(*) from analysis, ";
    $sql .= $table." where analysis.analysis_id=".$table.".analysis_id ".
      "group by logic_name\n";
    $query_numbers = $self->run_query($sql, $query_db, $query_numbers);
    $ref_numbers = $self->run_query($sql, $ref_db, $ref_numbers);
  }
  print "There are ".keys(%$query_numbers)." analysis in the query database ".
    $query_db->dbname." and ".keys(%$ref_numbers)." analysis in the reference ".
      "database ".$ref_db->dbname." all from ".@$table_names." tables\n";
  printf("%-15s %-15s %-15s\n", 'logic_name', 'query_count', 'ref_count');
  printf("%-15s %-15s %-15s\n", '----------', '-----------', '---------');
  foreach my $logic_name(keys(%$query_numbers)){
    printf("%-15s %-15s %-15s\n", $logic_name, $query_numbers->{$logic_name},
           $ref_numbers->{$logic_name});
  }
}


sub db_type_check{
  my ($self, $db) = @_;
  if(!$db->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
    throw("Can't run the RuleManager with $db you need a ".
          "Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
  }
  return;
}



sub run_query{
  my ($self, $sql, $db, $hash) = @_;
  my $sth = $db->prepare($sql);
  $sth->execute;
  while(my($logic_name, $count) = $sth->fetchrow){
    if(!$hash->{$logic_name}){
      $hash->{$logic_name} = $count;
    }else{
      warning($logic_name." already appears in ".$hash." with count ".
            $hash->{$logic_name}." from database ".$db->dbname);
    }
  }
  return $hash;
}

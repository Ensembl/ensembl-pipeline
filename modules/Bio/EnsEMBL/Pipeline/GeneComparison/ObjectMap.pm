#
# Written by Eduardo Eyras
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;

=head1 SYNOPSIS

    my $object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();

    for some $i, $j
    $object_map->match($object_list1[$i], $object_list2[$j], $score);

    my @list1 = $object_map->list1();
    my @list2 = $object_map->list2();
    
    foreach my $object ( @list1 ){
       my @partners = $object_map->partners( $object ) 
	   ...
	   }
    
=head1 DESCRIPTION

Class to hold a map between any two sets of objects

=head1 CONTACT

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap;

use diagnostics;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
@ISA = qw(Bio::EnsEMBL::Root);

######################################################################

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->{_list1} = [];
  $self->{_list2} = [];
  return $self;  
}

############################################################

sub match{
  my ($self, $obj1, $obj2 , $score) = @_;
  $self->add_to_list1( $obj1 );
  $self->add_to_list2( $obj2 );
  $self->partners( $obj1, $obj2 );
  $self->partners( $obj2, $obj1 );
  if ($score){
    $self->{_score}{$obj1}{$obj2} = $score;
    $self->{_score}{$obj2}{$obj1} = $score;
  }
}

############################################################

# stored the score of a match between two elements

sub score{
  my ($self, $obj1, $obj2 , $score) = @_;
  unless ($obj1 && $obj2 ){
    $self->throw("lacking arguments in score()");
  }
  #unless ( $self->{_score}{$obj1}{$obj2} && $self->{_score}{$obj2}{$obj1} ){
  #  $self->{_score}{$obj1}{$obj2} = 0;
  #  $self->{_score}{$obj2}{$obj1} = 0;
  #}
  if ($score){
    $self->{_score}{$obj1}{$obj2} = $score;
    $self->{_score}{$obj2}{$obj1} = $score;
  }
  return $self->{_score}{$obj1}{$obj2};
}    



############################################################

sub add_to_list1{
    my ($self,$obj) = @_;
    if (defined $obj){
	push( @{$self->{_list1}}, $obj );
    }
    return @{$self->{_list1}};
}

############################################################

sub add_to_list2{
    my ($self,$obj) = @_;
    if (defined $obj){
	push( @{$self->{_list2}}, $obj );
    }
    return @{$self->{_list2}};
}

############################################################

sub list1{
    my ($self) = @_;
    my %obj_hash;
    foreach my $obj ( @{$self->{_list1}} ){
	$obj_hash{$obj} = $obj;
    }
    my @list = values %obj_hash;
    return @list;
}

############################################################

sub list2{
    my ($self) = @_;
    my %obj_hash;
    foreach my $obj ( @{$self->{_list2}} ){
	$obj_hash{$obj} = $obj;
    }
    my @list = values %obj_hash;
    return @list;
}

############################################################

sub partners{
    my ($self,$obj1,$obj2) = @_;
    unless ($obj1){
	$self->warn("calling method without argument");
	return;
    }
    unless ( $self->{_partners}{$obj1} ){
	$self->{_partners}{$obj1} = [];
    }
    if ($obj2){
	push ( @{$self->{_partners}{$obj1}}, $obj2 );
    }
    return @{$self->{_partners}{$obj1}};
}    

############################################################
# this method takes all the pairs with the corresponding scores
# as held in the object and finds the best pairs using
# the 'stable-marriage' algorithm.

# it return 

# The Algorithm has the following condition for solution:
# there are no two elements that they're not paired-up to each other but
# they have better score with each other ( $score_matrix is higher ) 
# than with their current corresponding partners

# the main different between this optimization algorithm
# and a 'best-reciprocal-pairs' approach is that
# 'stable-marriage' produces rather 'best-available-pairs', so it keeps matches which
# not being maximal are still optimal. It warranties only on pair
# per element, and this is the best available one, i.e. ' you like
# C. Schiffer but she likes D. Copperfield more than she likes you so
# you have to stick to the next one in your priority list if available'.

sub stable_marriage{
  my ($self) = @_;
  my $verbose = 0;

  # get one list
  my @list1 = $self->list1;
  
  # to keep track of which ones have been married already:
  my %married_object1;
  my %married_object2;

  # to store the potential partners in list2 of each object in list1
  my %candidates_in_list2;

  # to store the chosen partner for each object
  my %partner;

  foreach my $e1 ( @list1 ){
    push( @{$candidates_in_list2{ $e1 }}, $self->partners( $e1 ) );
  }
  
  my @unmarried_ones = @list1;

 MARRIAGE:
  while ( @unmarried_ones ){

    my $object1 = shift @unmarried_ones;
 
    # sort the potential partners by score in descending order
    @{ $candidates_in_list2{ $object1 } } =  
      map  { $_->[1] }
    sort { $b->[0] <=> $a->[0] }
    map  { [ $self->score( $object1, $_ ), $_] } @{ $candidates_in_list2{ $object1 } };
    
    # go over the partners until you get married or run out of partners  
  PARTNER:
    while( @{ $candidates_in_list2{ $object1 }} && !defined($married_object1{$object1}) ){
      
      #print STDERR "checking partner list for $this_target\n";
      my $potential_partner_in_list2 = shift( @{ $candidates_in_list2{ $object1 } } );
      
      print STDERR "looking at $potential_partner_in_list2\n" if $verbose;
      # check whether it is already married
      if ( $married_object2{ $potential_partner_in_list2 } 
	   && $married_object2{ $potential_partner_in_list2 } == 1 ){
	
	# is it married to another target?
	if ( $partner{$potential_partner_in_list2} 
	     &&  !( $partner{ $potential_partner_in_list2 } eq $object1 ) ){
	  
	  # is it a 'worse' marriage?
	  if ( $self->score( $partner{ $potential_partner_in_list2 }, $potential_partner_in_list2) 
	       < $self->score( $object1, $potential_partner_in_list2 ) ){
	    
	    # put the divorced one back into the pool only if it has more potential partners
	    if ( @{ $candidates_in_list2{ $partner{ $potential_partner_in_list2} } } ){
	      push ( @unmarried_ones, $partner{ $potential_partner_in_list2} );
	    }

	    # divorce the 'worse partner'
	    print STDERR "divorcing ".$partner{ $potential_partner_in_list2}."\n" if $verbose;
	    delete $married_object1{ $partner{ $potential_partner_in_list2 } };
	    delete $partner{ $partner{ $potential_partner_in_list2 } };
	    delete $partner{ $potential_partner_in_list2 };
	    
	    # let be happier marriage
	    $married_object1{ $object1 } = 1;
	    $married_object2{ $potential_partner_in_list2 } = 1;
	    $partner{ $potential_partner_in_list2 } = $object1;
	    $partner{ $object1 } = $potential_partner_in_list2;
	    print STDERR "new marriage: $object1 -- $potential_partner_in_list2\n" if $verbose;
	    next MARRIAGE;
	    
	  }
	  else{
	    # look at the next potential partner in list2
	    next PARTNER;
	  }
	}
	# hmm, this object2 ( in list2)  is married, to whom?
	elsif ( $partner{ $potential_partner_in_list2 } eq $object1 ) {
	  # hey, we have already a happy couple
	  $partner{ $object1 } = $potential_partner_in_list2;
	  next MARRIAGE;
	}
	elsif ( !defined( $partner{ $potential_partner_in_list2 } ) ){
	  # we have a cheater!
	  $married_object2{ $potential_partner_in_list2 } = 0;
	  next PARTNER;
	}
      }
      else{
	
	# this object2 ( in list 2 ) is still single, let be marriage:
	$married_object1{ $object1 } = 1;
	$married_object2{ $potential_partner_in_list2 } = 1;
	$partner{ $potential_partner_in_list2 } = $object1;
	$partner{ $object1 } = $potential_partner_in_list2;
	#print STDERR "setting partner{ $object1 } = $potential_partner_in_list2\n";
	next MARRIAGE;
      }
      
    } # end of PARTNER
    
  }   # end of MARRIAGE
  
  my $new_object_map = Bio::EnsEMBL::Pipeline::GeneComparison::ObjectMap->new();
  foreach my $object1 ( keys %married_object1 ){
    $new_object_map->match( $partner{ $partner{$object1} }, $partner{$object1}, $self->score( $object1, $partner{$object1}) );
  }
  return $new_object_map;
  
}
############################################################

1;




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
    $object_map->match($object_list1[$i], $object_list2[$j]);

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
    my ($self, $obj1, $obj2 ) = @_;
    $self->add_to_list1( $obj1 );
    $self->add_to_list2( $obj2 );
    $self->partners( $obj1, $obj2 );
    $self->partners( $obj2, $obj1 );
}

sub add_to_list1{
    my ($self,$obj) = @_;
    if (defined $obj){
	push( @{$self->{_list1}}, $obj );
    }
    return @{$self->{_list1}};
}


sub add_to_list2{
    my ($self,$obj) = @_;
    if (defined $obj){
	push( @{$self->{_list2}}, $obj );
    }
    return @{$self->{_list2}};
}

sub list1{
    my ($self) = @_;
    my %obj_hash;
    foreach my $obj ( @{$self->{_list1}} ){
	$obj_hash{$obj} = $obj;
    }
    my @list = values %obj_hash;
    return @list;
}


sub list2{
    my ($self) = @_;
    my %obj_hash;
    foreach my $obj ( @{$self->{_list2}} ){
	$obj_hash{$obj} = $obj;
    }
    my @list = values %obj_hash;
    return @list;
}


sub partners{
    my ($self,$obj1,$obj2) = @_;
    unless ($obj1){
	$self->warn("caling method without argument");
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

1;




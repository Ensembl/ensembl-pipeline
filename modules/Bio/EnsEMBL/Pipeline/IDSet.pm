

package Bio::EnsEMBL::Pipeline::IDSet;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);


sub new{
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  my ($list) = $self->_rearrange([qw(ID_LIST)], @args);

  $self->{'_ID_list'} = []; 

  if($list){
    if(ref($list) eq 'ARRAY'){
      $self->{'_ID_list'} = $list;
    }else{
      $self->throw("ID_LIST need to be a listref not a $list $!");
    }
  }
  return $self;
}



sub add_element{
  my ($self, @elements) = @_;
  if(@elements){
    push(@{$self->{'_ID_list'}}, @elements);
  }
}


sub delete_element{
  my ($self, @elements) = @_;
  my @ids = @{$self->{'_ID_list'}};
  foreach my $e(@elements){
    my @new_ids = grep $_ != $e, @ids;
    @ids = @new_ids;
    @new_ids = ();
  }
  
  @{$self->{'_ID_list'}} = @ids;
}

sub ID_list{
  my ($self) = @_;
  return $self->{'_ID_list'};
}


sub calculate{
  my ($self, $idlist) = @_;
  
  $idlist->isa('Bio::EnsEMBL::Pipeline::IDSet') ||$self->throw("Need to passs these methods (and/or/not) an Bio::EnsEMBL::Pipeline::IDSet object otherwise won't work you passed in $idlist : $!");
  
  my %count;
  my @own_ids =  @{$self->{'_ID_list'}};
  my @comp_ids = @{$idlist->ID_list};
  foreach my $e(@own_ids, @comp_ids) {$count{$e}++}
  return \%count;
}

sub and{
  my ($self, $idlist) = @_;
  
  my %count = %{$self->calculate($idlist)};
  my @and;
  foreach my $e(keys(%count)){
    if($count{$e} == 2){
      push(@and, $e);
    }
  }
  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -ID_LIST => \@and,
						);
  return $idset;
}


sub not{
  my ($self, $idlist) = @_;
  my %count = %{$self->calculate($idlist)};
  my @and;
  my @own_ids =  @{$self->{'_ID_list'}};
  foreach my $e(@own_ids){
    if($count{$e} == 1){
      push(@and, $e);
    }
  }
  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -ID_LIST => \@and,
						);
  return $idset;
}


sub or{
  my ($self, $idlist) = @_;
  
  my %count = %{$self->calculate($idlist)};
  my @and = keys(%count);
  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -ID_LIST => \@and,
						);
  return $idset;
}

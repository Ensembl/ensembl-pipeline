use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::RDB::Slice;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::Task::RDB;
use Bio::EnsEMBL::Pipeline::IDSet;


@ISA = ('Bio::EnsEMBL::Pipeline::Task::RDB');


#the new is used from the base class as this constructor wouldn't need
#to do any additional work


=head2 generate_input_ids

  Arg [1]   : none
  Function  : fetchs a list of contig names and creates an IDSet using
  them
  Returntype: Bio::EnsEMBL::Pipeline::IDSet 
  Exceptions: none
  Caller    : 
  Example   : my $names = $self->generate_input_idsl

=cut



sub generate_input_ids{
  my ($self, $size, $overlap) = @_;

  if(!$size){
    $self->throw("To generate Slice based input_ids need a size\n");
  }

  if(!$overlap){
    $overlap = 0;
  }

  my @input_ids;

  my @chromosomes = @{$self->get_chromosomes};

  foreach my $chr(@chromosomes){
    my $length = $chr->length;
    my $count = 1;
    while ($count < $length) {
      my $start = $count;IDs 
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      push(@input_ids, $input_id);
      $count = $count + $size - $overlap;
    }
  }

  my $idset = Bio::EnsEMBL::Pipeline::IDSet(
					    -id_list => \@input_ids,
					   );

  return $idset;

  
}





=head2 get_input_ids

  Arg [1]   : none
  Function  : returns an IDSet of contig names
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: none
  Caller    : 
  Example   : my input_ids = $self->get_input_ids

=cut



sub get_input_ids{
  my ($self) = @_;

  if(!$self->{'input_ids'}){
    my $idset = $self->generate_input_ids;
    $self->{'input_ids'} = $idset;
  }
  
  return $self->{'input_ids'};
}



#warning this method assumes the length in the chromosome table is
#correct

sub get_chromosomes{
  my ($self) = @_;

  if(!$self->{'chromosomes'}){
    my $chr_adp = $self->db->get_ChromosomeAdaptor;
    my $chromosomes = $chr_adp->fetch_all;
    $self->{'chromsomes'} = $chromosomes;
  }

  return $self->{'chromosomes'};
}

use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;

use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::Root');


=head2 new

  Arg [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Function  : creates an InputIDFactory object
  Returntype: Bio::EnsEMBL::Pipeline::Utils::InputIDFactory
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  
  my $self = bless({}, $class);

  $self->{'db'} = undef;

  my ($db)=$self->_rearrange([qw(DB)], @_);

  $self->db($db) if($db);

  $self->throw("you need to pass at least a DBAdaptor to an InputIDFactory") unless($self->db);

  return $self;
}



sub db{
  my $self = shift;

  if(@_){
    $self->{'db'} = shift;
  }

  return $self->{'db'};
}

=head2 generate_input_ids

  Arg [1]   : none
  Function  : on the basis of whats in config decides which method to 
  call to generate the input_ids
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: throws if the type isn't recognised'
  Caller    : 
  Example   : 

=cut



sub generate_slice_input_ids {
    my ($self,$cs_name, $cs_version, $size,$overlap) = @_;

    $overlap = 0 if (!$overlap);
    $size = 0 if (!$size);
    if ($size < 0) {
       $self->throw("Slice size must be >= 0. Currently " . $size);
    }

    my @ids = @{$self->get_slice_names($cs_name, $cs_version,
                                       $size,$overlap)};

    return @ids;
}


sub get_slice_names{
  my ($self, $cs_name, $cs_version, $size, $overlap) = @_;

  $overlap = 0 if (!$overlap);
  $size = 0 if (!$size);
  $cs_version = '' unless($cs_version);
  my $csa = $self->db->get_CoordSystemAdaptor();
  my $sa = $self->db->get_SliceAdaptor();

  my @slices = @{$sa->fetch_all($cs_name, $cs_version)};
  my @ids;
  foreach my $slice(@slices){
    push(@ids, $slice->name);
  }

  return \@ids;
}





=head2 get_slice_names

  Arg [1]   : size, int
  Arg [2]   : overlap, int
  Function  : produces a set of slice names based on the size and overlap
  specified in the format chr_name.start-end
  Returntype:  Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: throws if it has no core db connection
  Caller    : 
  Example   : 

=cut


sub generate_non_redundant_input_ids{
  my ($self) = @_;

  my @ids = @{$self->get_non_redundant_ids()};

  return @ids;
}


sub get_non_redundant_ids{
  my ($self) = @_;

  my @slices = @{$self->db->get_SliceAdaptor->fetch_all_non_redundant()};
  my @ids;
  foreach my $slice(@slices){
    push(@ids, $slice->name);
  }

  return \@ids;
  
}


sub get_filenames{
  my ($self, $dir, $regex) = @_;
  if(!$dir){
    $self->throw("need a directory inorder to fetch the filenames to be used as input_ids $!");
  }

  my @input_ids;

  opendir(DIR, $dir);   
  my @allfiles = readdir DIR;
  closedir DIR;
	
  foreach my $f(@allfiles) {
    if($f eq '.' || $f eq '..'){
      next;
    }elsif(-d $f){
      next;
    }else{
      my $file;
      if($regex){
	if($f =~ m|$regex|){
	  $file = $f;
	}
      }else{
	$file = $f;
      }
      push(@input_ids, $file) if($file);
    }    
  }
  
  return \@input_ids;
}


sub generate_filename_input_ids{
   my ($self, $dir, $regex) = @_;

   my $ids = $self->get_filenames($dir, $regex);

   return @$ids;
}

sub get_translation_ids{
  my ($self) = @_;

  my $ids = $self->db->get_TranslationAdaptor->list_dbIDs;

  return $ids;
}


sub generate_translation_id_input_ids{
  my ($self) = @_;

  my $ids = $self->get_translation_ids;
  
  print STDERR "have ".@$ids." ids\n";
  return @$ids;
}




1;

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

sub generate_contig_input_ids{
  my ($self) = @_;

  my @ids  = @{$self->get_contig_names};

  print "Ids @ids\n";
  return @ids;

}

sub generate_slice_input_ids {
    my ($self,$size,$overlap) = @_;

    $overlap = 0 if (!defined($overlap));

    if ($size <= 0) {
       $self->throw("Slice size must be > 1. Currently " . $size);
    }

    my @ids = @{$self->get_slice_names($size,$overlap)};

    return @ids;
}


=head2 get_contig_names

  Arg [1]   : none
  Function  : uses the core dbconnection to get a list of contig names
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: throws if there is no db connection
  Caller    : 
  Example   : 

=cut

sub get_contig_names{
    my ($self) = @_;
    
    if(!$self->db){
	$self->throw("if you getting contig names InputIDFactory needs a dbconnection to a core db");
    }
    my $rawcontig_adaptor = $self->db->get_RawContigAdaptor;

    my $names = $rawcontig_adaptor->fetch_all_names;

    return $names;
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


sub get_slice_names{
  my ($self, $size, $overlap) = @_;

  if(!$self->db) {
    $self->throw("if you're getting slice names InputIDFactory needs a dbconnection to a core db");
  }

  my @input_ids;

  my @chromosomes = @{$self->db->get_ChromosomeAdaptor->fetch_all};

  foreach my $chr (@chromosomes){

    my $query = "select min(chr_start) from assembly where chromosome_id = " . $chr->dbID;
    my $sth   = $self->db->prepare($query);
    my $res   = $sth->execute;

    my $chrstart;
   
    while (my $ref = $sth->fetchrow_arrayref) {
      $chrstart = $ref->[0];
    }

    $query = "select max(chr_end) from assembly where chromosome_id = " . $chr->dbID;
    $sth   = $self->db->prepare($query);
    $res   = $sth->execute;

    my $length;

    while (my $ref = $sth->fetchrow_arrayref) {
      $length = $ref->[0];
    }

    my $count = $chrstart;

    while ($count < $length) {
      my $start = $count;
      my $end   = $count + $size - 1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr->chr_name . "." . $start . "-" .  $end;

      push(@input_ids, $input_id);

      $count = $count + $size - $overlap;
    }
  }

  return \@input_ids;

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

use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::Utils::InputIDFactory;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::IDSet;


@ISA = ('Bio::EnsEMBL::Root');



=head2 new

  Arg [1]   : Bio::EnsEMBL::Pipeline::Config
  Arg [2]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Arg [3]   : string
  Function  : creates an InputIDFactory object
  Returntype: Bio::EnsEMBL::Pipeline::Task::Utils::InputIDFactory
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  
  my $self = bless({}, $class);

  $self->{'config'} = undef;
  $self->{'db'} = undef;
  $self->{'taskname'} = undef;

  my ($config, $db, $taskname)=$self->_rearrange([qw(CONFIG DB TASKNAME)], @_);

  $self->config($config) if($config);
  $self->db($db) if($db);
  $self->taskname($taskname) if($taskname);
  $self->throw("you need to pass at least a config object and a taskname to an InputIDFactory $!") unless($self->config && $self->taskname);

  return $self;
}



=head2 getters/setters

  Arg [1]   : generally objecst or string
  Function  : will set an arg to a variable and return that arg 
  subsquently
  Returntype: objects or strings
  Exceptions: none
  Caller    : 
  Example   : 

=cut




sub config{
  my $self = shift;

  if(@_){
    $self->{'config'} = shift;
  }

  return $self->{'config'};
}

sub db{
  my $self = shift;

  if(@_){
    $self->{'db'} = shift;
  }

  return $self->{'db'};
}

sub taskname{
  my $self = shift;

  if(@_){
    $self->{'taskname'} = shift;
  }

  return $self->{'taskname'};
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



sub generate_input_ids{
  my ($self) = @_;

  my $type_string = $self->config->get_parameter($self->taskname, 'input_id');
  my (@other_info) = split /:/, $type_string;
  my $type = lc(shift @other_info);
  my $idset;
  if($type eq 'contig'){
    $idset = $self->get_contig_names;
  }elsif($type eq 'slice'){
    my ($size, $overlap) = @other_info;
    my $warn = 0;
    if(!$size){
      print STDERR "InputIDFactory:124 No size defined being set to 1MB\n";
      $size = 1000000;
    }
    if(!$overlap){
      print STDERR "InputIDFactory:128 No overlap defined being set to 0\n";
      $overlap = 0;
    }
    $idset = $self->get_slice_names($size, $overlap);
  }elsif($type eq 'chromosome'){
    $idset = $self->get_chromosome_names;
  }elsif($type eq 'file'){
    my ($dir) = @other_info;
    $idset = $self->get_file_names($dir);
  }else{
    $self->throw("don't recognise input_id type $type from string $type_string");
  }

  return $idset;
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
    $self->throw("if you getting contig names InputIDFactory needs a dbconnection to a core db $!");
  }
  my $rawcontig_adaptor = $self->db->get_RawContigAdaptor;

  my $names = $rawcontig_adaptor->fetch_all_names;

  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -id_list => $names,
						);

  return $idset;
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

  if(!$self->db){
    $self->throw("if you getting contig names InputIDFactory needs a dbconnection to a core db $!");
  }

  my @input_ids;

  my @chromosomes = @{$self->get_Chromosomes};

  foreach my $chr(@chromosomes){
    my $length = $chr->length;
    my $count = 1;
    while ($count < $length) {
      my $start = $count;
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr->chr_name . "." . $start . "-" .  $end;
      push(@input_ids, $input_id);
      $count = $count + $size - $overlap;
    }
  }

  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -id_list => \@input_ids,
						);

  return $idset;

}


=head2 get_chromosome_names

  Arg [1]   : none
  Function  : produces slice type input_ids in the format chr_name.1-length
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: throws if no core db connection
  Caller    : 
  Example   : 

=cut



sub get_chromosome_names{
  my ($self) = @_;

  if(!$self->db){
    $self->throw("if you getting contig names InputIDFactory needs a dbconnection to a core db $!");
  }

  my @input_ids;

  my @chromosomes = @{$self->get_Chromosomes};

  foreach my $chr(@chromosomes){
    my $length = $chr->length;
    my $input_id = $chr->chr_name . ".1-" .  $length;
    push(@input_ids, $input_id);
    
  }

  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -id_list => \@input_ids,
						);

  return $idset;

}
#this method assumes the lengths in the chromosome table are correct


=head2 get_Chromosomes

  Arg [1]   : none 
  Function  : gets chromosome objects from the database
  Returntype: array ref of Bio::EnsEMBL::Chromosomes
  Exceptions: none
  Caller    : 
  Example   : 

=cut


sub get_Chromosomes{
  my ($self) = @_;

  if(!$self->{'chromosomes'}){
    my $chr_adp = $self->db->get_ChromosomeAdaptor;
    my $chromosomes = $chr_adp->fetch_all;
    $self->{'chromosomes'} = $chromosomes;
  }
  return $self->{'chromosomes'};
}



=head2 get_file_names

  Arg [1]   : directory path (string)
  Function  : to get a list of filenames from the directory
  Returntype: Bio::EnsEMBL::Pipeline::IDSet
  Exceptions: throws if no directory is passed in 
  Caller    : 
  Example   : 

=cut


sub get_file_names{
  my ($self, $dir) = @_;
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
    }else{
      push(@input_ids, $f);
    }    
  }
  my $idset = Bio::EnsEMBL::Pipeline::IDSet->new(
						 -id_list => \@input_ids,
						);
  
  return $idset;
}


1;

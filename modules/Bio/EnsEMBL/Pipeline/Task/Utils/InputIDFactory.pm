use strict;
use warnings;
package Bio::EnsEMBL::Pipeline::Task::Utils::InputIDFactory;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::IDSet;


@ISA = ('Bio::EnsEMBL::Root');



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
}




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
      $warn = 1;
      $size = 1000000;
    }
    if(!$overlap){
      $warn = 1;
      $overlap = 0;
    }
    $self->warn("your type string $type_string specifed slice without giving a size or overlap assuming you want 1MB pieces with no overlap\n");
    $idset = $self->get_slice_names($size, $overlap);
  }elsif($type eq 'chromosome'){
    $idset = $self->get_chromosome_names;
  }elsif($type eq 'file'){
    my ($dir) = @other_info;
    $idset = $self->get_file_names($dir);
  }else{
    $self->throw("don't recognise input_id type $type from string $type_string");
  }
}



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


sub get_slice_names{
  my ($self, $size, $overlap) = @_;

  if(!$self->db){
    $self->throw("if you getting contig names InputIDFactory needs a dbconnection to a core db $!");
  }

  my @input_ids;

  my @chromosomes = @{$self->get_chromosomes};

  foreach my $chr(@chromosomes){
    my $length = $chr->length;
    my $count = 1;
    while ($count < $length) {
      my $start = $count;
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

sub get_chromosome_names{
  my ($self) = @_;

  if(!$self->db){
    $self->throw("if you getting contig names InputIDFactory needs a dbconnection to a core db $!");
  }

  my @input_ids;

  my @chromosomes = @{$self->get_chromosomes};

  foreach my $chr(@chromosomes){
    my $length = $chr->length;
    my $input_id = $chr . ".1-" .  $length;
    push(@input_ids, $input_id);
    
  }

  my $idset = Bio::EnsEMBL::Pipeline::IDSet(
					    -id_list => \@input_ids,
					   );

  return $idset;

}
#this method assumes the lengths in the chromosome table are correct
sub get_chromosomes{
  my ($self) = @_;

  if(!$self->{'chromosomes'}){
    my $chr_adp = $self->db->get_ChromosomeAdaptor;
    my $chromosomes = $chr_adp->fetch_all;
    $self->{'chromsomes'} = $chromosomes;
  }

  return $self->{'chromosomes'};
}


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
  my $idset = Bio::EnsEMBL::Pipeline::IDSet(
					    -id_list => \@input_ids,
					   );
  
  return $idset;
}


1;

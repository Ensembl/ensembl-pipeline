package Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Analysis;
use vars qw(@ISA);

@ISA = ('Bio::EnsEMBL::Root');


use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

=head2 new

  Arg [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Arg [2]   : int, toggle for slice based input_ids
  Arg [3]   : int toggle for single input_ids
  Arg [4]   : int toggle for filename based input_ids
  Arg [5]   : int 
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

  my ($db, 
      $slice, 
      $single, 
      $file,
      $translation_id, 
      $slice_size, 
      $slice_overlaps, 
      $seq_level,
      $top_level,
      $dir, 
      $regex, 
      $single_name,
      $verbose, 
      $logic_name, 
      $input_id_type,
      $insert_analysis, 
      $coord_system,
      $coord_system_version,
      $seq_region_name,
      $seq_region_start,
      $seq_region_end,
      $hap_pair,
      $include_non_reference)=rearrange([qw(DB 
                                            SLICE 
                                            SINGLE 
                                            FILE
                                            TRANSLATION_ID 
                                            SLICE_SIZE 
                                            SLICE_OVERLAPS 
                                            SEQ_LEVEL
                                            TOP_LEVEL 
                                            DIR 
                                            REGEX 
                                            SINGLE_NAME
                                            VERBOSE 
                                            LOGIC_NAME 
                                            INPUT_ID_TYPE
                                            INSERT_ANALYSIS 
                                            COORD_SYSTEM
                                            COORD_SYSTEM_VERSION
                                            SEQ_REGION_NAME
                                            SEQ_REGION_START
                                            SEQ_REGION_END
                                            HAP_PAIR
 					    INCLUDE_NON_REFERENCE
                                           )], @_);

  $slice = 0 unless ($slice);
  $single = 0 unless ($single);
  $file = 0 unless ($file);
  $hap_pair = 0 unless ($hap_pair);
  $include_non_reference = 0 unless ($include_non_reference);
  $translation_id = 0 unless($translation_id);

  if(!$db){
    throw("You can't create and store input_ids without a dbadaptor\n");
  }
  $self->db($db);
  if(!$slice && !$file && !$translation_id && !$single && !$hap_pair){
    throw("You must define one of these options SLICE, FILE, SINGLE ".
          "TRANSLATION_ID, HAP_PAIR for the input id factory to work");
  }
  print "slice: ",$slice,"\n";
print "file: ",$file,"\n";
print "translation_id: ",$translation_id,"\n";
print "single: ",$single,"\n";
print "hap_pair: ",$hap_pair,"\n";
  print "SUM: ",$slice+$file+$translation_id+$single+$hap_pair,"\n";

  if(($slice+$file+$translation_id+$single+$hap_pair) > 1){
    throw("You must only specify one of these options SLICE, FILE, SINGLE ".
          "TRANSLATION_ID, HAP_PAIR for the input id factory to work");
  }
  $self->slice($slice) if($slice);
  $self->top_level($top_level) if($top_level);
  $self->seq_level($seq_level) if($seq_level);
  $self->coord_system($coord_system) if($coord_system);
  if($slice && !$self->coord_system){
      throw("You must specify a coordinate system if you want slice ".
            "input ids created\n");
  }
#  if($hap_pair && !$self->coord_system){
#    throw("You must specify a coordinate system if you want slice ".
#          "input ids created\n");
#  }

  $self->coord_system_version($coord_system_version) 
      if($coord_system_version);
  $self->seq_region_name($seq_region_name)
      if defined $seq_region_name;
  $self->seq_region_start($seq_region_start)
      if defined $seq_region_start;
  $self->seq_region_end($seq_region_end)
      if defined $seq_region_end;

  $self->file($file) if($file);
  $self->single($single) if($single);
  $self->include_non_reference($include_non_reference) if($include_non_reference);
  $self->translation_id($translation_id) if($translation_id);
  $self->hap_pair($hap_pair) if($hap_pair);
  if(!$logic_name){
    throw("Must have a logic_name otherwise don't know which analysis to ".
          "store the input ids under");
  }
  if($insert_analysis && !$input_id_type){
    throw("if you want your analysis object to be inserted into the ".
          "database you must also provide an input_id_type");
  }
  my $analysis = $self->get_analysis($logic_name, $input_id_type, 
                                     $insert_analysis);
  $self->logic_name($logic_name);
  $self->slice_size($slice_size) if($slice_size);
  $self->slice_overlaps($slice_overlaps) if($slice_overlaps);
  $self->dir($dir) if($dir);
  $self->regex($regex) if($regex);
  $self->single_name($single_name) if($single_name);
  return $self;
}


#container methods

sub db{
  my ($self, $db) = @_;

  if($db){
    if(!$db->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
      throw("Can't run the RuleManager with $db you need a ".
            "Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{'dbadaptor'} = $db;
  }

  return $self->{'dbadaptor'};
}

sub stateinfocontainer{
  my ($self, $adaptor) = @_;

  if($adaptor){
    $self->{'stateinfocontainer'} = $adaptor;
  }
  if(!$self->{'stateinfocontainer'}){
    my $stateinfocontainer = $self->db->get_StateInfoContainer;
    $self->{'stateinfocontainer'} = $stateinfocontainer;
  }
  return $self->{'stateinfocontainer'};
}
sub analysis_adaptor{
  my ($self, $adaptor) = @_;

  $self->{'analysis_adaptor'} = $adaptor;

  if(!$self->{'analysis_adaptor'}){
    $self->{'analysis_adaptor'} = $self->db->get_AnalysisAdaptor;
  }
  return $self->{'analysis_adaptor'};
}

sub slice{
  my $self = shift;
  $self->{'slice'} = shift if(@_);
  return $self->{'slice'};
}

sub file{
  my $self = shift;
  $self->{'file'} = shift if(@_);
  return $self->{'file'};
}

sub translation_id{
  my $self = shift;
  $self->{'translation_id'} = shift if(@_);
  return $self->{'translation_id'};
}

sub single{
  my $self = shift;
  $self->{'single'} = shift if(@_);
  return $self->{'single'};
}

sub hap_pair{
  my $self = shift;
  $self->{'hap_pair'} = shift if(@_);
  return $self->{'hap_pair'};
}

sub include_non_reference{
  my $self = shift;
  $self->{'include_non_reference'} = shift if(@_);
  return $self->{'include_non_reference'};
}

sub slice_size{
  my $self = shift;
  $self->{'slice_size'} = shift if(@_);
  return $self->{'slice_size'};
}
sub slice_overlaps{
  my $self = shift;
  $self->{'slice_overlaps'} = shift if(@_);
  return $self->{'slice_overlaps'};
}
sub dir{
  my $self = shift;
  $self->{'dir'} = shift if(@_);
  return $self->{'dir'};
}
sub regex{
  my $self = shift;
  $self->{'regex'} = shift if(@_);
  return $self->{'regex'};
}
sub single_name{
  my $self = shift;
  $self->{'single_name'} = shift if(@_);
  return $self->{'single_name'} || 'genome';
}
sub coord_system{
  my $self = shift;
  $self->{'coord_system'} = shift if(@_);
  return $self->{'coord_system'};
}
sub coord_system_version{
  my $self = shift;
  $self->{'coord_system_version'} = shift if(@_);
  return $self->{'coord_system_version'};
}
sub seq_region_name {
  my $self = shift;
  $self->{'seq_region_name'} = shift if(@_);
  return $self->{'seq_region_name'};
}
sub seq_region_start {
  my $self = shift;
  $self->{'seq_region_start'} = shift if(@_);
  return $self->{'seq_region_start'};
}
sub seq_region_end {
  my $self = shift;
  $self->{'seq_region_end'} = shift if(@_);
  return $self->{'seq_region_end'};
}
sub top_level{
  my $self = shift;
  $self->{'top_level'} = shift if(@_);
  if($self->{'top_level'}){
    $self->coord_system('toplevel');
  }
  return $self->{'top_level'};
}
sub seq_level{
  my $self = shift;
  $self->{'seq_level'} = shift if(@_);
  if($self->{'seq_level'}){
    $self->coord_system('seqlevel');
  }
  return $self->{'seq_level'};
}


sub get_analysis{
  my ($self, $logic_name, $input_id_type, $insert) = @_;

  my $analysis;
  if($logic_name && $input_id_type && $insert){
    $analysis = Bio::EnsEMBL::Pipeline::Analysis->new;
    $analysis->logic_name($logic_name);
    $analysis->input_id_type($input_id_type);
    $self->analysis_adaptor->store($analysis);
  }elsif($logic_name && !$insert){
    $analysis = $self->analysis_adaptor->fetch_by_logic_name($logic_name);
  }
  if($analysis){
    $self->{'analysis'} = $analysis;
  }
  return $self->{'analysis'};
}
sub logic_name{
  my $self = shift;
  $self->{'logic_name'} = shift if(@_);
  if(!$self->{'logic_name'}){
    $self->{'logic_name'} = $self->get_analysis->logic_name;
  }
  return $self->{'logic_name'};
}

sub input_id_type{
  my $self = shift;
  $self->{'input_id_type'} = shift if(@_);
  if(!$self->{'input_id_type'}){
    $self->{'input_id_type'} = $self->get_analysis->input_id_type;
  }
  return $self->{'input_id_type'};
}

sub input_ids{
  my ($self, $input_ids) = @_;

  if($input_ids){
   throw("Must has an array ref of input_ids not a $input_ids ") 
      unless(ref($input_ids) eq 'ARRAY');
    $self->{'input_ids'} = $input_ids;
  }

  return $self->{'input_ids'};
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
  my $ids;
  if($self->slice){
    $ids = $self->get_slice_names;
  }elsif($self->file){
    $ids = $self->get_filenames;
  }elsif($self->single){
    $ids = $self->get_single;
  }elsif($self->translation_id){
    $ids = $self->get_translation_id;
  }elsif($self->hap_pair){
    $ids = $self->get_hap_pairs;
  }else{
    throw("Reaching this point means you haven't created InputIDFactory ".
          "without selecting what type of input_id to create this won't ".
          "work");
  }
  $self->input_ids($ids);
  return $ids;
}



sub get_slice_names{
  my ($self) = @_;

  $self->slice_size(0) if(!$self->slice_size);
  $self->slice_overlaps(0) if(!$self->slice_overlaps);
  $self->coord_system_version('') if(!$self->coord_system_version);
  if ($self->slice_size && $self->slice_size < 0) {
    throw("Slice size must be >= 0. Currently " . 
          $self->slice_size);
  }
  
  my $csa = $self->db->get_CoordSystemAdaptor();
  my $sa = $self->db->get_SliceAdaptor();
  
  my $slices;
  if ($self->seq_region_name) {
    my $sname = sprintf("%s:%s:%s:%s:%s:",
                        $self->coord_system,
                        $self->coord_system_version,
                        $self->seq_region_name,
                        $self->seq_region_start,
                        $self->seq_region_end);

    $slices = [$sa->fetch_by_name($sname)];
  } else {
    $slices = $sa->fetch_all($self->coord_system, 
                             $self->coord_system_version,
			     $self->include_non_reference);
  }  

  if($self->slice_size > 0){
    $slices = split_Slices($slices,$self->slice_size,$self->slice_overlaps);
  }

  my @ids;
  foreach my $slice(@$slices){
    push(@ids, $slice->name);
  }

  return \@ids;
}

sub get_hap_pairs{
  my ($self) = @_;

  my @ids;

  my $assembly_exception_feature_adaptor =$self->db->get_AssemblyExceptionFeatureAdaptor(); 
  my @assembly_exception_features = @{$assembly_exception_feature_adaptor->fetch_all};

  $self->coord_system_version('') if(!$self->coord_system_version);
  
  foreach my $aef(@assembly_exception_features){
    if ($aef->type eq 'HAP'){
       
      my $sname = sprintf("%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:",
                          $aef->slice->coord_system->name,
                          $aef->slice->coord_system->version,
                          $aef->slice->seq_region_name,
                          $aef->start,
                          $aef->end,
                          $aef->alternate_slice->coord_system->name,
                          $aef->alternate_slice->coord_system->version,
                          $aef->alternate_slice->seq_region_name,
                          $aef->alternate_slice->start,
                          $aef->alternate_slice->end,);
          
      push(@ids, $sname);
      
    }
  }
  return \@ids;
}

sub get_filenames{
  my ($self) = @_;

  if(!$self->dir){
    $self->throw("need a directory inorder to fetch the filenames to be used as input_ids $!");
  }

  my @input_ids;

  opendir(DIR, $self->dir);   
  my @allfiles = readdir DIR;
  closedir DIR;
  my $regexp = $self->regex();

  foreach my $f(@allfiles) {
    if($f eq '.' || $f eq '..'){
      next;
    }elsif(-d $f){
      next;
    }else{
      my $file;
      if($regexp){

        if($f =~ m|$regexp|){
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


sub get_translation_id{
  my ($self) = @_;

  my $ids = $self->db->get_TranslationAdaptor->list_dbIDs;

  return $ids;
}


sub get_single{
  my ($self) = @_;
  my @ids = ($self->single_name);
  return \@ids;
}



sub store_input_ids{
  my ($self) = @_;
  my $ids = $self->input_ids;
  foreach my $id(@$ids){
    eval{
      $self->stateinfocontainer->store_input_id_analysis($id, $self->get_analysis, '');
    };
    if($@){
      throw("Error storing input_id $id for analysis ".
            $self->get_analysis->logic_name."\n".$@);
    }
  }
  return 1;
}


sub get_id_hash{
  my ($self) = @_;
  my $ids = $self->input_ids;
  my $id_hash = {};
  $id_hash->{$self->input_id_type} = {};
  foreach my $id(@$ids){
    $id_hash->{$self->input_id_type}->{$id} = 1;
  }
  return $id_hash;
}


1;

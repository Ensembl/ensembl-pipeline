package Bio::EnsEMBL::Pipeline::Runnable::Protein_Annotation;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Copy qw(mv cp);
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::ProteinFeature;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);



sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    
    $self->{'_flist'}     = [];    # an array of Bio::SeqFeatures
    $self->{'_sequence'}  = undef; # location of Bio::Seq object
    $self->{'_workdir'}   = undef; # location of temp directory
    $self->{'_filename'}  = undef; # file to store Bio::Seq object
    $self->{'_results'}   = undef; # file to store results of Profile
    $self->{'_threshold'} = undef; # Value of the threshod
    $self->{'_protected'} = [];    # a list of files protected from 
    #deletion ???
    
    #print STDERR "@args\n";
    my( $query, $analysis, $parameters, $program, $database) = rearrange
      ([qw(QUERY
           ANALYSIS
           PARAMETERS
           PROGRAM
           DATABASE
          )], @args);
    
    
    $self->query ($query) if ($query);       
    $self->analysis ($analysis) if ($analysis);
    if(!$self->query){
      throw("need to have a query defined in order to run protein ".
            "annotation");
    }
    if(!$self->analysis){
      throw("need to have a analysis defined in order to run protein ".
            "annotation");
    }
    if(!$program){
      $program = $self->analysis->program;
    }
    $self->program($self->find_executable($program)) if($program);
    if(!$parameters){
      $parameters = $self->analysis->parameters;
    }
    $self->parameters($parameters);
    if(!$database){
      $database = $self->analysis->db_file;
    }
    $self->database($database);
    return $self; # success - we hope!
}

sub query{
    my ($self, $seq) = @_;
    
    if ($seq) {
      eval {
        $seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI")
	    };
      
      if (!$@) {
        $self->{'_sequence'} = $seq ;
      }else{
        if(-e $seq){
          $self->{'_sequence'} = $seq ;
        }else{
          throw("You must provide either a Bio::Seq or a file which ".
                " exists not $seq");
        }
      }
      my $filename = $self->get_tmp_file('/tmp/', 'annotation.'.$$, 
                                         'results');
      $self->results($filename);
    }
    return $self->{'_sequence'};
  }


sub program {
    my $self = shift;
    my $location = shift;
    if ($location) {
      unless (-e $location) {
        throw($self->program." not found at $location");
      }
      $self->{'_program'} = $location;
    }
    return $self->{'_program'};
}

sub analysis {
    my $self = shift;
    if (@_) {
        $self->{'_analysis'} = shift;
    }
    return $self->{'_analysis'};
} 

sub parameters{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'parameters'} = $value;
    }
    return $obj->{'parameters'};

}

sub database {
    my $self = shift;
    my $database = shift;
    if(!$self->{'_database'}){
      $self->{'_database'} = $database;
    }
    return $self->{'_database'};
} 


sub protein_files{
  my ($self) = @_;

  if(!$self->{'_protein_files'}){
    $self->{'_protein_files'} = {};
  }
  if($self->multiprotein){
    #print STDERR "Don't need to separate file\n";
    my $filename = $self->get_tmp_file('/tmp/', $self->analysis->logic_name.
                                       ".$$", 'seq');
    cp($self->query, $filename);
    $self->file($filename);
    $self->{'_protein_files'}{$filename} = 1;
  }else{
    my $in  = Bio::SeqIO->new(-file => $self->query, '-format' =>'Fasta');
    while(my $tmpseq = $in->next_seq()){
      #print STDERR "Producing file for ".$tmpseq->display_id."\n";
      my $stub = $self->analysis->logic_name.".".$tmpseq->display_id.".$$";
      my $filename = $self->get_tmp_file($self->workdir, $stub, "seq");
      $filename = $self->write_protein_file($filename, $tmpseq);
      $self->{'_protein_files'}{$filename} = $tmpseq->display_id;
      $self->file($filename);
    }
  }
  return $self->{'_protein_files'};
}

sub write_protein_file{
  my ($self, $filename, $seq) = @_;

  eval {
      my $out = Bio::SeqIO->new(-file => ">$filename", -format => 'fasta');
      $out->write_seq($seq);
      $out->close();
  };
  $@ and $self->throw("Could not write temporary sequence file $filename [$@]");

  return $filename;
}

sub multiprotein{
  my ($self) = @_;

  throw($self->program. "'s module must implement this method to define ".
        " if it can handle multi sequence fasta files or needs to take ".
        " one protein at once");
}

sub run{
  my ($self, $dir) = @_;

  $self->workdir ('/tmp') unless ($self->workdir($dir));
  $self->checkdir;

  if(-s $self->query){
    #the input is a sequence file 
    my %files = %{$self->protein_files};
    foreach my $file(keys(%files)){
      my $id = $files{$file};
      #print STDERR "Running with ".$file."\n";
      $self->filename($file);
      $self->run_analysis();
      $self->parse_results($id);
    }
  }else{
    eval {
      $self->query->isa ("Bio::PrimarySeqI") || 
        $self->query->isa ("Bio::SeqI")
    };
    if (!$@) {
      #The input is a sequence object
      # write sequence to file
      my $filename = $self->get_tmp_file('/tmp/', $self->analysis->logic_name.
                                         ".$$", 'seq');
      $self->filename($filename);
      $self->file($filename);
      $self->writefile;
      # run program
      $self->run_analysis;
      # parse output
      $self->parse_results;
    }else{
      throw("Can't run if ".$self->query." isn't either a Bio::Seq or a ".
            " file which has a size greater than 0");
    }
  }
  $self->deletefiles;
}


sub create_protein_feature{
  my ($self, $start, $end, $score, $seqname, $hstart, $hend, $hseqname,
     $analysis, $p_value, $percent_id) = @_;


  print STDERR "Create:have id ".$seqname."\n";
  my $fp = Bio::EnsEMBL::ProteinFeature->new(
                                             -start    => $start,
                                             -end      => $end,
                                             -hstart   => $hstart,
                                             -hend     => $hend,
                                             -percent_id => $percent_id,
                                             -score    => $score,
                                             -p_value  => $p_value,
                                             -hseqname => $hseqname,
                                             -seqname => $seqname,
                                             -analysis => $analysis,
                                            );

  return $fp;
}



sub add_to_output{
  my ($self, $f) = @_;
  if($f){
    #print STDERR "Adding $f to list\n";
    push(@{$self->{'_flist'}}, $f);
  }
  #print STDERR "have ". @{$self->{'_flist'}}." results\n";
}

sub output {
    my ($self) = @_;
    return @{$self->{'_flist'}};
}

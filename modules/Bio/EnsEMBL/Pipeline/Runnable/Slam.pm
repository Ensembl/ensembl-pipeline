#
# Written by Jan-Hinnerk Vogel
#
# Copyright GRL/EBI 2004
#
# You may distribute this module under the same terms as perl itself

## run
## output


package Bio::EnsEMBL::Pipeline::Runnable::Slam;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;



# Inherits from Interface Bio::EnsEMBL::Pipeline::RunnableI
 @ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
#  my ($class,@args) = @_;
#  my $self = $class->SUPER::new(@args);

  my ($class,@args) = @_;
  my $self={};  #construct empty hash
  bless $self,$class;


  my (
      $first_slice,
      $first_slice_name,
      $second_slice,
      $second_slice_name,
      $slam_bin,
      $options
      ) = $self->_rearrange([qw(
				FIRST_SLICE
				FIRST_SLICE_NAME
				SECOND_SLICE
				SECOND_SLICE_NAME				
				SLAM_BIN
				OPTIONS
			       )
			    ], @args
			   );

      $self->{_output} = [];


  ################################################################
  # Setting $first_slice

  if( $first_slice ){
    unless ($first_slice->isa("Bio::PrimarySeqI") ) {
      $self->throw("Query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->set_first_slice( $first_slice );
  }
  else{
    $self->throw("Slam needs a first  slice-object: $first_slice ");
  }


  ################################################################
  # Setting $second_slice

  if( $second_slice ){
    unless ($second_slice->isa("Bio::PrimarySeqI") ) {
      $self->throw("Query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->set_second_slice( $second_slice );
  }
  else{
    $self->throw("Slam needs a second  slice-object: $second_slice ");
  }

  ################################################################
  # Name of $first_slice

  if ( $first_slice_name) {
    $self->setSliceName("first_slice",$first_slice_name);
  }else {
    $self->setSliceName("first_slice");
  }

  ################################################################
  # Name of $second_slice
  if ( $second_slice_name) {
    $self->setSliceName("second_slice",$second_slice_name);
  }else {
    $self->setSliceName("second_slice");
  }

  ################################################################
  # Path to used slam-binary if not provided by call

  if ($slam_bin) {
      $self->set_slam_bin($slam_bin);
  } else {
      $self->set_slam_bin('/nfs/acari/jhv/slam/prog/slam.pl');
  }

  #################################################################
  # Passing additional options to basic options for calling slam.pl

  my $basic_options = " ";
  if ($options){    $basic_options .= " ".$options; }
  $self->options( $basic_options);

  return $self;
}



sub DESTROY {
  my $self = shift;
#  unlink $self->_query_file;
}



####runn
####################################################################################runn
################################################################################ 

sub run {
  my ($self) = @_;

  # name of results file
  $self->results($self->workdir . "/results.$$");

  $self->_write_first_seq;
  $self->_write_second_seq;

    my $command =  $self->{_slam_bin} .
                    " ".$self->{_first_filename} .
		      " ".$self->{_second_filename} .
			" ".$self->options;

  print STDERR "Exonerate command : $command\n"     if $self->_verbose;


  $command =  $self->{_slam_bin} ." query_first_slice.11548 query_second_slice.11548";


  open( EXO, "$command |" ) || $self->throw("Error running Slam $!");
#  $self->parse_results(\*EXO);
  close(EXO);





  return 1
}






############################################################
# get/set methods  setti
############################################################

sub set_first_slice {
  my ($self,$fslice) = @_;
  if ($fslice) {
    $self->{_first_slice} = $fslice;
  }
}

sub set_second_slice {
  my ($self,$sslice) = @_;
  if ($sslice) {
    $self->{_second_slice} = $sslice;
  }
}

sub set_slam_bin {
  my ($self, $location) = @_;
  if ($location) {
    $self->throw("Slam not found at $location: $!\n") unless (-e $location);
    $self->{_slam_bin} = $location ;
  }
  return $self->{_slam_bin};
}


sub _write_first_seq {
  my ($self)     = shift;
  my $name       = $self->{_first_slice_name};
  my $slice      = $self->{_first_slice};
  my $filename   = $self->workdir."/query_" . "$name.$$";
     $self->{_first_filename}=$filename;
  my $seqout     = Bio::SeqIO->new('-file'   => ">$filename" , '-format' => 'Fasta' );
     $seqout->write_seq($slice);
}

sub _write_second_seq {
  my ($self)     = shift;
  my $name       = $self->{_second_slice_name};
  my $slice      = $self->{_second_slice};
  my $filename   = $self->workdir."/query_" . "$name.$$";
     $self->{_second_filename}=$filename;
  my $seqout     = Bio::SeqIO->new('-file'   => ">$filename" , '-format' => 'Fasta' );
     $seqout->write_seq($slice);
}

sub output {
  my ($self, @output) = @_;

  if (@output) {
    unless( $self->{_output} ){
      $self->{_output} = [];
    }
    push( @{$self->{_output}}, @output );
  }
  return @{$self->{_output}};
}






############################################################





sub _verbose {
  my $self = shift;
  if (@_){
    $self->{_verbose} = shift;
  }
  return $self->{_verbose}
}


sub setSliceName {
  my ($self,$target,$name) = @_;
  if ($target eq "first_slice"){
    if (!$name){
      $self->{_first_slice_name}="first_slice";
    }else{
      $self->{_first_slice_name}=$name;
    }
  }elsif ($target eq "second_slice") {
    if (!$name){
      $self->{_second_slice_name}="second_slice";
    }else{
      $self->{_second_slice_name}=$name;
    }
  }
}

sub getSliceNames {
  my $self = shift;
  my @sliceNames = ($self->_second_slice_name,$self->_first_slice_name);
  return @sliceNames;
}

sub pH {
  my $self = shift;
  my %tmp = %{$self};

  print "Key:\t\t\tValue:\n";
  print "------------------------------\n";
  
  foreach (keys %tmp) {
    print "-$_-\t\t\t-$tmp{$_}-\n";
  }
}

1;


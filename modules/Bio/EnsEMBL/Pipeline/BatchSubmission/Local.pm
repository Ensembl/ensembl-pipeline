package Bio::EnsEMBL::Pipeline::BatchSubmission::Local;


use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Pipeline::BatchSubmission;
use vars qw(@ISA);
use strict;
use File::Copy;

@ISA = qw(Bio::EnsEMBL::Pipeline::BatchSubmission);


sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'commandline'} = undef;
  return $self;
}





sub construct_command_line{
  my($self, $command) = @_; 
  
  if(!$command){
    throw("Can't create command line if not command is passed in");
  }
  $self->command($command);
  
  my $command_line = $command;
  $command .= " 1> ".$self->temp_stdout." ";
  $command .= "2> ".$self->temp_stderr." ";
  $self->command($command);
}


sub open_command_line{
  my ($self, $verbose) = @_;
  if($verbose){
    print STDERR "Opening ".$self->command."\n";
  }
  eval{
    system($self->command);
  };
  if($@){
    throw("FAILED to open commandline locally $@");
  }
  eval{
    if($self->stderr_file){
      $self->copy_output($self->temp_stderr, $self->stderr_file);
    }
    if($self->stdout_file){
      $self->copy_output($self->temp_stdout, $self->stdout_file);
    }
  };
  if($@){
    throw("FAILED to copy output to defined location $@ keeping temporary".
          "output")
  }else{
    $self->delete_temp_output;
  }
}


sub copy_output{
  my ($self, $file_to_copy, $new_location) = @_;
  eval{
    copy($file_to_copy, $new_location);
  };
  if($@){
    throw("Copy of temporary output ".$file_to_copy." to new location ".
          $new_location." failed $@");
  }
  return $new_location;
}

sub delete_temp_output{
  my ($self) = @_;
  if(-e $self->temp_stderr){
    unlink $self->temp_stderr;
  }
  if(-e $self->temp_stdout){
    unlink $self->temp_stdout;
  }
}


sub temp_stderr{
  my ($self, $file) = @_;
  if($file){
    $self->{'temp_stderr'} = $file;
  }
  if(!$self->{'temp_stderr'}){
    my $num = int(rand(100000));
    $self->{'temp_stderr'} = "/tmp/temp_local_out.".$num.".$$.err";
  }
  return $self->{'temp_stderr'};
}

sub temp_stdout{
  my ($self, $file) = @_;
  if($file){
    $self->{'temp_stdout'} = $file;
  }
  if(!$self->{'temp_stdout'}){
    my $num = int(rand(100000));
    $self->{'temp_stdout'} = "/tmp/temp_local_out.".$num.".$$.out";
  }
  return $self->{'temp_stdout'};
}


sub job_stats{
  my ($self) = @_;
  my $hash = {};
  $hash->{1} = 1;
  return $hash;
}

sub temp_filename{
  my ($self) = @_;
  return '/tmp/';
}

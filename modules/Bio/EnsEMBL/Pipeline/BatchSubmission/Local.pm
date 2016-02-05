=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Pipeline::BatchSubmission::Local - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Pipeline::BatchSubmission::Local;


use warnings ;
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
# This is not the real job_id of the command but it's better than nothing...
  $self->id($$);
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

sub memstring_to_resource {
    return '';
}

1;

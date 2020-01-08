=head1 LICENSE


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

=head1 AUTHORS

Laura Clarke 

=head1 NAME

Bio::EnsEMBL::Pipeline::BatchSubmission - 

=head1 SYNOPSIS

This is just a base class: objects should be created via child classes
Bio::EnsEMBL::Pipeline::BatchSubmission::*

my $batchjob = Bio::EnsEMBL::Pipeline::BatchSubmission::XXX->new(
             -STDOUT     => $stdout_file,
             -STDERR     => $stderr_file,
             -PARAMETERS => @args,
             -PRE_EXEC   => $pre_exec,
             -QUEUE      => $queue,
             -JOBNAME    => $jobname,
             -NODES      => $nodes,
             -RESOURCE   => $resource,
             -MEMORY     => $memory_string
             );

$batch_job->construct_command_line('test.pl');
$batch_job->open_command_line();


=head1 DESCRIPTION

Generic base class for specific BatchSubmission modules found in 
BatchSubmission directory for handling pipeline jobs in a distributed 
environment with load sharing software such as LSF, PBS,etc.

All get/sets and generic methods found here, while specific methods such as
construct_command_line have to be implemented in the specific child 
classes Bio::EnsEMBL::Pipeline::BatchSubmission::*

the specific methods which must be implemented are

construct_command_line, a method to build the submission statement
open_command_line, a method to open the submission statement 
copy_output ,a method copy the output of a job at to the defined 
             destination if desired
delete_output, a method to delete the output from the temporary location
               that held it

If the MEMORY option is to be used, memstring_to_resource() needs
to be implmented (see LSF.pm for an example).  It is supposed to
allow the user to specify a memory requirement using a simple syntax
rather than using a more complicated resource specification.  The
memstring_to_resource() method does the work of translating the simple
syntax into a resource requirement specification.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::BatchSubmission;


use warnings ;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw();

sub new {
  my ( $class, @args ) = @_;

  my $self = bless {}, $class;

  $self->{'stdout'}     = undef;
  $self->{'stderr'}     = undef;
  $self->{'id'}         = undef;
  $self->{'parameters'} = undef;
  $self->{'pre_exec'}   = undef;
  $self->{'command'}    = undef;
  $self->{'queue'}      = undef;
  $self->{'jobname'}    = undef;
  $self->{'nodes'}    = undef;  #must by a space delimited line of nodes
  $self->{'resource'} = undef;

  my ( $stdout, $stderr,  $parameters, $pre_exec, $command,
       $queue,  $jobname, $resource,   $nodes,    $memstring )
    = rearrange( [ 'STDOUT',  'STDERR', 'PARAMETERS', 'PRE_EXEC',
                   'COMMAND', 'QUEUE',  'JOBNAME',    'RESOURCE',
                   'NODES',   'MEMORY',
                 ],
                 @args );

  if ($stdout)     { $self->stdout_file($stdout) }
  if ($stderr)     { $self->stderr_file($stderr) }
  if ($parameters) { $self->parameters($parameters) }
  if ($pre_exec)   { $self->pre_exec($pre_exec) }
  if ($command)    { $self->command($command) }
  if ($queue)      { $self->queue($queue) }
  if ($jobname)    { $self->jobname($jobname) }
  if ($resource)   { $self->resource($resource) }
  if ($nodes)      { $self->nodes($nodes) }
  if ($memstring)  { $self->memstring_to_resource($memstring) }

  return $self;
} ## end sub new

# The memory_conversion() utility method converts a string like
# "350MB" into an integer representing a memory size of a given
# unit, e.g. "0.35" (if $unit is 'GB'), "350" (if $unit is 'MB') or
# "350000" (if $unit is 'KB').  Supported units are 'B' (bytes), 'KB'
# (kilo-bytes), 'MB', (mega-bytes), and 'GB' (giga-bytes).  This method
# is used by the memstring_to_resource() method.
sub memory_conversion {
  my ( $self, $memstring, $unit ) = @_;

  # This is using 1000 bytes per kilobytes etc., not 1024.
  my %byte_multiplier =
    ( 'B' => 1, 'KB' => 1000, 'MB' => 1000000, 'GB' => 1000000000 );

  my ( $number, $suffix ) = ( $memstring =~ /([0-9.]+)\s*([KMG]?B)/i );

  if ( !exists( $byte_multiplier{ uc($suffix) } ) ) {
    die("Unknown unit in memory string: $suffix");
  }
  elsif ( !exists( $byte_multiplier{ uc($unit) } ) ) {
    die("Unknown unit in argument: $unit");
  }

  # Convert from whatever unit is in the supplied string to bytes.
  my $bytes = $number*$byte_multiplier{ uc($suffix) };

  # Convert from bytes to whatever unit the user asked for.
  return $bytes/$byte_multiplier{ uc($unit) };
}


##################
#accessor methods#
##################

sub stdout_file{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'stdout'} = $arg;
   }

   return $self->{'stdout'};
}



sub stderr_file{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'stderr'} = $arg;
   }

   return $self->{'stderr'};
}


sub id{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'id'} = $arg;
   }

   return $self->{'id'};
}

sub parameters{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'parameters'} = $arg;
   }

   return $self->{'parameters'};
}

sub pre_exec{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'pre_exec'} = $arg;
   }

   return $self->{'pre_exec'};
}

sub command{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'command'} = $arg;
   }

   return $self->{'command'};
}


sub queue{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'queue'} = $arg;
   }

   return $self->{'queue'};
}

sub jobname{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'jobname'} = $arg;
   }

   return $self->{'jobname'};
}

sub nodes{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'nodes'} = $arg;
   }

   return $self->{'nodes'};
}

sub resource{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'resource'} = $arg;
   }

   return $self->{'resource'};
}

sub temp_filename {
   my ($self, $arg) = @_;

   if($arg){
     $self->{'tmp_jobfilename'} = $arg;
   }

   return $self->{'tmp_jobfilename'};
}

#############
#run methods#
#############

sub construct_command_line{
  my ($self) = @_;

  throw("this method construct command line must be implemented");
}


sub open_command_line{
  my ($self)= @_;
  
  throw("open_command_line isn't implemented yet");

}

sub copy_output{
  my ($self) = @_;

  throw("this method copy_output must be implemented");
}

sub delete_output{
  my ($self) = @_;

  throw("this method delete_output must be implemented");
}

sub job_stats {
  my ($self) = @_;

  throw("this method job_stats must be implemented");
}


###############
#check methods#
###############

sub is_db_overloaded {
    my ($self) = @_;

    warning("If you need to take care of the load of your DB, is_db_overloaded must be implemented");
    return 0;
}

1;

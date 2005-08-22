# Base Class for handling job submission via Load Sharing software
#
# Cared for by Laura Clarke 
#
# Copyright Laura Clarke 
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::BatchSubmission

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
             -RESOURCE   => $resource
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

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::BatchSubmission;


use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw();

sub new{
  my ($class, @args) = @_;
  my $self = bless {},$class;   
  $self->{'stdout'} = undef;
  $self->{'stderr'} = undef;
  $self->{'id'} = undef;
  $self->{'parameters'} = undef;
  $self->{'pre_exec'} = undef;
  $self->{'command'} = undef;
  $self->{'queue'} = undef;
  $self->{'jobname'} = undef;
  $self->{'nodes'} = undef; #must by a space delimited line of nodes
  $self->{'resource'} = undef;

  my($stdout, $stderr, $parameters, $pre_exec, $command,$queue, $jobname,
     $resource, $nodes) = rearrange([qw(STDOUT
                                        STDERR
                                        PARAMETERS
                                        PRE_EXEC 
                                        COMMAND 
                                        QUEUE 
                                        JOBNAME
                                        RESOURCE
                                        NODES)],@args);
  
  if($stdout){
    $self->stdout_file($stdout);
  }
  if($stderr){
    $self->stderr_file($stderr);
  }
  if($parameters){
    $self->parameters($parameters);
  }
  if($pre_exec){
    $self->pre_exec($pre_exec);
  }
  if($command){
    $self->command($command);
  }
  if($queue){
    $self->queue($queue);
  }
  if($jobname){
    $self->jobname($jobname);
  }
  if($resource){
    $self->resource($resource);
  }
  if($nodes){
    $self->nodes($nodes);
  }

  return $self;
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

1;

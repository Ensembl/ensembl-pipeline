package Bio::EnsEMBL::Pipeline::BatchSubmission;

BEGIN {
  require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Root);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);    
  $self->{'stdout'} = undef;
  $self->{'stderr'} = undef;
  $self->{'id'} = undef;
  $self->{'parameters'} = undef;
  $self->{'pre_exec'} = undef;
  $self->{'command'} = undef;
  $self->{'queue'} = undef;
  $self->{'jobname'} = undef;
  $self->{'nodes'} = undef; #must by a space delimited line of nodes

  my($stdout, 
     $stderr, 
     $parameters, 
     $pre_exec, 
     $command, 
     $queue, 
     $jobname) = $self->_rearrange([qw(STDOUT 
				       STDERR 
				       PARAMETERS 
				       PRE_EXEC 
				       COMMAND 
				       QUEUE 
				       JOBNAME
				       NODES)],@args);

  if(defined($stdout)){
    $self->stdout_file($stdout);
  }
  if(defined($stderr)){
    $self->stderr_file($stderr);
  }
  if(defined($parameters)){
    $self->parameters($parameters);
  }
  if(defined($pre_exec)){
    $self->pre_exec($pre_exec);
  }
  if(defined($command)){
    $self->command($command);
  }
  if(defined($queue)){
    $self->queue($queue);
  }
  if(defined($jobname)){
    $self->jobname($jobname);
  }

  return $self;
}

##################
#accessor methods#
##################

sub stdout_file{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'stdout'} = $arg;
   }

   return $self->{'stdout'};
}



sub stderr_file{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'stderr'} = $arg;
   }

   return $self->{'stderr'};
}


sub id{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'id'} = $arg;
   }

   return $self->{'id'};
}

sub parameters{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'parameters'} = $arg;
   }

   return $self->{'parameters'};
}

sub pre_exec{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'pre_exec'} = $arg;
   }

   return $self->{'pre_exec'};
}

sub command{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'command'} = $arg;
   }

   return $self->{'command'};
}


sub queue{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'queue'} = $arg;
   }

   return $self->{'queue'};
}

sub jobname{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'jobname'} = $arg;
   }

   return $self->{'jobname'};
}

sub nodes{
   my ($self, $arg) = @_;

   if(defined($arg)){
     $self->{'nodes'} = $arg;
   }

   return $self->{'nodes'};
}


#############
#run methods#
#############

sub construct_command_line{
  my($self) = @_;

  $self->throw("this method construct command line must be implemented!\n");
}


sub open_command_line{
  my ($self)= @_;
  
  $self->throw("open_command_line isn't implemented yet\n");

}

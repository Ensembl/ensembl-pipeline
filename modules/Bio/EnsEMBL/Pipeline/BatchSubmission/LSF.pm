# Ensembl Pipeline module for handling job submission via Platform LSF 
# load sharing software
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

Bio::EnsEMBL::Pipeline::BatchSubmission::LSF

=head1 SYNOPSIS

my $batchjob = Bio::EnsEMBL::Pipeline::BatchSubmission::LSF->new(
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

This module provides an interface to the Platform LSF load sharing software and
its commands. It implements the method construct_command_line which is not 
defined in the base class and which enables the pipeline to submit jobs in a 
distributed environment using LSF.

See base class Bio::EnsEMBL::Pipeline::BatchSubmission for more info

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::BatchSubmission::LSF;


use Bio::EnsEMBL::Pipeline::BatchSubmission;
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Pipeline::BatchSubmission);


sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'bsub'} = undef;
  
  return $self;
 
}


##################
#accessor methods#
##################

sub bsub{
  my($self, $arg) = @_;

  if(defined($arg)){
    $self->{'bsub'} = $arg;
  }

  return $self->{'bsub'};

}

##other accessor are in base class##

######################
#command line methods#
######################

sub construct_command_line{
  my($self, $command, $stdout, $stderr) = @_; 
#command must be the first argument then if stdout or stderr aren't definhed the objects own can be used
  
  if(!$command){
    $self->throw("cannot create bsub if nothing to submit to it : $!\n");
  }
  my $bsub_line;
  $self->command($command);
  if($stdout){
    $bsub_line = "bsub -o ".$stdout;
  }else{
    $bsub_line = "bsub -o ".$self->stdout_file;
  }
  if($self->nodes){
    my $nodes = $self->nodes;
    # $nodes needs to be a space-delimited list
    $nodes =~ s/,/ /;
    $nodes =~ s/ +/ /;
    # undef $nodes unless $nodes =~ m{(\w+\ )*\w};
    $bsub_line .= " -m '".$nodes."' ";
  }
  if(my $res = $self->resource){
    $res = qq{-R '$res'} unless $res =~ /^-R/;
    $bsub_line .= " $res ";
  }
  $bsub_line .= " -q ".$self->queue    if $self->queue;
  $bsub_line .= " -J ".$self->jobname  if $self->jobname;
  $bsub_line .= " ".$self->parameters." "  if defined $self->parameters;
  if($stderr){
    $bsub_line .= " -e ".$stderr;
  }else{
    $bsub_line .= " -e ".$self->stderr_file;
  }
  $bsub_line .= " -E \"".$self->pre_exec."\"" if defined $self->pre_exec; 
  ## must ensure the prexec is in quotes ##
  $bsub_line .= " ".$command;
  $self->bsub($bsub_line);
  
}



sub open_command_line{
  my ($self, $verbose)= @_;

  open(SUB, $self->bsub." 2>&1 |");
  my $lsf;
  while(<SUB>){
    print STDERR if($verbose);
    if (/Job <(\d+)>/) {
      $lsf = $1;
    }
  }
  $self->id($lsf);
  close(SUB);
}


sub get_pending_jobs {
  my($self, %args) = @_;

  my ($user)  = $args{'-user'}  || $args{'-USER'}  || undef;
  my ($queue) = $args{'-queue'} || $args{'-QUEUE'} || undef;

  my $cmd = "bjobs -p";
  $cmd   .= " -q $queue" if $queue;
  $cmd   .= " -u $user"  if $user;
  $cmd   .= " | grep PEND";

  open CMD, "$cmd 2>&1 |" or do {
    return undef;
  };

  my @lines = <CMD>;

  close CMD or do {
    return undef;
  };

  return 0 if $lines[0] =~ /No pending job found/;
  return scalar @lines;
}


sub get_job_time{
  my ($self) = @_;
  my $command = "bjobs -l";
  #print $command."\n";
  my %id_times;
  open(BJOB, "$command |") or $self->throw("couldn't open pipe to bjobs");
  my $job_id;
  while(<BJOB>){
    chomp;
    if(/Job\s+\<(\d+)\>/){
      $job_id = $1;
    }elsif(/The CPU time used/){
      my ($time) = $_ =~ /The CPU time used is (\d+)/;
      $id_times{$job_id} = $time;
    }
  }
  close(BJOB);
  #or $self->throw("couldn't close pipe to bjobs");
  return \%id_times;
}



sub check_existance{
  my ($self, $id, $verbose) = @_;
  if(!$id){
    die("Can't run without an LSF id");
  }
  my $command = "bjobs ".$id."\n";
  #print STDERR "Running ".$command."\n";
  my $flag = 0; 
  open(BJOB, "$command 2>&1 |") or $self->throw("couldn't open pipe to bjobs");
  while(<BJOB>){
    print STDERR if($verbose);
    chomp;
    if ($_ =~ /No unfinished job found/) {
      #print "Set flag\n";
      $flag = 1;
    } 
    my @values = split;
    if($values[0] =~ /\d+/){
      return $values[0];
    }
  }
  close(BJOB);
  print STDERR "Have lost ".$id."\n" if($verbose);
  return undef;
}


sub kill_job{
  my ($self, $job_id) = @_;

  my $command = "bkill ".$job_id;
  system($command);
}

sub stdout_file{
   my ($self, $arg) = @_;

   if($arg){
     $self->{'stdout'} = $arg;
   }

   if(!$self->{'stdout'}){
     $self->{'stdout'} ='/dev/zero'
   }
   return $self->{'stdout'};
}



sub stderr_file{
   my ($self, $arg) = @_;

   if ($arg){
     $self->{'stderr'} = $arg;
   }
   if(!$self->{'stderr'}){
     $self->{'stderr'} ='/dev/zero'
   }
   return $self->{'stderr'};
}



sub temp_filename{
  my ($self) = @_;

  $self->{'lsf_jobfilename'} = $ENV{'LSB_JOBFILENAME'};
  return $self->{'lsf_jobfilename'};
}


sub temp_outfile{
  my ($self) = @_;

  $self->{'_temp_outfile'} = $self->temp_filename.".out";

  return $self->{'_temp_outfile'};
}

sub temp_errfile{
  my ($self) = @_;

  $self->{'_temp_errfile'} = $self->temp_filename.".err";
  

  return $self->{'_temp_errfile'};
}


sub submission_host{
  my ($self) = @_;

  $self->{'_submission_host'} = $ENV{'LSB_SUB_HOST'};
  

  return $self->{'_submission_host'};
}

sub lsf_user{
  my ($self) = @_;

 
  $self->{'_lsf_user'} = $ENV{'LSFUSER'};
  

  return $self->{'_lsf_user'};
}


sub copy_output{
  my ($self, $stderr_file, $stdout_file) = @_;

  $stderr_file = $self->stderr_file if(!$stderr_file);
  $stdout_file = $self->stdout_file if(!$stdout_file);
  my $err_file = $self->temp_errfile;
  my $out_file  = $self->temp_outfile;

  if(!$self->temp_filename){
    my ($p, $f, $l) = caller;
    $self->warn("The lsf environment variable LSB_JOBFILENAME is not defined".
                " we can't copy the output files which don't exist $f:$l");
    return;
  }
  my $command = $self->copy_command;
  if(-e $err_file){
    
    my $err_copy = $command." ".$err_file." ".$self->lsf_user."@".$self->submission_host.":".$stderr_file." 2>&1 ";
   

    if(system($err_copy)){
      $self->throw("Couldn't execute ".$err_copy);
    }
  }
  if(-e $out_file){
    my $out_copy = $command." ".$out_file." ".$self->lsf_user."@".$self->submission_host.":".$stdout_file." 2>&1";
   
    if(system($out_copy)){
      $self->throw("Couldn't execute ".$out_copy);
    }
   
  }

}


sub delete_output{
  my ($self) = @_;
  
  unlink $self->temp_errfile if(-e $self->temp_errfile);
  unlink $self->temp_outfile if(-e $self->temp_outfile);
}

sub copy_command{
  my ($self, $arg) = @_;

  if($arg){
    $self->{'_copy_command'} = $arg;
  }

  return $self->{'_copy_command'} || 'lsrcp';
}


1;

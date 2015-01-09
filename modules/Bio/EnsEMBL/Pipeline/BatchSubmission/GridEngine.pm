=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Steve Searle

=head1 NAME

Bio::EnsEMBL::Pipeline::BatchSubmission::GridEngine - 

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

This module provides an interface to the Sun Grid Engine load sharing softwa
re and its commands. It implements the method construct_command_line which is 
not defined in the base class and which enables the pipeline to submit jobs 
in a distributed environment using Sun Grid Engine.

See base class Bio::EnsEMBL::Pipeline::BatchSubmission for more info


=head1 APPENDIX

The rest of the documentation details each of the object methods. Interna
l
methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Pipeline::BatchSubmission::GridEngine;

use warnings ;
use Bio::EnsEMBL::Pipeline::BatchSubmission;
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Pipeline::BatchSubmission);


sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  return $self;
 
}


######################
#command line methods#
######################

sub construct_command_line{
  my($self, $command, $stdout, $stderr) = @_; 
  #print STDERR "creating the command line\n";
#command must be the first argument then if stdout or stderr aren't definhed the objects own can be used
  if(!$command){
    $self->throw("cannot create qsub if nothing to submit to it : $!\n");
  }
  my $qsub_line;
  $self->command($command);
  if($stdout){
    $qsub_line = "qsub -V -cwd -v FINAL_STDOUT=".$stdout;
  }else{
    $qsub_line = "qsub -V -cwd -v FINAL_STDOUT=".$self->stdout_file;
  }
  if($stderr){
    $qsub_line .= " -v FINAL_STDERR=".$stderr;
  }else{
    $qsub_line .= " -v FINAL_STDERR=".$self->stderr_file;
  }
  $qsub_line .= " -o /tmp -e /tmp";

# Depends on queues being made for each node with name node.q
  if($self->nodes){
    my $nodes = $self->nodes;
    # $nodes needs to be a space-delimited list
    $nodes =~ s/,/.q,/;
    $qsub_line .= " -q ".$nodes." ";
  }

  if (defined($self->queue) && $self->queue ne "") {$qsub_line .= " -l ".$self->queue;}

  $qsub_line .= " -N ".$self->jobname  if defined $self->jobname;

  $qsub_line .= " ".$self->parameters." "  if defined $self->parameters;

  $qsub_line .= " -v PREEXEC=\"".$self->pre_exec."\"" if defined $self->pre_exec; 

  ## must ensure the prexec is in quotes ##
  my $ge_wrapper = "ge_wrapper.pl";
  unless (-x $ge_wrapper) {
    $ge_wrapper = __FILE__;
    $ge_wrapper =~ s:GridEngine.pm:../../../../../scripts/ge_wrapper.pl:;
    print $ge_wrapper . "\n";
    my $caller = caller(0);
    $self->throw("ge_wrapper not found - needs to be set in $caller\n") unless -x $ge_wrapper;
  }

  $qsub_line .= " $ge_wrapper \"".$command . "\"";
  $self->command($qsub_line);
  #print "have command line\n";
}



sub open_command_line{
  my ($self)= @_;

  print STDERR $self->command."\n";
  print STDERR "opening command line\n";
  open(SUB, $self->command." 2>&1 |");
  my $geid;
  while(<SUB>){
    if (/your job (\d+)/) {
      $geid = $1;
    }
  }
  print STDERR "have opened ".$self->command."\n";
  print STDERR "geid ".$geid."\n";
  $self->id($geid);
  close(SUB);
}

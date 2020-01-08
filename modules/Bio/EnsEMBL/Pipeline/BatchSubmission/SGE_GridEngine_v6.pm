
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
# Ensembl Pipeline module for handling job submission via Sun Grid Engine
# load sharing software
#
# Cared for by Steve Searle
#
# Copyright Steve Searle
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::BatchSubmission::GridEngine

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

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Interna
l
methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::BatchSubmission::SGE_GridEngine_v6;

use warnings;
use Bio::EnsEMBL::Pipeline::BatchSubmission;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Pipeline::Config::General;
use vars qw(@ISA);

# command out by GT
#use strict;

@ISA = qw(Bio::EnsEMBL::Pipeline::BatchSubmission);

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    return $self;

}

##################
#accessor methods#
##################

sub qsub {
    my ( $self, $qsub_line ) = @_;

    if ( defined($qsub_line) ) {
        $self->{_qsub} = $qsub_line;
    }
    return $self->{_qsub};
}

##other accessor are in base class##

######################
#command line methods#
######################

sub construct_command_line {
    my ( $self, $command, $stdout, $stderr ) = @_;

    #print STDERR "creating the command line\n";
    #command must be the first argument then if stdout or stderr aren't definhed the objects own can be used
    if ( !$command ) {
        $self->throw("cannot create qsub if nothing to submit to it : $!\n");
    }

    my $qsub_line = "";

    $test = $self->stdout_file if defined $self->stdout_file;

    $test =~ s/\:/\_/g;    #SGE does not like : in job name so change to _

    # this section writes a little wrapper script for the job itself which will be executed by
    # the sun grid engine.

    open( QSUB, '>', $test . '.sh' );

    #print QSUB "#!/usr/bin/tcsh\n";
    print QSUB "#!/bin/bash\n";

    if ( !defined $SGE_PERL5LIB_ENV_SCRIPT || length($SGE_PERL5LIB_ENV_SCRIPT) == 0 ) {
        throw(  "Missing config variable : SGE_PERL5LIB_ENV_SCRIPT !!! \n"
              . "You have to define a setup script in Pipeline/Config/General.pm to set up your PERL5LIB env" );
    }
    if ( -e $SGE_PERL5LIB_ENV_SCRIPT && -r $SGE_PERL5LIB_ENV_SCRIPT ) {
        print "SGE_GridEngine_v6.pm : using setup script $SGE_PERL5LIB_ENV_SCRIPT\n";
    }
    else {
        throw(
"Your setup script SGE_PERL5LIB_ENV_SCRIPT specified in Pipeline/Config/General.pm is not readable or does not exist\n"
        );
    }
    my $env_setup_script = $SGE_PERL5LIB_ENV_SCRIPT;
    print QSUB "source $env_setup_script\n";

    $self->command($command);

    if ($stdout) {

        #$qsub_line = "qsub -V -cwd -v FINAL_STDOUT=".$stdout;
        # add by GT
        #print QSUB "setenv FINAL_STDOUT $stdout\n";
        print QSUB "export FINAL_STDOUT=$stdout\n";
        $stdout =~ s/:+/_/g;    # outout file can't include ::
        $qsub_line .= " -o $stdout ";    # jhv
    }
    else {

        #$qsub_line = "qsub -V -cwd -v FINAL_STDOUT=".$self->stdout_file;
        $ffstdout = $self->stdout_file;
        $ffstdout =~ s/:+/_/g;           # outout file can't include ::
                                         #print QSUB "setenv FINAL_STDOUT $ffstdout\n";
        print QSUB "export FINAL_STDOUT=$ffstdout\n";
        print "XXX SETTING STDOUT $ffstdout\n";
        $qsub_line .= " -o $ffstdout ";    # jhv
    }

    if ($stderr) {

        #$qsub_line .= " -v FINAL_STDERR=".$stderr;
        #add by GT
        #print QSUB "setenv FINAL_STDERR $stderr\n";
        print QSUB "export FINAL_STDERR=$stderr\n";
        $stderr =~ s/:+/_/g;               # outout file can't include ::
        $qsub_line .= " -e $stderr ";      # jhv
    }
    else {
        $ffstderr = $self->stderr_file;
        $ffstderr =~ s/:+/_/g;             # outout file can't include ::
                                           # $qsub_line .= " -v FINAL_STDERR=".$self->stderr_file;
                                           # add by GT
                                           #$stderrprint = .$self->stderr_file;
                                           #print QSUB "setenv FINAL_STDERR $ffstdout.err\n";
        print QSUB "export FINAL_STDERR=$ffstdout.err\n";
        print "XXX SETTING STDERR $ffstdout \n";
        $qsub_line .= " -e $ffstderr ";    # jhv
    }

    #$qsub_line .= " -o /tmp -e /tmp";

    #add by GT
    #print QSUB "setenv SGE_STDOUT_PATH /tmp\n";
    print QSUB "export SGE_STDOUT_PATH=/tmp\n";

    #print QSUB "setenv SGE_STDERR_PATH /tmp\n";
    print QSUB "export SGE_STDERR_PATH=/tmp\n";

    # Depends on queues being made for each node with name node.q
    if ( $self->nodes ) {
        my $nodes = $self->nodes;

        # $nodes needs to be a space-delimited list
        $nodes =~ s/,/.q,/;

        #$qsub_line .= " -q ".$nodes." ";
        #print "qsub _line 1 is  $qsub_line GridEngine.pm\n";
        print " nodes = $nodes \n";
    }

    # remove by GT
    #  if (defined($self->queue) && $self->queue ne "") {$qsub_line .= " -l ".$self->queue;}

    #$qsub_line .= " -N ".$self->jobname  if defined $self->jobname;

    #print "qsub _line 2 is  $qsub_line in GridEngine.pm\n";
    #print "self = $self \n";

    #$qsub_line .= " ".$self->parameters." "  if defined $self->parameters;
    #print "qsub _line 3 is  $qsub_line in GridEngine.pm\n";

    #$qsub_line .= " -v PREEXEC=\"".$self->pre_exec."\"" if defined $self->pre_exec;

    #print "qsub _line 4 is  $qsub_line in GridEngine.pm\n";

    #add by GT
    $ffpreexec1 = $self->pre_exec;

    # print QSUB "/usr/local/ensembl/bin/"."$ffpreexec1\n";
    print QSUB "$ffpreexec1\n";

    ## must ensure the prexec is in quotes ##
    my $ge_wrapper = "ge_wrapper.pl";
    unless ( -x $ge_wrapper ) {
        $ge_wrapper = __FILE__;

        #$ge_wrapper =~ s:GridEngine.pm:../../../../../scripts/ge_wrapper.pl:;
        $ge_wrapper =~ s:SGE_GridEngine_v6.pm:../../../../../scripts/ge_wrapper.pl:;
        print "Setting ge_wrapper to [$ge_wrapper]\n";
        my $caller = caller(0);
        $self->throw("ge_wrapper not found - needs to be set in $caller\n") unless -x $ge_wrapper;
    }

    #$qsub_line .= " $ge_wrapper \"".$command . "\"";
    #  $self->qsub($qsub_line);

    print "ge_wrapper is [$ge_wrapper]\n";
    print "command is\n";
    print "$command\n";

    # add by GT
    $gewrappercommand = "$ge_wrapper \"" . $command . "\"";
    print QSUB "$gewrappercommand\n";

    # change qsub_line by GT
    #$qsub_line = "qsub -l slot_limit=1  -V" . $qsub_line . " $test.sh\n";
    $qsub_line = "qsub  -V" . $qsub_line . " $test.sh\n";

    $self->qsub($qsub_line);

    print "SGE: have command line\n";

    #  print "CMD : $qsub_line\n\n";

    close(QSUB);
}

sub open_command_line {
    my ($self) = @_;

    print "run open_command_line in GridEngine.pm\n";

    print STDERR $self->qsub . " 2>&1 \n";
    print STDERR "opening command line\n";
    open( SUB, $self->qsub . " 2>&1 |" );
    my $geid;
    while (<SUB>) {
        if (/Your job (\d+)/) {
            $geid = $1;
        }
    }
    print STDERR "have opened " . $self->qsub . "\n";
    print STDERR "geid " . $geid . "\n";
    $self->id($geid);
    close(SUB);
}

#add by GT
sub temp_filename {
    print "run temp_filename in GridEngine.pm\n";
    my ($self) = @_;

    #for ( keys %ENV ) {
    #   print STDERR "key $_ $ENV{$_}\n";
    #}
    $self->{'sge_jobfilename'} = $ENV{'JOB_NAME'};

    print "ENV{JOB_NAME}=$ENV{'JOB_NAME'} in GridEngine.pm\n";

    print " JOBNAME " . $self->{'sge_jobfilename'} . "\n";
    return $self->{'sge_jobfilename'};
}

#add by GT
sub job_stats {
    print "run job_stats in GridEngine.pm\n";
    my ( $self, $verbose ) = @_;
    my $command = "qstat";

    # Need to sleep to make sure SGE is displaying correct info
    #print "sleeping sge \n";
    sleep(20);

    local *QSTAT;
    open( QSTAT, "$command 2>&1 |" ) or throw("couldn't open pipe to qstat");

    my %jobs;
  LINE:
    while ( my $line = <QSTAT> ) {

        #print "QSTAT $line";
        chomp($line);
        if ( $line =~ /^job-ID/ ) {
            next LINE;    # skip the header returned by qstat
        }
        if ( $line =~ /--------------/ ) {
            next LINE;    # skip the 2nd line
        }
        $line =~ s/^\s+//;

        #if ($_ =~ / /) {
        #  last LINE;
        #}
        my @values = split /\s+/, $line;

        # 422 0.55500 BLOOD=chun ensembl      qw    09/10/2010 16:35:46
        # record status of job; needs also to be configured in BatchQUeue STATUS_TO_COUNT
        $jobs{ $values[0] } = $values[4];

        #print "recording job $values[0] as $values[4]\n";
    }

    #sleep(10);
    return \%jobs;
}

sub memstring_to_resource {
    return '';
}
sub resource {
    return '';
}

1;

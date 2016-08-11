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

    Bio::EnsEMBL::Pipeline::Config::BatchQueue

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
    use Bio::EnsEMBL::Pipeline::Config::BatchQueue qw();

=head1 DESCRIPTION

    Configuration for pipeline batch queues. Specifies per-analysis
    resources and configuration, e.g. so that certain jobs are run only
    on certain nodes.

    It imports and sets a number of standard global variables into the
    calling package. Without arguments all the standard variables are
    set, and with a list, only those variables whose names are provided
    are set. The module will die if a variable which doesn't appear in
    its C<%Config> hash is asked to be set.

    The variables can also be references to arrays or hashes.

    Edit C<%Config> to add or alter variables.

    All the variables are in capitals, so that they resemble environment
    variables.

    To run a job only on a certain host, you have to add specific
    resource-requirements. This can be useful if you have special
    memory-requirements, for example if you like to run the job only
    on linux 64bit machines or if you want to run the job only on a
    specific host group. The commands bmgroup and lsinfo show you
    information about certain host-types / host-groups.

    Here are some example resource-statements / sub_args statements:

        sub_args => '-m bc_hosts',  # only use hosts of host-group 'bc_hosts'
                                    # (see bmgroup)
        sub_args => '-m bc1_1',     # only use hosts of host-group 'bc1_1'

        resource => 'select[type==X86_64]',     # use Linux 64 bit machines only
        resource => 'select[model==X86_64]',    # only run on X86_64 hosts

        resource => 'alpha',    # only run on DEC alpha
        resource => 'linux',    # run on any machine capable of running
                                # 32-bit X86 Linux apps

        # Note: On the Sanger farm, all machines are X86_64 Linux hosts.

=head2 Database throttling

    Do to find tokens for your MySQL servers 
    bhosts -s | grep tok

    This runs a job on a linux host, throttles genebuild2:3306 to not have
    more than 2000 active connections, 10 connections per job in the
    duration of the first 10 minutes when the job is running (means 200
    hosts * 10 connections = 2000 connections):

      resource => 'select[linux] rusage[myens_build2tok=10:duration=10]',

    Running on 'linux' hosts with not more than 2000 active connections
    for genebuild3 and genebuild7, 10 connections per job to each db-instance
    for the first 10 minutes:

      resource => 'select[linux] rusage[myens_build3tok=10:myens_build7tok=10:duration=10]',

    Running on hosts of model 'X86_64' hosts with not more than 200
    active connections to genebuild3, 10 connections per job for the first
    10 minutes:

      resource => 'select[model==X86_64] rusage[myens_build3tok=10:duration=10]',

    Running on hosts of host_group bc_hosts with not more than 200
    active connections to genebuild3, 10 connections per job for the first
    10 minutes:

      resource => 'rusage[myens_build3tok=10:duration=10]',

      sub_args =>'-m bc_hosts',

=head2 Memory requirements

    There are three ways of specifying LSF memory requirements.
    
=head3 Full resource string

    To allocate 4Gb of memory:

      resource  => 'select[mem>4000] rusage[mem=4000]',
      sub_args  => '-M 4000'

=head3 Short memory specification string

    To allocate 4Gb of memory, all of the below settings means the same
    thing (the case of the unit specification is ignored):

      memory => '4Gb'

      memory => '4000Mb'

      memory => '4000000Kb'

    Also see the section on memory arrays when using retries. 
      
=head2 Retrying on failure

    When a job fails due to insufficient memory, or due to run-time
    constraints, it may be retried.  The number of retries and the
    settings to use when retrying may be specified in two different
    ways.

    The maximum number of retries is set with the 'retries' option:

      retries => 3 # Run the job a maximum of four times (three retries).

=head3 Using 'retry_'

    When retrying, the pipeline submission code will look for options
    prefixed with 'retry_' and use these.  For example:

      memory        => '500Mb',     # Use 500 Mb for first run
      retry_memory  => '1Gb'        # Use 1 Gb for the retries

    The options that may be prefixed in this way are:

      queue,
      sub_args,
      resource,
      memory

=head3 Using arrays

    Instead of using the 'retry_' prefix, the original option may
    instead hold an array, like this:

      # Use 0.5Gb for the first run, 1Gb for the second, and 1.5Gb for
      # the third (and any later) run:
      memory => ['500Mb', '1Gb', '1500Mb']

    If the 'retries' value is larger than the length of the array, the
    last value of the array will be re-used.

=cut


package Bio::EnsEMBL::Pipeline::Config::BatchQueue;

use strict;
use vars qw(%Config);

%Config = (
  DEFAULT_BATCH_QUEUE => 'normal',
  DEFAULT_BATCH_SIZE => 10,
  DEFAULT_CLEANUP => 'no',
  DEFAULT_LSF_PERL => '/software/ensembl/central/bin/perl',
  DEFAULT_LSF_PRE_EXEC_PERL => '/software/ensembl/central/bin/perl',
  DEFAULT_OUTPUT_DIR => '',
  DEFAULT_RESOURCE => 'rusage[myens_build1tok=10,myens_build2tok=10,myens_build3tok=10,myens_build4tok=10,myens_build5tok=10,myens_build6tok=10,myens_build7tok=10,myens_build8tok=10]',
  DEFAULT_RETRIES => 2,
  DEFAULT_RETRY_QUEUE => 'long',
  DEFAULT_RETRY_RESOURCE => 'rusage[myens_build1tok=10,myens_build2tok=10,myens_build3tok=10,myens_build4tok=10,myens_build5tok=10,myens_build6tok=10,myens_build7tok=10,myens_build8tok=10]',
  DEFAULT_RETRY_SUB_ARGS => '',
  DEFAULT_RUNNABLEDB_PATH => 'Bio/EnsEMBL/Analysis/RunnableDB',
  DEFAULT_RUNNER => '',
  DEFAULT_SUB_ARGS => '',
  DEFAULT_VERBOSITY => 'WARNING',
  JOB_LIMIT => 10000,
  JOB_STATUSES_TO_COUNT => [
    'PEND'
  ],
  MARK_AWOL_JOBS => 1,
  MAX_JOB_SLEEP => 3600,
  MIN_JOB_SLEEP => 120,
  QUEUE_MANAGER => 'LSF',
  SLEEP_PER_JOB => 30,
  QUEUE_CONFIG => [
    {
      logic_name => 'cdna_update',
      batch_size => 10,
      output_dir => '/lustre/scratch109/sanger/cgg/cdna_update_human85/output',
      resource => '-R"select[mem>2000] rusage[mem=2000, myens_build12tok=20:duration=5, myens_build12tok=20]" -M 2000',
      retries => 1,
      retry_resource => '-R"select[mem>4000] rusage[mem=4000, myens_build12tok=20:duration=5, myens_build12tok=20]" -M 4000'
    },
    {
      logic_name => 'cdna_update_2',
      batch_size => 10,
      output_dir => '/lustre/scratch109/sanger/cgg/cdna_update_human85/output2',
      resource => '-R"select[mem>4000] rusage[mem=4000, myens_build12tok=20:duration=5, myens_build12tok=20]" -M 4000',
      retries => 1,
      retry_resource => '-R"select[mem>4000] rusage[mem=4000, myens_build12tok=20:duration=5, myens_build12tok=20]" -M 4000'
    },
    {
      logic_name => 'cdna_update_3',
      batch_size => 1,
      output_dir => '/lustre/scratch109/sanger/cgg/cdna_update_human85/output3',
      resource => '-R"select[mem>4000] rusage[mem=4000, myens_build12tok=20:duration=5, myens_build12tok=20]" -M 4000',
      retries => 1,
      retry_resource => '-R"select[mem>4000] rusage[mem=4000, myens_build12tok=20:duration=5, myens_build12tok=20]" -M 4000'
    }
  ]
);

sub import {
  my ($callpack) = caller(0);    # Name of the calling package
  my $pack = shift;              # Need to move package off @_

  # Get list of variables supplied, or else all
  my @vars = @_ ? @_ : keys(%Config);
  if ( !@vars ) { return }

  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw(" .
    join( ' ', map { '$' . $_ } @vars ) . ")";
  if ($@) { die $@ }

  foreach (@vars) {
    if ( defined( $Config{$_} ) ) {
      no strict 'refs';
      # Exporter does a similar job to the following
      # statement, but for function names, not
      # scalar variables
      *{"${callpack}::$_"} = \$Config{$_};
    }
    else {
      die("Error: Config: $_ not known\n");
    }
  }
} ## end sub import

1;

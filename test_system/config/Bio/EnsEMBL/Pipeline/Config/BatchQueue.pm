=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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

      resource =>
        'select[linux && myens_build2tok>2000] ' .
        'rusage[myens_build2tok=10:duration=10]',

    Running on 'linux' hosts with not more than 2000 active connections
    for genebuild3 and genebuild7, 10 connections per job to each db-instance
    for the first 10 minutes:

      resource =>
        'select[linux && myens_build3tok>200 && myens_build7tok>2000] ' .
        'rusage[myens_build3tok=10:myens_build7tok=10:duration=10]',

    Running on hosts of model 'X86_64' hosts with not more than 200
    active connections to genebuild3, 10 connections per job for the first
    10 minutes:

      resource =>
        'select[model==X86_64 && myens_build3tok>200] ' .
        'rusage[myens_build3tok=10:duration=10]',

    Running on hosts of host_group bc_hosts with not more than 200
    active connections to genebuild3, 10 connections per job for the first
    10 minutes:

      resource =>
        'select[myens_build3tok>200] ' .
        'rusage[myens_build3tok=10:duration=10]',

      sub_args =>'-m bc_hosts',

=head2 Memory requirements

    There are two ways of specifying LSF memory requirements.

=head3 Full resource string

    To allocate 4Gb of memory:

      resource  => 'select[mem>4000] rusage[mem=4000]',
      sub_args  => '-M 4000000'

=head3 Short memory specification string

    To allocate 4Gb of memory, all of the below settings means the same
    thing (the case of the unit specification is ignored):

      memory => '4Gb'

      memory => '4000Mb'

      memory => '4000000Kb'

=head2 Retrying on failure

    When a job fails due to insufficient memory, or due to run-time
    constraints, it may be retried.  The number of retries and the
    settings to use when retrying may be specified in two different
    ways.

    The maximum number of retries is set with the 'retries' option:

      reties => 3 # Run the job a maximum of four times (three retries).

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

  # Depending on the job-submission-system you're using, use LSF, you
  # can also use 'Local'.
  #
  # For mor info look into:
  # /ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/BatchSubmission

  QUEUE_MANAGER => 'LSF', # use "SGE_GridEngine_v6" for running in
                          # ensembl cloud evironment

  DEFAULT_BATCH_SIZE  => 10,
  DEFAULT_RETRIES     => 3,
  DEFAULT_BATCH_QUEUE => '',  # Put in the queue of your choice, e.g. 'normal'
  DEFAULT_RESOURCE    => 'select[myens_build1tok>800] rusage[myens_build1tok=10]',
  DEFAULT_SUB_ARGS    => '',
  DEFAULT_OUTPUT_DIR =>
    ( defined( $ENV{'TESTROOT'} ) ? $ENV{'TESTROOT'} : '.' ) .
    '/test_system_output',
  DEFAULT_CLEANUP     => 'no',
  DEFAULT_VERBOSITY   => 'WARNING',

  # The two variables below are to overcome a bug in LSF.  We're
  # currently running the pre-exec with a different perl.  LSF currently
  # unsets the LD_LIBRARY_PATH which we need for certain 64bit libraries
  # in pre-exec commands (more info see LSF_LD_SECURITY variable).

  DEFAULT_LSF_PRE_EXEC_PERL =>'/usr/local/ensembl32/bin/perl',
  # ONLY use 32bit perl for lsf -pre-exec jobs
  DEFAULT_LSF_PERL =>'/usr/local/ensembl32/bin/perl',
  # ONLY use ensembl64/bin/perl for memory jobs > 4 gb

  # SANGER farm: Don't forget to source source
  # /software/intel_cce_80/bin/iccvars.csh for big mem jobs.

  # At <this number of jobs> RuleManager will sleep for a certain period
  # of time.  If you effectively want this never to run set the value
  # to something very high, e.g. 100000.  This is important for queue
  # managers which cannot cope with large numbers of pending jobs (e.g.
  # early LSF versions and SGE).
  JOB_LIMIT => 10000,

  JOB_STATUSES_TO_COUNT => ['PEND'],    # These are the jobs which will
                                        # be counted. valid statuses
                                        # for this array are RUN, PEND,
                                        # SSUSP, EXIT, DONE ; use 'qw'
                                        # for Sun Grid Engine

  MARK_AWOL_JOBS => 1,
  MAX_JOB_SLEEP  => 3600,   # The maximun time to sleep for when job limit
                            # reached
  MIN_JOB_SLEEP => 120, # The minimum time to sleep for when job limit reached
  SLEEP_PER_JOB => 30,  # The amount of time to sleep per job when job limit
                        # reached

  DEFAULT_RUNNABLEDB_PATH => 'Bio/EnsEMBL/Analysis/RunnableDB',

  DEFAULT_RUNNER         => '',
  DEFAULT_RETRY_QUEUE    => 'long',
  DEFAULT_RETRY_SUB_ARGS => '',
  DEFAULT_RETRY_RESOURCE => '',

  QUEUE_CONFIG => [
    { logic_name      => 'repeatmask',
      batch_size      => 2,
      retries         => 5,
      runner          => '',
      queue           => 'normal',
      verbosity       => 'INFO',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
      lsf_perl        => '/usr/local/ensembl32/bin/perl',

      # Most RepeatMasker jobs need not more than 500MB.
      # Some jobs might fail (unlikely with 1M slices), but they will
      # defintely pass with 2GB.
      memory => [ '500Mb', '2Gb', '4Mb' ],

      retry_queue      => '',
      retry_batch_size => 5,
    },
    { logic_name => 'genscan',
      batch_size => 500,
      retries    => 3,
      runner     => '',

      memory => [ '400Mb', '700Mb' ],

      retry_queue    => '',
    },
    { logic_name => 'firstef',
      batch_size => 1000,
      retries    => 3,
      runner     => '',

      memory => [ '300Mb', '500Mb' ],

      retry_queue    => '',
    },
    { logic_name => 'cpg',
      batch_size => 100,
      retries    => 3,
      runner     => '',

      # cpg generally uses little memory
      # This is unlikely to be used as most cpg jobs uses less than
      # 100MB (default memory limit on Sanger LSF).
      memory => [ '100Mb', '200Mb' ],

      retry_queue => '',
    },
    { logic_name => 'job_using_more_than_4_gig_of_memory',
      batch_size => 10,
      retries    => 3,
      runner     => '',
      lsf_perl   => '/usr/local/ensembl64/bin/perl',

      resource => '',
      sub_args => '',

      retry_resource => '',
      retry_sub_args => '',
      retry_queue    => '',
    },
    { logic_name => 'eponine',
      batch_size => 1000,
      retries    => 3,
      runner     => '',

      # eponine will probably need about 1GB of memory, it seems to use
      # anywhere between 500MB and 700MB.  Use 2Gb if that's not enough.
      memory => [ '1Gb', '2Gb' ],
    },
    { logic_name  => 'trf',
      batch_size  => 1000,
      retries     => 3,
      runner      => '',
      retry_queue => '',

      # trf is a borderline case for the 100MB limit, give it 200MB.
      # For really big things, it might need more, give it 1GB
      memory => [ '200Mb', '1Gb' ],
    },
    { logic_name  => 'trnascan',
      batch_size  => 1000,
      retries     => 3,
      runner      => '',
      retry_queue => '',

      # trnascan has a similar memory footprint to trf
      memory => [ '200Mb', '1Gb' ],
    },
    { logic_name => 'uniprot',
      batch_size => 10,
      retries    => 5,
      runner     => '',

      memory => [ '300Mb', '500Mb', '1Gb' ],

      sub_args => '-sp 100',  # Run before vertrna and unigene
    },
    { logic_name => 'vertrna',
      batch_size => 50,
      retries    => 3,
      runner     => '',

      memory => '500Mb',

      sub_args => '-sp 95',    # Run after uniprot
    },
    { logic_name => 'unigene',
      batch_size => 50,
      retries    => 3,
      runner     => '',

      memory => '500Mb',

      sub_args => '-sp 90',    # Run after vertrna
    },
    {
      logic_name => 'sim_consensus',
      batch_size => 100,
      resource       => 'select[mem>2000] rusage[mem=2000]',
      retries        => 3,
      sub_args       => '-M 2000000',
      runner         => '',
      retry_queue    => '',
      retry_resource => '',
      retry_sub_args => '',
    },
    {
      logic_name => 'utr_addition',
      batch_size => 100,
      resource       => 'select[mem>1500] rusage[mem=1500]',
      retries        => 3,
      sub_args       => '-M 1500',
      runner         => '',
      retry_queue    => '',
      retry_resource => '',
      retry_sub_args => '',
    },
    {
      logic_name => 'LayerAnnotation',
      batch_size => 100,
      resource       => 'select[mem>1500] rusage[mem=1500]',
      retries        => 3,
      sub_args       => '-M 1500000',
      runner         => '',
      retry_queue    => '',
      retry_resource => '',
      retry_sub_args => '',
    },
    {
      logic_name => 'ensembl',
      batch_size => 100,
      resource       => 'select[mem>1000] rusage[mem=1000]',
      retries        => 3,
      sub_args       => '-M 1000000',
      runner         => '',
      retry_queue    => '',
      retry_resource => '',
      retry_sub_args => '',
    },
    {
      logic_name => 'prints',
      batch_size => 500,
      queue => 'normal',
      retries => 2,
      output_dir => '',
      resource => 'select[mem>2000] rusage[mem=2000]',
      sub_args => '-M 2000000',
    },
    {
      logic_name => 'tmhmm',
      batch_size => 1000,
      queue => 'normal',
      retries => 2,
      output_dir => '',
    },
    {
      logic_name => 'ncoils',
      batch_size => 1000,
      queue => 'small',
      retries => 2,
      output_dir => '',
    },
    {
      logic_name => 'signalp',
      batch_size => 20,
      queue => 'normal',
      retries => 2,
      output_dir => '',
    },
    {
      logic_name => 'seg',
      batch_size => 1,
      queue => 'small',
      retries => 2,
      output_dir => '',
      resource => 'select[mem>400] rusage[mem=400]',
      sub_args => '-M 400000',
    },
    {
      logic_name => 'pirsf',
      batch_size => 20,
      queue => 'normal',
      retries => 2,
      output_dir => '',
      resource => 'select[mem>2000] rusage[mem=2000]',
      sub_args => '-M 2000000',
    },
    {
      logic_name => 'smart',
      batch_size => 200,
      queue => 'normal',
      retries => 2,
      output_dir => '',
    },
    {
      logic_name => 'superfamily',
      batch_size => 500,
      queue => 'normal',
      retries => 2,
      output_dir => '',
      resource => 'select[mem>300] rusage[mem=300]',
      sub_args => '-M 300000',
    },
    {
      logic_name => 'tigrfam',
      batch_size => 1500,
      queue => 'normal',
      retries => 2,
      output_dir => '',
      resource => 'select[mem>500] rusage[mem=500]',
      sub_args => '-M 500000',
    },
    {
      logic_name => 'pfam',
      batch_size => 200,
      queue => 'normal',
      retries => 2,
      output_dir => '',
      resource => 'select[mem>500] rusage[mem=500]',
      sub_args => '-M 500000',
    },
    {
      logic_name => 'pfscan',
      batch_size => 100,
      queue => 'normal',
      retries => 2,
      output_dir => '',
    },


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

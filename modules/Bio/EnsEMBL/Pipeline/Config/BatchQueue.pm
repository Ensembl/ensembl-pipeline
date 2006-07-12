# EnsEMBL module for Bio::EnsEMBL::Pipeline::Config::BatchQueue;
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::BatchQueue

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::Config::BatchQueue;
    use Bio::EnsEMBL::Pipeline::Config::BatchQueue qw();

=head1 DESCRIPTION

Configuration for pipeline batch queues. Specifies per-analysis
resources and configuration, e.g. so that certain jobs are run
only on certain nodes.

It imports and sets a number of standard global variables into the
calling package. Without arguments all the standard variables are set,
and with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Config> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%Config> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Config::BatchQueue;

use strict;
use vars qw(%Config);

%Config = (
	QUEUE_MANAGER       => 'LSF',
	DEFAULT_BATCH_SIZE  => 1,
	DEFAULT_RETRIES     => 1,
	DEFAULT_BATCH_QUEUE =>
	  'normal',    # put in the queue  of your choice, eg. 'normal'
	DEFAULT_OUTPUT_DIR      => '/ecs4/scratch4/ml6/pipeline_default',
	DEFAULT_CLEANUP         => 'yes',
	DEFAULT_VERBOSITY       => 'WARNING',
	MAX_PENDING_JOBS        => 500,
	AUTO_JOB_UPDATE         => 1,
	SAVE_RUNTIME_INFO       => 1,
	DEFAULT_RUNNABLEDB_PATH => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
	DEFAULT_RUNNER          =>
'/ecs4/work5/finished/production_pipe/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Finished/runner.pl',
	JOB_LIMIT => 1000,    # the maximum number of farm jobs with the status below
	JOB_STATUSES_TO_COUNT => ['PEND'],    # these are the jobs which will be
	                                      # counted
	     # valid statuses for this array are RUN, PEND, SSUSP, EXIT, DONE
	MARK_AWOL_JOBS => 1,
	MAX_JOB_SLEEP  => 3600,    # the maximun time to sleep for when job limit
	                           # reached
	MIN_JOB_SLEEP => 120,  # the minium time to sleep for when job limit reached
	SLEEP_PER_JOB => 30,   # the amount of time to sleep per job when job limit
	                       # reached
	## MYSQL QUEUE PARAMETERS
	QUEUE_HOST => 'otterpipe1',
	QUEUE_NAME => 'pipe_queue',
	
	## BIG MEMORY JOB PARAMETERS
	BIG_MEM_PRIORITY =>
	  90,  # job priority value used for the out_of_memory jobs (must be unique)
	BIG_MEM_QUEUE    => 'bigmem',    # Big memory farm queue
	BIG_MEM_RESOURCE =>
	  'select[mem>2500] rusage[mem=1500]',    # Big memory resource requirement
	
	## URGENT JOB PARAMETERS
	URGENT_JOB_PRIORITY => 99,
	URGENT_INPUTID_FILE =>
	  '/ecs4/work5/finished/production_pipe/tmp/urgent_input_id',

	QUEUE_CONFIG => [
		{
			logic_name => 'RepeatMasker',
			batch_size => 2,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			verbosity  => 'INFO',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
			priority        => [80],
		},
		{
			logic_name => 'trf',
			batch_size => 100,
			resource   =>
			  'select[otterp1<500 && linux] rusage[otterp1=10:duration=2]',
			retries         => 1,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [75],
		},
		{
			logic_name => 'CpG',
			batch_size => 100,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
			priority        => [75],

		},
		{
			logic_name => 'Eponine',
			batch_size => 100,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
			priority        => [70],
		},
		{
			logic_name => 'marker',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [70],
		},
		{
			logic_name => 'Genscan',
			batch_size => 10,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
			priority        => [70],
		},
		{
			logic_name      => 'Fgenesh',
			batch_size      => 10,
			resource        => 'select[type==ALPHA5]',
			retries         => 1,
			sub_args        => '',
			nodes           => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
			priority        => [70],
		},
		{
			logic_name => 'tRNAscan',
			batch_size => 20,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
			priority        => [70],
		},
		{
			logic_name => 'Uniprot_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [52,32],
		},
		{
			logic_name => 'Uniprot',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [66,46],
		},
		{
			logic_name => 'Est2genome_human_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [54,34],
		},
		{
			logic_name => 'Est2genome_human',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [64,44],
		},
		{
			logic_name => 'Est2genome_mouse_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [54,34],
		},
		{
			logic_name => 'Est2genome_mouse',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [64,44],
		},
		{
			logic_name => 'Est2genome_other_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [54,34],
		},
		{
			logic_name => 'Est2genome_other',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [64,44],
		},
		{
			logic_name => 'Est2genome_fish_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [54,34],
		},
		{
			logic_name => 'Est2genome_fish',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [64,44],
		},
		{
			logic_name => 'Est2genome_clusters_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [54,34],
		},
		{
			logic_name => 'Est2genome_clusters',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [64,44],
		},
		{
			logic_name => 'EST2genome_xtrop_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [54,34],
		},
		{
			logic_name => 'EST2genome_xtrop',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [64,44],
		},
		{
			logic_name => 'Halfwise',
			batch_size => 1,
			resource   =>
			  'select[mypfam<100] rusage[mypfam=50:duration=10:decay=1]',
			retries         => 3,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [58],

		},
		{
			logic_name => 'Halfwise_new',
			batch_size => 1,
			resource   =>
			  'select[mypfam<100] rusage[mypfam=50:duration=10:decay=1]',
			retries         => 3,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [58],

		},
		{
			logic_name => 'refseq_human',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [60,40],
		},
		{
			logic_name => 'refseq_zebrafish',
			batch_size => 2,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [60,40],
		},
		{
			logic_name => 'refseq_mouse',
			batch_size => 2,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [60,40],
		},
		{
			logic_name => 'refseq_rat',
			batch_size => 2,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [60,40],
		},
		{
			logic_name => 'vertrna_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [56,36],
		},
		{
			logic_name => 'vertrna',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
			priority        => [68,48],
		},

		#									#
		#	ncRNA Pipeline configuration	#
		#									#
		{
			logic_name      => 'BlastmiRNA',
			batch_size      => 500,
			resource        => 'linux',
			retries         => 1,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
		},
		{
			logic_name      => 'miRNA',
			batch_size      => 1,
			resource        => 'linux',
			retries         => 1,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
		},
		{
			logic_name => 'BlastRfam',
			batch_size => 50,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'normal',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
		},
		{
			logic_name      => 'Rfam',
			batch_size      => 1,
			resource        => 'linux',
			retries         => 1,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
		},
		{
			logic_name      => 'BlastWait',
			batch_size      => 1,
			resource        => 'linux',
			retries         => 1,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
		},
		{
			logic_name      => 'SubmitmiRNA',
			batch_size      => 1,
			resource        => 'linux',
			retries         => 1,
			sub_args        => '',
			runner          => '',
			queue           => 'normal',
			cleanup         => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
		},

		#									#
		# Exonerate analysis configuration	#
		#									#
		{
			logic_name => 'Uniprot_exonerate_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
		},
		{
			logic_name => 'Est_exonerate_human_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
		},
		{
			logic_name => 'Est_exonerate_mouse_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
		},
		{
			logic_name => 'Est_exonerate_other_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 1,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
		},
		{
			logic_name => 'vertrna_exonerate_raw',
			batch_size => 1,
			resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
			retries    => 3,
			sub_args   => '',
			runner     => '',
			queue      => 'long',
			cleanup    => 'yes',
			runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB/Finished',
		},

	]
);

sub import {
	my ($callpack) = caller(0);    # Name of the calling package
	my $pack       = shift;        # Need to move package off @_
	                               # Get list of variables supplied, or else all
	my @vars = @_ ? @_ : keys(%Config);
	return unless @vars;

	# Predeclare global variables in calling package
	eval "package $callpack; use vars qw("
	  . join( ' ', map { '$' . $_ } @vars ) . ")";
	die $@ if $@;

	foreach (@vars) {
		if ( defined $Config{$_} ) {
			no strict 'refs';

			# Exporter does a similar job to the following
			# statement, but for function names, not
			# scalar variables:
			*{"${callpack}::$_"} = \$Config{$_};
		}
		else {
			die "Error: Config : $_ not known\n";
		}
	}
}

1;

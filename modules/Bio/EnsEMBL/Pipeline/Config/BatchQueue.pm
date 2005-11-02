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
	   DEFAULT_RETRIES     => 3,
	   DEFAULT_BATCH_QUEUE => 'normal', # put in the queue  of your choice, eg. 'normal'
	   DEFAULT_OUTPUT_DIR  => '/ecs4/scratch4/ck1/errors',
	   DEFAULT_CLEANUP     => 'yes',
	   MAX_PENDING_JOBS    => 500,
	   AUTO_JOB_UPDATE     => 1,
	   SAVE_RUNTIME_INFO   => 1,
	   QUEUE_CONFIG        => [

				   {
				    logic_name => 'RepeatMask',
				    batch_size => 2,
					#nodes      => 'ecs4_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes',
				   },
				   {
				    logic_name => 'GC',
				    batch_size => 100,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',

				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes',
				   },
				   {
				    logic_name => 'Full_dbSTS',
				    batch_size => 10,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Eponine',
				    batch_size => 100,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Genscan',
				    batch_size => 10,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'CpG',
				    batch_size => 100,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {		
				    logic_name => 'tRNAscan',
				    batch_size => 20,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'trf',
				    batch_size => 100,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Fgenesh',
				    batch_size => 10,
				    resource   => 'select[type==ALPHA5]',
				    #retries    => 3,
				    sub_args   => '',
				    nodes      => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Uniprot',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Est2genome_human',
				    batch_size => 1,
				    #nodes      => 'bc_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Est2genome_mouse',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    #	 nodes      => 'ecs4_nodes',
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Est2genome_other',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    #	 nodes      => 'ecs4_nodes',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Est2genome_clusters',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    #	 nodes      => 'ecs4_nodes',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Est2genome_human_old',
				    batch_queue => 'acarilong',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Est2genome_mouse_old',
				    batch_queue => 'acarilong',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Est2genome_other_old',
				    batch_queue => 'acarilong',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },

				   {
				    logic_name => 'Full_dbGSS',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Markers1',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Markers2',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'Halfwise',
				    batch_size => 1,
				    resource   => 'select[mypfam<100] rusage[mypfam=50:duration=10:decay=1]',
				    #retries    => 3,
				    sub_args   => '',
				    #	 nodes      => 'ecs_nodes',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'refseq_human',
				    batch_size => 1,
				    #	 nodes      => 'ecs_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'refseq_zebrafish',
				    batch_size => 2,
				    #	 nodes      => 'ecs_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 2,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'refseq_mouse',
				    batch_size => 2,
				    #	 nodes      => 'ecs_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 2,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'refseq_rat',
				    batch_size => 2,
				    #	 nodes      => 'ecs_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 2,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'refseq_chicken',
				    batch_size => 2,
				    #	 nodes      => 'ecs_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 2,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'vertrna',
				    batch_size => 1,
				    #	 nodes      => 'ecs_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'zfishEST',
				    batch_size => 1,
				    #	 nodes      => 'ecs4_nodes',
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'contig_vs_contig.exp',
				    batch_size => 1,
				    resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    #retries    => 3,
				    sub_args   => '',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'SubmitTranscript',
				    batch_size => 1,
				    # pfhmm is a multithreaded code so need to specifically specify number of CPU to use
				    resourece  => 'select[otterp1<500 && ncpus<4] rusage[otterp1=10:duration=2] span [hosts=1]';
				    #retries    => 3,
				    sub_args   => '-n 2',
				    cleanup    => 'yes',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'SubmitTranscriptChunk',
				    batch_size => 1,
				    #resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    resourece  => 'select[otterp1<500 && ncpus<4] rusage[otterp1=10:duration=2] span [hosts=1]';
				    #retries    => 3,
				    sub_args   => '-n 2',
				    cleanup    => 'yes',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				   {
				    logic_name => 'SubmitProteome',
				    batch_size => 1,
				    #resource   => 'select[otterp1<500] rusage[otterp1=10:duration=2]',
				    resourece  => 'select[otterp1<500 && ncpus<4] rusage[otterp1=10:duration=2] span [hosts=1]';
				    #retries    => 3,
				    sub_args   => '-n 2',
				    cleanup    => 'yes',
				    runner     => '',
				    queue      => 'normal',
				    cleanup    => 'yes'
				   },
				  ]
	  );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else all
    my @vars = @_ ? @_ : keys(%Config);
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;
    
    foreach (@vars) {
	if (defined $Config{ $_ }) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Config{ $_ };
	} else {
	    die "Error: Config: $_ not known\n";
	}
    }
}

1;

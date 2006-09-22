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
use LowCoverageGeneBuildConf;
use vars qw(%Config);

%Config = (
  QUEUE_MANAGER       => 'LSF',
  DEFAULT_BATCH_SIZE  => 10,
  DEFAULT_RETRIES     => 3,
  DEFAULT_BATCH_QUEUE => '', # put in the queue  of your choice, eg. 'acari'
  DEFAULT_OUTPUT_DIR  => $LC_scratchDIR."raw_computes",
  DEFAULT_CLEANUP     => 'no',	
  DEFAULT_VERBOSITY   => 'WARNING',

  AUTO_JOB_UPDATE     => 1,
  JOB_LIMIT => 100000, # at this number of jobs RuleManager will sleep for 
                      # a certain period of time if you effectively want this never to run set 
                      # the value to very high ie 100000 for a certain period of time
  JOB_STATUSES_TO_COUNT => ['PEND'], # these are the jobs which will be
                                     # counted
                                     # valid statuses for this array are RUN, PEND, SSUSP, EXIT, DONE   
  MARK_AWOL_JOBS      => 0,
  MAX_JOB_SLEEP       => 3600, # the maximun time to sleep for when job limit 
                               # reached
  MIN_JOB_SLEEP      => 120, # the minimum time to sleep for when job limit reached
  SLEEP_PER_JOB      => 30, # the amount of time to sleep per job when job limit 
                            # reached
  DEFAULT_RUNNABLEDB_PATH => 'Bio/EnsEMBL/Pipeline/RunnableDB',      

  DEFAULT_RUNNER => '',

  QUEUE_CONFIG       => [
    {
      logic_name => 'RepeatMask',
      batch_size => 100,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',   
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Supp_RepeatMask',
      batch_size => 100,
      # For running the tail
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',   
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Ab_initio_RepeatMask',
      batch_size => 100,
      # For running the tail      
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',   
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Genscan',        
      batch_size => 200,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',    
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Uniprot',        
      batch_size => 200,
      # For finishing off      
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Vertrna',        
      batch_size => 200,
      # For finishing off      
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Unigene',        
      batch_size => 200,
      # For finishing off      
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal', 
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'CpG',
      batch_size => 500,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Dust',
      batch_size => 500,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'tRNAscan',
      batch_size => 500,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'TRF',
      batch_size => 500,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Eponine',        
      batch_size => 200,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',    
    },
    {
      logic_name => 'marker',        
      batch_size => 10,
      resource   => '',
      retries    => 4,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',    
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'cdna_wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',     
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'exonerate_wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',   
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'hum_cdna_wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',     
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'hum_exonerate_wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',   
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Best_Wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',    
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Blast_Wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Pmatch_Wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',    
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Final_Wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',     
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'EST_Wait',        
      batch_size => 1,
      resource   => '',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',   
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Pmatch',        
      batch_size => 1,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',
      output_dir => $LC_scratchDIR."pmatch",
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'BestPmatch',        
      batch_size => 1,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."bestpmatch",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'TargettedGenewise',        
      batch_size => 4,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."targettedgenewise",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'SimilarityGenewise',        
      batch_size => 1,
      resource   => '',
      retries    => 0,
      sub_args   => '',
      runner     => '',
      queue      => 'long',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."sim_gw_85_nomini_rerun",
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'combined_gw_e2g',        
      batch_size => 100,
     # resource   => '',
     # resource   => 'select[myia64f<800 && ecs4my3350<2000 && my$LC_DBHOST<2000] rusage[myia64f=10:duration=10:decay=1:ecs4my3350=10:duration=10:decay=1:my$LC_DBHOST=10:duration=10:decay=1]',
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."humutrs",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'ensembl',        
      batch_size => 100,
      #resource   => 'select[myia64f<800 && ecs4my3351<400] rusage[myia64f=10:duration=30:decay=1:ecs4my3351=10:duration=30:decay=1]',
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."genebuild",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'pseudogene',        
      batch_size => 50,
     # resource   => '',
      #resource   => 'select[myia64f<800 && ecs4my3351<500 && my$LC_DBHOST<800] rusage[myia64f=10:duration=20:decay=1:ecs4my3351=10:duration=20:decay=1:my$LC_DBHOST=10:duration=20:decay=1]',
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."pseudo",     
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'GeneCombiner',        
      batch_size => 1,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."genecombiner",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    }, 
    {
      logic_name => 'cdna_exonerate',        
      batch_size => 1,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."cdna_exonerate",     
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'hum_cdna_exonerate',        
      batch_size => 1,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."hum_cdna_exonerate",     
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'est_exonerate',        
      batch_size => 1,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."est_exonerate",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'cdna_genebuilder',        
      batch_size => 200,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."cdna_genebuild",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'hum_cdna_genebuilder',        
      batch_size => 200,
      #resource   => 'select[myia64f<2000 && ecs4my3350<2000] rusage[myia64f=10:duration=10:decay=1:ecs4my3350=10:duration=10:decay=1]',
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."hum_cdna_genebuild",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'est_genebuilder',        
      batch_size => 1,
      resource   => '',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."est_genebuild",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'tmhmm',        
      batch_size => 10,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'scanprosite',        
      batch_size => 1,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'ncoils',        
      batch_size => 10,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'pfscan',        
      batch_size => 1,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'Pfam',        
      batch_size => 1,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 0,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'Signalp',        
      batch_size => 10,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'Seg',        
      batch_size => 1,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'Prints',        
      batch_size => 10,
      resource   => 'select[my$LC_DBHOST<800] rusage[my$LC_DBHOST=10:duration=10:decay=1]',
      retries    => 3,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      output_dir => $LC_scratchDIR."protein_annotation",     
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
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

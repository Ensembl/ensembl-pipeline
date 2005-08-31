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
  QUEUE_MANAGER       => 'Local',
  DEFAULT_BATCH_SIZE  => 10,
  DEFAULT_RETRIES     => 3,
  DEFAULT_BATCH_QUEUE => '', # put in the queue  of your choice, eg. 'acari'
  DEFAULT_OUTPUT_DIR  => '',	
  DEFAULT_CLEANUP     => 'yes',	
  JOB_LIMIT => 10000, # at this number of jobs RuleManager will sleep for 
                      # a certain period of time if you effectively want this never to run set 
                      # the value to very high ie 100000 for a certain period of time
  JOB_STATUSES_TO_COUNT => ['PEND'], # these are the jobs which will be
                                     # counted
                                     # valid statuses for this array are RUN, PEND, SSUSP, EXIT, DONE   
  MARK_AWOL_JOBS      => 1,
  MAX_JOB_SLEEP       => 3600, # the maximun time to sleep for when job limit 
                               # reached
  MIN_JOB_SLEEP      => 120, # the minium time to sleep for when job limit reached
  SLEEP_PER_JOB      => 30, # the amount of time to sleep per job when job limit 
                            # reached
  DEFAULT_RUNNABLEDB_PATH => 'Bio/EnsEMBL/Analysis/RunnableDB',      

  DEFAULT_RUNNER => '',

  QUEUE_CONFIG       => [
    {
      logic_name => 'RepeatMask',
      batch_size => 10,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Genscan',        
      batch_size => 10,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Genefinder',        
      batch_size => 10,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Snap',        
      batch_size => 5,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Fgenesh',        
      batch_size => 10,
      resource   => 'select[model==ES451250]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Swall',        
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'long',
      cleanup    => 'no', 
     runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Vertrna',        
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'long',
      cleanup    => 'no',  
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Unigene',        
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'long',
      cleanup    => 'no',     
     runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'CpG',
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',
 runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Dust',
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',
     runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'tRNAscan',
      batch_size => 200,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',
     runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'TRF',
      batch_size => 200,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Eponine',        
      batch_size => 10,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no', 
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',    
    },
    {
      logic_name => 'marker',        
      batch_size => 10,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'FirstEF',        
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',        
      runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },
    {
      logic_name => 'Uniprot',        
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'long',
     runnabledb_path => 'Bio/EnsEMBL/Analysis/RunnableDB',
    },          
     {
      logic_name => 'Pmatch',        
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',
      output_dir => '',
      runnabledb_path => 'Bio/EnsEMBL/Pipeline/RunnableDB',
    },
    {
      logic_name => 'BestPmatch',        
      batch_size => 1,
      resource   => 'select[model==IBMBC2800]',
      retries    => 1,
      sub_args   => '',
      runner     => '',
      queue      => 'normal',
      cleanup    => 'no',    
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

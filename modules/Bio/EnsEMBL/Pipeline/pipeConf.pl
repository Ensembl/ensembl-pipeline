# Copyright EMBL-EBI 2000
# Author: Arne Stabenau
# Creation: 11.07.2000
# Last modified SCP 12.04.2001

# configuration information
# give useful keynames to things

# if states are involved maybe good to have statenames in key.
# parameters should be avail after compile time so
# that you can make changes at runtime.

# some of these options can be specified on the command line (e.g. to
# the RuleManager script) and will override these defaults.
# it may also be possible to specify environment variables like
# ENS_<OPT> (these will be overridden by the settings below).


BEGIN {
package main;

%pipeConf = ( 
    'nfstmp.dir' => '/work1/scp/out',
                               # working directory for err/outfiles
    'pep_file'   => '/scratch1/ensembl/mongin/pep_tmp.fa',
    'DBI.driver' => 'mysql',
    'dbhost'     => 'ecs1e',
    'dbname'     => 'simon_test',
    'dbuser'     => 'ensadmin',
    'queue'      => 'acari',   # farm queue
    'batchsize'  => 1,         # no of jobs to send to LSF together
    'bindir'     => '/usr/local/ensembl/bin',
    'datadir'    => '/usr/local/ensembl/data',
    'usenodes'   => 'ecsnodes',        # farm nodes to use (default all)
    'autoupdate' => 1,          # true->update InputIdAnalysis via Job
    'runner'     => '',          # path to runner.pl, needed by Job.pm
);
}

1;

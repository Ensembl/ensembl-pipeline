# EnsEMBL module for Bio::EnsEMBL::Pipeline::Config::General;
#
# You may distribute this module under the same terms as perl itself


=head1 NAME

Bio::EnsEMBL::Pipeline::Config::General

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::Config::General;
    use Bio::EnsEMBL::Pipeline::Config::General qw();

=head1 DESCRIPTION

General pipeline configuration.

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


package Bio::EnsEMBL::Pipeline::Config::General;

use strict;
use vars qw(%Config);

%Config = (

  # binaries, libraries and data files
  BIN_DIR  => '/usr/local/ensembl/bin',
  DATA_DIR => '/usr/local/ensembl/data',
  LIB_DIR  => '/usr/local/ensembl/lib',


  # temporary working space (e.g. /tmp)
  PIPELINE_WORK_DIR   => '',

  # default input_file_directory
  PIPELINE_INPUT_DIR => '',

  # default target_file_directory
  PIPELINE_TARGET_DIR => '',

  #regex for spiting up slice input_ids
  SLICE_INPUT_ID_REGEX => '(\S+)\.(\d+)-(\d+):?([^:]*)',	   
    
  PIPELINE_REPEAT_MASKING => ['RepeatMask'],	

  SNAP_MASKING => [],

  MAX_JOB_TIME => 86400, # the max number of seconds a job should be 
                         # spending in CPU time in LSF before being killed 86400 is 24 hours

  KILLED_INPUT_IDS => '', # a path to produse a file of killed
                          # input_ids, good to have it local to the machine you running
                          # rulemanager on to prevent nfs stress
  RENAME_ON_RETRY => 1, # toggle to see if you want the stdout/err
                        #files renamed when a job is retried
    
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

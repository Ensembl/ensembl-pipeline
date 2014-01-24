=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

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


=cut


package Bio::EnsEMBL::Pipeline::Config::General;

use strict;
use LowCoverageGeneBuildConf;
use vars qw(%Config);

%Config = (

  RENAME_ON_RETRY => 1, # toggle to see if you want the stdout/err
                        #files renamed when a job is retried
  SGE_PERL5LIB_ENV_SCRIPT => "", # path to script to set up PERL5LIB enviromnet for Sun Grid Engine
                                 # you find an example in ensembl-pipeline/scripts/set_PERL5LIB_SGE.sh 

    
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

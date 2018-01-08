=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Config::General

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::General;
    use Bio::EnsEMBL::Analysis::Config::General qw();

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

package Bio::EnsEMBL::Analysis::Config::General;

use strict;
use LowCoverageGeneBuildConf;
use vars qw(%Config);

%Config = (

    # binaries, libraries and data files
    BIN_DIR  => '/usr/local/ensembl/bin',
    DATA_DIR => '/usr/local/ensembl/data',
    LIB_DIR  => '/usr/local/ensembl/lib',

    # Path where the parser and parameter files for FirstEF program are allocated
    PARAMETERS_DIR => '/vol/software/linux-i386/farm/lib/firstef/parameters/',
    PARSE_SCRIPT => '/vol/software/linux-i386/farm/lib/firstef/FirstEF_parser.pl',

    ANALYSIS_WORK_DIR => '/tmp',
    # LC_REPMASK_CHOICE will be something like [ 'RepeatMask' ]
    ANALYSIS_REPEAT_MASKING => $LC_REPMASK_CHOICE,

    CORE_VERBOSITY => 'WARNING',
    LOGGER_VERBOSITY => 'OFF',


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

package Bio::EnsEMBL::Analysis::Config::General;

use strict;
use vars qw(%Config);

%Config = (

           # binaries, libraries and data files
           BIN_DIR  => '/usr/local/ensembl/bin',
           DATA_DIR => '/usr/local/ensembl/data',
           LIB_DIR  => '/usr/local/ensembl/lib',

           ANALYSIS_WORK_DIR => '/tmp',
           ANALYSIS_REPEAT_MASKING => ['RepeatMask'],      
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

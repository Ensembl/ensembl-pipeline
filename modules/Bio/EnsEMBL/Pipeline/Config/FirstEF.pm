# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::FirstEF - imports global variables used by the FirstEF Runnable/RunnableDB

=head1 SYNOPSIS
    use Bio::EnsEMBL::Pipeline::Config::FirstEF;
    use Bio::EnsEMBL::Pipeline::Config::FirstEF qw(  );

=head1 DESCRIPTION

General is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%General> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%General> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::Config::FirstEF;

use strict;
use vars qw( %FirstEF );

# Hash containing config info
%FirstEF = (
	    # Input Id regular expression (to parse something like chr.start-end)
	    FEF_INPUTID_REGEX => '(\S+)\.(\d+)-(\d+)',

	    # Temp directory for LSF output
	    FEF_TMPDIR   => '',

	    # The name of the file that will contain all of the LSF bsub commands.
	    FEF_BSUB_FILE => '',

	    # Name of the batch queue to use.
	    FEF_QUEUE => 'acari',

	    # Stipulate the size of the slices that firstef should run on.
	    FEF_CHUNKSIZE => ''.

	    # Name of the script that actually runs the analysis for each slice.
	    FEF_RUN_SCRIPT => '',

	    # Reference database where assemble, chromosome, dna, etc information is stored.
	    FEF_REFDBHOST => '',
	    FEF_REFDBUSER => '',
	    FEF_REFDBNAME => '',

	    # The database where the firstef features will be written
	    FEF_WRITEDBHOST => '',
	    FEF_WRITEDBUSER => '',
	    FEF_WRITEDBNAME => '',
	    FEF_WRITEDBPASS => '',

	    # Directory where firstef.* binaries and FirstEF_parser.pl 
	    # are located.  At runtime, the runnable determines which
	    # platform specific binary should be used.
	    FEF_APPLICATION_DIR => '/usr/local/ensembl/firstef/',

	    # Usually, this directory is a sub-directory of the firstef
	    # application directory called 'parameters'.  This directory
	    # contains all of the training data that firstef uses for
	    # identifying first exons.
	    FEF_PARAMETER_DIR   => '/usr/local/ensembl/firstef/parameters',
	   );

sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_
  
  # Get list of variables supplied, or else
  # all of General:
  my @vars = @_ ? @_ : keys( %General );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $General{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$General{ $_ };
	} else {
	    die "Error: General: $_ not known\n";
	}
    }
}

1;

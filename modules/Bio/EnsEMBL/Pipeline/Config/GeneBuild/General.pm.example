# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::GeneBuild::General - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS
    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General;
    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General qw(  );

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


package Bio::EnsEMBL::Pipeline::Config::GeneBuild::General;

use strict;
use vars qw( %General );

# Hash containing config info
%General = (
	    

	    # This is a regexp for parsing input ids - eg
	    # chrname.start-end would be matched by (\S+)\.(\d+)-(\d+)
	    GB_INPUTID_REGEX => '(^\S+)\.(\d+)-(\d+)',

	    # skip blastminigenewise stage? Only relevant for Contig_BMG
	    GB_SKIP_BMG   =>  0,
	    GB_BMG_FILTER => 0,  #whether to filter blast hits from blast on score
	    GB_BMG_SCORE_CUTOFF => 0,	  #the score to cut off for the filter we use 150

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

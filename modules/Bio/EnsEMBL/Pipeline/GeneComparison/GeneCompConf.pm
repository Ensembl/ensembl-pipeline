# configuration file for run_GeneComparison
# based on the config files for the pipeline
#
# Written by Eduardo Eyras
# eae@sanger.ac.uk

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneConf - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf;
    use Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf qw(  );

=head1 DESCRIPTION

GeneCompConf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%GeneConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%GeneConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Pipeline::GeneComparison::GeneCompConf;

use strict;
use vars qw( %GeneCompConf );

# Hash containing config info
%GeneCompConf = (
		 
		 # database with the annotation/benchmark genes to compare to
		 DBHOST1    => "",
	         DBNAME1    => "",
		 PATH1      => "",    
		 DBUSER1    => "",

		 # genetypes is an arrayref so that we can include more than one type
		 GENETYPES1 => [""],
		 

		 # database with the prediction/test genes to compare with the annotation/benchmark
		 DBHOST2    => "",                    
		 DBNAME2    => "",
		 PATH2      => "",
		 DBUSER2    => "",
		 
		 # genetypes is an arrayref so that we can include more than one type
		 GENETYPES2 => [""], 
	     
);

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %GeneCompConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw(". join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;

    foreach (@vars) {
	if ( defined $GeneCompConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$GeneCompConf{ $_ };
	} else {
	    die "Error: GeneCompConf: $_ not known\n";
	}
    }
}

1;

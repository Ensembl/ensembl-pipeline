#
# BioPerl module for Bio::EnsEMBL::Analysis::BioperlDBConf;
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::BioperlDBConf - imports global variables needed for using bioperlDB

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::BioperlDBConf;
    use Bio::EnsEMBL::Pipeline::BioperlDBConf qw(  );

=head1 DESCRIPTION

BioperlDBConf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%BioperlDBConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%BioperlDBConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::BioperlDBConf;

use strict;
use vars qw( %BioperlDBConf );

# Hash containing config info
%BioperlDBConf = (
		  BIOPERLDB               => 0,
		  BPNAME                  => '',	     
		  BPUSER                  => '',
		  BPBIODB_ID              => '',
		  BP_SUPPORTING_DATABASES => '',
		  BP_DBI_DRIVER           => '',
		  
	    );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %BioperlDBConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $BioperlDBConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$BioperlDBConf{ $_ };
	} else {
	    die "Error: BioperlDBConf: $_ not known\n";
	}
    }
}

1;

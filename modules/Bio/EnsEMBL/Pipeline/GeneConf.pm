#
# BioPerl module for Bio::EnsEMBL::Analysis::GeneConf;
#
# Cared for by Tim Hubbard <th@sanger.ac.uk>
#
# Copyright Tim Hubbard, James Gilbert
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::GeneConf - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::GeneConf;
    use Bio::EnsEMBL::Pipeline::GeneConf qw(  );

=head1 DESCRIPTION

GeneConf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn't appear in its
C<%GeneConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%GeneConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

B<Tim Hubbard> email th@sanger.ac.uk
B<James Gilbert> email jgrg@sanger.ac.uk

=cut

#'

package Bio::EnsEMBL::Pipeline::GeneConf;

use strict;
use vars qw( %GeneConf );

my $prefix='SEPT20';

# Hash containing config info
%GeneConf = (
	    EXON_ID_SUBSCRIPT => $prefix.'E',
	    EXON_ID_DIGITS => 11,
	    TRANSCRIPT_ID_SUBSCRIPT => $prefix.'T',
	    TRANSCRIPT_ID_DIGITS => 11,
	    GENE_ID_SUBSCRIPT => $prefix.'G',
	    GENE_ID_DIGITS => 11,
	    PROTEIN_ID_SUBSCRIPT => $prefix.'P',
	    PROTEIN_ID_DIGITS => 11,
	    );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %GeneConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $GeneConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$GeneConf{ $_ };
	} else {
	    die "Error: GeneConf: $_ not known\n";
	}
    }
}

1;

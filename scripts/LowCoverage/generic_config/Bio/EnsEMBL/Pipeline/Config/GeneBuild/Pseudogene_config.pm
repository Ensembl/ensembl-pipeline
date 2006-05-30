# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pseudogene_config - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS
    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pseudogene_config;
    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pseudogene_config qw(  );

=head1 DESCRIPTION

Pseudogene_config is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Pseudogene_config> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%Pseudogene_config> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pseudogene_config;

use strict;
use vars qw( %Pseudogene_config );

# Hash containing config info
%Pseudogene_config = (	

		# configs for the introns in repeats test

		# total length of introns
		PS_MAX_INTRON_LENGTH   => '5000',

		# %coverd  repeats
		PS_MAX_INTRON_COVERAGE => '80',

		# %coverd by repeats
		PS_MAX_EXON_COVERAGE   => '20',	
		
		PS_NUM_FRAMESHIFT_INTRONS  => 1,

		PS_NUM_REAL_INTRONS  => 1,

		# configs for the spliced elsewhere tests

		# %ID of a tbalstx of the (presumed) retrotransposed query sequence to its 
		# homolog that is spliced elsewhere in the genome. hits falling below 
		# this cutoff are ignored (80%) is suggested
		PS_PERCENT_ID_CUTOFF   => 80,

		# ratio of the spans of the retrotransposed gene vs its spliced homologue
		# spliced / retrotransposed
		# ie: 1 is the same length genes 
		# many retrotransposed genes have a ratio > 10
		# used to make retrotransposition decision
		PS_SPAN_RATIO          => 1.5,

		# mimimum number of exons for the spliced gene to have
		PS_MIN_EXONS           => 2,
		# limit the number of hits from the align feature table
		PS_SQL_LIMIT           => 100,
		# minimum score of the align feature required
		PS_FEATURE_SCORE       => 300,

	       );	

sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_
  
  # Get list of variables supplied, or else
  # all of Pseudogene_config:
  my @vars = @_ ? @_ : keys( %Pseudogene_config );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $Pseudogene_config{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Pseudogene_config{ $_ };
	} else {
	    die "Error: Pseudogene_config: $_ not known\n";
	}
    }
}

1;

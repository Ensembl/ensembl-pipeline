# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch - imports global variables used by EnsEMBL gene building

=head1 SYNOPSIS
    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch;
    use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch qw(  );

=head1 DESCRIPTION

Pmatch is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Pmatch> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%Pmatch> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

for more information on how we run out pmatch analysis see 
running_the_genebuild.txt in the cvs module ensembl-docs

=cut


package Bio::EnsEMBL::Pipeline::Config::GeneBuild::Pmatch;

use strict;
use vars qw( %Pmatch );

# Hash containing config info
%Pmatch = (   
	   GB_PFASTA => $ENV{PMATCH_FASTA}, #location of fasta file to be written by new_prepare_proteome.pl and used in pmatch analysis
	    
	   # path to pmatch executable
	   GB_PMATCH      => '/usr/local/ensembl/bin/pmatch',
	    
	   # maximum allowable intron length for stitching individual pmatch hits back into a gene region
	   GB_PMATCH_MAX_INTRON => '50000',  
	   GB_PMATCH_MASKING => [], 
	   #the above setting will lead to no repeat being masked out, if [''] is used all repeats in the table will be masked out
	   #if only some of the sets of repeats wants to be masked out put logic names of the analyses in the array e.g ['RepeatMask', 'Dust']
	   GB_PMATCH_SOFTMASK => 0, #set to one if you want the repeats softmasked ie lower case bases rather than uppercase N's
	   GB_INITIAL_PMATCH_LOGICNAME => 'Pmatch', #logic_name of first pmatch analysis  
	   GB_FINAL_PMATCH_LOGICNAME => 'BestPmatch', #logicname of best in genome analysis
	  );

sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_
  
  # Get list of variables supplied, or else
  # all of Pmatch:
  my @vars = @_ ? @_ : keys( %Pmatch );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $Pmatch{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$Pmatch{ $_ };
	} else {
	    die "Error: Pmatch: $_ not known\n";
	}
    }
}

1;

package Bio::EnsEMBL::Analysis::Config::GeneBuild::Pmatch;

use strict;
use vars qw( %Config );

# Hash containing config info

%Config = (


            #Under each the two keys "PMATCH_BY_LOGIC" and "BESTPMATCH_BY_LOGIC",
            #the DEFAULT hash must be present and should contain default settings.
            #Non-standard settings should be provided in a second hash that follows
            #the DEFAULT one with the same structure, using the logic_name as the
            #second hash's key, e.g. "PMATCH".

            # ***** WATCH OUT ***** #

            #The second hash (mainly its key) MUST be present EVEN IF there are no
            #non-standard settings to be defined. In that case, just leave the second
            #hash empty with no values, e.g.
            #
            #PMATCH => {
            #           },
            #
            #Missing out the DEFAULT and/or the second hash will result in the config
            #not being read properly and the analysis will fail right from the start.

           PMATCH_BY_LOGIC =>
           {
            DEFAULT =>{
                       # Before testing Pmatch analysis in pipeline test system, modify the protein file path, which should look like
                       # something like this:
                       #
                       # PROTEIN_FILE => '/your/working/directory/for/running/the/test/homo_sapiens/data/pmatch_proteins.fa.',
                       #
                       # PROTEIN_FILE must always be specified. The "working directory" will always be in ensembl-pipeline/test_system as you'll 
                       # have to be in that directory to run ANY test in the first place.
                       
                       PROTEIN_FILE => '/your/cvs/checkout/dir/ensembl-pipeline/test_system/homo_sapiens/data/pmatch_proteins.fa',
                       MIN_COVERAGE => 25,
                       BINARY_LOCATION => '/usr/local/ensembl/bin/pmatch',
                       REPEAT_MASKING => [],
                       MAX_INTRON_LENGTH => 50000,
                       OUTPUT_DB => 'REFERENCE_DB',
                      },
            PMATCH => {
                       },

            #If you wish any of the settings to be non standard generate a hash
            #with the same structure but change DEFAULT to be your analysis
            #logic_name
           },
           BESTPMATCH_BY_LOGIC =>
           {
            DEFAULT =>{ 
                        PMATCH_LOGIC_NAME => 'Pmatch',
                        MIN_COVERAGE => 25,
                        INPUT_DB => 'REFERENCE_DB',
                        OUTPUT_DB => 'REFERENCE_DB',
                       },
            BESTPMATCH => {
                       },
           },
);

sub import {
  my ($callpack) = caller(0); # Name of the calling package
  my $pack = shift; # Need to move package off @_

  # Get list of variables supplied, or else everything
  my @vars = @_ ? @_ : keys( %Config );
  return unless @vars;
  
  # Predeclare global variables in calling package
  eval "package $callpack; use vars qw("
    . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $Config{$_} ) {
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

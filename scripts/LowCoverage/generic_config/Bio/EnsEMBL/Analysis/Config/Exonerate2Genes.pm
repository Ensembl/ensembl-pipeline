#
# package Bio::EnsEMBL::Pipeline::Config::ExonerateTranscript
# 
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Config::Exonerate2Genes

=head1 SYNOPSIS

    use Bio::EnsEMBL::Pipeline::Config::Exonerate2Genes;

=head1 DESCRIPTION

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first declared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Exonerate> hash is asked to be set.

Since The RunnableDB that this config controls can be used for 
inferring transcript structures from (different sets of) EST, cDNA 
and proteins, and several uses may be required in the same pipeline, 
this Config contains one primary config variable, EXONERATE_TRANSCRIPT_CONFIG.
This is hash keyed off logic name, each entry of which is a hash
containing the variable that affect the behaviour of the RunnableDB.
When the RunnableDB instance is created, the correct entry is identified
by logic name and value for a corresponding set of local variables are
set.

=head1 CONTACT

=cut


package Bio::EnsEMBL::Analysis::Config::Exonerate2Genes;

use strict;
use LowCoverageGeneBuildConf;
use vars qw( %Config );

# Hash containing config info
%Config = (
  EXONERATE_CONFIG_BY_LOGIC => {
    DEFAULT => {
      GENOMICSEQS => $LC_GENOMICSEQS,
      QUERYTYPE   => undef,
      QUERYSEQS   => undef,
      IIDREGEXP   => undef,
      OUTDB       => undef,
      FILTER      => undef,
      COVERAGE_BY_ALIGNED => undef,
#      OPTIONS             => '--exhaustive FALSE --model est2genome --softmasktarget --score 500 --fsmmemory 800  --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14',
      OPTIONS             => '--exhaustive FALSE --model est2genome --softmasktarget --score 500 --fsmmemory 800  --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14 --bestn 10',
    },

    cdna_exonerate => {
      # GENOMICSEQS obtained from DEFAULT
      QUERYTYPE  => 'dna',
      QUERYSEQS  => $LC_workDIR."seqdata/cDNAs/",
      # IIDREGEXP not set; input ids are file names
      OUTDB => { -dbname => $LC_DBprefix."cDNA",
                 -host   => $LC_DBHOST,
                 -port   => $LC_DBPORT,
                 -user   => $LC_DBUSER,
                 -pass   => $LC_DBPASS,
      },               
      FILTER => { OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                  PARAMETERS => {
                          -coverage => 90,
                          -percent_id => 97,
                          -best_in_genome => 1,
                          -reject_processed_pseudos => 1,
                  },              
      },
    },
    hum_cdna_exonerate => {
      # GENOMICSEQS obtained from DEFAULT
      QUERYTYPE  => 'dna',
      QUERYSEQS  => $LC_workDIR."seqdata/human_cDNAs/",
      # IIDREGEXP not set; input ids are file names
      OUTDB => { -dbname => $LC_DBprefix."hum_cDNA",
                 -host   => $LC_DBHOST,
                 -port   => $LC_DBPORT,
                 -user   => $LC_DBUSER,
                 -pass   => $LC_DBPASS,
      },               
      FILTER => { OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                  PARAMETERS => {
                          -coverage 		    => 40,
                          -percent_id 		    => 40,
                          -best_in_genome 	    => 1,
                          -reject_processed_pseudos => 1,
                  },              
      },
    },
  }
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

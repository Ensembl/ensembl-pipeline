# Ensembl module for Bio::EnsEMBL::Analysis::Config::Blast
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Config::Blast

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Config::Blast;
  
  use Bio::EnsEMBL::Analysis::Config::Blast qw(BLAST_CONTIG);

=head1 DESCRIPTION

This is a module needed to provide configuration for the
blast RunnableDBs. This informs the various blast runnabledbs what
Parser and Filter objects they should instantiate and also of any
constructor arguments which should go to the Blast Runnable note any
Blast constructor arguments will be overridden by the same key in
analysis table parameters column

BLAST_CONFIG is an hash of hashes which contains analysis specific
settings and is keyed on logic_name

Important values are the key of the hash which should be the same as the
the logic name of the analysis you are configuring
and blast parser which should be a perl path to a object to parser the 
blast report. 

the blast_parser object should fit the standard interface which is to 
have a method called parse_file which accepts a filename as an argument and
returns a set of results all other constructor arguements are optional but 
the two parser objects which currently exist both need a regex, query type
and database type the types needed for different flavours of blast
can be found below and the two parser objects BPliteWrapper and 
FilterBPlite both live in Bio/EnsEMBL/Analysis/Tools

blast_filter should be a perl path to a filter object The filter module 
isn't obligatory, if none is specified then the blast results won't be 
filters after parsing (note some parsers may do filtering) The only filter
object which currently exists is 
Bio::EnsEMBL::Analysis::Tools::FeatureFilter. All filter objects must
have a filter_results method but aside from that there is no other 
requirements. Any constructor args should be specified in the hash. 
FeatureFilter can take min_score, max_pvalue coverage prune and hardprune

blast_params is any constructor parameters for which ever blast module
you are using note these will be overridden by values got from the 
parameters column of the analysis table

AB_INITIO_LOGICNAME if for BlastGenscanPep/DNA runnabledbs for which
ab initio predictions to use when running the blast and they currently 
default to Genscan

this example file contains a default setting which is what we use for out
BlastGenscanPep against Swall 
=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Config::Blast;

use strict;
use vars qw(%Config);


# PARSER_PARAMS-table for Blast-Configuration (below):
#
# Analysis  query_type    database_type
#-------------------------------------
#  blastp     pep            pep
#  blastn     dna            dna
#  blastx     dna            pep         (i.e.Uniprot-wublastx)
# tblastn     pep            dna
# tblastx     dna            dna
#
#
# For running an unfiltered blastx against Uniprot with the logic_name "Uniprot" in the 
# analysis-table, use query_type='dna' and database-type 'pep'.
#
#    ########  Example-configuration for blastx against UNIPROT  #########
#    Uniprot =>
#            {
#             BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
#             PARSER_PARAMS => {
#                               -regex => '^\w+\s+(\w+)',
#                               -query_type => 'dna',       # see PARSER_PARAMS-table
#                               -database_type => 'pep',    # see PARSER_PARAMS-table
#                               -threshold_type => 'PVALUE',
#                               -threshold => 0.01,
#                              },
#             BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
#             FILTER_PARAMS => { },
#             BLAST_PARAMS => {
#                              -unknown_error_string => 'FAILED',
#                              -type => 'wu',
#                             },



%Config = (
  BLAST_CONFIG => {
    DEFAULT => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::BPliteWrapper',
      PARSER_PARAMS => {
        -regex => '^(\w+)',
        -query_type => undef,
        -database_type => undef,
      },
      BLAST_FILTER => undef,
      FILTER_PARAMS => {},
      BLAST_PARAMS => {
        -unknown_error_string => 'FAILED',
        -type => 'wu',
      }
    }
  },
  BLAST_AB_INITIO_LOGICNAME => 'Genscan',
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

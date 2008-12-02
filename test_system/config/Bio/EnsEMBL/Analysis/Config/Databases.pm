# Ensembl module fohr Bio::EnsEMBL::Analysis::Config::Databases
#
# Copyright (c) 2006 Ensembl
#
# POD documentation - main docs before the code
#

=head1 NAME

Bio::EnsEMBL::Analysis::Config::Databases 

=head1 SYNOPSIS

    use Bio::EnsEMBL::Analysis::Config::Databases ; 

=head1 DESCRIPTION

Databases.pm is the main configuration file which holds the different 
parameters (usernames, hosts, passwords, ports, database-names) to 
connect to different databases used in the Ensembl-Analysis pipeline. 

It imports and sets a number of standard global variables into the
calling package. Without arguments all the standard variables are set,
and with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn\'t appear in its
C<%Config> hash is asked to be set.

 
A common way to get an DBAdaptor in a module is 
 
 print "Loading database : ".  $$DATABASES{REFERENCE_DB}{"-dbname"} . "\n";  

 my $ref_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( %{ $$DATABASES{REFERENCE_DB} }) ;

 
 OR if you write a RunnableDB:  


   use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
   @ISA = qw ( Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild ) 

   my $genes_db = $self->get_dbadaptor("GENEBUILD_DB");

 
 OR for normal scripts :

   use Bio::EnsEMBL::Analysis::Tools::Utilities qw( get_db_adaptor_by_string );
   get_db_adaptor_by_string("GENEBUILD_DB")

  
 
The variables can also be references to arrays or hashes.

Edit C<%Config> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables. 

Databases is a pure ripoff of humConf written by James Gilbert.
humConf is based upon ideas from the standard perl Env environment
module.


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Config::Databases ; 

use strict;
use vars qw(%Config);

%Config= (


  DATABASES => { 
                  
                 # The REFERENCE_DB (formely known as GB_DB) holds sequence + repeats + features 
                 # from raw computes (e.g. ab-inito predictions, dna- or protein alignments ) 
                  
                 REFERENCE_DB => 
                                 { 
                                   -dbname => 'at6_system_test',
                                   -host => 'genebuild2',
                                   -port => '3306',
                                   -user => 'ensadmin',
                                   -pass => 'ensembl',
                                  },
  
 
                 # The GENEWISE_DB holds genes made by FPC_TargettedGenewise or FPC_BlastMiniGenewise 
                 # (TGE_gw or similarity_genewise - genes ) ( formerly GB_GW_DB ) 

                 GENEWISE_DB => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },


                 # The EXONERATE_DB ( formerly GB_cDNA ) holds alignments to cDNA's or 
                 # EST's gene-structtures made by exonerate (Exonerate2Genes.pm) 
                 
                 EXONERATE_DB => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },


                 # The BLESSED_DB (formerly GB_BLESSED) holds the 'blessed' gene-set ( if there is one ) 

                 BLESSED_DB => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
                                 
                   
                 # The UTR_DB (formerly GB_COMB) holds genes made by the UTR-addtion-run 
                 # Combine_Genewises_and_E2Gs.pm writes to UTR_DB 
                   
                 UTR_DB  => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
                                 
                 # 
                 # GENEBUILD_DB (formerly GB_FINALDB) is the Database where 
                 # GeneBuilder.pm writes it's results to this database and 
                 # The Pseudogene-code READS from this database 
                 # 

                 GENEBUILD_DB =>
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
         
                        
                 # PSEUDO_DB holds the pseudo-genes
           
                 PSEUDO_DB => 
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
                                 
                   
                 # COALESCER_DB is the DB where TranscriptCoalescer writes it's results to
                   
                 COALESCER_DB =>
                                 { 
                                   -dbname => '',
                                   -host => '',
                                   -port => '',
                                   -user => '',
                                   -pass => '',
                                  },
                                 
                 # OrthologueEvalutor, FindMissingGenes, etc. write it's results into ORTHOLOGUE_DB 

                  ORTHOLOGUE_DB => {
                                     -dbname => '',
                                     -host =>   '',
                                     -port =>   '',
                                     -user =>   '',
                                     -pass =>   '',
                                   },


                # This database is use mainly with human an mouse for merging the 
                # havana get set with the ensembl genebuild 

                  HAVANA_DB  => {
                                     -dbname => '',
                                     -host =>   '',
                                     -port =>   '',
                                     -user =>   '',
                                     -pass =>   '',
                                   },

                # This database is use mainly with human variation genewise alignments
                  VARIATION_DB  => {
                                     -dbname => '',
                                     -host =>   '',
                                     -port =>   '',
                                     -user =>   '',
                                     -pass =>   '',
                                   },
 
 
                 
                # Database which stores the ditags...
 
                  DITAG_DB => {
                                     -dbname => '',
                                     -host =>   '',
                                     -port =>   '',
                                     -user =>   '',
                                     -pass =>   '',
                                   },
 
                 #############
                 ### Databases for RunnableDB/CopyGenes.pm
                 #############

                 COPY_SOURCE_DB => {
                     -dbname => '',
                     -host => '',
                     -port => '',
                     -user => '',
                     -pass => '',
                 },

                 COPY_TARGET_DB => {
                     -dbname => '',
                     -host => '',
                     -port => '',
                     -user => '',
                     -pass => '',
                },
             },
             DNA_DBNAME => "REFERENCE_DB",      

                   
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

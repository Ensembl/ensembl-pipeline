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
		 # to run exonerate to compare two sets of cDNAs:
		 
		 # module to be called, wrapped around Exonerate.pm
		 COMP_EXONERATE_RUNNABLE     => 'Bio::EnsEMBL::Pipeline::GeneComparison::cDNA_Comparison',
		 
		 # where the cDNA chunks are, and how many there are
		 COMP_cDNA_CHUNKDIR               => '/ecs2/work1/eae/japan_cdnas/chunks',
		 COMP_cDNA_CHUNKNUMBER            => 3000,
		 

		 # transcript set to compare to the cdnas
		 COMP_TRANSCRIPTS            => '/ecs2/work1/eae/japan_cdnas/ensembl_genes/ensembl_transcripts',

		 # file with the cdnas you want to compare with
		 COMP_cDNA_FILE              => '/ecs2/work1/eae/japan_cdnas/japan_human_cDNAs.fa',
		 # (japan_human_cDNAs.fa contains refseqs as well)
                 
		 # where the results will go
		 COMP_TMPDIR                 => '/ecs2/work1/eae/NCBI_30/cdnas/jobs_dir/cdna_ensembl_mapping_TEST',
		 # where the job-file with the bsub lines will go:
		 COMP_EXONERATE_BSUBS         => '/ecs2/work1/eae/NCBI_30/cdnas/jobs_dir/cdna_ensembl_mapping_TEST/compare_exonerate_cDNA.jobs',

		 # queue where the jobs will be run
		 COMP_QUEUE                   => 'acari',
		 
		 # path to ensembl-pipeline/scripts/GeneComparison/
		 COMP_SCRIPTDIR               => '/nfs/acari/eae/ensembl-branch-121/ensembl-pipeline/scripts/GeneComparison/',
		 
		 

		 # mouse_5.3.1
		 # annotation/benchmark genes
		 #DBHOST1    => "ecs1b",
	         #DBNAME1    => "mouse_whitehead_0401_rikens",
		 #PATH1      => "CHR",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["genomewise"],
		 
		 # prediction genes
		 #DBHOST2    => 'ecs1e',                    
		 #DBNAME2    => 'mouse_whitehead_0401_denormalised',
		 #PATH2      => 'CHR',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["ensembl"], 
	    
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST => 'ecs1e',                    
		 #REF_DBNAME => 'mouse_whitehead_0401_denormalised',
		 #REF_PATH   => 'CHR',
		 #REF_DBUSER => 'ensro',

		 # chr20 NCBI_28

		 # annotation/benchmark genes
		 #DBHOST1    => "ecs1d",
	         #DBNAME1    => "homo_sapiens_sanger_4_28",
		 #PATH1      => "NCBI_28",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["HUMACE-Known","HUMACE-Novel_CDS"],
		 
		 # prediction genes
		 #DBHOST2    => 'ecs1d',                    
		 #DBNAME2    => 'homo_sapiens_core_4_28',
		 #PATH2      => 'NCBI_28',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["ensembl"], 

		 # prediction genes
		 #DBHOST2    => 'ecs1f',                    
		 #DBNAME2    => 'ens_NCBI_28_est',
		 #PATH2      => 'NCBI_28',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["genomewise"], 
		 
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST => 'ecs1d',                    
		 #REF_DBNAME => 'homo_sapiens_core_4_28',
		 #REF_PATH   => 'NCBI_28',
		 #REF_DBUSER => 'ensro',
		 
		 ############################################
		 # comparison between mouse+riken and mouse #
		 ############################################
		 # annotation/benchmark genes
		 #DBHOST1    => 'ecs1e',
	         #DBNAME1    => 'mouse_5_3_riken_genebuild',
		 #PATH1      => "CHR",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["Riken_ensembl"],

		 # prediction genes
		 #DBHOST2    => 'ecs1e',                    
		 #DBNAME2    => 'mouse_whitehead_0401_denormalised',
		 #PATH2      => 'CHR',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["ensembl"], 

		 #####################################################
		 # comparison between mouse+riken and riken-estgenes #
		 #####################################################

		 # annotation/benchmark genes
		 #DBHOST1    => 'ecs1e',
		 #DBNAME1    => 'mouse_5_3_riken_genebuild',
		 #PATH1      => "CHR",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["Riken_ensembl"],

		 # prediction genes
		 #DBHOST2    => 'ecs1b',                    
		 #DBNAME2    => 'mouse_whitehead_0401_rikens',
		 #PATH2      => 'CHR',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["genomewise"], 
 
		 #####################################################
		 # comparison between mouse and riken-estgenes #
		 #####################################################

		 # annotation/benchmark genes
		 #DBHOST1    => 'ecs1e',
	         #DBNAME1    => 'mouse_whitehead_0401_denormalised',
		 #PATH1      => "CHR",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["ensembl"],

		 # prediction genes
		 #DBHOST2    => 'ecs1b',                    
		 #DBNAME2    => 'mouse_whitehead_0401_rikens',
		 #PATH2      => 'CHR',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["genomewise"], 
		 
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST => 'ecs1e',                    
		 #REF_DBNAME => 'mouse_whitehead_0401_denormalised',
		 #REF_PATH   => 'CHR',
		 #REF_DBUSER => 'ensro',
		 
		 #################################################################
		 # comparison between ensembl genes NCBI_28 and est_genes NCBI_28
		 #################################################################

		 # annotation/benchmark genes
		 #DBHOST1    => 'ecs1d',
	         #DBNAME1    => 'homo_sapiens_sanger_6_28',
		 #PATH1      => "NCBI_28",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["HUMACE-Known","HUMACE-Novel_CDS"],

		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST => 'ecs1d',                    
		 #REF_DBNAME => 'homo_sapiens_sanger_6_28',
		 #REF_PATH   => 'NCBI_28',
		 #REF_DBUSER => 'ensro',

		 # prediction genes
		 #DBHOST2    => 'ecs1e',                    
		 #DBNAME2    => 'ens_NCBI_28',
		 #PATH2      => 'NCBI_28',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["ensembl"], 
		 
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST2 => 'ecs1d',                    
		 #REF_DBNAME2 => 'homo_sapiens_core_6_28',
		 #REF_PATH2   => 'NCBI_28',
		 #REF_DBUSER2 => 'ensro',
		 
		 


		 #################################################################
		 # comparison between ensembl genes NCBI_29 and human curation
		 #################################################################

		 # annotation/benchmark genes
		 #DBHOST1    => 'ecs2d',
	         #DBNAME1    => 'homo_sapiens_sanger_6_29',
		 #PATH1      => "NCBI_29",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["HUMACE-Known","HUMACE-Novel_CDS"],

		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST1 => 'ecs2d',                    
		 #REF_DBNAME1 => 'homo_sapiens_sanger_6_29',
		 #REF_PATH1   => 'NCBI_29',
		 #REF_DBUSER1 => 'ensro',

		 # prediction genes
		 #DBHOST2    => 'ecs1a',                    
		 #DBNAME2    => 'ens_NCBI_29',
		 #PATH2      => 'NCBI_29',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["ensembl"], 
		 
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST2 => 'ecs1e',                    
		 #REF_DBNAME2 => 'NCBI29_raw',
		 #REF_PATH2   => 'NCBI_29',
		 #REF_DBUSER2 => 'ensro',


		 ###############################################
		 # comparison between cDNA+EST and annotation  #
		 ###############################################

		 # annotation/benchmark genes
		 #DBHOST1    => 'ecs1a',                    
		 #DBHOST1    => 'ecs1d',                    
		 #DBNAME1    => 'homo_sapiens_sanger_6_29',
		 #DBNAME1    => 'ens_NCBI_29',
		 #PATH1      => 'NCBI_29',    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["ensembl"], 
		 #GENETYPES1 => ["HUMACE-Known","HUMACE-Novel_CDS"],
		 
		 # prediction genes
		 #DBHOST2    => 'ecs1e',                    
		 #DBNAME2    => 'ens_NCBI_29_est2',
		 #PATH2      => 'NCBI_29',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["genomewise_MASS"], 
		 
		 
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST1 => 'ecs1a',                    
		 #REF_DBNAME1 => 'ens_NCBI_29',
		 #REF_PATH1   => 'NCBI_29',
		 #REF_DBUSER1 => 'ensro',

		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST2 => 'ecs1a',                    
		 #REF_DBNAME2 => 'ens_NCBI_29',
		 #REF_PATH2   => 'NCBI_29',
		 #REF_DBUSER2 => 'ensro',

		 #################################################################
		 # comparison between ensembl genes NCBI_29 and est_genes NCBI_29
		 #################################################################

		 # annotation/benchmark genes
		 #DBHOST1    => 'ecs1a',
	         #DBNAME1    => 'ens_NCBI_29',
		 #PATH1      => "NCBI_29",    
		 #DBUSER1    => "ensro",
		 #GENETYPES1 => ["ensembl"],

		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST1 => 'ecs1a',                    
		 #REF_DBNAME1 => 'ens_NCBI_29',
		 #REF_PATH1   => 'NCBI_29',
		 #REF_DBUSER1 => 'ensro',

		 # prediction genes
		 #DBHOST2    => 'ecs1a',                    
		 #DBNAME2    => 'ens_NCBI_29_cdna',
		 #PATH2      => 'NCBI_29',
		 #DBUSER2    => 'ensro',
		 #GENETYPES2 => ["genomewise"], 
		 
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 #REF_DBHOST2 => 'ecs1e',                    
		 #REF_DBNAME2 => 'NCBI29_raw',
		 #REF_PATH2   => 'NCBI_29',
		 #REF_DBUSER2 => 'ensro',
		 
		 #################################################################
		 # comparison between ensembl genes NCBI_30 and human annotation #
		 #################################################################

		 # annotation/benchmark genes
		 DBHOST1    => 'ecs2a',
	         DBNAME1    => 'otter_merged_chrs_with_anal',
		 PATH1      => 'SANGER',    
                 PORT1      => 3306,
		 DBUSER1    => 'ensro',
		 GENETYPES1 => ["Known","Novel_CDS"],

		 # reference database ( one with common sequence, golden path, contig, etc... )
		 REF_DBHOST1 => 'ecs2a',                    
		 REF_DBNAME1 => 'otter_merged_chrs_with_anal',
		 REF_PATH1   => 'SANGER',
                 REF_PORT1   => 3306,
		 REF_DBUSER1 => 'ensro',

		 # prediction genes
		 DBHOST2    => 'ecs2c',                    
		 DBNAME2    => 'otter_chr20_old_annotation',
		 PATH2      => 'VEGA',
                 PORT2      => 19322,
		 DBUSER2    => 'ensro',
		 GENETYPES2 => ["HUMACE-Known","HUMACE-Novel_CDS"],
		 
		 # reference database ( one with common sequence, golden path, contig, etc... )
		 REF_DBHOST2 => 'ecs2c',                    
		 REF_DBNAME2 => 'otter_chr20_old_annotation',
                 REF_PORT2   => 19322,
		 REF_PATH2   => 'VEGA',
		 REF_DBUSER2 => 'ensro',

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

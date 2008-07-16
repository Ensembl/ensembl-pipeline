package MitConf;

use strict;
use vars qw( %MitConf );


%MitConf = (
	    #location of genbank file
	    MIT_GENBANK_FILE => '/nfs/acari/sw4/cvs_checkout/ensembl-personal/sw4/Builds/Guineapig/Data/NC_000884.gb',
	    # database to put sequnece and genes into
	    MIT_DBNAME => 'sw4_cavia_core',
	    MIT_DBHOST => 'genebuild5',
	    MIT_DBUSER => 'ensadmin',
	    MIT_DBPASS => 'ensembl',
	    MIT_DBPORT => '3306',
	    # logic name of analysis object to be assigned to coding genes
	    # Non coding genes are assigned ncRNA
	    MIT_LOGIC_NAME => 'ensembl',
	    # if want the verbose output and sequence output after the load
	    MIT_DEBUG => '',
	    # Name of the mitochondrial chromosome
	    MIT_NAME => 'MT',
	    # Names of the sequences for the coord systems
	    # Will try and parse this information out of the 
	    # Genbank file if left blank
	    MIT_CONTIG_SEQNAME => '',
	    MIT_CHROMOSOME_SEQNAME => '',
	    MIT_SUPERCONTIG_SEQNAME => '',
	    MIT_CLONE_SEQNAME => '',
	    # Name of top level in coord system
	    MIT_TOPLEVEL => 'scaffold',
	    # Different oganisms use different mitochondial codons:
	    # Vertebrate Mitochondrial (2)
	    # Yeast Mitochondrial (3)
	    # Mold, Protozoan, Coelenterate Mito. and Myco/Spiroplasma (4)
	    # Invertebrate Mitochondrial (5)
	    # Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear (6)
	    # Echinoderm Mitochondrial (9)
	    # Euploid Nuclear (10)
	    # Bacterial (11)
	    # Alternative Yeast Nuclear (12)
	    # Ascidian Mitochondrial (13)
	    # Flatworm Mitochondrial (14)
	    # Blepharisma Macronuclear (15)
	    # Chlorophycean Mitochondrial (16)
	    # Trematode Mitochondrial (21)
	    MIT_CODON_TABLE  =>  '2',
	    # Types
	    MIT_GENE_TYPE => 'protein_coding',
	    MIT_TRNA_TYPE => 'Mt_tRNA',
	    MIT_RRNA_TYPE => 'Mt_rRNA',
	   );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %MitConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $MitConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$MitConf{ $_ };
	} else {
	    die "Error: MitConf: $_ not known\n";
	}
    }
}

1;



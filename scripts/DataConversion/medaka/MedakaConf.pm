package MedakaConf;

use strict;
use vars qw( %MedakaConf );


%MedakaConf = (
	    #location of gff file
	    MED_GFF_FILE => '/nfs/acari/ba1/Projects/medaka1/play.gff',
	    #'/ecs4/work6/searle/medaka1/HdrR_200510/annotation/medaka200510_predictedgene_01.gff',
	    
	    # database to put genes into
	    MED_DBNAME => '', #ba1_medaka_1_gff
	    MED_DBHOST => '', #ecs2
	    MED_DBUSER => '', #ensro
	    MED_DBPASS => '',
	    MED_DBPORT => '', #3362
	    
	    MED_LOGIC_NAME => 'ensembl', 
	    MED_GENE_TYPE => 'protein_coding',
	   );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %MedakaConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $MedakaConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$MedakaConf{ $_ };
	} else {
	    die "Error: MitConf: $_ not known\n";
	}
    }
}

1;



package WormBaseConf;

use strict;
use vars qw( %WormBaseConf );


%WormBaseConf = (
		 #location of agp and gff files
		 WB_CHR_INFO => [
				 { 
				  chr_name => 'I',
				  agp_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_I.agp',
				  length => '15080471',
				  gff_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_I.gff',
				 },
				 { 
				  chr_name => 'II',
				  agp_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_II.agp',
				  length => '15279300',
				  gff_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_II.gff',
				 },
				 { 
				  chr_name => 'III',
				  agp_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_III.agp',
				  length => '13783262',
				  gff_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_III.gff',
				 },
				 { 
				  chr_name => 'IV',
				  agp_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_IV.agp',
				  length => '17493790',
				  gff_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_IV.gff',
	#			  gff_file => '/ecs2/work1/lec/code/elegans_conversion/test.gff',
				 },
				 { 
				  chr_name => 'V',
				  agp_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_V.agp',
				  length => '20916335',
				  gff_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_V.gff',
				 },
				 { 
				  chr_name => 'X',
				  agp_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_X.agp',
				  length => '17705013',
				  gff_file => '/ecs2/work1/lec/code/elegans_conversion/CHROMOSOME_X.gff',
				 }
			       ],
		 WB_AGP_TYPE => 'elegans_91',
		 
		 # database to put sequnece and genes into
		 WB_DBNAME => 'elegans_maintrunk',
		 WB_DBHOST => 'ecs1d',
		 WB_DBUSER => 'ecs1dadmin',
		 WB_DBPASS => 'TyhRv',
		 #path to index based on clone name
		 WB_CLONE_INDEX => ['/ecs2/work1/lec/code/elegans_conversion/wormbase_clones'],
		 WB_GAPS => [{
			      -chr_name => '',
			      -start => '',
			      -end => '',
			     }
			    ],
		 WB_LOGIC_NAME => 'wormbase',
		 WB_DEBUG => 1,
		);

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of GeneConf:
    my @vars = @_ ? @_ : keys( %WormBaseConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $WormBaseConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$WormBaseConf{ $_ };
	} else {
	    die "Error: WormBaseConf: $_ not known\n";
	}
    }
}

1;

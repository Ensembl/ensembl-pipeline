package WormBaseConf.pm

use strict;
use vars qw( %WormBaseConf );


%WormBaseConf = (
		 #location of agp and gff files
		 WB_CHR_INFO => [
				 { 
				  chr_name => 'I',
				  agp_file => '',
				  length => '15080471',
				  gff_file => '',
				 }
				 { 
				  chr_name => 'II',
				  agp_file => '',
				  length => '15279300',
				  gff_file => '',
				 }
				 { 
				  chr_name => 'III',
				  agp_file => '',
				  length => '13783262',
				  gff_file => '',
				 }
				 { 
				  chr_name => 'IV',
				  agp_file => '',
				  length => '17493790',
				  gff_file => '',
				 }
				 { 
				  chr_name => 'V',
				  agp_file => '',
				  length => '20916335',
				  gff_file => '',
				 }
				 { 
				  chr_name => 'X',
				  agp_file => '',
				  length => '17705013',
				  gff_file => '',
				 }
			       ],
		 WB_AGP_TYPE => '',
		 
		 # database to put sequnece and genes into
		 WB_DBNAME => '',
		 WB_DBHOST => '',
		 WB_DBUSER => '',
		 WB_DBPASS => '',
		 #path to index based on clone name
		 WB_CLONE_INDEX => [''];
		 
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
	if ( defined $GeneConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$GeneConf{ $_ };
	} else {
	    die "Error: GeneConf: $_ not known\n";
	}
    }
}

1;

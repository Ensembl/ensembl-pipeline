# Author: Emmanuel Mongin
# Creation: 03.19.2001


=head1 Run_protein_RunnableDB

=head2 Description

This script will submit run_protein_RunnableDB

=cut

use strict;

BEGIN {
    require "Bio/EnsEMBL/Pipeline/Prot_analysis_Conf.pl";
}

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;

my %db_conf =  %::db_conf;

my %scripts_conf = %::scripts_conf;



my $prosite_chunk     = $scripts_conf{'prosite_chunk'};
my $profile_chunk     = $scripts_conf{'profile_chunk'};
my $pfam_chunk        = $scripts_conf{'pfam_chunk'};
my $prints_chunk      = $scripts_conf{'prints_chunk'};
my $scanprosite_chunk = $scripts_conf{'scanprosite_chunk'};
my $tmhmm_chunk       = $scripts_conf{'tmhmm_chunk'};
my $coils_chunk       = $scripts_conf{'coils_chunk'};
my $signalp_chunk     = $scripts_conf{'signalp_chunk'};
my $seg_chunk         = $scripts_conf{'seg_chunk'};

push(@2berun,$prosite_chunk);
push(@2berun,$profile_chunk);
push(@2berun,$pfam_chunk);
push(@2berun,$prints_chunk);
push(@2berun,$scanprosite_chunk);
push(@2berun,$tmhmm_chunk);
push(@2berun,$coils_chunk);
push(@2berun,$signalp_chunk);
push(@2berun,$seg_chunk);



my @2berun;

foreach my $2run(@2berun) {
    if ($2run) {
	#This 

    }
}

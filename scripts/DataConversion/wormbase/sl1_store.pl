#!/usr/local/ensembl/bin/perl  -w

use strict;
use WormBase;
use WormBaseConf;
use Clone2Acc;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Chromosome;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::RawContig;

$| = 1;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $WB_DBHOST,
					    -user => $WB_DBUSER,
					    -dbname => $WB_DBNAME,
					    -pass  => $WB_DBPASS,
					   );


# adding assembly type to meta table


my $analysis_adaptor = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($WB_LOGIC_NAME);
my $operon_analysis = $analysis_adaptor->fetch_by_logic_name($WB_SL1_LOGIC_NAME);

foreach my $chromosome_info(@{$WB_CHR_INFO}) {

  print "handling ".$chromosome_info->{'chr_name'}." with files ".$chromosome_info->{'agp_file'}." and ".$chromosome_info->{'gff_file'}."\n" if($WB_DEBUG);
  

 my $chr = $db->get_SliceAdaptor->fetch_by_chr_start_end($chromosome_info->{'chr_name'}, 1, ($chromosome_info->{'length'}+1));
 
  
  my @operons = @{&parse_SL1($chromosome_info->{'gff_file'}, $chr, $operon_analysis)};
  my $non_transforming = &write_simple_features(\@operons, $db);


  
}



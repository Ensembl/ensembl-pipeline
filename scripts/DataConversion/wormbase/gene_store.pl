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

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host => $WB_DBHOST,
   -user => $WB_DBUSER,
   -dbname => $WB_DBNAME,
   -pass  => $WB_DBPASS,
   -port  => $WB_DBPORT,
  );


# adding assembly type to meta table


my $analysis_adaptor = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($WB_LOGIC_NAME);

foreach my $chromosome_info(@{$WB_CHR_INFO}) {

  print "handling ".$chromosome_info->{'chr_name'}." with files ".$chromosome_info->{'agp_file'}." and ".$chromosome_info->{'gff_file'}."\n" if($WB_DEBUG);
  

  my $chr = $db->get_SliceAdaptor->fetch_by_region('Chromosome', $chromosome_info->{'chr_name'}, 1, ($chromosome_info->{'length'}, 1, $WB_NEW_COORD_SYSTEM_VERSION));
 
  my $genes = &parse_gff($chromosome_info->{'gff_file'}, $chr, $analysis);
  &write_genes($genes, $db);

 
  
  

  my @genes = @{$chr->get_all_Genes};
  open(TRANSLATE, "+>>".$WB_NON_TRANSLATE) or die "couldn't open ".$WB_NON_TRANSLATE." $!";
  TRANSLATION: foreach my $gene(@genes){
      my $translation = &translation_check($gene,$db);
      if($translation){
	next TRANSLATION;
      }else{
 	print TRANSLATE $gene->stable_id." from ".$chromosome_info->{'chr_name'}." doesn't translate\n";
	next TRANSLATION;
      }
    }
  close(TRANSLATE);

}

#!/usr/local/ensembl/bin/perl  -w

use strict;
use WormBase;
use WormBaseConf;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


$| = 1;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $WB_DBHOST,
					    -user => $WB_DBUSER,
					    -dbname => $WB_DBNAME,
					    -pass  => $WB_DBPASS,
					    -port => $WB_DBPORT,
					   );
my $analysis_adaptor = $db->get_AnalysisAdaptor;
my $analysis = $analysis_adaptor->fetch_by_logic_name($WB_PSEUDO_LOGIC_NAME);

foreach my $chromosome_info(@{$WB_CHR_INFO}) {

  print "handling ".$chromosome_info->{'chr_name'}." with files ".$chromosome_info->{'agp_file'}." and ".$chromosome_info->{'gff_file'}."\n" if($WB_DEBUG);
  

 my $chr = $db->get_SliceAdaptor->fetch_by_region('Chromosome', $chromosome_info->{'chr_name'}, 1, ($chromosome_info->{'length'}, 1, $WB_NEW_COORD_SYSTEM_VERSION));

  my $genes = &parse_pseudo_gff($chromosome_info->{'gff_file'}, $chr, $analysis);
  my $non_transforming =  &write_genes($genes, $db, 1);

 

  my @genes = @{$chr->get_all_Genes_by_type($analysis->logic_name)};
  print STDERR "have ".@genes." pseudogenes\n";
}

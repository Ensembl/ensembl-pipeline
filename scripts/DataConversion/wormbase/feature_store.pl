#!/usr/local/ensembl/bin/perl  -w

use strict;
use WormBase;
use WormBaseConf;
use Clone2Acc;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
$| = 1;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $WB_DBHOST,
					    -user => $WB_DBUSER,
					    -dbname => $WB_DBNAME,
					    -pass  => $WB_DBPASS,
					    -port => $WB_DBPORT,
					   );

my $analysis_logic_name = shift;
my $analysis_adaptor = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($analysis_logic_name);
if(!$analysis){
  die "COULD NOT GET ANALYSIS ADAPTOR FOR $analysis_logic_name";
}

foreach my $chromosome_info(@{$WB_CHR_INFO}) {

  if($chromosome_info->{'chr_name'} && $chromosome_info->{'agp_file'} && $chromosome_info->{'gff_file'}){
    print "handling ".$chromosome_info->{'chr_name'}." with files ".$chromosome_info->{'agp_file'}.
          " and ".$chromosome_info->{'gff_file'}."\n" if($WB_DEBUG);

    my $chr = $db->get_SliceAdaptor->fetch_by_region('Chromosome', $chromosome_info->{'chr_name'}, 1,
						     ($chromosome_info->{'length'}, 1, $WB_NEW_COORD_SYSTEM_VERSION));
    
    #MtDNA does not have Expr_profile, RNAi or Operon so may cause hassles for these
    next if (($chromosome_info->{'chr_name'} eq 'MtDNA') 
            && ($analysis_logic_name eq "Expression_profile" 
	    or $analysis_logic_name eq "RNAi" 
	    or $analysis_logic_name eq "Operon"));    
    
    #read the features from gff file
    if($analysis_logic_name eq "Expression_profile" or $analysis_logic_name eq "RNAi" or
       $analysis_logic_name eq "Operon"){

      my $features;
      if($analysis_logic_name eq "Expression_profile"){
	$features = parse_expr($WB_workDIR."".$chromosome_info->{'gff_file'}, $chr, $analysis);
      }
      elsif($analysis_logic_name eq "RNAi"){
	$features = parse_rnai($WB_workDIR."".$chromosome_info->{'gff_file'}, $chr, $analysis);
      }
      elsif($analysis_logic_name eq "Operon"){
	$features = parse_operons($WB_workDIR."".$chromosome_info->{'gff_file'}, $chr, $analysis);
      }
      print "have ".scalar @$features." features ($analysis_logic_name).\n" if($WB_DEBUG);

      #save the features to db
      &write_simple_features($features, $db);
    }
    elsif($analysis_logic_name eq "Pseudogene"){
      #parse out pseudogenes
      my $pseudogenes = parse_pseudo_gff($WB_workDIR."".$chromosome_info->{'gff_file'}, $chr, $analysis);
      print "have ".scalar @$pseudogenes." pseudogenes / tRNA genes.\n" if($WB_DEBUG);

      #store pseudogenes
      &write_genes($pseudogenes, $db, 1);
    }
    elsif($analysis_logic_name eq "tRNA"){
      #parse out tRNA-genes
      #there only seem to be tRNA in MtDNA (tRNAscan in CHROMOSOMES)
      my $tRNAgenes = parse_tRNA_genes($WB_workDIR."".$chromosome_info->{'gff_file'}, $chr, $analysis);
      print "have ".scalar @$tRNAgenes." tRNA genes.\n" if($WB_DEBUG);

      #store tRNA-genes
      &write_genes($tRNAgenes, $db);
    }
    elsif($analysis_logic_name eq "rRNA"){
      #parse out rRNA-genes
      my $rRNAgenes = parse_rRNA_genes($WB_workDIR."".$chromosome_info->{'gff_file'}, $chr, $analysis);
      print "have ".scalar @$rRNAgenes." rRNA genes.\n" if($WB_DEBUG);

      #store tRNA-genes
      &write_genes($rRNAgenes, $db);
    }       
    elsif($analysis_logic_name eq "wormbase"){
      #parse out real worm genes
      my $genes = parse_gff($WB_workDIR."".$chromosome_info->{'gff_file'}, $chr, $analysis);
      print "have ".scalar @$genes." genes.\n" if($WB_DEBUG);

      #store genes
      &write_genes($genes, $db);

      #check translations
      my @genes = @{$chr->get_all_Genes};
      open(TRANSLATE, "+>>".$WB_NON_TRANSLATE) or die "couldn't open ".$WB_NON_TRANSLATE." $!";
      TRANSLATION: foreach my $gene(@genes){
	my $translation = &translation_check($gene);
	if($translation){
	  next TRANSLATION;
	}else{
	  print TRANSLATE $gene->stable_id." from ".$chromosome_info->{'chr_name'}." doesn't translate\n";
	  next TRANSLATION;
	}
      }
      close(TRANSLATE);
    }
    else{
      print "\n ANALYSIS NOT DEFINED ($analysis_logic_name)!\n";
      exit 1;
    }
  }
  else{
    print "skipping chromosome.\n" if($WB_DEBUG);
  }

}

exit 0;

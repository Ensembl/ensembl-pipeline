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
my $meta_sql = "insert into meta(meta_key, meta_value) values(\'assembly.default\', ?)";

my $meta_sth = $db->prepare($meta_sql);
$meta_sth->execute($WB_AGP_TYPE); 

my $analysis_adaptor = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($WB_LOGIC_NAME);
my $operon_analysis = $analysis_adaptor->fetch_by_logic_name($WB_OPERON_LOGIC_NAME) if($WB_OPERON_LOGIC_NAME);
my $rnai_analysis = $analysis_adaptor->fetch_by_logic_name($WB_RNAI_LOGIC_NAME) if($WB_RNAI_LOGIC_NAME);
my $expr_analysis = $analysis_adaptor->fetch_by_logic_name($WB_EXPR_LOGIC_NAME) if($WB_EXPR_LOGIC_NAME);
my $sl1_analysis = $analysis_adaptor->fetch_by_logic_name($WB_SL1_LOGIC_NAME);
my $sl2_analysis = $analysis_adaptor->fetch_by_logic_name($WB_SL2_LOGIC_NAME);
foreach my $chromosome_info(@{$WB_CHR_INFO}) {

  print "handling ".$chromosome_info->{'chr_name'}." with files ".$chromosome_info->{'agp_file'}." and ".$chromosome_info->{'gff_file'}."\n" if($WB_DEBUG);

  my $chromosome = Bio::EnsEMBL::Chromosome->new(-chr_name => $chromosome_info->{'chr_name'},
						 -length => $chromosome_info->{'length'},
						 -adaptor =>$db->get_ChromosomeAdaptor, 
						);

  $db->get_ChromosomeAdaptor->store($chromosome); # adding chromosome to database
  open(FH, $chromosome_info->{'agp_file'}) or die "couldn't open ".$chromosome_info->{'agp_file'}." $!";
  my $fh = \*FH;
  my ($seq_ids, $non_ids) = &get_seq_ids($fh); # getting embl accs from agp
  open(NON, "+>>".$WB_SEQ_IDS) or die "couldn't open ".$WB_SEQ_IDS." $!";
    foreach my $id(@$non_ids){
      print NON $id." doesn't fit format on chromosome ".$chromosome_info->{'chr_name'}."\n";
    }
  close(NON);
  seek($fh, 0, 0); #resetting the fh to the start of the agp file 
  my $pfetch = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();
  
  my $obda = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(-db => $WB_CLONE_INDEX);

  my %seqs = %{&get_sequences_pfetch($seq_ids, $pfetch)};
  my %chr_hash = %{&agp_parse($fh, $chromosome->dbID, $WB_AGP_TYPE)};
  
  foreach my $id(keys(%seqs)){
    my $seq = $seqs{$id};
    if($seq->length < $chr_hash{$id}->{contig_end}){
      my ($acc) = $id =~ /(\S+)\.\d+/;
      my $clone_name = $WB_ACC_2_CLONE->{$acc};
      my $clone_seq = $obda->get_Seq_by_acc($clone_name);
      $clone_seq->id($id);
      $clone_seq->desc('');
      if($clone_seq->length > $chr_hash{$id}->{contig_end}){
	warn "sequence ".$clone_name." ".$id." id the wrong length in both the embl accession and the wormbase data not much i can do\n";
	delete($seqs{$id});
      }else{
	$seqs{$id} = $clone_seq;
      }
    }
  }
  my %contig_id;
  foreach my $id(keys(%seqs)){
    my $seq = $seqs{$id};
    my ($acc, $version) = $id =~ /(\S+)\.(\d+)/;
    my $clone_name = $WB_ACC_2_CLONE->{$acc};
    my $contig_id = $id.".1.".$seq->length;
    my $time = time;
    my $contig = &make_Contig($contig_id, $seq->seq, $seq->length);
    my $clone = &make_Clone($clone_name, 1, $acc, $version, 3, $contig, $time, $time);
    eval{
      $db->get_CloneAdaptor->store($clone);
    };
    if($@){
      die("couldn't store ".$clone->id." ".$clone->embl_id." $@");
    }
    $contig_id{$id} = $contig->dbID;
    if(!$WB_RAW_CONTIGS){
      my $sql = "insert into input_id_analysis(input_id, analysis_id,  created) values('$contig_id', $WB_SUBMIT_CONTIG_ID, now())";
      my $sth = $db->prepare($sql);
      $sth->execute($contig_id, $WB_SUBMIT_CONTIG_ID);
    }
  }

  foreach my $id(keys(%chr_hash)){
    my $agp = $chr_hash{$id};
    my $contig = $contig_id{$id};

    &insert_agp_line($agp->{'chromosome_id'}, $agp->{'chr_start'}, $agp->{'chr_end'}, $agp->{'superctg_name'}, $agp->{'superctg_start'}, $agp->{'superctg_end'}, $agp->{'superctg_ori'}, $contig, $agp->{'contig_start'}, $agp->{'contig_end'}, $agp->{'contig_ori'}, $agp->{'type'}, $db);
  }

  my $chr = $db->get_SliceAdaptor->fetch_by_chr_start_end($chromosome_info->{chr_name}, 1, ($chromosome_info->{length} - 1));
  my $genes = &parse_gff($chromosome_info->{'gff_file'}, $chr, $analysis);
  my $non_transforming =  &write_genes($genes, $db);

  open(TRANSFORM, "+>>".$WB_NON_TRANSFORM) or die "couldn't open ".$WB_NON_TRANSFORM." $!";
  foreach my $id(keys(%$non_transforming)){
    print TRANSFORM $id." gene wouldn't transform on chromsome ".$chromosome_info->{'chr_name'}."\n";
    
  }
  my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end($chromosome_info->{chr_name}, 1, ($chromosome_info->{length} - 1));

  my @genes = @{$slice->get_all_Genes};
  open(TRANSLATE, "+>>".$WB_NON_TRANSLATE) or die "couldn't open ".$WB_NON_TRANSLATE." $!";
  TRANSLATION: foreach my $gene(@genes){
      my $translation = &translation_check($gene);
      if($translation){
	next TRANSLATION;
      }else{
	#my ($clone_name) = $gene->stable_id =~ /(\S+)\.\d+/;
#	$gene->transform;
#	my @exons = @{$gene->get_all_Exons};
#	my $old_contig = $exons[0]->contig;
#	my $clone_seq = $obda->get_Seq_by_acc($clone_name);
#	my ($embl_acc) = $old_contig->name =~ /(\S+\.\d+)\.\d+\.\d+/;
#	my ($version) = $embl_acc =~ /\S+\.(\d+)/;
#	my $contig_id = $clone_name.".".$version.".1.".$clone_seq->length;
#	my $time = time;
#	my $contig = &make_Contig($contig_id, $clone_seq->seq, $clone_seq->length);
#	my $clone = &make_Clone($clone_name, 2, $embl_acc, $version, 3, $contig, $time, $time);
#	eval{
#	  $db->get_CloneAdaptor->store($clone);
#	};
#	if($@){
#	  die("couldn't store ".$clone->id." ".$clone->embl_id." $!");
#	}
#	my $assembly_sql = "update assembly set contig_id = ? where contig_id = ?";
#	my $assembly_sth = $db->prepare($assembly_sql);
#	$assembly_sth->execute($contig->dbID, $old_contig->dbID);
#	my $exon_sql = "update exon set contig_id = ? where contig_id = ?";
#	my $exon_sth = $db->prepare($exon_sql);
#	$exon_sth->execute($contig->dbID, $old_contig->dbID);
#	my $new_seq_gene = $db->get_GeneAdaptor->fetch_by_stable_id($gene->stable_id);
#	my $checked_gene = &translation_check($new_seq_gene);
#	if(!$checked_gene){
#	  print TRANSLATE "gene ".$new_seq_gene->stable_id." doesn't translate on either embl sequence ".$old_contig->name." or wormbase sequence ".$contig->name." on chromosome ".$chromosome_info->{'chr_name'}."\n";
      
      
	print TRANSLATE $gene->stable_id." from ".$chromosome_info->{'chr_name'}." doesn't translate\n";
	next TRANSLATION;
      }
    }
  close(TRANSLATE);
  
  if($operon_analysis){
    my @operons = @{&parse_operons($chromosome_info->{'gff_file'}, $chr, $operon_analysis)};
    $non_transforming = &write_operons(\@operons, $db);
  }
  if($rnai_analysis){
    my @operons = @{&parse_rnai($chromosome_info->{'gff_file'}, $chr, $rnai_analysis)};
    $non_transforming = &write_simple_features(\@operons, $db);
  }
  if($expr_analysis){
    my @operons = @{&parse_rnai($chromosome_info->{'gff_file'}, $chr, $expr_analysis)};
    $non_transforming = &write_simple_features(\@operons, $db);
  }
  if($sl1_analysis){
    my @operons = @{&parse_rnai($chromosome_info->{'gff_file'}, $chr, $sl1_analysis)};
    $non_transforming = &write_simple_features(\@operons, $db);
  }
  if($sl2_analysis){
    my @operons = @{&parse_rnai($chromosome_info->{'gff_file'}, $chr, $sl2_analysis)};
    $non_transforming = &write_simple_features(\@operons, $db);
  }

  close($fh);
}



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


my $meta_sql = "insert into meta(meta_key, meta_value) values(\'assembly.default\', ?)";

my $meta_sth = $db->prepare($meta_sql);
$meta_sth->execute($WB_AGP_TYPE); 

my $analysis_adaptor = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name($WB_LOGIC_NAME);

foreach my $chromosome_info(@{$WB_CHR_INFO}){

  my $chromosome = Bio::EnsEMBL::Chromosome->new(-chr_name => $chromosome_info->{'chr_name'},
						 -length => $chromosome_info->{'length'},
						 -adaptor =>$db->get_ChromosomeAdaptor, 
						);

  $db->get_ChromosomeAdaptor->store($chromosome);
  open(FH, $chromosome_info->{'agp_file'});
  my $fh = \*FH;
  my $seq_ids = &get_seq_ids($fh);
  
  seek($fh, 0, 0);
  my $pfetch = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();
  
  my $obda = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(-db => $WB_CLONE_INDEX);

  my %seqs = %{&get_sequences_pfetch($seq_ids, $pfetch)};
  my %chr_hash = %{&agp_parse($fh, $chromosome->dbID, $WB_AGP_TYPE)};
  
  foreach my $id(keys(%seqs)){
    my $seq = $seqs{$id};
   # print STDERR "have sequence ".$seq->length." and contig end ".$chr_hash{$id}->{contig_end}."\n";
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
    my $clone     = new Bio::EnsEMBL::Clone;
    my $contig    = new Bio::EnsEMBL::RawContig; 
    my $clone_name = $WB_ACC_2_CLONE->{$acc};
    $clone->htg_phase(3);
    $clone->id($clone_name); 
    $clone->embl_id($seq->id);
    $clone->version(1);
    $clone->embl_version($version);
    my $contig_id = $id.".1.".$seq->length;
    $contig->name($contig_id);
    $contig->seq($seq->seq);
    $contig->length($seq->length);
    $contig->adaptor($db->get_RawContigAdaptor);
    $clone->add_Contig($contig);
    my $time = time;
    $clone->created($time);
    $clone->modified($time);
    eval{
      $db->get_CloneAdaptor->store($clone);
    };
    if($@){
      die("couldn't store ".$clone->id." ".$clone->embl_id." $@");
    }
    $contig_id{$id} = $contig->dbID;
  }

  foreach my $id(keys(%chr_hash)){
    my $agp = $chr_hash{$id};
    my $contig = $contig_id{$id};
    my $chr_id = $agp->{'chromosome_id'};
    my $chr_start = $agp->{'chr_start'};
    my $chr_end = $agp->{'chr_end'};
    my $superctg_name = $agp->{'superctg_name'};
    my $superctg_start = $agp->{'superctg_start'};
    my $superctg_end = $agp->{'superctg_end'};
    my $superctg_ori = $agp->{'superctg_ori'};
    my $contig_start = $agp->{'contig_start'};
    my $contig_end = $agp->{'contig_end'};
    my $contig_ori = $agp->{'contig_ori'};
    my $type = $agp->{'type'};

    my $sql = "insert into assembly(chromosome_id, chr_start, chr_end, superctg_name, superctg_start, superctg_end, superctg_ori, contig_id, contig_start, contig_end, contig_ori, type) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
    my $sth = $db->prepare($sql);
    $sth->execute($chr_id, $chr_start, $chr_end, $superctg_name, $superctg_start, $superctg_end, $superctg_ori, $contig, $contig_start, $contig_end, $contig_ori, $type); 
  }

  my $chr = $db->get_SliceAdaptor->fetch_by_chr_start_end($chromosome_info->{chr_name}, 1, ($chromosome_info->{length} - 1));
  my $genes = &parse_gff($chromosome_info->{'gff_file'}, $chr);
  #my ($non_translating, $non_transforming) = &write_genes($genes, $db);
  &write_genes($genes, $db);
  #print "there are ".keys(%$non_translating)." clones with non translating genes\n";
#  print "there are ".keys(%$non_transforming)." clones with non translating genes\n";

#  foreach my $clone_name(keys(%$non_translating)){
#    print "remaking ".$clone_name." sequence\n";
#    my $wormbase_seq = $obda->get_Seq_by_acc($clone_name);
#    my $clone_adaptor = $db->get_CloneAdaptor;
#    my $old_clone = $clone_adaptor->fetch_by_name($clone_name);
#    my ($old_contig) = @{$old_clone->get_all_Contigs};
#    my $contig_id = $old_contig->dbID;
#    my $embl_acc = $old_clone->embl_id;
#    my $version = $old_clone->embl_version;
#    $clone_adaptor->remove($old_clone);
#    my $clone     = new Bio::EnsEMBL::Clone;
#    my $contig    = new Bio::EnsEMBL::RawContig; 
#    $clone->htg_phase(3);
#    $clone->id($clone_name); 
#    $clone->embl_id($embl_acc);
#    $clone->version(1);
#    $clone->embl_version($version);
#    my $contig_name = $embl_acc.".1.".$wormbase_seq->length;
#    $contig->name($contig_name);
#    $contig->seq($wormbase_seq->seq);
#    $contig->length($wormbase_seq->length);
#    $contig->adaptor($db->get_RawContigAdaptor);
#    $clone->add_Contig($contig);
#    my $time = time;
#    $clone->created($time);
#    $clone->modified($time);
#    eval{
#      $db->get_CloneAdaptor->store($clone);
#    };
#    if($@){
#      die("couldn't store ".$clone->id." ".$clone->embl_id." $!");
#    }
#    my $assembly_sql = "update assembly set contig_id = ? where contig_id = ?";
#    my $assembly_sth = $db->prepare($assembly_sql);
#    $assembly_sth->execute($contig->dbID, $contig_id);
#    my $gene_adaptor = $db->get_GeneAdaptor;
#    foreach my $gene(@{$non_translating->{$clone_name}}){
#      foreach my $exon(@{$gene->get_all_Exons}){
#	$exon->contig($contig);
#      }
#      print "checking translation of ".$gene->stable_id."\n";
#      my $translated_gene = &translation_check($gene);
#      if($translated_gene){
#	$gene_adaptor->store($translated_gene);
#      }else{
#	print STDERR "gene ".$gene->stable_id." doesn't translate at all something odd going on\n";
#      }
#    }
#  }

  my $slice = $db->get_SliceAdaptor->fetch_by_chr_start_end($chromosome_info->{chr_name}, 1, ($chromosome_info->{length} - 1));

  my @genes = @{$slice->get_all_Genes};

  TRANSLATION: foreach my $gene(@genes){
    my $translation = &translation_check($gene);
    if($translation){
      next TRANSLATION;
    }else{
      #print "gene ".$gene->stable_id." doesn't translate\n";
      #foreach my $transcript(@{$gene->get_all_Transcripts}){
      #	&display_exons(@{$transcript->get_all_Exons});
      #	&non_translate($transcript);
      #      }
      my ($clone_name) = $gene->stable_id =~ /(\S+)\.\d+/;
      #print STDERR "have clone ".$clone_name."\n";
      $gene->transform;
      my @exons = @{$gene->get_all_Exons};
      my $old_contig = $exons[0]->contig;
      my $clone_seq = $obda->get_Seq_by_acc($clone_name);
      my $clone     = new Bio::EnsEMBL::Clone;
      my $contig    = new Bio::EnsEMBL::RawContig; 
      $clone->htg_phase(3);
      $clone->id($clone_name);
      #print "old contig name ".$old_contig->name."\n";
      my ($embl_acc) = $old_contig->name =~ /(\S+\.\d+)\.\d+\.\d+/;
      my ($version) = $embl_acc =~ /\S+\.(\d+)/;
      #print "embl_acc = ".$embl_acc." version ".$version."\n";
      $clone->embl_id($embl_acc);
      $clone->version(1);
      $clone->embl_version($version);
      my $contig_name = $clone_name.".".$version.".1.".$clone_seq->length;
      $contig->name($contig_name);
      $contig->seq($clone_seq->seq);
      $contig->length($clone_seq->length);
      $contig->adaptor($db->get_RawContigAdaptor);
      $clone->add_Contig($contig);
      my $time = time;
      $clone->created($time);
      $clone->modified($time);
      eval{
	$db->get_CloneAdaptor->store($clone);
      };
      if($@){
	die("couldn't store ".$clone->id." ".$clone->embl_id." $!");
      }
      #print "changing contig dbid from ".$old_contig->dbID." to ".$contig->dbID."\n";
      my $assembly_sql = "update assembly set contig_id = ? where contig_id = ?";
      my $assembly_sth = $db->prepare($assembly_sql);
      $assembly_sth->execute($contig->dbID, $old_contig->dbID);
      my $exon_sql = "update exon set contig_id = ? where contig_id = ?";
      my $exon_sth = $db->prepare($exon_sql);
      $exon_sth->execute($contig->dbID, $old_contig->dbID);
      my $new_seq_gene = $db->get_GeneAdaptor->fetch_by_stable_id($gene->stable_id);
      my $checked_gene = &translation_check($new_seq_gene);
      if($checked_gene){
	#print "gene translates hurrah\n";
      }else{
	print $new_seq_gene->stable_id." still doesn't translate bugger\n";
	#$gene->transform;
	#my @exons = @{$new_seq_gene->get_all_Exons};
	#print "have exon ".$exons[0]."\n";
	#print "gene lies on contig ".$exons[0]->contig->name."\n";
	#&display_exons(@exons);
	#my ($transcript) = @{$new_seq_gene->get_all_Transcripts};
	#&non_translate($transcript);
      }
    }
  }
}


sub display_exons{
  my (@exons) = @_;

  @exons = sort{$a->start <=> $b->start || $a->end <=> $b->end} @exons;

  
  foreach my $e(@exons){
       print $e->stable_id."\t ".$e->start."\t ".$e->end."\t ".$e->strand."\t ".$e->phase."\t ".$e->end_phase."\n";
    }
  
}

sub non_translate{
  my (@transcripts) = @_;
  
  foreach my $t(@transcripts){
    
    my @exons = @{$t->get_all_Exons};
#    print "transcript sequence :\n".$t->seq."\n";
    foreach my $e(@exons){
      my $seq = $e->seq;
      my $pep0 = $seq->translate('*', 'X', 0);
      my $pep1 = $seq->translate('*', 'X', 1);
      my $pep2 = $seq->translate('*', 'X', 2);
      print "exon sequence :\n".$e->seq->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 0 frame\n ".$pep0->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 1 phase\n ".$pep2->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 2 phase\n ".$pep1->seq."\n\n";
      print "\n\n";
      
    }
    
  }
}



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
      die("couldn't store ".$clone->id." ".$clone->embl_id." $!");
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
  my ($non_translating, $non_transforming) = &write_genes($genes, $db);

  print "there are ".keys(%$non_translating)." clones with non translating genes\n";
  print "there are ".keys(%$non_transforming)." clones with non translating genes\n";

  foreach my $clone_name(keys(%$non_translating)){
    my $wormbase_seq = $obda->get_Seq_by_acc($clone_name);
    my $clone_adaptor = $db->get_CloneAdaptor;
    my $clone = $clone_adaptor->fetch_by_name($clone_name);
    my ($contig) = @{$clone-get_all_Contigs};
    my $contig_id = $contig->dbID;
    my $embl_acc = $clone->embl_id;
    my $version = $clone->embl_version;
    $clone_adaptor->remove($clone);
    my $clone     = new Bio::EnsEMBL::Clone;
    my $contig    = new Bio::EnsEMBL::RawContig; 
    my $clone_name = $WB_ACC_2_CLONE->{$acc};
    $clone->htg_phase(3);
    $clone->id($clone_name); 
    $clone->embl_id($acc);
    $clone->version(1);
    $clone->embl_version($version);
    my $contig_name = $id.".1.".$wormbase_seq->length;
    $contig->name($contig_name);
    $contig->seq($wormbase_seq->seq);
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
      die("couldn't store ".$clone->id." ".$clone->embl_id." $!");
    }
  }

}

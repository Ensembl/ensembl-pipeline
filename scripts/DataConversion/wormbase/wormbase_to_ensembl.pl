#!/usr/local/ensembl/bin -w

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
foreach my $chromosome_info(@{$WB_CHR_INFO}){

  my $chromosome = Bio::EnsEMBL::Chromosome->new(-chr_name => $chromosome_info->{'chr_name'},
						 -length => $chromosome_info->{'length'},
						);
  $db->get_ChromosomeAdaptor->store($chromosome);
  open(FH, $chromosome_info->{'agp_file'});
  my $fh = \*FH;
  my $seq_ids = &get_seq_ids($fh);
  
  my $pfetch = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new(options => '-A');
  my $obda = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(-db => $WB_CLONE_INDEX);

  my %seqs = %{&get_sequences($seq_ids, $pfetch)};
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

  foreach my $id(keys(%seqs)){
    my $seq = $seqs{$id};
    my $clone     = new Bio::EnsEMBL::Clone;
    my $contig    = new Bio::EnsEMBL::RawContig; 
    my ($version) = $clone_ids[0] =~ /\S+\.(\d+)/;
    my $clone_name = $WB_ACC_2_CLONE->{$acc};
    $clone->htg_phase(3);
    $clone->id($clone_name); 
    $clone->embl_id($seq->id);
    $clone->version(1);
    $clone->embl_version($version);
    my $contig_id = $clone_ids[0].".1.".$seq->length;
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
  }
}

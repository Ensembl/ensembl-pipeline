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
					    -port => $WB_DBPORT,
					   );


# adding assembly type to meta table
my $meta_sql = "insert into meta(meta_key, meta_value) values(\'assembly.default\', ?)";

my $meta_sth = $db->prepare($meta_sql);
$meta_sth->execute($WB_AGP_TYPE); 



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
    if(!$WB_RAW_COMPUTES){
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


  
}

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

my $clone_cs = &store_coord_system($db, $WB_CLONE_SYSTEM_NAME, 
                                   '', 0, 1, 1);
my $chromosome_cs = &store_coord_system($db, $WB_CHROMOSOME_SYSTEM_NAME, 
                                        $WB_AGP_TYPE, 1, 0, 1);


foreach my $chromosome_info(@{$WB_CHR_INFO}) {

  print "handling ".$chromosome_info->{'chr_name'}." with files ".$chromosome_info->{'agp_file'}." and ".$chromosome_info->{'gff_file'}."\n" if($WB_DEBUG);
  
  my $chromosome = &store_slice($db, $chromosome_info->{'chr_name'}, 1, $chromosome_info->{'length'}, 1, $chromosome_cs);

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
  my %chr_hash = %{&agp_parse($fh, 
                              $chromosome->adaptor->get_seq_region_id($chromosome), 
                              $WB_AGP_TYPE)};
  my %contig_id;
  foreach my $id(keys(%seqs)){
    my $seq = $seqs{$id};
    my $strand = $chr_hash{$seq->id}->{'contig_ori'};
    my $contig = &store_slice($db, $id, 1, $seq->length, $strand, $clone_cs, $seq->seq);
    $contig_id{$id} = $contig->adaptor->get_seq_region_id($contig);
  }
  
   foreach my $id(keys(%chr_hash)){
    my $agp = $chr_hash{$id};
    my $contig = $contig_id{$id};

    &insert_agp_line($agp->{'chromosome_id'}, $agp->{'chr_start'}, 
                     $agp->{'chr_end'}, $contig, $agp->{'contig_start'},
                     $agp->{'contig_end'}, $agp->{'contig_ori'}, 
                     $db);
  }
}

my $mc = $db->get_MetaContainer();
$mc->store_key_value('assembly.mapping', 
                     $chromosome_cs->name.":".$chromosome_cs->version."|".
                     $clone_cs->name);


#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use WormBaseConf;


my $old_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ecs1d',
					    -user => 'ensro',
					    -dbname => 'alistair_elegans_newschema',
					    -pass  => '',
					   );

my $new_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ecs1d',
					    -user => 'ecs1dadmin',
					    -dbname => 'elegans_maintrunk',
					    -pass  => 'TyhRv',
					   );



my %new_contigs;
my $sql = 'select name from contig';
my $sth = $new_db->prepare($sql);
$sth->execute;
my $old_rca = $old_db->get_RawContigAdaptor();
my $new_rca = $new_db->get_RawContigAdaptor();
RC: while (my($name) = $sth->fetchrow){
 
  my $old_rc = $old_rca->fetch_by_name($name);

  if(!$old_rc){
    $new_contigs{$name} = 1;
    next RC;
  }

  my $new_rc = $new_rca->fetch_by_name($name);

 my @similarity_features = @{$old_rc->get_all_SimilarityFeatures()};
  my $pfa = $new_db->get_ProteinAlignFeatureAdaptor;
  my $dfa = $new_db->get_DnaAlignFeatureAdaptor;
  foreach my $align(@similarity_features){
    $align->contig($new_rc);
    $align->dbID('');
    if($align->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
      $align->adaptor($pfa);
      $pfa->store($align);
    }elsif($align->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
      $align->adaptor($dfa);
      $dfa->store($align);
    }
  }

  my @prediction_transcripts = @{$old_rc->get_all_PredictionTranscripts};
  my $pta = $new_db->get_PredictionTranscriptAdaptor;
  foreach my $pt(@prediction_transcripts){
    my @exons = @{$pt->get_all_Exons};
    foreach my $e(@exons){
      $e->contig($new_rc);
      $e->adaptor('');
      $e->dbID('');
    }
    $pt->dbID('');
    $pt->adaptor($pta);
    $pta->store($pt);
  }

  my @simple_features = @{$old_rc->get_all_SimpleFeatures};
  my $sfa = $new_db->get_SimpleFeatureAdaptor;
  foreach my $sf(@simple_features){
    $sf->contig($new_rc);
    $sf->dbID('');
    $sf->adaptor($sfa);
    $sfa->store($sf);
  }
  my $count = 0;
  my @repeat_features = @{$old_rc->get_all_RepeatFeatures};
 
  my $rfa = $new_db->get_RepeatFeatureAdaptor;
  foreach my $rf(@repeat_features){
    $rf->contig($new_rc);
    $rf->dbID('');
    $rf->adaptor($rfa);
    $rf->repeat_consensus->dbID('');
    $rfa->store($rf);
    $count++;
  }
 
}


open(FH, ">".$WB_NEW_CONTIGS);
my $time = time;
foreach my $name(keys(%new_contigs)){
  print FH $name."\t".$WB_SUBMIT_CONTIG_ID."\t".$time."\n";
}

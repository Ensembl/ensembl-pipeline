#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use WormBaseConf;


my $old_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ecs1f',
					    -user => 'ensro',
					    -dbname => 'elegans_maintrunk',
					    -pass  => '',
					   );

my $new_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ecs1b',
					    -user => 'ensadmin',
					    -dbname => 'elegans_94',
					    -pass  => 'ensembl',
					   );



my %new_contigs;
my $sql = 'select name from contig';
my $sth = $new_db->prepare($sql);
$sth->execute;
my $old_rca = $old_db->get_RawContigAdaptor();
my $new_rca = $new_db->get_RawContigAdaptor();
RC: while (my($name) = $sth->fetchrow){
  print "sorting contig ".$name."\n";
  my $old_rc = $old_rca->fetch_by_name($name);

  if(!$old_rc){
    $new_contigs{$name} = 1;
    next RC;
  }

  my $new_rc = $new_rca->fetch_by_name($name);
  my $analysis_adaptor = $new_db->get_AnalysisAdaptor();
  my @similarity_features = @{$old_rc->get_all_SimilarityFeatures()};
  my $pfa = $new_db->get_ProteinAlignFeatureAdaptor;
  my $dfa = $new_db->get_DnaAlignFeatureAdaptor;
  my %analysis_hash;
  my @dnafs;
  my @pepfs;
  foreach my $align(@similarity_features){
    #print STDERR "have align feature ".$align." of type ".$align->analysis->logic_name." matched to ".$align->hseqname." on contig ".$align->contig->name."\n";
    if(!$analysis_hash{$align->analysis->logic_name}){
      my $analysis = $analysis_adaptor->fetch_by_logic_name($align->analysis->logic_name);
      if(!$analysis){
	warn "haven't got analysis of type ".$align->analysis->logic_name." is ".$new_db->name." can't store this\n";
      }else{
	$analysis_hash{$analysis->logic_name} = $analysis;
      }
    }
    $align->contig($new_rc);
    $align->dbID('');
    $align->analysis($analysis_hash{$align->analysis->logic_name});
    if($align->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
      $align->adaptor($pfa);
      push(@pepfs, $align);
    }elsif($align->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
      $align->adaptor($dfa);
      push(@dnafs, $align);
    }
  }
  $pfa->store(@pepfs);
  $dfa->store(@dnafs);

  my @new_pts;
  my @prediction_transcripts = @{$old_rc->get_all_PredictionTranscripts};
  my $pta = $new_db->get_PredictionTranscriptAdaptor;
  foreach my $pt(@prediction_transcripts){
    if(!$analysis_hash{$pt->analysis->logic_name}){
      my $analysis = $analysis_adaptor->fetch_by_logic_name($pt->analysis->logic_name);
      if(!$analysis){
	warn "haven't got analysis of type ".$pt->analysis->logic_name." is ".$new_db->name." can't store this\n";
      }else{
	$analysis_hash{$analysis->logic_name} = $analysis;
      }
    }
    my @exons = @{$pt->get_all_Exons};
    foreach my $e(@exons){
      $e->contig($new_rc);
      $e->adaptor('');
      $e->dbID('');
    }
    $pt->dbID('');
    $pt->adaptor($pta);
    $pt->analysis($analysis_hash{$pt->analysis->logic_name});
    push(@new_pts, $pt);
  }
  $pta->store(@new_pts);
  my @new_sfs;
  my @simple_features = @{$old_rc->get_all_SimpleFeatures};
  my $sfa = $new_db->get_SimpleFeatureAdaptor;
  foreach my $sf(@simple_features){
    if(!$analysis_hash{$sf->analysis->logic_name}){
      my $analysis = $analysis_adaptor->fetch_by_logic_name($sf->analysis->logic_name);
      if(!$analysis){
	warn "haven't got analysis of type ".$sf->analysis->logic_name." is ".$new_db->dbname." can't store this\n";
      }else{
	$analysis_hash{$analysis->logic_name} = $analysis;
      }
    }
    $sf->contig($new_rc);
    $sf->dbID('');
    $sf->adaptor($sfa);
    $sf->analysis($analysis_hash{$sf->analysis->logic_name});
    push(@new_sfs);
  }
  $sfa->store(@new_sfs);
  my $count = 0;
  my @new_rfs;
  my @repeat_features = @{$old_rc->get_all_RepeatFeatures};
 
  my $rfa = $new_db->get_RepeatFeatureAdaptor;
  foreach my $rf(@repeat_features){
    if(!$analysis_hash{$rf->analysis->logic_name}){
      my $analysis = $analysis_adaptor->fetch_by_logic_name($rf->analysis->logic_name);
      if(!$analysis){
	warn "haven't got analysis of type ".$rf->analysis->logic_name." is ".$new_db->name." can't store this\n";
      }else{
	$analysis_hash{$analysis->logic_name} = $analysis;
      }
    }
    $rf->contig($new_rc);
    $rf->dbID('');
    $rf->adaptor($rfa);
    $rf->repeat_consensus->dbID('');
    $rf->analysis($analysis_hash{$rf->analysis->logic_name});
    push(@new_rfs, $rf);
    $count++;
  }
  $rfa->store(@new_rfs);
 
}


open(FH, ">".$WB_NEW_CONTIGS);
my $time = time;
foreach my $name(keys(%new_contigs)){
  print FH $name."\t".$WB_SUBMIT_CONTIG_ID."\t".$time."\n";
}

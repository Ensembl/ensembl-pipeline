use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 7;}

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Runnable::ExonerateArray;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::MiscFeature;
use Bio::SeqIO;

ok(1);

ok(my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host =>'ecs2',
                                             -user  =>'ensro',
                                             -port => 3361,
                                             -dbname=>'yuan_human_core_23'));


ok(my $seqio = new Bio::SeqIO(-file   => 'data/array_test.fa', 
                                 -format => 'fasta'));

my @seqs;

while (my $seq = $seqio->next_seq) {
  push(@seqs,$seq);
}

ok(my $exonerate = Bio::EnsEMBL::Pipeline::Runnable::ExonerateArray->new(
   -db          => $db, 	
   -query_seqs  => \@seqs,
   -query_type  => 'DNA',
   -database    => 'data/array_db.fa',
   -target_type => 'DNA',
   -exonerate   => '/usr/local/ensembl/bin/exonerate-0.8.2',
   -verbose     => 1,
   #-options     => '--exhaustive FALSE --model est2genome --softmasktarget --score 500 --fsmmemory 800  --saturatethreshold 100 --dnahspthreshold 60 --dnawordlen 14 --showalignment FALSE --showvulgar FALSE --percent 90'
   ));

ok($exonerate->run);

ok(my @misc_features = $exonerate->output());

ok(display(@misc_features));

sub display {
  my @misc_features = @_;

  foreach my $misc_feat (@misc_features) {
    my @all_attribs = @{$misc_feat->get_all_Attributes()};
    my @all_sets = @{$misc_feat->get_all_MiscSets()};

    foreach my $attrib (@all_attribs) {
      print $attrib->code,"\t",$attrib->value,"\n";
    }
    foreach my $set (@all_sets) {
      print $set->code,"\t",$set->name,"\n";
    }
    print $misc_feat->slice->name,"\t",$misc_feat->start,"\t",$misc_feat->end,"\t",$misc_feat->strand,"\n";
  }
  return 1
}


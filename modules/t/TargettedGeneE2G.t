use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 3;} 
use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneE2G;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);

ok(my $ens_test = EnsTestDB->new());
ok($ens_test->do_sql_file("t/TargettedGene.dump"));
ok($ens_test->module("Bio::EnsEMBL::DBSQL::DBAdaptor"));
ok(my $db = $ens_test->get_DBSQL_Obj);

my $id='ctg22fin4:17398358,17400418:O60755:AF073799';

my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneE2G';
ok(my $tge = "$runnable"->new(-dbobj    => $db,
                              -input_id => $id));	

ok9$tge->fetch_input());
ok($tge->run());

ok(my @genes = $tge->output());
print STDERR "found " . scalar(@genes) . " genes\n";

# check those genes
my $count=0;
foreach my $gene(@genes){
   $count++;
   my $ecount = 0;
   print STDERR "gene $count\n";
     foreach my $exon ( @{$gene->get_all_Exons} ) {
       $ecount++;
       print STDERR 	"exon $ecount\t" .
			$exon->contig_id . "\t" . 
			$exon->start     . "\t" . 
			$exon->end       . "\t" . 
			$exon->strand    . "\n";
       foreach my $sf (@{$exon->get_all_supporting_features}) {
	  print STDERR "\tsupporting_features: " . 
                       $sf->seqname     . "\t" .
		       $sf->start       . "\t" .
		       $sf->end         . "\t" .
		       $sf->strand      . "\t" .
		       $sf->score    . "\t" .
		       $sf->hseqname    . "\t" .
		       $sf->hstart      . "\t" .
		       $sf->hend        . "\n";
	}
     }
}

ok($tge->write_output());














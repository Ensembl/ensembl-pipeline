use lib 't';
use Test;
use strict;

BEGIN { $| = 1; plan test => 3;} 

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1);
ok(my $ens_test = EnsTestDB->new());

ok($ens_test->do_sql_file("t/dumps/analysis.dump"));
ok($ens_test->do_sql_file("t/dumps/dna.dump"));
ok($ens_test->do_sql_file("t/dumps/repeat.dump"));
ok($ens_test->do_sql_file("t/dumps/dna_align_feature.dump"));
ok($ens_test->do_sql_file("t/dumps/assembly.dump"));

ok($ens_test->module("Bio::EnsEMBL::DBSQL::DBAdaptor"));

ok(my $db = $ens_test->get_DBSQL_Obj);

my $id = 'chr1.75000-100000';

my $fetcher  = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::FPC_BlastMiniEst2Genome';
my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('est2genome');

my $runner = "$runnable"->new(-db         => $db,
			    -input_id   => $id,
                            -analysis   => $analysis);	


ok($runner);
ok($runner->fetch_input());
ok($runner->run());

ok(my @genes = $runner->output());

print STDERR "found " . scalar(@genes) . " genes\n";

# check those genes
my $count=0;
foreach my $gene(@genes){
   $count++;
   my $ecount = 0;
   print STDERR "gene $count\n";
     foreach my $exon ( @{$gene->get_all_Exons()} ) {
       $ecount++;
       print STDERR 	"exon $ecount\t" .
			$exon->contig_id . "\t" . 
			$exon->start     . "\t" . 
			$exon->end       . "\t" . 
			$exon->strand    . "\n";
       foreach my $sf(@{$exon->get_all_supporting_features}) {
	  print STDERR "\tsupporting_features: " . 
                       $sf->seqname     . "\t" .
		       $sf->start       . "\t" .
		       $sf->end         . "\t" .
		       $sf->strand      . "\t" .
		       $sf->hseqname    . "\t" .
		       $sf->hstart      . "\t" .
		       $sf->hend        . "\n";
	}
     }
}

ok($runner->write_output());














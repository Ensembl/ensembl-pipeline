use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 11;}

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise;
use Bio::PrimarySeq;
use Bio::Seq;
use Bio::SeqIO;

ok(1); 

ok(my $ens_test = EnsTestDB->new());
ok($ens_test->do_sql_file("t/dumps/TargettedGene.dump"));
ok(my $db = $ens_test->get_DBSQL_Obj);

my @id = ('16:13850772,13850951:HBA_HUMAN:');


ok(my $fetcher  = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch);

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise');
ok(my $ana_adaptor = $db->get_AnalysisAdaptor);             
ok(my $ana = $ana_adaptor->fetch_by_logic_name('human_swall_protein'));


my $tgw = "$runnable"->new(-db         => $db,	
                           -input_id   => @id,
			   -analysis => $ana,	
			   -seqfetcher => $fetcher);		 

ok($tgw);

$tgw->fetch_input();

ok(1);

$tgw->run();

ok(1);

ok(my @genes = $tgw->output());

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

# no gene_externla table as yet ...
$tgw->write_output();

ok(1);














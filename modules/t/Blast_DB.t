use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 17;}

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::Blast;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Analysis;

ok(1);

ok(my $ens_test = EnsTestDB->new());

ok($ens_test->do_sql_file("t/dumps/analysis.dump"));

ok($ens_test->do_sql_file("t/dumps/dna.dump"));

ok(my $db = $ens_test->get_DBSQL_Obj);

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::Blast');

ok(my $ana_adaptor = $db->get_AnalysisAdaptor);

ok($ana_adaptor);

ok(my $ana = $ana_adaptor->fetch_by_logic_name('wublastx'));

ok($ana_adaptor->exists( $ana ));

my $id = 'test:1:AC099340:1:1000000:1';
my $contig_id = 'AC099340';
my $pwd = `pwd`; chomp $pwd;
my $dbfile = "$pwd/t/data/testpep.fa";

my $blastdb = new Bio::EnsEMBL::Pipeline::Runnable::BlastDB(-dbfile => $dbfile,
							    -type    => 'PROTEIN');

$blastdb->run;

$ana->db_file($dbfile);

$::fasta_header_re{'testpep.fa'} = '^(\w+)\s+';

ok(my $runobj = "$runnable"->new(-db         => $db,
				 -input_id   => $id,
				 -analysis   => $ana));


ok($runobj->fetch_input);
my @runnables = $runobj->runnable;
$runnables[0]->add_regex("$pwd/t/data/testpep.fa",'^(\w+)\s+');

ok($runobj->run);

ok(my @out = $runobj->output);

ok(display(@out));

ok($runobj->write_output());

ok(my @features = @{$db->get_RawContigAdaptor()->fetch_by_name($contig_id)->get_all_SimilarityFeatures()});

ok(display(@features));

sub display {
  my @results = @_;

  foreach my $obj (@results) {
    print $obj->gffstring . "\n";

    if ($obj->sub_SeqFeature) {

      foreach my $exon ($obj->sub_SeqFeature) {
        print STDERR "Sub: ".$exon->gffstring."\n";
      }

    }
  }
  return 1;
}

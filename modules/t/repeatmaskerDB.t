use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 17;}

use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker;
use Bio::EnsEMBL::Analysis;

ok(1);
    
ok(my $ens_test = EnsTestDB->new());
ok($ens_test->do_sql_file("t/runnabledb.dump"));
ok(my $db = $ens_test->get_DBSQL_Obj);

ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::RepeatMasker');
ok(my $ana_adaptor = $db->get_AnalysisAdaptor);
ok($ana_adaptor);
ok(my $ana = $ana_adaptor->fetch_by_logic_name('RepeatMask'));

ok($ana_adaptor->exists( $ana ));

my $id = 'small.contig';

my $runobj = "$runnable"->new(-db         => $db,
															-input_id   => $id,
                              -analysis   => $ana );
ok($runobj);

$runobj->fetch_input;;

ok(1);

$runobj->run;

ok(1);

ok(my @out = $runobj->output);

display(@out);

ok(1);
$runobj->write_output();

ok(1);


ok(my @repeats = @{$db->get_RawContigAdaptor()->fetch_by_name($id)->get_all_RepeatFeatures()});

display(@repeats);

ok(1);

sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results)
    {
       print STDERR ($obj->seqname . "\t" . $obj->start . "\t" . $obj->end . "\t" . $obj->repeat_consensus_id . "\t" . $obj->hstart . "\t" . $obj->hend . "\n");
       if ($obj->sub_SeqFeature)
       {
            foreach my $exon ($obj->sub_SeqFeature)
            {
                print STDERR "Sub: ".$exon->gffstring."\n";
            }
       }
    }
}

use lib 't';
use strict;
use Test;

BEGIN { $| = 1; plan test => 16; }

use EnsTestDB;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep;
use Bio::EnsEMBL::Exon;

ok(1);
ok(my $ens_test = EnsTestDB->new());
ok($ens_test->do_sql_file("t/dumps/blastgenscanpepDB.dump"));
ok(my $db = $ens_test->get_DBSQL_Obj);

ok(my $runnable    = 'Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep');
ok(my $ana_adaptor = $db->get_AnalysisAdaptor());
ok(my $ana         = $ana_adaptor->fetch_by_logic_name('blastgenscanPEP'));

my $id = 'Z84721.1.1.43058';

my $pwd = `pwd`;  chomp($pwd);
my $database = $pwd . "/t/data/mini_protein.fa";

$ana->db_file($database);

my $runobj = "$runnable"->new(  -db		=> $db,
				-input_id	=> $id,
                                -analysis	=> $ana );

$runobj->threshold(1);

ok(1);

$runobj->fetch_input;;

ok(1);

$runobj->run;

ok(1);

ok(my @out = $runobj->output);

display(@out);

$runobj->write_output();

ok(1);

ok(my $contig = $db->get_RawContigAdaptor()->fetch_by_name($id));

ok(my @features = @{$contig->get_all_SimilarityFeatures()});

display(@features);

ok(my @genscan_peptides =  @{$contig->get_all_PredictionTranscripts()});

#check re-translation of features can be matched within peptides
my $all_features_found = 1;
foreach my $feature (@features) {
    next unless ($feature->isa("Bio::EnsEMBL::DnaPepAlignFeature"));

    if ($all_features_found == 1) {
	 $all_features_found = 0;
	 my $exon = Bio::EnsEMBL::Exon->new();
	 $exon->temporary_id           ($feature->seqname);
	 $exon->start        ($feature->start);
	 $exon->end          ($feature->end);
	 $exon->strand       ($feature->strand);
	 $exon->contig_id    ($id);
	 $exon->attach_seq   ($contig);
	 my $count = 1;
	 foreach my $pep (@genscan_peptides)
	 {
	     my $full_pep = $pep->translate;
	     my $exon_pep = $exon->translate->seq;

	     #$exon_pep =~ s/^\M//i; #remove leading M's
	     $exon_pep =~ s/\*$//; 
	     $exon_pep =~ s/X$//; 
	     
	     if (!defined($exon_pep) || index ($full_pep, $exon_pep) > -1) {
		 $all_features_found = 1;
		 last;
	     }
	 }
    }
}

ok($all_features_found);

##############################################################################
sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results)
    {
       print ($obj->gffstring."\n");
       if ($obj->sub_SeqFeature)
       {
            foreach my $exon ($obj->sub_SeqFeature)
            {
                print "Sub: ".$exon->gffstring."\n";
            }
       }
    }
}







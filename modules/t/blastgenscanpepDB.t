## Bioperl Test Harness Script for Modules
##
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------



## We start with some black magic to print on failure.
BEGIN { 
        $| = 1; print "1..8\n"; 
	    use vars qw($loaded); 
      }

END {   print "not ok 1\n" unless $loaded;  }

use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::Utils;
use Bio::EnsEMBL::Exon;
use Data::Dumper;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
# Load some data into the db
#$ens_test->do_sql_file("t/blastgenscanpepDB.dump");
    
# Get an EnsEMBL db object for the test db
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensrv4',
                                      -dbname => 'ensembl_freeze17_michele',
                                      -user => 'ensadmin');
print "ok 2\n";    

my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanPep';
my $parameters = '-THRESHOLD => 1, -ARGS => -hspmax 1000 nogap';
my $ana = Bio::EnsEMBL::Pipeline::Analysis->new (   -db             => 'swir',
                                                    -db_file        => 'swir',
                                                    -db_version     => '__NONE__',
                                                    -program        => 'wublastp',
                                                    -program_file   => 'wublastp',
                                                    -module         => $runnable,
                                                    -module_version => 1,
                                                    -gff_source     => 'blastp',
                                                    -gff_feature    => 'similarity',
                                                    -parameters     => $parameters );

unless ($ana)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }
my $id = 'AP000074.1.1.100000';

my $contig = $db->get_Contig($id);
my @genscan_peptides = $contig->get_genscan_peptides;

my $runobj = "$runnable"->new(  -dbobj      => $db,
			                    -input_id   => $id,
                                -analysis   => $ana );
unless ($runobj)
{ print "not ok 4\n"; }
else
{ print "ok 4\n"; }

$runobj->fetch_input;;
$runobj->run;

my @out = $runobj->output;
unless (@out)
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }
#display(@out);

$runobj->write_output();
my $contig = $db->get_Contig($id);
my @features = $contig->get_all_SimilarityFeatures();
#display(@features);

unless (@features)
{ print "not ok 6\n"; }
else
{ print "ok 6\n"; }

my @genscan_peptides =  $contig->get_genscan_peptides;

unless (@genscan_peptides)
{ print "not ok 7: (Data error or bug in RawContig)\n"; }
else
{ print "ok 7\n"; }

#check re-translation of features can be matched within peptides
my $all_features_found = 1;
foreach my $feature (@features)
{
    next unless ($feature->isa("Bio::EnsEMBL::FeaturePair"));
    if ($all_features_found == 1)
    {
        $all_features_found = 0;
        my $exon = Bio::EnsEMBL::Exon->new();
        $exon->id           ($feature->hseqname);
        $exon->start        ($feature->start);
        $exon->end          ($feature->end);
        $exon->strand       ($feature->strand);
        $exon->phase        ($feature->phase);
        $exon->contig_id    ($id);
        #$exon->end_phase($feat->end_phase);
        $exon->attach_seq   ($contig->primary_seq);
        my $count = 1;
        foreach my $pep (@genscan_peptides)
        {
            my $full_pep = $pep->translate->seq;
            my $exon_pep = $exon->translate->seq;
            $exon_pep =~ s/^\M//i; #remove leading M's
	    $exon_pep =~ s/\*$//; 
	    $exon_pep =~ s/X$//; 
            print STDERR "PEP: $full_pep\n";
            print STDERR "SEQ: $exon_pep\n";
            
            if (!defined($exon_pep) || index ($full_pep, $exon_pep) > -1)
            {
                print STDERR "$count ---------------MATCHED-----------------------\n";
                print STDERR "$count: Feature: ".$feature->hseqname
                ." matched in ".$pep->id."\n";
                print STDERR "$count: PHASE: ".$feature->phase." End ".$feature->end_phase
                                      ." Start ".$feature->start." - ".$feature->end
                                      ."(".$feature->length
                                      .") Pep ".$feature->hstart." - ".$feature->hend
                                      ." Strand ".$feature->strand."\n";
                print STDERR "$count --------------------------------------------\n";
                
                $all_features_found = 1;
                last;
            }
            #<DELETE>
            else
            {
                #print STDERR "$count: NOT FOUND\n";
                #print STDERR "$count: Feature ".$feature->hseqname."\n";
            }
            #</DELETE>
            $count ++;
        }
        #<DELETE>
        unless ($all_features_found)
        {
            print STDERR "$count *********************FAILED******************\n";
            print STDERR "Feature ".$feature->hseqname."\n";
            print STDERR "$count: PHASE: ".$feature->phase." End ".$feature->end_phase
                                      ." Start ".$feature->start." - ".$feature->end
                                      ."(".$feature->length
                                      .") Pep ".$feature->hstart." - ".$feature->hend
                                      ." Strand ".$feature->strand."\n";
            print STDERR "$count **********************************************\n";
            
            
            for (my $phase = 0; $phase < 3; $phase++)
            {
                $exon->phase($phase);
                my $ex_pep = $exon->translate->seq();
                $ex_pep =~ s/^M//i;
                print STDERR "EXON: $ex_pep\n";
                foreach my $pep (@genscan_peptides)
                {
                    my $peptide = $pep->translate->seq();
                    #print STDERR "GENSCAN: \n";
                    print STDERR "SUCCESS\n" if (index($peptide, $ex_pep) > -1);
                }
                
            }
            $all_features_found = 1;
        }
        #</DELETE>
    }
}

unless ($all_features_found)
{ print "not ok 8\n"; }
else
{ print "ok 8\n"; }

##############################################################################
sub display {
    my @results = @_;
    #Display output
    print STDERR "RESULTS FROM RUN\n";
    foreach my $obj (@results)
    {
       print STDERR ($obj->gffstring."\n");
       if ($obj->sub_SeqFeature)
       {
            foreach my $exon ($obj->sub_SeqFeature)
            {
                print STDERR "Sub: ".$exon->gffstring."\n";
            }
       }
    }
}

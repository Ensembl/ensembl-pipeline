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
        $| = 1; print "1..6\n"; 
	    use vars qw($loaded); 
      }

END {   print "not ok 1\n" unless $loaded;  }

use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanDNA

$loaded = 1;

print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();

$ens_test->do_sql_file("t/dumps/blastgenscanpepDB.dump");
    
my $db = $ens_test->get_DBSQL_Obj;

print "ok 2\n";    

###########################
# Build the analysis object
###########################

my $runnable    = 'Bio::EnsEMBL::Pipeline::RunnableDB::BlastGenscanDNA';

my $ana_adaptor = $db->get_AnalysisAdaptor();
my $ana = $ana_adaptor->fetch_by_logic_name('blastgenscanDNA');

my $database = `pwd`;
chomp($database);

$database .= "/t/data/mini_mrna.fa";
$ana->db_file($database);
$ana->parameters('-B=10');

unless ($ana) {
  print "not ok 3\n"; 
} else {
  print "ok 3\n"; 
}


$ana_adaptor->exists( $ana );

#####################
# Create the runnable
#####################

my $id = 'Z84721.1.1.43058';

my $runobj = "$runnable"->new(-db         => $db,
			      -input_id   => $id,	
			      -analysis   => $ana );
unless ($runobj) {
 print "not ok 4\n"; 
} else { 
  print "ok 4\n"; 
}

##################
# Run the runnable
##################

$runobj->fetch_input;;
$runobj->run;

my @out = $runobj->output;

unless (@out) {
  print "not ok 5\n"; 
} else {
  print "ok 5\n"; 
}

display(@out);

##################
# Write the output
##################

$runobj->write_output();

#########################################
# Retrieve the features from the database
#########################################

my $contig   = $db->get_RawContigAdaptor()->fetch_by_name($id);
my @features = @{$contig->get_all_SimilarityFeatures()};

display(@features);

unless (@features) {
  print "not ok 6\n"; 
} else {
  print "ok 6\n"; 
}


##############################################################################
sub display {
    my @results = @_;
    #Display output
    foreach my $obj (@results)
    {
       print ($obj->gffstring."\n");
       print ("PHASE: ".$obj->phase."\n") if (defined($obj->phase));
       
       if ($obj->sub_SeqFeature)
       {
            foreach my $exon ($obj->sub_SeqFeature)
            {
                print "Sub: ".$exon->gffstring."\n";
                print ("PHASE: ".$exon->phase."\n");
            }
       }
    }
}

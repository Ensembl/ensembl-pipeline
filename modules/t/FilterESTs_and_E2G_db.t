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
	    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
	    use vars qw($loaded); 
      }

END {   print "not ok 1\n" unless $loaded;  }

use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G;
use Bio::EnsEMBL::Analysis;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new(
   {schema_sql => [ '../sql/table.sql', '../../ensembl/sql/table.sql', '../sql/est.sql' ]}
   ); # Must remember to create the correst est tables - THIS OVERRIDES THE CONFIG FILE

# Load some data into the db
$ens_test->do_sql_file("t/FilterESTs_and_E2G_db.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    

my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G';
my $ana_adaptor = $db->get_AnalysisAdaptor;
my $ana = $ana_adaptor->fetch_by_logic_name('exonerate_e2g');

unless ($ana)
{ print "not ok 3\n"; }
else
{ print "ok 3\n"; }

my $id ='cb25.fpc0829.1800001-1850000';
$ana_adaptor->exists( $ana );
my $runobj = "$runnable"->new(-db      => $db,
			      -input_id   => $id,
			      -analysis   => $ana );
unless ($runobj)
{ print "not ok 4\n"; }
else
{ print "ok 4\n"; }

$runobj->fetch_input;
$runobj->run;

my @out = $runobj->output;
unless (@out)
{ print "not ok 5\n"; }
else
{ print "ok 5\n"; }


$runobj->write_output();

print "ok 6\n";

#Hmm, maybe should try retrieving things.

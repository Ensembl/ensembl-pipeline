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

use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G;
use Bio::EnsEMBL::Analysis;
use Test;

# Use the ESTConf here and then later reset vars to internal test dbs.  Unless 
# this is done first, when a new FilterESTs_and_E2G runnable is made it will 
# read/reread global variables from ESTconf.pm.  This would mean that the 
# test then starts using whatever dbs are configured external to this test.
# If this were to happen, innocently running this test will write junk into 
# user dbs that may well be configured for a production run - erk.

use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_REFDBHOST
					EST_REFDBNAME
					EST_REFDBUSER
					EST_DBNAME
					EST_DBHOST
					EST_DBUSER 
					EST_DBPASS
					EST_SOURCE
					EST_INDEX
					EST_MIN_PERCENT_ID
					EST_MIN_COVERAGE
					EST_INPUTID_REGEX
					EST_GENETYPE
					EST_FEATFILT_COVERAGE
					EST_FEATFILT_MINSCORE
				       );

BEGIN { $| = 1; plan test => 9;}

ok(1);   #1
    
my $ens_test = EnsTestDB->new(
   {schema_sql => [ '../sql/table.sql', '../../ensembl/sql/table.sql', '../sql/est.sql' ]}
   ); # Must remember to create the correst est tables - THIS OVERRIDES THE CONFIG FILE

# Load some data into the db
$ens_test->do_sql_file("t/dumps/FilterESTs_and_E2G_db.dump");
    
# Get an EnsEMBL db object for the test db
ok(my $db = $ens_test->get_DBSQL_Obj); #2

# Set global variable to our own potted database.
$EST_REFDBHOST	       = $ens_test->{'host'};
$EST_REFDBNAME	       = $db->{'_dbname'};
$EST_REFDBUSER	       = $ens_test->{'readonly_user'};
$EST_DBNAME	       = $db->{'_dbname'};
$EST_DBHOST	       = $ens_test->{'host'};
$EST_DBUSER	       = $ens_test->{'user'};
$EST_DBPASS	       = $ens_test->{'pass'};
$EST_SOURCE	       = 'worm_ests';
$EST_INDEX	       = 't/data/ests';
$EST_MIN_PERCENT_ID    = 70; 
$EST_MIN_COVERAGE      = 90;
$EST_INPUTID_REGEX     = '(^\S+\.\S+)\.(\d+)-(\d+)';
$EST_GENETYPE	       = 'exonerate_e2g';
$EST_FEATFILT_COVERAGE = 500;
$EST_FEATFILT_MINSCORE = 10;



ok(my $runnable = 'Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G');  #3
ok(my $ana_adaptor = $db->get_AnalysisAdaptor);                               #4
ok(my $ana = $ana_adaptor->fetch_by_logic_name('exonerate_e2g'));             #5

my $id ='cb25.fpc0829.1800001-1850000';

ok($ana_adaptor->exists( $ana ));       #6

ok(my $runobj = "$runnable"->new(-db       => $db,
	      		         -input_id => $id,
			         -analysis => $ana ));       #7

ok($runobj->fetch_input);    #8
ok($runobj->run);            #9

$runobj->write_output;

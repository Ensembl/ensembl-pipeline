#! /usr/local/ensembl/bin/perl -w


use strict;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Pipeline::Config::FirstEF qw (
						FEF_WRITEDBNAME
						FEF_WRITEDBHOST
						FEF_WRITEDBUSER
						FEF_WRITEDBPASS
						FEF_REFDBNAME 
						FEF_REFDBHOST
						FEF_REFDBUSER
					       );

my $input_id;
my $write = 0;

&GetOptions( 
	     'input_id:s'  => \$input_id,
             'write'       => \$write
	     );




my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $FEF_WRITEDBNAME,
					     -host   => $FEF_WRITEDBHOST,
					     -user   => $FEF_WRITEDBUSER,
					     -pass   => $FEF_WRITEDBPASS);

my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $FEF_REFDBNAME, 
						-host   => $FEF_REFDBHOST,
						-user   => $FEF_REFDBUSER);

$db->dnadb($dnadb);

my $analysis_adaptor = $db->get_AnalysisAdaptor;
my $analysis = $analysis_adaptor->fetch_by_logic_name('firstef');

my $firstef = Bio::EnsEMBL::Pipeline::RunnableDB::FirstEF->new(-input_id  => $input_id,
							       -db        => $db,
							       -analysis  => $analysis );

$firstef->fetch_input;
$firstef->run;

my @output = $firstef->output;

print "Found " . scalar @output . " features\n";

if ($write) {
  $firstef->write_output;
}

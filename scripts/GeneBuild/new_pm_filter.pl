#!/usr/local/bin/perl 

use Bio::EnsEMBL::Pipeline::RunnableDB::Pmatch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

$| = 1;

my ($dbhost, $dbname, $dbpass, $dbuser, $id, $logic, $write);

&GetOptions( 
	    'dbname:s' => \$dbname,
	    'dbhost:s' => \$dbhost,
	    'dbuser:s' => \$dbuser,
	    'dbpass:s' => \$dbpass,
	    'input_id:s' => \$id,
	    'analysis:s' => \$logic,
	    'write' => $write,
	   );


if(!$dbname || !$dbhost || !$dbuser || !$dbpass || !$id){
  print STDERR "Usage = -dbname -dbuser -dbhost -dbpass -input_id -analysis -write\n";
  print STDERR "input id must be in the format specified in GB_INPUTID_REGEX in Bio::EnsEMBL::Pipeline::Config::GeneBuild::General\n";
  exit(0)
}


my $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(-host => $dbhost,
						      -user => $dbuser,
						      -pass => $dbpass,
						      -dbname => $dbname,
						     );


if(!$logic){
  print STDERR "You have provided no logic name assuming it is Pmatch\n";
  $logic = 'Pmatch';
}

my $analysis_adaptor = $db->get_AnalysisAdaptor;
my $ana = $analysis_adaptor->fetch_by_logic_name($logic);

print STDERR "input_id = ".$id."\n dbobj = ".$db."\n anaysis = ".$ana."\n";

my $runnable = Bio::EnsEMBL::Pipeline::RunnableDB::Pmatch->new (
								-db    => $db,
								-input_id => $id,
								-analysis => $ana,								
							       );

$runnable->fetch_input();
$runnable->run();

if($write){
print STDERR "Writing output to database\n";
$runnable->write_output;
}
print STDERR "have ".$runnable->output." results\n";

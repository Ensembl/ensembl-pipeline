#!/usr/local/bin/perl 

use Bio::EnsEMBL::Pipeline::RunnableDB::BestPmatch;
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
	    'write!' => \$write,
	   );
							

if(!$dbname || !$dbhost || !$dbuser || !$dbpass){
  print STDERR "Usage = -dbname -dbuser -dbhost -dbpass -input_id -analysis -write\n";
  exit(0)
}

if(!$logic){
  print STDERR "You have provided no logic name assuming it is BestPmatch\n";
  $logic = 'BestPmatch';
}


my $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(-host => $dbhost,
						      -user => $dbuser,
						      -pass => $dbpass,
						      -dbname => $dbname,
						     );

my $analysis_adaptor = $db->get_AnalysisAdaptor;

my $ana = $analysis_adaptor->fetch_by_logic_name($logic);
			     


if(!$id){
  $id = 'genome';
}
print "input_id = ".$id."\n dbobj = ".$db."\n anaysis = ".$ana."\n";

my $runnable = Bio::EnsEMBL::Pipeline::RunnableDB::BestPmatch->new (
								-db    => $db,
								-input_id => $id,
								-analysis => $ana,								
							       );
#print "created runnable\n";
$runnable->fetch_input();
$runnable->run();
if($write){
print STDERR "Writing output to database\n";
$runnable->write_output;
}
print STDERR "have ".$runnable->output." results\n";

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Pipeline::Config;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::PipelineManager;
# set a sgnal handler 
$SIG{INT} = \&_shut_down;

# Parse command line
my @files  = ();
# should there be defauts for these?
my $host = '';
my $user = '';
my $pass = '';
my $port = '';
my $dbname = '';
GetOptions ('file=s'   => \@files,
	    'host=s'   => \$host,
	    'port=s'   => \$port,
	    'user=s'   => \$user,
	    'pass=s'   => \$pass,
	    'dbname=s' => \$dbname);

# create config with file or db
my $config;

if (@files > 0) {

  $config = new Bio::EnsEMBL::Pipeline::Config(-FILES => @files);

} else {

  # create DBAdaptor, then Config
  my $dbobj = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor('-host'   => $host,
							   '-port'   => $port,
							   '-dbname' => $dbname,
							   '-user'   => $user,
							   '-pass'   => $pass);

  $config = new Bio::EnsEMBL::Pipeline::Config(-DB => $dbobj);

}

# instantiate PipelineManager
my $pipeline_manager = new Bio::EnsEMBL::Pipeline::PipelineManager($config);

# call run
$pipeline_manager->run();

# Called when signal is recieved
sub _shut_down {
  print STDERR "\n\nReceived INTERRUPT, will shutdown gracefully...\n\n";
  $pipeline_manager->{'stop'} = 1;

}

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Pipeline::Config;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::PipelineManager;
# set a signal handler 
$SIG{INT}  = \&_shut_down;
$SIG{TERM} = \&_shut_down;

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


sub usage() {
  print STDERR <<EOF
startPipeline (file_options | database_options)

  file_options: -file <file1> [-file <file2> [...]]
  file1, file2, etc. are the names of configuration files.  One of the
  configuration files must specify a [PIPELINE_DATABASE] that the configuration
  and other pipeline control information will be written to.

  database_options: -dbname <dbname> [-host <host>] [-port <port>] [-user user]
                    [-pass <pass>]
  Information necessary to connect to a mysql database where an existing 
  configuration is already stored.

EOF
;

}

# create config with file or db
my $config;

if (@files > 0) {

  $config = new Bio::EnsEMBL::Pipeline::Config(-FILES => \@files);

} else {

  die(usage()) if(!$dbname);

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





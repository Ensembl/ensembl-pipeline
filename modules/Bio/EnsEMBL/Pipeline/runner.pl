#!/usr/local/ensembl/bin/perl

# default pipeline runner script
# this script is passed as part of the batch submission request
# and is therefore the first script to be executed on the
# remote side

# note that you may need to alter the #! line above to suit your
# local perl installation


use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Config::General;
use Sys::Hostname;
use Getopt::Long;


#parameters for Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

my $host;
my $dbname;
my $dbuser;

my $port            = '3306';
my $pass            = undef;
my $job_id;


GetOptions(
    'host=s'       => \$host,
    'port=n'       => \$port,
    'dbname=s'     => \$dbname,
    'dbuser=s'     => \$dbuser,
    'pass=s'       => \$pass,
    'check!'       => \$check,
    'output_dir=s' => \$output_dir
)
or die ("Couldn't get options");

if( defined $check ) {
  my $host = hostname();
  if ( ! -e $output_dir ) {
    die "output dir $output_dir doesn't exist according to host $host";
  }
  my $deadhostfile = "$output_dir/deadhosts";
  open( FILE, $deadhostfile ) or exit 0;
  while( <FILE> ) {
    chomp;
    if( $host eq $_ ) {
      die "Cant use this host";
    }
  }
  exit 0;

}

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $host,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $pass,
    -port   => $port
)
or die ("Failed to create Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor to db $dbname \n");

my $job_adaptor = $db->get_JobAdaptor();

# Print a list of jobs in this 'batch' - this is perhaps useful later
# to parse these output files for CPU/Mem usage etc.
print STDOUT "LSF Batch summary\n";
print STDOUT "Time ", scalar localtime time, " (", time, ")\n";
print STDOUT "Job ID\tinput ID\tanalysis ID\n";
foreach my $id (@ARGV) {
    my $job = $job_adaptor->fetch_by_dbID($id);
    print STDOUT join ("\t", $id, $job->input_id, $job->analysis->logic_name), "\n";
}

while( $job_id = shift ) {

  print STDERR "Fetching job " . $job_id . "\n";

  my $job         = $job_adaptor->fetch_by_dbID($job_id);

  if( !defined $job) {
    print STDERR ( "Couldnt recreate job $job_id\n" );
    next;
  }
  print STDERR "Running job $job_id\n";
  print STDERR "Module is " . $job->analysis->module . "\n";
  print STDERR "Input id is " . $job->input_id . "\n";
  print STDERR "Files are " . $job->stdout_file . " " . $job->stderr_file . "\n";

  eval {
    $job->run_module;
  };
  $pants = $@;

  if ($pants) {
    print STDERR "Job $job_id failed: [$pants]";
  }

  print STDERR "Finished job $job_id\n";
}


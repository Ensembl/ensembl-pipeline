#!/usr/local/ensembl/bin/perl

BEGIN {
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";
}


use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Sys::Hostname;

use Getopt::Long;


#parameters for Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor

my $host;
my $dbname;
my $dbuser;

my $port            = '3306';
my $pass            = undef;
my $module          = undef;

my $test_env; #not used at present but could be used to test LSF machine

my $object_file;
my $job_id;

print STDERR join( " ", @ARGV ),"\n";

GetOptions(
    'host=s'    => \$host,
    'port=n'    => \$port,
    'dbname=s'  => \$dbname,
    'dbuser=s'  => \$dbuser,
    'pass=s'    => \$pass,
    'check!'    => \$check,
    'objfile=s' => \$object_file,
    'module=s'  => \$module
)
or die ("Couldn't get options");

if( defined $check ) {
  my $host = hostname();
  if ( ! -e $::pipeConf{'nfstmp.dir'} ) {
    die "no nfs connection";
  }
  my $deadhostfile = $::pipeConf{'nfstmp.dir'}."/deadhosts";
  open( FILE, $deadhostfile ) or exit 0;
  while( <FILE> ) {
    chomp;
    if( $host eq $_ ) {
      die "Cant use this host";
    }
  }

  # tests for DB existence - these probably shouldn't be hard-wired in ...
  if (defined (my $dir = $ENV{"BLASTDB"})) {
    -e "$dir/unigene.seq" or warn "Not found unigene";
    -e "$dir/sptr" or warn "Not found sptr";
  }
  exit 0;

}

print STDERR "In runner\n";
my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
    -host   => $host,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $pass,
    -port   => $port,
    -perlonlyfeatures  => 1,
    -perlonlysequences => 1
)
or die ("Failed to create Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor to db $dbname \n");

print STDERR "Connected to database\n";


print STDERR "Getting job adaptor\n";

my $job_adaptor = $db->get_JobAdaptor();

# Print a list of jobs in this 'batch' - this is perhaps useful later
# to parse these output files for CPU/Mem usage etc.
print STDOUT "LSF Batch summary\n";
print STDOUT "Time ", scalar localtime time, " (", time, ")\n";
# print STDOUT "LSF ID: ", $job_adaptor->fetch_by_dbID($ARGV[0])->LSF_id, "\n";
print STDOUT "Job ID\tinput ID\tanalysis ID\n";
foreach my $id (@ARGV) {
    my $job = $job_adaptor->fetch_by_dbID($id);
    print STDOUT join ("\t", $id, $job->input_id, $job->analysis->logic_name), "\n";
    # print STDOUT $id, "\n";
}

# print STDERR "Got job adapter\n";

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
    $job->runLocally;
  };
  $pants = $@;

  if ($pants) {
    print STDERR "Job $job_id failed: [$pants]";
  }

  print STDERR "Finished job $job_id\n";
}
print STDERR "Leaving runner\n";

$db->{'_db_handle'}->DESTROY();


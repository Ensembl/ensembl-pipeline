#!/usr/local/bin/perl 

BEGIN {
#    shift (@INC,"/nfs/disk100/humpub/mq2/ensembl-pipeline/modules");
#    unshift (@INC,"/nfs/disk100/humpub/mq2/ensembl/modules");
#    unshift (@INC,"/nfs/disk100/humpub/birney/ensembl-pipeline/modules");
#    unshift (@INC,"/nfs/disk100/humpub/birney/ensembl/modules");
#    unshift (@INC,"/nfs/disk100/humpub/birney/bioperl-live");
#    unshift (@INC,"/nfs/disk100/humpub/michele/bioperl-06");
#    unshift (@INC,"/nfs/disk100/humpub/mq2");
#    unshift (@INC,"/nfs/disk65/mq2/ensembl/modules");
#    unshift (@INC,"/nfs/disk65/mq2/ensembl-pipeline/modules");
    require "Bio/EnsEMBL/Pipeline/pipeConf.pl";

}


use Bio::EnsEMBL::Pipeline::DBSQL::Obj;
use Sys::Hostname;

use Getopt::Long;


#parameters for Bio::EnsEMBL::Pipeline::DBSQL::Obj
my $host            = 'ensrv3.sanger.ac.uk';
my $port            = '3306';
my $dbname          = 'arne_anaTest';
my $dbuser          = 'ensadmin';
my $pass            = undef;
my $module          = undef;
my $test_env; #not used at present but could be used to test LSF machine

my $object_file;
my $job_id;

GetOptions(     'host=s'     => \$host,
                'port=n'     => \$port,
                'dbname=s'   => \$dbname,
                'dbuser=s'   => \$dbuser,
                'pass=s'     => \$pass,
                'job=n'      => \$job_id,
		'check!'    => \$check,
                'objfile=s'  => \$object_file,
                'module=s'   => \$module ) or die ("Couldn't get options");

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
  exit 0;
}


if( defined $job_id) {
my $db = Bio::EnsEMBL::Pipeline::DBSQL::Obj->new 
    (  -host   => $host,
       -user   => $dbuser,
       -dbname => $dbname,
       -pass   => $pass,
       -port   => $port,
       -perlonlyfeatures  => 1,
       -perlonlysequences => 1 )
    or die ("Failed to create Bio::EnsEMBL::Pipeline::Obj to db $dbname \n");


    
  my $job_adaptor = $db->get_JobAdaptor();
  my $job         = $job_adaptor->fetch_by_dbID($job_id);

  if( !defined $job) {
    die( "Couldnt recreate job $job_id" );
  }

  $job->runInLSF;
} else {
  die( "Called runner without job_id" );
}

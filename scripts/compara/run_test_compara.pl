#!/usr/local/bin/perl


use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer;



use Getopt::Long;

my $host   = 'ecs1b';
my $port   = undef;
my $dbname = 'abel_mouse_human';
my $dbuser = 'ensadmin';
my $pass   = 'ensembl';




&GetOptions( 
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname,
	     'dbuser:s'  => \$dbuser,
	     'pass:s'    => \$pass,
	     );


my $input_id = shift;

if( !defined $input_id ) {
   die "Must call run_test_compara with an input id!";
 }


$db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -host => $host,
						    -user => $dbuser,
						    -pass => $pass,
						    -dbname => $dbname );

my $tag1;
my $contig1;

my $tag2;
my $contig2;

if( $input_id =~ /(\S+):(\S+)::(\S+):(\S+)/ ) {
  $tag1 = $1;
  $contig1 = $2;
  $tag2 = $3;
  $contig2 = $4;
} else {
  die "Input id should be yadda:contig_id::yadda:contig_id";
}


$rundb = Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer->new(
								-dbobj => $db,
								-input_id => $input_id );


# run the runnabledb

$rundb->fetch_input();
$rundb->run();

$rundb->write_output();


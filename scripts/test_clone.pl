#!usr/local/bin/perl

use strict;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Pipeline::DBSQL::Clone;
use Getopt::Long;
use Bio::Seq;

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'timtest';
my $dbuser = 'elia';
my $dbpass = '';
my $module = 'Bio::EnsEMBL::Pipeline::DBSQL::Obj';
my $help;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'h|help'     => \$help
	     );

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass;debug=10";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
print "Testing create_Clone....\n";
$db->create_Clone('dummy','SU','22');
print "Created clone with id 'dummy', clone_group 'SU' and chromosome '22'\n";
print "Testing get_Clone using 'dummy' as an id....\n";
my $clone=$db->get_Clone('dummy');
print "Ok, got clone...\n";

print "Testing clone contents...\n";
print "Clone ".$clone->disk_id.":\n"; 
print "      clone_group ".$clone->clone_group."\n";
print "      chromosome ".$clone->chromosome."\n";      
print "      last_check ".localtime($clone->last_check)."\n";
print "      created ".localtime($clone->created)."\n";
print "      modified ".localtime($clone->modified)."\n";
print "      dna_update_state ".$clone->dna_update_state."\n";
print "      update_state ".$clone->update_state."\n";
print "      update_label ".$clone->update_label."\n";
print "      update_date ".$clone->update_date."\n";
print "      internal_lock ".$clone->internal_lock."\n";
print "      external_lock ".$clone->external_lock."\n";

print "Deleting dummy clone...\n";
$db->delete_Clone('dummy');
print "Test successful!\n";
$db->DESTROY;




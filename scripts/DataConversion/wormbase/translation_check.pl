#!/usr/local/ensembl/bin/perl

use strict;

use WormBaseConf;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


$| = 1;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $WB_DBHOST,
					    -user => $WB_DBUSER,
					    -dbname => $WB_DBNAME,
					    -pass  => $WB_DBPASS,
					   );


my $sql = "select transcript_id from transcript";
my $transcript_adaptor = $db->get_TranscriptAdaptor;
my $sth = $db->prepare($sql);
$sth->execute;

while (my ($dbID) = $sth->fetchrow){
  print "checking ".$dbID."\n";
  my $transcript = $transcript_adaptor->fetch_by_dbID($dbID);
  
  my $pep = $transcript->translate->seq;
  if($pep =~ /\*/){
    print "transcript ".$transcript->stable_id." dbid ".$transcript->dbID." doesn't translate\n";
  }
}


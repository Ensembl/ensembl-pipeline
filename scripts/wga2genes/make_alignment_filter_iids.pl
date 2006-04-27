#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

my (
    $dbname,
    $dbhost,
    $dbuser,
    $dbport,
    $dbpass,
    $seq_dbname,
    $seq_dbuser,
    $seq_dbpass,
    $seq_dbport,
    $seq_dbhost,
    $target_logic,
    $write,
);

$dbuser = 'ensro';
$dbport = 3306;
$seq_dbuser = "ensro";
$seq_dbport = 3306;
$write = 0;

&GetOptions(
            'dbname=s' => \$dbname,
            'dbuser=s' => \$dbuser,
            'dbhost=s' => \$dbhost,
            'dbport=s' => \$dbport,
            'dbpass=s' => \$dbpass,
            'seqdbname=s' => \$seq_dbname,
            'seqdbhost=s' => \$seq_dbhost,
            'seqdbport=s' => \$seq_dbport,
            'seqdbuser=s' => \$seq_dbuser,
            'seqdbpass=s' => \$seq_dbpass,
            'logic=s' => \$target_logic,
            'write' => \$write,
            );

my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new(
	'-dbname' => $dbname,
	'-host' => $dbhost,
	'-user' => $dbuser,
	'-port' => $dbport,
);

my $seq_db;
if (defined $seq_dbname) {
  $seq_db = Bio::EnsEMBL::DBSQL::DBAdaptor->
    new(
	'-dbname' => $seq_dbname,
	'-host' => $seq_dbhost,
	'-port' => $seq_dbport,
	'-user' => $seq_dbuser,
        '-pass' => $seq_dbpass,
        );
} else {
  $seq_db = $db;
}

my $ana = $db->get_AnalysisAdaptor->fetch_by_logic_name($target_logic);
if (not defined $ana) {
  die "Could not find analysis with logic name '$target_logic' in pipe db\n";
}


my @iids;
foreach my $sl (@{$seq_db->get_SliceAdaptor->fetch_all('toplevel')}) {
  my $iid = $sl->seq_region_name;
 
  push @iids, $sl->seq_region_name;
}


if ($write) {
  my $s_inf_cont = $db->get_StateInfoContainer;
   
  foreach my $iid (@iids) {
    eval {
      $s_inf_cont->store_input_id_analysis( $iid, $ana, '' );
    };
    if ($@) {
      print STDERR "Input id $iid already present\n";
    } else {
      print STDERR "Stored input id $iid\n";
    }
  } 
} else {
  foreach my $iid (@iids) {
    printf("INSERT into input_id_analysis(input_id, input_id_type, analysis_id " . 
           " values('%s', '%s', %d);\n", 
           $iid, 
           $ana->input_id_type,
           $ana->dbID);
    
  }
}



#!/usr/local/ensembl/bin/perl 
# this is a script to check pipeline sanity before you run your pipeline

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks;
use Bio::EnsEMBL::Pipeline::Config::General;
use Bio::EnsEMBL::Pipeline::Config::BatchQueue;

my $host     = '';
my $dbname   = '';
my $dbuser   = '';
my $pass     = '';
my $port     = '';
my $info; #prints out more info if true

&GetOptions(
            'dbhost:s'          => \$host,
            'dbname:s'          => \$dbname,
            'dbuser:s'          => \$dbuser,
            'dbpass:s'          => \$pass,	    
            'dbport:s'          => \$port,
            'info!' => \$info,
           );




my $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new
  ( -host   => $host,
    -dbname => $dbname,
    -user   => $dbuser,
    -pass   => $pass,
    -port   => $port
  );

my $sanity = Bio::EnsEMBL::Pipeline::Utils::PipelineSanityChecks->new
  (
   -DB => $db,
  );


#checks general config values are setup as expected
if(!&config_sanity_check($info)){
  print ("Your pipeline configuration isn't set up correctly");
}

#checks the batchqueue file is setup as exepcted
if(!&batchqueue_sanity_check($db, $info)){
  print ("Your BatchQueue.pm isn't setup correctly\n");
}

#checks all the modules in the analysis table compile
if(!&module_compile($db, $info)){
  print ("Your modules don't all compile properly\n");
}

my @rules = $db->get_RuleAdaptor->fetch_all;

if(!&accumulator_sanity($sanity, \@rules)){
  print ("Your accumulators aren't setup correctly\n");
}

&db_sanity_check($sanity);

$sanity->rule_type_sanity(\@rules, $info);



sub accumulator_sanity{
  my $sanity = shift;
  my $accumulators = 1;
  my $rules = shift;
  return $sanity->accumulator_sanity_check($rules, $accumulators);
}

sub module_compile{
  my $db = shift;
  my $info = shift;
  my @analyses = @{$db->get_AnalysisAdaptor->fetch_all};
  my $ok = 1;
  ANALYSIS:foreach my $analysis(@analyses){
    my $module = $analysis->module;
    if(!$module ){
      print STDERR "analysis ".$analysis->logic_name." has no module\n";
      next ANALYSIS;
    }
    if($module eq  'Dummy'){
      next ANALYSIS;
    }
    my $perl_path;
    my $runnable_db_path = 'Bio/EnsEMBL/Pipeline/RunnableDB';
    if($module =~ /::/){
      print STDERR "Module contains path info already\n" if($info);
      $module =~ s/::/\//g;
      $perl_path = $module;
    }else{
      $perl_path = $runnable_db_path."/".$module;
    }
    eval {
      require $perl_path.".pm";
    };
    if($@){
      $ok = 0;
      print "module ".$perl_path." doesn't compile for ".
        $analysis->logic_name." [$@]\n";
    }
  }
  return $ok;
}



sub batchqueue_sanity_check{
  my $db = shift;
  my $info = shift;
  my $ok = 1;
  if($DEFAULT_BATCH_QUEUE){
    print "Your default queue is ".$DEFAULT_BATCH_QUEUE."\n" if($info);
  }else{
    print "You have no DEFAULT_BATCH_QUEUE\n";
    $ok = 0;
  }
  if(-d $DEFAULT_OUTPUT_DIR){
    print "Your default output dir is ".$DEFAULT_OUTPUT_DIR."\n" if($info);
  }else{
    print "Your DEFAULT_OUTPUT_DIR $DEFAULT_OUTPUT_DIR is either not ".
      "defined or doesn't exist\n";
    $ok = 0;
  }
  my $number_of_analyses = scalar(@$QUEUE_CONFIG);
  my @analyses = @{$db->get_AnalysisAdaptor->fetch_all};
  if($number_of_analyses != scalar(@analyses)){
    print "You have written configuration for $number_of_analyses but you ".
      "have ".@analyses." analyses in your database\n";
    $ok = 0;
    if($info){
      my %configured;
      foreach my $q (@$QUEUE_CONFIG){
        my $ln = $q->{'logic_name'};
        $configured{$ln} = 1;
      }
      foreach my $a(@analyses){
        if(!$configured{$a->logic_name}){
          print STDERR $a->logic_name." doesn't have its own configuration\n";
        }
      }
    }
  }
  return $ok;
}


sub config_sanity_check {
  my $ok = 1;
  my $info = shift;
  print STDERR "checking config sanity\n" if($info);
  unless ($QUEUE_MANAGER) {
    print "Need to specify QUEUE_MANAGER in Config/BatchQueue.pm\n";
    $ok = 0;
  }
  unless ($LIB_DIR) {
    print "Need to specify LIB_DIR in Config/General.pm\n";
    $ok = 0;
  }
  unless ($DATA_DIR) {
    print "Need to specify DATA_DIR in Config/General.pm\n";
    $ok = 0;
  }
  unless ($BIN_DIR) {
    print "Need to specify BIN_DIR in Config/General.pm\n";
    $ok = 0;
  }
  
  return $ok;
}



sub db_sanity_check{
  my ($sanity) = @_;
  eval{
    $sanity->db_sanity_check;
  };
  if($@){
    print STDERR "Your database doesn't pass all the sanity checks\n".
      "$@\n";
  }
}

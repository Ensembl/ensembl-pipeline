#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (GB_DBNAME
							     GB_DBHOST
							     GB_DBUSER
							     GB_DBPASS
							    );
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts qw (	    
							   GB_SIZE
							   GB_GENEBUILD_PIPELINE
							   GB_SUBMIT_GENEBUILD_LOGICNAME
							   GB_PMATCH_PIPELINE
							   GB_SUBMIT_PMATCH_LOGICNAME
							   GB_PMATCH_SIZE
							   GB_PMATCH_CHUNKS
							   
							  );

use Bio::EnsEMBL::DBSQL::DBAdaptor;

if($GB_DBUSER eq 'ensadmin' && $GB_DBPASS eq ''){
  print "You cannot have dbuser set to ensadmin with no dbpass set!\nPlease correct the entries in GeneBuild config files\n";
  exit(1);
}


foreach my $arg($GB_DBNAME, $GB_DBHOST, $GB_DBUSER){
  if ($arg eq '' ){
    print "You need to set various parameters in GeneBuild config files\n" .  
      "Here are your current values for required settings: \n" .
      "dbname      => $GB_DBNAME\n" .
      "dbhost      => $GB_DBHOST\n" .
      "dbuser      => $GB_DBUSER\n" .
      "dbpass      => $GB_DBPASS\n";
    
    exit(1);
  }
}

my %chrhash;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $GB_DBHOST,
					      -user   => $GB_DBUSER,
					      -pass   => $GB_DBPASS,
					      -dbname => $GB_DBNAME,
					     );
&get_chrlengths($db);



my $analysis_adaptor = $db->get_AnalysisAdaptor;


if($GB_GENEBUILD_PIPELINE){
  my @input_ids;
  my $analysis = $analysis_adaptor->fetch_by_logic_name($GB_SUBMIT_GENEBUILD_LOGICNAME);
  foreach my $chr(keys %chrhash) {
    #print STDERR "making input ids for ".$chr."\n";
    my $length = $chrhash{$chr};
    my $count = 1;
    my $size = $GB_SIZE;
    while ($count < $length) {
      my $start = $count;
      my $end   = $count + $size -1;
      
      if ($end > $length) {
	$end = $length;
      }
      
      my $input_id = $chr . "." . $start . "-" .  $end;
      push(@input_ids, $input_id);
      $count = $count + $size;
    }
  }
  &insert_input_ids($analysis, $db, @input_ids);
}

if($GB_PMATCH_PIPELINE){
  my @input_ids;
  my $analysis = $analysis_adaptor->fetch_by_logic_name($GB_SUBMIT_PMATCH_LOGICNAME);
  my $size = $GB_PMATCH_SIZE;
  foreach my $chr(keys %chrhash) {
    #print STDERR "making input ids for ".$chr."\n";
    my $length = $chrhash{$chr};
    if(!$GB_PMATCH_CHUNKS){
      my $input_id = $chr . ".1-" .  $length;
      push(@input_ids, $input_id);
    }else{
      foreach my $chr(keys %chrhash) {
	#print STDERR "making input ids for ".$chr."\n";
	my $length = $chrhash{$chr};
	my $count = 1;
	while ($count < $length) {
	  my $start = $count;
	  my $end   = $count + $size -1;
	  
	  if ($end > $length) {
	    $end = $length;
	  }
	  
	  my $input_id = $chr . "." . $start . "-" .  $end;
	  push(@input_ids, $input_id);
	  $count = $count + $size;
	}
  }
    }
    
  }
  &insert_input_ids($analysis, $db, @input_ids);
}

sub insert_input_ids{
  my ($analysis, $db, @input_ids) = @_;
  
  foreach my $input_id(@input_ids){
    
    my $time = time;
    my $sql = "insert into input_id_analysis(input_id, analysis_id, created) values('$input_id', '".$analysis->dbID."', now())";
    #print STDERR $sql."\n";
    my $sth = $db->prepare($sql);
    $sth->execute();
  }
}

sub get_chrlengths{
  my($db) = @_;

  my $q = "SELECT c.name, max(a.chr_end) 
           FROM   chromosome c, assembly a
           WHERE  c.chromosome_id = a.chromosome_id
           GROUP BY c.name";
  #print STDERR "trying query ".$q."\n";
  my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
  
  while( my ($chr, $length) = $sth->fetchrow_array) {
    $chrhash{$chr} = $length;
  }
  #print "have ".keys(%chrhash)." chromosomes\n";
}

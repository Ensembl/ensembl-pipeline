#!/usr/local/bin/perl -w

=head1 NAME

  populate_analysis_table.pl

=head1 SYNOPSIS
 
  populate_analysis_table.pl 
  populate_analysis_table.pl -insert

=head1 DESCRIPTION
  
  This script generates analysis table entries for the various genebuild databases. It reads
  the entries currently in the reference database GB_DB and generates the other analyses 
  required for all stages of the genebuild by reading from the genebuild configuration files 
  (Config::GeneBuild::* and GeneCombinerConf)

  It prints the sql required to insert these entries to STDOUT, and optionally also inserts 
  them into the various databases. Checks to make sure only one entry inserted per analysis 
  type if the same db is being used for more than one stage!

=head1 OPTIONS
  
  -insert - if this is included on the command line, the script will attempt to insert the 
  lines into the analysis tables of all the gene build databases (ie not GB_DB) specified
  in Config::GeneBuild::Databases

  Options are to be set in GeneBuild config files
  The important ones for this script are:

     GeneBuild::Databases::GB_DBNAME   
     GeneBuild::Databases::GB_DBHOST   
     GeneBuild::Databases::GB_DBUSER   
     GeneBuild::Databases::GB_DBPASS   

     GeneBuild::Databases::GB_GW_DBNAME   
     GeneBuild::Databases::GB_GW_DBHOST   
     GeneBuild::Databases::GB_GW_DBUSER   
     GeneBuild::Databases::GB_GW_DBPASS      

     GeneBuild::Databases::GB_COMB_DBNAME   
     GeneBuild::Databases::GB_COMB_DBHOST   
     GeneBuild::Databases::GB_COMB_DBUSER   
     GeneBuild::Databases::GB_COMB_DBPASS  

     GeneBuild::Databases::GB_cDNA_DBNAME   
     GeneBuild::Databases::GB_cDNA_DBHOST   
     GeneBuild::Databases::GB_cDNA_DBUSER   
     GeneBuild::Databases::GB_cDNA_DBPASS  

     GeneBuild::Databases::GB_FINALDBNAME   
     GeneBuild::Databases::GB_FINALDBHOST   
     GeneBuild::Databases::GB_FINALDBUSER   
     GeneBuild::Databases::GB_FINALDBPASS  

     GeneCombinerConf::FINAL_DBNAME
     GeneCombinerConf::FINAL_DBHOST
     GeneCombinerConf::FINAL_DBUSER
     GeneCombinerConf::FINAL_DBPASS

     GeneBuild::Scripts::GB_LENGTH_RUNNABLES
							   
     GeneCombinerConf::GENECOMBINER_RUNNABLES

=cut

use strict;

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_DBNAME
							     GB_DBHOST
							     GB_DBUSER
							     GB_DBPASS
							     GB_GW_DBNAME
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							     GB_COMB_DBNAME
							     GB_COMB_DBHOST
							     GB_COMB_DBUSER
							     GB_COMB_DBPASS
							     GB_cDNA_DBNAME
							     GB_cDNA_DBHOST
							     GB_cDNA_DBUSER
							     GB_cDNA_DBPASS
							     GB_FINALDBNAME
							     GB_FINALDBHOST
							     GB_FINALDBUSER
							     GB_FINALDBPASS
							    );

use Bio::EnsEMBL::Pipeline::GeneCombinerConf             qw (
							     FINAL_DBNAME
							     FINAL_DBHOST
							     FINAL_DBUSER
							     FINAL_DBPASS
							     GENECOMBINER_RUNNABLES
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts   qw (
							     GB_LENGTH_RUNNABLES
							     GB_PMATCH_RUNNABLES
							    );


use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use English;

$OUTPUT_FIELD_SEPARATOR = "\n";

# global vars
my $insert = 0; # 1 to auto insert
my @analyses = ();
my $insert_prefix = "insert into analysis (analysis_id, created, logic_name, db, db_version, db_file, program, program_version, program_file, parameters, module, module_version, gff_source, gff_feature ) values (";
my $insert_postfix = ");";

# get input options
&GetOptions( 
	    'insert' => \$insert,
	   );

# do something useful

&analyses_from_pipeline;
&analyses_from_config;
print "@analyses\n";;

if($insert){ &insert; }

sub analyses_from_pipeline {
  # get old analyses from GB_DB
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					      '-host'   => $GB_DBHOST,
					      '-user'   => $GB_DBUSER,
					      '-pass'   => $GB_DBPASS,
					      '-dbname' => $GB_DBNAME,
					     );

  my $query = "select * from analysis";

  my $sth = $db->prepare($query);
  my $res = $sth->execute;
  my @fields = qw (logic_name
		   db
		   db_version
		   db_file
		   program
		   program_version
		   program_file
		   parameters
		   module
		   module_version
		   gff_source
		   gff_feature
		  );
  my $counter = 0;
  while (my $row = $sth->fetchrow_hashref){
    $counter++;    

    my $analysis = $insert_prefix;
    $analysis .= $row->{analysis_id}           . ",";
    $analysis .= "'" . $row->{created}         . "',";
    
    foreach my $field(@fields){
      my $string = (defined($row->{$field}) && $row->{$field} ne "") ? $row->{$field} : "NULL";
      $analysis .= "'" . $string . "',";
    }
    
    $analysis =~ s/,$//; # remove trailing comma
    $analysis .= $insert_postfix;
    #print STDERR $analysis."\n";
    push (@analyses, $analysis);
  }
  
}

sub analyses_from_config{

  # track logic names to make sure we don't try to insert with the same one more than once
  my %logic_names;
  
  # from genebuild configs
  
    #print STDERR "genebuild analysis\n";
    foreach my $runnable_hash(@{$GB_LENGTH_RUNNABLES}){
      my $entry = $insert_prefix;
      my $analysis = $runnable_hash->{analysis};
      my $runnable = $runnable_hash->{runnable};
      #print STDERR "running with ".$analysis." ".$runnable."\n";
      if(defined $logic_names{$analysis}){
	my ($p, $f, $l) = caller;
	#print STDERR $f.":".$l."\n";
	print STDERR "Already have an entry for $analysis (" . $logic_names{$analysis} . "); please fix this in the config file and rerun the script!\n";
	exit(0);
      }
      #print STDERR "adding ".$runnable." with key ".$analysis."\n";
      $logic_names{$analysis} = $runnable;
      
      
      $entry .= "'\\N', now(), '$analysis', 'NULL', '1', 'NULL', '$analysis', '1', 'NULL', 'NULL', '$runnable', 'NULL', '$analysis', 'gene'";
      $entry .= $insert_postfix;
      
      push (@analyses, $entry);
    }
 
    
 
    #print STDERR "genecombiner analysis\n";
    # from genecombiner config
    foreach my $runnable_hash(@{$GENECOMBINER_RUNNABLES}){
      my $entry = $insert_prefix;
      my $analysis = $runnable_hash->{analysis};
      my $runnable = $runnable_hash->{runnable};
      #print STDERR "running with ".$analysis." ".$runnable."\n";
      if(defined $logic_names{$analysis}){
	my ($p, $f, $l) = caller;
	#print STDERR $f.":".$l."\n";
	print STDERR "Already have an entry for $analysis (" . $logic_names{$analysis} ."); please fix this in the config file and rerun the script!\n";
	exit(0);
      }
      $logic_names{$analysis} = $runnable;
      
      $entry .= "'\\N', now(), '$analysis', 'NULL', '1', 'NULL', '$analysis', '1', 'NULL', 'NULL', '$runnable', 'NULL', '$analysis', 'gene'";
      $entry .= $insert_postfix;
      
      push (@analyses, $entry);
      
    }
  
    #print STDERR "pmatch analysis\n";
    foreach my $runnable_hash(@{$GB_PMATCH_RUNNABLES}){
      my $entry = $insert_prefix;
      my $analysis = $runnable_hash->{analysis};
      my $runnable = $runnable_hash->{runnable};
      #insert into analysis (analysis_id, created, logic_name, db, db_version, db_file, program, program_version, program_file, parameters, module, module_version, gff_source, gff_feature ) values (
      if(defined $logic_names{$analysis}){
	my ($p, $f, $l) = caller;
	#print STDERR $f.":".$l."\n";
	print STDERR "Already have an entry for $analysis (" . $logic_names{$analysis} . "); please fix this in the config file and rerun the script!\n";
	exit(0);
      }
      $logic_names{$analysis} = $runnable;
      
      $entry .= "'\\N', now(), '$analysis', 'NULL', '1', 'NULL', '$analysis', '1', 'NULL', 'NULL', '$runnable', 'NULL', '$analysis', 'feature'";
      $entry .= $insert_postfix;
      
      push (@analyses, $entry);
      
    }
  
}

sub insert {
  my %dbs; # keep track of host/name pairs to make sure we don't insert into the same db more than once

  # first GeneBuild::Databases
  # I'm sure there's a more automatic way to generate the db details and convert them
  # and avoid all this hard coding. Humph.

  # GB_DB
  $dbs{$GB_DBNAME} = $GB_DBHOST;
  
  # GB_GWDB
  #print STDERR "inserting into genewise db\n";
  if(   (!defined $dbs{$GB_GW_DBNAME}) || 
	(defined $dbs{$GB_GW_DBNAME} && $dbs{$GB_GW_DBNAME} ne $GB_GW_DBHOST)){
    $dbs{$GB_GW_DBNAME} = $GB_GW_DBHOST;
    &insert_into_db($GB_GW_DBNAME, $GB_GW_DBHOST, $GB_GW_DBUSER, $GB_GW_DBPASS);
  }
  
  # GB_COMB_DB
  #print STDERR "inserting into gb comb db\n";
  if(   (!defined $dbs{$GB_COMB_DBNAME}) || 
	(defined $dbs{GB_COMB_DBNAME} && $dbs{$GB_COMB_DBNAME} ne $GB_COMB_DBHOST)){
    $dbs{$GB_COMB_DBNAME} = $GB_COMB_DBHOST;
    &insert_into_db($GB_COMB_DBNAME, $GB_COMB_DBHOST, $GB_COMB_DBUSER, $GB_COMB_DBPASS);
  }
  #print STDERR "inserting into cdna db\n";
  # GB_cDNA_DB - don;t think we need to populate this analysis table ...
  if(   (!defined $dbs{$GB_cDNA_DBNAME}) || 
	(defined $dbs{$GB_cDNA_DBNAME} && $dbs{$GB_cDNA_DBNAME} ne $GB_cDNA_DBHOST)){
    $dbs{$GB_cDNA_DBNAME} = $GB_cDNA_DBHOST;
    &insert_into_db($GB_cDNA_DBNAME, $GB_cDNA_DBHOST, $GB_cDNA_DBUSER, $GB_cDNA_DBPASS);
  }

  # GB_FINALDB
  #print STDERR "inserting into finaldb\n";
  if(   (!defined $dbs{$GB_FINALDBNAME}) || 
	(defined $dbs{$GB_FINALDBNAME} && $dbs{$GB_FINALDBNAME} ne $GB_FINALDBHOST)){
    $dbs{$GB_FINALDBNAME} = $GB_FINALDBHOST;
    &insert_into_db($GB_FINALDBNAME, $GB_FINALDBHOST, $GB_FINALDBUSER, $GB_FINALDBPASS);
  }
  
  # Now GeneCombiner
  #print STDERR "Inserting into GeneCombiner db\n";
  if(   (!defined $dbs{$FINAL_DBNAME}) || 
	(defined $dbs{$FINAL_DBNAME} && $dbs{$FINAL_DBNAME} ne $FINAL_DBHOST)){
    $dbs{$FINAL_DBNAME} = $FINAL_DBHOST;
    &insert_into_db($FINAL_DBNAME, $FINAL_DBHOST, $FINAL_DBUSER, $FINAL_DBPASS);
  }
}

sub insert_into_db{
  my ($name, $host, $user, $pass) = @_;

  if($name eq '' || $host eq '' || $user eq ''){
    print "can't insert into db unless we have name[$name], host[$host] & user[$user] details\n";
    return;
  }
  
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					      '-host'   => $host,
					      '-user'   => $user,
					      '-pass'   => $pass,
					      '-dbname' => $name,
					     );

  #print "inserting into $name @ $host\n";  
  foreach my $analysis(@analyses){
    my $sth = $db->prepare($analysis);
    $sth->execute;
  }
}

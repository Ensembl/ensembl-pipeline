#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::Config::Protein_Annotation::Analysis;
use Bio::EnsEMBL::Pipeline::Config::Protein_Annotation::General;

my ($dbhost, $dbname, $dbpass, $dbuser, $input_id_analysis, $rules, $analysis);

&GetOptions( 
	    'dbname:s' => \$dbname,
	    'dbhost:s' => \$dbhost,
	    'dbuser:s' => \$dbuser,
	    'dbpass:s' => \$dbpass,
	    'input_id_analysis' => \$input_id_analysis,
	    'rules' => \$rules,
	    'analysis' => \$analysis,
	   );


if(!$dbname || !$dbhost || !$dbuser || !$dbpass){
  print STDERR "Usage = -dbname -dbuser -dbhost -dbpass\n";
  exit(0)
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host  => $dbhost,
					    -user  => $dbuser,
					    -dbname=> $dbname,
					    -pass => $dbpass,
					   );

my @analyses = @$PA_ANALYSIS_TYPE;
#print STDERR $PA_SINGLE_DUMMY." ".$PA_CHUNK_DUMMY." ".$PA_PROTEOME_DUMMY."\n";
my $analysis_adaptor = $db->get_AnalysisAdaptor;
if($analysis){
  #print STDERR "have ".@analyses." to store\n";
  foreach my $a_info(@analyses){
    #print "analysis ".$a_info->{'logic_name'}."\n";
    my $a = Bio::EnsEMBL::Analysis->new();
    $a->logic_name($a_info->{'logic_name'});
    $a->module($a_info->{'module'});
    $a->program($a_info->{'program_file'});
    $a->program_file($a_info->{'program_file'});
    $a->db($a_info->{'db'});
    $a->db_file($a_info->{'db_file'});
    $a->gff_source($a_info->{'gff_source'});
    $a->gff_feature($a_info->{'gff_feature'});
    $a->parameters($a_info->{'parameters'});
    
    $analysis_adaptor->store($a);
  }

  my $sql = "insert into analysis(logic_name, created) values(?, now())";

  my $sth = $db->prepare($sql);
  $sth->execute($PA_SINGLE_DUMMY);
  $sth->execute($PA_CHUNK_DUMMY);
  $sth->execute($PA_PROTEOME_DUMMY);
}

if($rules){
  foreach my $a_type(@analyses){
    my $ana_obj = $analysis_adaptor->fetch_by_logic_name($a_type->{'logic_name'});
    if(!$ana_obj){
      print STDERR "the analysis object with this logic_name ".$a_type->{'logic_name'}." doesn't exist you need to specify the -analysis flag if you want the script to import analysis objects as well\n";
      exit(0);
    }
    
    my $dummy_obj = $analysis_adaptor->fetch_by_logic_name($a_type->{'chunk_size'});
    if(!$dummy_obj){
      print STDERR "have ".$ana_obj->logic_name." but no ".$a_type->{'chunk_size'}." can't really create rules\n";
      exit(0);
    }
    my $goal_sql = "insert into rule_goal(goal) values(?)";
    my $condition_sql = "insert into rule_conditions(rule_id, condition) values(?, ?)";

    my $goal_sth = $db->prepare($goal_sql);
    $goal_sth->execute($ana_obj->dbID);
    my $goal_dbID = $goal_sth->{'mysql_insertid'};
    my $condition_sth = $db->prepare($condition_sql);
    $condition_sth->execute($goal_dbID, $dummy_obj->logic_name);
  }
}

my %dummy_hash;

if($input_id_analysis){
  
  
 TYPE:foreach my $a_type(@analyses){
    if($dummy_hash{$a_type->{'chunk_size'}}){
      next TYPE;
    }else{
      my $single;
      if($a_type->{'chunk_size'} eq $PA_SINGLE_DUMMY){
      	if(!$dummy_hash{$a_type->{'chunk_size'}}){
	  $single = $analysis_adaptor->fetch_by_logic_name($PA_SINGLE_DUMMY);
	  $dummy_hash{$single->logic_name} = $single;
	}else{
	  $single = $dummy_hash{$a_type->{'chunk_size'}};
	}
	if(!$single){
	  print STDERR "Dont' have an analysis object for ".$a_type->{'chunk_size'}." can't put in input_id_analysis\n";
	  exit(0);
	}
	my $sql = "insert into input_id_analysis(input_id, analysis_id, created) select translation_id, ".$single->dbID.", now() from translation";
	my $sth = $db->prepare($sql);
	$sth->execute;
      }elsif($a_type->{'chunk_size'} eq $PA_CHUNK_DUMMY){
	if(!$dummy_hash{$a_type->{'chunk_size'}}){
	  $single = $analysis_adaptor->fetch_by_logic_name($PA_CHUNK_DUMMY);
	  $dummy_hash{$single->logic_name} = $single;
	}else{
	  $single = $dummy_hash{$a_type->{'chunk_size'}};
	}
	if(!$single){
	  print STDERR "Dont' have an analysis object for ".$a_type->{'chunk_size'}." can't put in input_id_analysis\n";
	  exit(0);
	  
	}
	my $sql = "insert into input_id_analysis(input_id, analysis_id, created) values(?, ?, now());";
	my $sth = $db->prepare($sql);
	my $dir = $PA_CHUNKS_DIR;


	opendir(DIR, $dir);
	print "opened ".$dir."\n";
	my @allfiles = readdir DIR;
	closedir DIR;
	
	foreach my $f(@allfiles) {
	  if($f eq '.' || $f eq '..'){
	    next;
	  }
	  $sth->execute($f, $single->dbID);
	}
      }elsif($a_type->{'chunk_size'} eq $PA_PROTEOME_DUMMY){
	if(!$dummy_hash{$a_type->{'chunk_size'}}){
	  $single = $analysis_adaptor->fetch_by_logic_name($PA_PROTEOME_DUMMY);
	  $dummy_hash{$single->logic_name} = $single;
	}else{
	  $single = $dummy_hash{$a_type->{'chunk_size'}};
	}
	if(!$single){
	  print STDERR "Dont' have an analysis object for ".$a_type->{'chunk_size'}." can't put in input_id_analysis\n";
	  exit(0);
	}
	my $sql = "insert into input_id_analysis(input_id, analysis_id, created) values(?, ?, now());";
	my $sth = $db->prepare($sql);
	$sth->execute('proteome', $single->dbID);
      }else{
	print STDERR "Don't recognise ".$a_type->{'chunk_size'}." logic name, are you sure you have filled in the config correctly\n";
      }
    }
  }
}

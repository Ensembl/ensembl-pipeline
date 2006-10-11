#!/usr/local/ensembl/bin/perl -w

=pod

=head1 NAME

LowCoverageGeneBuild.pl will:
-----------------------------
-Help to partially automate the Genebuild process of Low-coverage (ie. 2x) genomes.


LowCoverageGeneBuild.pl will not:
---------------------------------
- Tell you which RepeatMasks to use.
  It is up to you to decide which RepeatMasking analyses (RepeatMask, Ab_initio_RepeatMask, Supp_RepeatMask) 
  should be used for the RepeatMask-dependent analyses (Genscan, Unigene, Uniprot, Vertrna).
  Usually RepeatMask and Supp_RepeatMask are sufficient. 
This is an adapted version of Felix's wormbase_to_ensembl.pl 



=head1 SYNOPSIS / OPTIONS

Fill in the necessary variables in the script and in the config file, LowCoverageGeneBuildConf.pm.
Call the script with one of the following options, in the following order:
   setup	 
       get the data files, create the database,load sequence, create input-ids, etc.
   check_seqs   
       checks that seqs in db are like those from BROAD
   run_RepeatMask_and_independent
       use the rulemanager to start running the 3 RepeatMasks and CpG, Dust, Eponine, TRF, tRNAscan
   compare
       compare the three RepeatMasks and allows you to decide which you want to use
   run_RepeatMask_dependent
       use the rulemanager to start running Genscan, Unigene, Uniprot, Vertrna
   dump_RepeatMasked_sequence
   	   produce a RepeatMasked sequence file 
   healthcheck  }
   backup       } -> for your own convenience, you could do this by hand
   clean        }
 

=head1 DESCRIPTION

The script sets up an ensembl database for a low-coverage genome in an automated fashion.
  1 - Assembly data are downloaded manually from BROAD and repeat data are downloaded manually from EBI.
  2 - The appropriate config files should be added to ensembl-config by this script.
  3 - The empty db is created, tables created, rules and analyses enetered, etc
  4 - The rulemanager runs RepeatMasking and RepeatMask-independent analyses
  5 - MANUALLY decide which RepeatMasking to use and fill in $LC_REPMASK_CHOICE (in LowCoverageGeneBuildConf.pm)
  6 - The rulemanager runs RepeatMask-dependent analyses
  7 - Back-up of db
  8 - Healthchecks
  9 - Clean up
  
Before starting the GeneBuild you will need to manually download the assembly from BROAD
and the repeat libraries from EBI, as outlined in README.


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut


#this program uses scripts from ensembl-personal - before cvs commiting this need to copy those scripts somewhere else and change the paths!

#have added an automatic taxonomy update bit - have not tested or committed it yet!

use strict;
use warnings;
use LowCoverageGeneBuildConf;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
$| = 1;

my $cmd;
my $status;
my $input;
my $comm = $ARGV[0];
my $perllib;
my $path;

if(!$comm or ($comm ne "setup" && $comm ne "check_seqs" && $comm ne "run_RepeatMask_and_independent" && $comm ne "compare"
              && $comm ne "run_RepeatMask_dependent" && $comm ne "dump_RepeatMasked_sequence" && $comm ne "healthcheck" && $comm ne "backup" && $comm ne "clean")){
  exec('perldoc', $0);
  exit 1;

} 
	 
set_env(1);

if ($comm eq "setup"){
  print STDERR "\nSet-up will make the necessary ensembl-config files, create your database, load sequence, set toplevel, create input-ids, ".
               "load_taxonomy, and load_analysis_descriptions.\n\n";
 
  #Get files  
  print STDERR "\nPreparing for setup... Please check that you have:\n".
        "- updated your CVS co\n".
        "- Manually downloaded assembly and repeat data (see README)".
        "- modified LowCoverageGeneBuildConf.pm ***\n".
	"- ssh bc-dev-64\n";


  #Make config dir in ensembl-config:
  #---------------------------------
  print STDERR "\n>> Making config dir in $LC_cvsDIR...\n";
  config_setup();  
  print STDERR "   ...config setup completed.\n";  

  #Creating new db
  #---------------
  print STDERR "\n** Would you like to DROP DATABASE -D $LC_DBNAME -h $LC_DBHOST -P $LC_DBPORT if it exists? \nyes or no:";
  my $drop = <STDIN>;  
  if ($drop =~ /yes/i){
    createDB();
  } else {
    print STDERR "You chose not to drop the database\n";
    exit 1;
  }
  
   
  #Loading assembly
  #----------------
  #Firstly, load the sequence in the given file (assembly.bases) into the database under a coord system called contig
  print STDERR "\n>> Loading sequence in $LC_ASSEMBLY into $LC_DBNAME...\n";
  print STDERR "  perl $LC_cvsDIR/ensembl-pipeline/scripts/load_seq_region.pl -dbhost $LC_DBHOST ".
               "-dbuser $LC_DBUSER -dbport $LC_DBPORT -dbpass $LC_DBPASS -dbname $LC_DBNAME ".
               "-coord_system_name contig -rank 2 -sequence_level -fasta_file $LC_ASSEMBLY\n";
  $cmd = "perl ".$LC_cvsDIR."/ensembl-pipeline/scripts/load_seq_region.pl ".
         "-dbhost $LC_DBHOST -dbuser $LC_DBUSER -dbport $LC_DBPORT -dbpass $LC_DBPASS -dbname $LC_DBNAME ". 
         "-coord_system_name contig -rank 2 -sequence_level -fasta_file $LC_ASSEMBLY";  
  if(system($cmd)){
    die "\n**Error with first load_seq_region**\n";
    exit 1;
  } else {
    print STDERR "   ...First load_seq_region successful\n";
  }
  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBro -P$LC_DBPORT -D$LC_DBNAME -e 'SELECT count(*) FROM seq_region'\n";
  system("mysql -h$LC_DBHOST -u$LC_DBro -P$LC_DBPORT -D$LC_DBNAME -e 'SELECT count(*) FROM seq_region'") 
         && warn "\nCan't count 1st seq_regions\n";
    
  
  #Secondly, load the assembled pieces from the agp file into the seq_region table
  #-------------------------------------------------------------------------------
  print STDERR "\n>> Loading assembled pieces from agp file into seq_region table...\n";
  print STDERR "  perl $LC_cvsDIR/ensembl-pipeline/scripts/load_seq_region.pl -dbhost $LC_DBHOST ".
                "-dbuser $LC_DBUSER -dbport $LC_DBPORT -dbpass $LC_DBPASS -dbname $LC_DBNAME ".
                "-coord_system_name scaffold -coord_system_version $LC_DEFAULT -rank 1 -agp_file $LC_ASSEMBLY_AGP";
  $cmd = "perl ".$LC_cvsDIR."/ensembl-pipeline/scripts/load_seq_region.pl ".
         "-dbhost $LC_DBHOST -dbuser $LC_DBUSER -dbport $LC_DBPORT -dbpass $LC_DBPASS -dbname $LC_DBNAME ". 
         "-coord_system_name scaffold -coord_system_version $LC_DEFAULT -rank 1 -agp_file $LC_ASSEMBLY_AGP";  
  if(system($cmd)){
    die "\n**Error with second load_seq_region**\n";
    exit 1;
  } else {
    print STDERR "   ...Second load_seq_region successful\n";
  }  
  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBro -P$LC_DBPORT -D$LC_DBNAME -e 'SELECT count(*) FROM seq_region'\n";
  system("mysql -h$LC_DBHOST -u$LC_DBro -P$LC_DBPORT -D$LC_DBNAME -e 'SELECT count(*) FROM seq_region'") 
         && warn "\nCan't count 2nd seq_regions\n";  
 
 
  #A couple of checks and changes:
  #------------------------------
  print STDERR "  mysql -D $LC_DBNAME -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT  -e 'SELECT * FROM coord_system'\n";
  system("mysql -D $LC_DBNAME -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT  -e 'SELECT * FROM coord_system'") 
         && warn "\nError with SELECT * from coord_system\n- -D $LC_DBNAME -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT\n";    

  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
               "'UPDATE coord_system SET attrib 'default_version,sequence_level' WHERE coord_system_id=1\n";
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e 'UPDATE coord_system ".
         "SET attrib=\"default_version,sequence_level\" WHERE coord_system_id=1'") 
         && warn "\nError with UPDATE coord_system\n";  

  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
               "'UPDATE coord_system SET attrib= 'default_version' WHERE coord_system_id=2'\n";
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e 'UPDATE coord_system ".
         "SET attrib=\"default_version\" WHERE coord_system_id=2'") 
         && warn "\nError with UPDATE coord_system\n";    

  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
               "'UPDATE coord_system SET version=NULL WHERE coord_system_id=1'\n";
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e 'UPDATE coord_system ".
         "SET version=NULL WHERE coord_system_id=1'") 
         && warn "\nError with UPDATE coord_system\n";  

  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
               "'INSERT INTO meta(meta_key,meta_value) VALUES ('assembly.default','$LC_DEFAULT')'\n";
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
         "\"INSERT INTO meta(meta_key,meta_value) VALUES ('assembly.default','$LC_DEFAULT')\"") 
           && warn "\nError with INSERT INTO meta\n";
 
   print STDERR "  mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
     "'INSERT INTO meta(meta_key,meta_value) VALUES ('assembly.date','$LC_ASSEMBLY_DATE')'\n";
   system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
     "\"INSERT INTO meta(meta_key,meta_value) VALUES ('assembly.date','$LC_ASSEMBLY_DATE')\"")
     && warn "\nError with INSERT INTO meta\n";

  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
    "'INSERT INTO meta(meta_key,meta_value) VALUES ('assembly.name','$LC_ASSEMBLY_NAME')'\n";
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e ".
    "\"INSERT INTO meta(meta_key,meta_value) VALUES ('assembly.name','$LC_ASSEMBLY_NAME')\"")
    && warn "\nError with INSERT INTO meta\n";

  #load genebuild info into meta table:
  #--------------------------------------
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e \"INSERT INTO ".
         "meta (meta_key, meta_value) VALUES ('genebuild.version','$LC_DATE')\"")
         && warn "\nProblem with loading genebuild version\n";
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME -e \"INSERT INTO ".
         "meta (meta_key, meta_value) VALUES ('genebuild.id',$LC_GeneBuilderID)\"")
         && warn "\nProblem with loading genebuild id\n";


  system("mkdir -p $LC_workDIR/assembly");

  #Loading a standard agp file into the assembly table
  #---------------------------------------------------
  print STDERR "\n>> Loading agp file into the assembly table...\n";
  print STDERR "$LC_cvsDIR/ensembl-pipeline/scripts/load_agp.pl ".
               "-dbhost $LC_DBHOST -dbuser $LC_DBUSER -dbport $LC_DBPORT -dbname $LC_DBNAME -dbpass $LC_DBPASS ".
               "-assembled_name scaffold -assembled_version $LC_DEFAULT -component_name contig -agp_file $LC_ASSEMBLY_AGP >&  $LC_workDIR/assembly/load_agp.log\n";
  $cmd = "$LC_cvsDIR/ensembl-pipeline/scripts/load_agp.pl -dbhost $LC_DBHOST -dbuser ".
         "$LC_DBUSER -dbport $LC_DBPORT -dbname $LC_DBNAME -dbpass $LC_DBPASS ".
	 "-assembled_name scaffold -assembled_version $LC_DEFAULT -component_name contig -agp_file $LC_ASSEMBLY_AGP >&  $LC_workDIR/assembly/load_agp.log"; 
  if(system($cmd)){
    die "\nError with load_agp.log\n";
    exit 1;
  } else {
    print STDERR "Submitted to bsub: load_agp.log\n\n";
  }  


 #insert a break or sleep here to wait for bsub
 #--------------------------------------------
  print STDERR "You need to wait for the job on the farm to be completed before you can continue.\n".
                "To check if the job is finished, type 'bjobs' in another window.\n".
		"WHEN THE JOB IS FINISHED, type 'yes':";
  $input = <STDIN>;  
  while ($input !~ /yes/i){
    print STDERR "\n** Type yes when you are ready to continue...  : ";
    $input = <STDIN>;
  }		


  print STDERR "  mysql -h$LC_DBHOST -u$LC_DBro -P$LC_DBPORT -D$LC_DBNAME -e 'SELECT count(*) FROM assembly'\n";
  system("mysql -h$LC_DBHOST -u$LC_DBro -P$LC_DBPORT -D$LC_DBNAME -e 'SELECT count(*) FROM assembly'") 
         && warn "\nCan't count assembly\n";  
    	
	  
  #Flag the set of non-redundant seq_regions as 'toplevel'.  
  #------------------------------------------------------
  print STDERR "\n>> Setting toplevel (the quick way)...\n";
#  print STDERR "  bsub -q normal -o $LC_workDIR/assembly/set_toplevel.log perl $LC_cvsDIR/ensembl-pipeline/scripts/set_toplevel.pl". 
#               " -dbhost $LC_DBHOST -dbuser $LC_DBUSER -dbport $LC_DBPORT -dbname $LC_DBNAME -dbpass $LC_DBPASS \n";
#  $cmd = "bsub -q normal -o $LC_workDIR/assembly/set_toplevel.log perl  ".
#         "$LC_cvsDIR/ensembl-pipeline/scripts/set_toplevel.pl -dbhost $LC_DBHOST -dbuser $LC_DBUSER ".
#         "-dbport $LC_DBPORT -dbname $LC_DBNAME -dbpass $LC_DBPASS";
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -p$LC_DBPASS -P$LC_DBPORT -D$LC_DBNAME -e 'INSERT IGNORE INTO attrib_type ( code, name, description ) VALUES ( \"toplevel\",\"Top Level\",\"Top Level Non-Redundant Sequence Region\")'") && warn "\nCan't load attrib_type for toplevel\n";  
  system("mysql -h$LC_DBHOST -u$LC_DBUSER -p$LC_DBPASS -P$LC_DBPORT -D$LC_DBNAME ".
         "-e 'INSERT INTO seq_region_attrib " .
             "SELECT distinct( sr.seq_region_id ), at.attrib_type_id, 1 " .
             "FROM   seq_region sr LEFT JOIN assembly a " .
             "ON sr.seq_region_id=a.cmp_seq_region_id, attrib_type at " .
             "WHERE  cmp_seq_region_id is NULL " .
             "AND    at.code=\"toplevel\"'"
        ) && die "Failed setting toplevel";
#  if(system($cmd)){
#    die "\nError submitting set_toplevel\n";
#    exit 1;
#  } else {
#    print STDERR "Submitted to bsub: set_toplevel";
#  }  
#
#
#  #insert a break or sleep here to wait for bsub
#  #---------------------------------------------
#  print STDERR "You need to wait for the 2nd job on the farm to be completed before you can continue.\n".
#                "To check if the job is finished, type 'bjobs' in another window.\n".
#		"WHEN THE JOB IS FINISHED, type 'yes':";
#  $input = <STDIN>;  
#  while ($input !~ /yes/i){
#    print STDERR "\n** Type yes when you are ready to continue...  : ";
#    $input = <STDIN>;
#  }
  
  #Setting up the analysis and rule tables
  #---------------------------------------
  print STDERR "\n>> Setting up the analysis and rule tables...\n";
  $cmd = "perl $LC_cvsDIR/ensembl-pipeline/scripts/analysis_setup.pl -dbname $LC_DBNAME -dbhost ".
         "$LC_DBHOST -dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT -read -file ".
	 "$LC_cvsDIR/ensembl-config/$LC_SPECIES/$LC_BUILD_VERSION/pipe_conf/analysis.conf";
  if(system($cmd)){
    die "\nError with analysis_setup\n\n";
    exit 1;
  }    
  $cmd = "perl $LC_cvsDIR/ensembl-pipeline/scripts/rule_setup.pl -dbname $LC_DBNAME -dbhost ".
         "$LC_DBHOST -dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT -read -file ".
	 "$LC_cvsDIR/ensembl-config/$LC_SPECIES/$LC_BUILD_VERSION/pipe_conf/rule.conf";
  if(system($cmd)){
    die "\n**Error with rule_setup**\n\n";
    exit 1;
  }  
  
  
  #checking the setup of analysis, rule, rule_goal, rule_condition and meta tables
  #---------------------------------------------------------------------------------
  print STDERR "\n>> Checking the setup of analysis, rule, rule_goal, rule_condition and meta tables...\n";
  print STDERR "  mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'SELECT ".
               "analysis_id, logic_name, program_file FROM analysis'\n";
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'SELECT analysis_id, logic_name, program_file FROM analysis'") 
         && warn "\nProblem with MySQL query (analysis table)\n";

  print STDERR "  mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'SELECT analysis.analysis_id, ".
               "input_id_type_analysis.input_id_type, analysis.logic_name, analysis.program_file FROM analysis, input_id_type_analysis ".
	       "WHERE analysis.analysis_id=input_id_type_analysis.analysis_id'\n" ;
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'SELECT analysis.analysis_id, ".
         "input_id_type_analysis.input_id_type, analysis.logic_name, analysis.program_file FROM analysis, ".
	 "input_id_type_analysis WHERE analysis.analysis_id=input_id_type_analysis.analysis_id'")
         && warn "\nProblem with MySQL query (analysis & input_id_analysis table)\n";
 
  print STDERR "  $LC_cvsDIR/ensembl-pipeline/scripts/RuleHandler.pl -dbname $LC_DBNAME -dbhost ".
               "$LC_DBHOST -dbuser $LC_DBro -dbport $LC_DBPORT -rules\n";
  system("$LC_cvsDIR/ensembl-pipeline/scripts/RuleHandler.pl -dbname $LC_DBNAME -dbhost $LC_DBHOST -dbuser $LC_DBro -dbport $LC_DBPORT -rules")
         && warn "\nProblem with RuleHandler\n";
 
  print STDERR "  mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from rule_goal'\n";
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from rule_goal'") 
         && warn "\nProblem with MySQL query (rule_goal table)\n";

  print STDERR "  mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from rule_conditions'\n";	 
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from rule_conditions'") 
         && warn "\nProblem with MySQL query (rule_condition)\n";
	 
  print STDERR "  mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from meta'\n";	 
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from meta'") 
         && warn "\nProblem with MySQL query (meta table)\n";

  print STDERR "  mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e ".
               "'select distinct(analysis.logic_name),analysis.analysis_id from repeat_feature, analysis ".
	       "where repeat_feature.analysis_id=analysis.analysis_id'\n";
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e ".
      "'select distinct(analysis.logic_name),analysis.analysis_id from repeat_feature, ".
      "analysis where repeat_feature.analysis_id=analysis.analysis_id'");


  #Load analysis descriptions
  #-------------------------
  print STDERR ">> Now loading the analysis descriptions table...\n";			  
  system("perl $LC_cvsDIR/ensembl-personal/lec/rhesus/scripts/load_analysis_descriptions.pl ".
	 "-dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbhost $LC_DBHOST -dbport $LC_DBPORT ".
	 "-dbname $LC_DBNAME $LC_cvsDIR/ensembl-personal/lec/rhesus/scripts/analysis.descriptions");

  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select logic_name, ".
         "description from analysis left join analysis_description on analysis.analysis_id = analysis_description.analysis_id'") 
         && warn "\nProblem with MySQL query (analysis_description table)\n";
  
  print STDERR "   ...Done loading the analysis descriptions.\n";
  print STDERR "\nAll of your analyses should have descriptions (except for SubmitContig)-\n".
  		"if not, insert descriptions into $LC_cvsDIR/ensembl-personal/lec/rhesus/scripts/analysis.descriptions\n".
		"and run $LC_cvsDIR/ensembl-personal/lec/rhesus/scripts/load_analysis_descriptions.pl\n"; #copy above command  


  #Load taxonomy and other info into meta table
  #---------------------------------------------
  print STDERR ">> Loading taxonomy and meta table information...";
  system("perl $LC_cvsDIR/ensembl-pipeline/scripts/load_taxonomy.pl -name '$LC_NAME' ".
         "-taxondbhost $TAXON_DBHOST -taxondbport $TAXON_DBPORT -taxondbname $TAXON_DBNAME ".
	 "-lcdbhost $LC_DBHOST -lcdbport $LC_DBPORT -lcdbname $LC_DBNAME -lcdbuser $LC_DBUSER -lcdbpass $LC_DBPASS");
  print STDERR "   ...Done loading taxonomy and meta table information.\n";


  #Make input_ids 
  #-------------
  print STDERR "\n>> Making input_ids...\n"; 
  print STDERR "  perl $LC_cvsDIR/ensembl-pipeline/scripts/make_input_ids -dbhost $LC_DBHOST ".
               "-dbuser $LC_DBUSER -dbport $LC_DBPORT -dbname $LC_DBNAME -dbpass $LC_DBPASS ".
	       "-seq_level -logic_name SubmitContig\n";
  $cmd = "perl $LC_cvsDIR/ensembl-pipeline/scripts/make_input_ids -dbhost $LC_DBHOST ".
         "-dbuser $LC_DBUSER -dbport $LC_DBPORT -dbname $LC_DBNAME -dbpass $LC_DBPASS -seq_level -logic_name SubmitContig";
  if(system($cmd)){
    die "\nError with making input ids\n";
    exit 1;
  } else {
    print STDERR "   ...Done making input_ids.\n";
  }
        
  print STDERR "\n\nSETUP IS FINISHED. (Now run 'check_seqs')\n\n";
  
  #Change config files to be user-specific and species-specific:
  #The config files import LowCoverageGeneBuildConf.pm and so access all
  #user-specific and species-specific information from this file.
  print STDERR  "Please run the following (manually on the command-line) to make the required output directories:\n\n";
  print STDERR  ">cd $LC_cvsDIR"."/ensembl-config/$LC_SPECIES/$LC_BUILD_VERSION/Bio/EnsEMBL/Pipeline/Config/ \n";
  print STDERR	">foreach dir ( `grep -i -e 'output_dir' BatchQueue.pm | gawk '{print \$3}'  | sed -e s/,//  | sed -e s/\\\"//g | sed -e s':\\.:\/:' | sed -e s':\$LC_scratchDIR:$LC_scratchDIR:'| sort -u` )\n"
				."	mkdir -p \$dir\n"
				."	end\n"
				."Don't forget to move back to this directory when you're ready to run the next step.\n";


} elsif ($comm eq "check_seqs"){
  print STDERR "mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME  -e 'SELECT * FROM seq_region WHERE coord_system_id=2 limit 10'\n";
  $cmd = "mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME  -e 'SELECT * FROM seq_region WHERE coord_system_id=2 limit 10'";
  if(system($cmd)){
    warn "\nError with checking seq_regions\n";
    exit 1;
  }
  

  #Choose 2 random scaffolds, in this case 1 and 1054.
  #--------------------------------------------------
  #Ideally these scaffolds should have more than one contig each.
  #To look for a scaffold that has more than one contig: grep -w scaffold_1054 assembly.agp
  print STDERR "\n  You may need to adjust the parsing of headers $LC_workDIR/assembly/assembly.agp.\n".
               "  Headers of the format '>scaffold_100.1-945' should work\n";
  print STDERR "\n\n** Please check manually for 2 scaffolds that have more than one contig each.\n".
               "eg. To see if scaffold 1054 has more than one contig, do 'grep -w scaffold_1054 $LC_workDIR/assembly/assembly.agp | wc -l' in a separate window\n".
	       "* Now please enter the two chosen scaffold numbers...\n";
  print STDERR "Enter first scaffold number:  ";
  my $number1 = <STDIN>;
  $number1 =~ s/\s//g;
  $number1 =~ s/\n//g;
  chomp($number1);
  print STDERR "\nEnter second scaffold number:  ";
  my $number2 = <STDIN>;
  $number2 =~ s/\s//g;
  $number2 =~ s/\n//g;
  chomp($number2);

  system("mkdir -p $LC_workDIR/seqdump");
    
  $status = 0;
  print STDERR "\n>> Dumping scaffold $number1 from db into $LC_workDIR/seqdump/toplevel::scaffold_$number1...\n";
  $cmd = "$LC_cvsDIR/ensembl-personal/searle/scripts/fetch_slice_seq.pl -host $LC_DBHOST -user $LC_DBro ".
         "-port $LC_DBPORT -dbname $LC_DBNAME -path $LC_DEFAULT toplevel::scaffold_".$number1." -outfile ".
	 "$LC_workDIR/seqdump/toplevel::scaffold_".$number1;
  print STDERR "$cmd\n";
  $status += system($cmd);
  if($status){ warn("\nError with fetch_slice_seq.pl, scaffold $number1") }
  
  print STDERR "\n>> Dumping scaffold $number2 from db into $LC_workDIR/seqdump/toplevel::scaffold_$number2...\n";
  $status = 0;
  $cmd = "$LC_cvsDIR/ensembl-personal/searle/scripts/fetch_slice_seq.pl -host $LC_DBHOST -user ".
         "$LC_DBro -port $LC_DBPORT -dbname $LC_DBNAME -path $LC_DEFAULT toplevel::scaffold_".$number2.
	 " -outfile $LC_workDIR/seqdump/toplevel::scaffold_".$number2;
  print STDERR "$cmd\n";
  $status += system($cmd);
  if($status){ warn("\nError with fetch_slice_seq.pl, scaffold $number2") }

  
  #Reformat sequence to 60 characters per row

  print STDERR ">> Reformatting $LC_SCAFFOLDS to 60 characters per row...\n";
  $status = 0;
  $status += system("$LC_cvsDIR/ensembl-personal/searle/scripts/reformat_fasta.pl $LC_SCAFFOLDS > $LC_workDIR/assembly/assembly.agp.fasta.60");
  if($status){ warn("\nError with reformatting $LC_SCAFFOLDS to 60 characters per row\n") }  
 
 
  #Find the id for scaffold:
  #-------------------------
  my %scaff_hash;
  @{$scaff_hash{$number1}} = @{find_matching_headers($number1)}; 
  @{$scaff_hash{$number2}} = @{find_matching_headers($number2)};
 

  #Now find the seqs:   
  #-----------------
  if (@{$scaff_hash{$number1}} != 1) {
    print "Problem with scaffold $number1 - no hits or more than 1 hit\n";
    exit 1;
  } 
  print STDERR ">> Finding fasta $scaff_hash{$number1}[0] in $LC_workDIR/assembly/assembly.agp.fasta.60...\n";
  $status = 0;
  $status += system("$LC_cvsDIR/ensembl-personal/searle//scripts/find_seq_in_fasta.pl -id $scaff_hash{$number1}[0] ".
                    "-file $LC_workDIR/assembly/assembly.agp.fasta.60 > $LC_workDIR/seqdump/scaffold_".$number1.".fa");
  if($status){ warn("\nError with finding $number1\n") }

  if (@{$scaff_hash{$number2}} != 1) {
        print "Problem with scaffold $number2 - no hits or more than 1 hit\n";
	exit 1;
    }
  print STDERR ">> Finding fasta $scaff_hash{$number2}[0] in $LC_workDIR/assembly/assembly.agp.fasta.60...\n";
  $status = 0;
  $status += system("$LC_cvsDIR/ensembl-personal/searle/scripts/find_seq_in_fasta.pl -id $scaff_hash{$number2}[0] ".
                    "-file $LC_workDIR/assembly/assembly.agp.fasta.60 > $LC_workDIR/seqdump/scaffold_".$number2.".fa");
  if($status){ warn("\nError with finding $number2\n") }


  #Check whether the seq in the db and the seq in the fasta file are identical:
  #---------------------------------------------------------------------------
  print STDERR "\n>> Diffing scaffold $number1...\n";
  system("diff $LC_workDIR/seqdump/scaffold_$number1.fa $LC_workDIR/seqdump/toplevel::scaffold_$number1");
  
  print STDERR "\n>> Diffing scaffold $number2...\n";
  system("diff $LC_workDIR/seqdump/scaffold_$number2.fa $LC_workDIR/seqdump/toplevel::scaffold_$number2");

  print STDERR "\nThe scaffolds should only differ on the header lines, the sequences should be identical.\n\n".
               "CHECK_SEQS IS FINISHED. (Now check for enough diskspace on $LC_scratchDIR, ssh into bc-dev-64, ".
	       "and run 'run_RepeatMask_and_independent' from there)\n";
  
  
} elsif ($comm eq "run_RepeatMask_and_independent"){

  print STDERR "It is recommended that you test each of the following analyses once to check for failures:\n".
               "RepeatMask, Ab_initio_RepeatMask, Supp_RepeatMask, CpG, Dust, Eponine, TRF, tRNAscan\n\n".
	       "An example of the command to use is:\n".
	       "$LC_cvsDIR/ensembl-analysis/scripts/test_RunnableDB  -dbname $LC_DBNAME -dbhost $LC_DBHOST ".
	       "-dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT ".
	       "-logic Supp_RepeatMask -input_id contig::contig_10030:1:5096:1\n";

  print STDERR "\n** Are you ready to start the rulemanager for RepeatMask_and_Independent analyses? \nyes or no:";
  my $ready = <STDIN>;  
  if ($ready =~ /yes/i){
	
    print STDERR "\nRulemanager will now run RepeatMask, Ab_initio_RepeatMask, Supp_RepeatMask, CpG, Dust, Eponine, TRF, tRNAscan...\n".
               "See file ".$LC_scratchDIR."/".$LC_BUILD_VERSION."_repeatmasking.out for errors.\n".
               "Don't forget to monitor your scripts.\n";
    system("perl $LC_cvsDIR/ensembl-pipeline/scripts/rulemanager.pl -dbname $LC_DBNAME -dbhost ".
         "$LC_DBHOST -dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT ".
         "-analysis RepeatMask -analysis Supp_RepeatMask -analysis Ab_initio_RepeatMask ".
         "-analysis Eponine -analysis CpG -analysis Dust -analysis TRF -analysis tRNAscan -shuffle ".
	 ">& ".$LC_scratchDIR."/".$LC_BUILD_VERSION."_repeatmasking.out ");

    print STDERR "\n\nrun_RepeatMask_and_independent appears to be finished\nHOWEVER it is a good idea to check the jobs table "
	      ."to ensure that all of the jobs have run\nand check ".$LC_scratchDIR."/".$LC_BUILD_VERSION."_repeatmasking.out too\n\n"
	      ."If there are still jobs, *wait until those which are running are finished*\nthen empty the job and job_status tables," 
	      ." remove the pipeline lock\nand rerun 'run_RepeatMask_and_independent'\n"
              ."Once all of the jobs have finished you are ready to run 'compare'\n";
	  #print STDERR "\n\nrun_RepeatMask_and_independent IS FINISHED. (Now run 'compare')\n";
  }else{
    print STDERR "RepeatMask_and_independent analyses have not been run\n\n";
    exit;
  }
     
} elsif ($comm eq "compare"){  

  print STDERR "\n>> Checking repeat types...\n";
  print STDERR "perl $LC_cvsDIR/ensembl/misc-scripts/repeats/repeat-types.pl -user ".
               "$LC_DBUSER -pass $LC_DBPASS -host $LC_DBHOST -port $LC_DBPORT -dbpattern $LC_DBNAME\n";
  system("perl $LC_cvsDIR/ensembl/misc-scripts/repeats/repeat-types.pl -user ".
         "$LC_DBUSER -pass $LC_DBPASS -host $LC_DBHOST -port $LC_DBPORT -dbpattern $LC_DBNAME");
  
  print STDERR "mysql -uensro -hLC_DBHOST -D$LC_DBNAME -e 'select count(*), repeat_type from repeat_consensus group by repeat_type'";
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select count(*), ".
         "repeat_type from repeat_consensus group by repeat_type' ");

  print STDERR "mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select analysis_id, ".
               "logic_name, program_file from analysis'\n";
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select analysis_id, logic_name, program_file from analysis'");

  print STDERR "mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from meta_coord'\n";
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select * from meta_coord' ");

  
  system("mkdir -p $LC_workDIR/repeat_libraries");

  print STDERR ">> Submitting jobs to check repeat coverage (output at $LC_workDIR/repeat_libraries/*.out)...\n";
  system("bsub -q normal -R linux -o $LC_workDIR/repeat_libraries/RepeatMask.out ".
         "$LC_cvsDIR/ensembl-personal/searle/monodelphis1/scripts/repeat_coverage.pl ".
         "-dbname $LC_DBNAME -host $LC_DBHOST -port $LC_DBPORT -repeat RepeatMask -path $LC_DEFAULT");
  
  system("bsub -q normal -R linux -o $LC_workDIR/repeat_libraries/Ab_initio_RepeatMask.out ".
         "$LC_cvsDIR/ensembl-personal/searle/monodelphis1/scripts/repeat_coverage.pl ".
         "-dbname $LC_DBNAME -host $LC_DBHOST -port $LC_DBPORT -repeat Ab_initio_RepeatMask -path $LC_DEFAULT");
  
  system("bsub -q normal -R linux -o $LC_workDIR/repeat_libraries/Supp_RepeatMask.out ".
         "$LC_cvsDIR/ensembl-personal/searle/monodelphis1/scripts/repeat_coverage.pl ".
         "-dbname $LC_DBNAME -host $LC_DBHOST -port $LC_DBPORT -repeat Supp_RepeatMask -path $LC_DEFAULT");

  system("bsub -q normal -R linux -o $LC_workDIR/repeat_libraries/Ab_and_supp_RepeatMask.out ".
         "$LC_cvsDIR/ensembl-personal/searle/monodelphis1/scripts/repeat_coverage.pl ".
         "-dbname $LC_DBNAME -host $LC_DBHOST -port $LC_DBPORT -repeat Ab_initio_RepeatMask -repeat Supp_RepeatMask -path $LC_DEFAULT");
  
  system("bsub -q normal -R linux -o $LC_workDIR/repeat_libraries/Supp_and_orig_RepeatMask.out ".
         "$LC_cvsDIR/ensembl-personal/searle/monodelphis1/scripts/repeat_coverage.pl ".
         "-dbname $LC_DBNAME -host $LC_DBHOST -port $LC_DBPORT -repeat RepeatMask -repeat Supp_RepeatMask -path $LC_DEFAULT");
  
  system("bsub -q normal -R linux -o $LC_workDIR/repeat_libraries/Ab_and_orig_RepeatMask.out ".
         "$LC_cvsDIR/ensembl-personal/searle/monodelphis1/scripts/repeat_coverage.pl ".
         "-dbname $LC_DBNAME -host $LC_DBHOST -port $LC_DBPORT -repeat RepeatMask -repeat Ab_initio_RepeatMask -path $LC_DEFAULT");
  
  system("bsub -q normal -R linux -o $LC_workDIR/repeat_libraries/all_three_RepeatMask.out ".
         "$LC_cvsDIR/ensembl-personal/searle/monodelphis1/scripts/repeat_coverage.pl ".
         "-dbname $LC_DBNAME -host $LC_DBHOST -port $LC_DBPORT -repeat RepeatMask -repeat Ab_initio_RepeatMask -repeat Supp_RepeatMask -path $LC_DEFAULT");
    
  system("$LC_cvsDIR/ensembl-personal/searle/scripts/countlc.pl $LC_SCAFFOLDS > $LC_workDIR/repeat_libraries/scaffolds.lowercase.out");

  print STDERR "You need to wait for the 7 jobs on the farm to be completed before you can continue.\n".
                "To check if the jobs are finished, type 'bjobs' in another window.\n".
		"** WHEN THE JOBS ARE FINISHED, type 'yes':";
  $input = <STDIN>;  
  while ($input !~ /yes/i){
    print STDERR "\n** Type yes when you are ready to continue...  : ";
    $input = <STDIN>;
  }
  
  print STDERR "\n\nCOMPARE IS FINISHED. \nYou need to now: \n".
               "(1) analyse results\n ". 
	       "    This this done by seeing how many bases are masked out by the ".
	       "different combinations of repeat masking compared to the total number of bases.\n".
	       "    To see the total number of bases:\n    more $LC_workDIR/repeat_libraries/scaffolds.lowercase.out\n".
	       "    To see how many bases are masked:\n    foreach file ($LC_workDIR/repeat_libraries/*.out)\n".
	       "      echo \$file\n      grep \^\"Total masked\" \$file\n    end\n".
               "(2) decide which combination of the 3 RepeatMasks to use (look at the 'total masked' number at end of *.out files)\n".
	       "(3) fill in LC_REPMASK_CHOICE in LowCoverageGeneBuildConf.pm".
	       "(4) and then run 'run_RepeatMask_dependent'\n";
  
} elsif ($comm eq "run_RepeatMask_dependent"){
  print STDERR "\n\nBefore you run the rulemanager you should test each of the analyses ".
               "(Genscan, Unigene, Uniprot and Vertrna) individually\n".
	       "now, this is not straightforward as you can only test Uniprot, Unigene and Vertrna once some Genscans have finished\n".
	       "so, test Genscan first like this:\n".
               "Use this command another window - make sure your Perl5Lib is pointing to the right place - \n".
	       "$LC_cvsDIR/ensembl-analysis/scripts/test_RunnableDB -dbname $LC_DBNAME -dbhost $LC_DBHOST ".
	       "-dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT -logic Genscan -input_id contig::contig_10030:1:5096:1\n".
               "\nThen set the rulemamager off like this *but kill it once some jobs have finished*:\n".
               "perl $LC_cvsDIR/ensembl-pipeline/scripts/rulemanager.pl -dbname $LC_DBNAME -dbhost ".
	       "$LC_DBHOST -dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT ".
               "-analysis Genscan >& $LC_scratchDIR/".$LC_BUILD_VERSION."_repeatmask_dependents.out &\n".
               "\nThen you can test the other analyses in turn using a contig which has been run through Genscan eg:\n".
               "$LC_cvsDIR/ensembl-analysis/scripts/test_RunnableDB -dbname $LC_DBNAME -dbhost $LC_DBHOST ".
               "-dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT -logic Uniprot -input_id contig::contig_10030:1:5096:1\n".
               "\nWhen they have all been tested you can run the rulemanager...\n";
  print STDERR "\nHere is a list of some contigs on which genscan has already been run (maybe empty if this is your first run of this step)\n";
  print STDERR "mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select distinct(sr.name),count(*) as npred,sr.length from prediction_transcript pt, seq_region sr where sr.seq_region_id=pt.seq_region_id limit 10'\n";
  system("mysql -h $LC_DBHOST -u $LC_DBro -P $LC_DBPORT -D $LC_DBNAME -e 'select sr.name,sr.length from prediction_transcript pt, seq_region sr where sr.seq_region_id=pt.seq_region_id limit 10'");

  print STDERR "\n** Are you ready to run the rulemanager on all analyses? \n yes or no:";  
  if (<STDIN>=~/yes/i){
    print STDERR "Rulemanager will now run Genscan, Uniprot, Unigene and Vertrna... don't forget to monitor your jobs\n".
    print STDERR "perl $LC_cvsDIR/ensembl-pipeline/scripts/rulemanager.pl -dbname $LC_DBNAME -dbhost ".
           "$LC_DBHOST -dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT ".
            "-analysis Genscan -analysis Uniprot -analysis Unigene -analysis Vertrna >& $LC_scratchDIR/".
	    $LC_BUILD_VERSION."_repeatmask_dependents.out\n";
    system("perl $LC_cvsDIR/ensembl-pipeline/scripts/rulemanager.pl -dbname $LC_DBNAME -dbhost ".
           "$LC_DBHOST -dbuser $LC_DBUSER -dbpass $LC_DBPASS -dbport $LC_DBPORT ".
            "-analysis Genscan -analysis Uniprot -analysis Unigene -analysis Vertrna >& $LC_scratchDIR/".
	    $LC_BUILD_VERSION."_repeatmask_dependents.out ");

    print STDERR "\n\nrun_RepeatMask_dependent appears to be finished\nHOWEVER it is a good idea to check the jobs table "
                ."to ensure that all of the jobs have run\nand check $LC_scratchDIR/".$LC_BUILD_VERSION."_repeatmask_dependents.out too\n\n"
                ."If there are still jobs, *wait until those which are running are finished*\nthen empty the job and job_status tables," 
                ." remove the pipeline lock\nand rerun 'run_RepeatMask_dependent'\n"
                ."Once all of the jobs have finished you are ready to run 'dump_RepeatMasked_sequence'\n";
				  
  }else{
  	  print STDERR "RepeatMask_dependent analyses have not been run\n";
  }  
  	
} elsif ($comm eq "dump_RepeatMasked_sequence"){
  print STDERR "\n>> Now going to dump RepeatMasked sequences...\nYou should have already filled in LC_REPMASK_CHOICE variable\n";
  print STDERR "\nYou have the following amount of space on $LC_workDIR :\n";
  system("df -h $LC_workDIR | gawk '{print \$4}' | grep G");
  
  print STDERR "\n** Are you ready to dump the repeat-masked sequence used in the raw compute? yes or no:";
  unless (<STDIN> =~/yes/i){
    exit 1;    
  }
  
  print STDERR ">> Dumping repeat-masked sequence used in the raw compute...\n";
  my $choice = '';	       
  foreach my $r (@$LC_REPMASK_CHOICE){
   print STDERR "$r\n";
   $choice = $choice."_".$r;
  }
  system("mkdir $LC_workDIR/seqdata");
  
  # used bigmem because unsure of how much memory will need
  # Could probably go with 'normal' instead of 'bigmem'
  # MEM: 234 Mbytes;  SWAP: 245 Mbytes;  NTHREAD: 3 (for armadillo)
  $cmd = "bsub -q bigmem -R 'select[mem>3399] rusage[mem=3399]' -o $LC_workDIR/seqdata/all_".
         $LC_DBNAME."".$choice.".out ".
         "$LC_cvsDIR/ensembl-personal/searle/scripts/fetch_slice_seq.pl -host $LC_DBHOST ".
	 "-path $LC_DEFAULT -port $LC_DBPORT -pass $LC_DBPASS -user $LC_DBUSER ".
         "-dbname $LC_DBNAME -softmask";
  
  foreach my $r(@$LC_REPMASK_CHOICE){
    $cmd = $cmd." -softrepeattype ".$r;
  }	 
  $cmd = $cmd." -softrepeattype Dust -all -outfile $LC_workDIR/seqdata/all_".$LC_DBNAME."".$choice.".fa"; 
  print STDERR "\n\ncommand = ".$cmd."\n";
  
  if(system($cmd)){
    warn "\nError submitting job to the farm.\n";
    exit 1;
  } 
  print STDERR "Use bjobs to see when this job is completed.\n";


} elsif ($comm eq "healthcheck") {

  system("mkdir $LC_workDIR/healthchecks");
  print STDERR "\nMake sure your species been included in Species.java and database.properties is pointing to the correct database";
  print STDERR "\nPost_genebuild healthcheck results are written to: $LC_workDIR/healthchecks/firstrun.out\n";

  print STDERR "\nTo run the healthcheck use these commands:\n";
  print STDERR "cd $LC_cvsDIR/ensj-healthcheck\n".
               "./compile-healthcheck.sh\n".
               "./run-healthcheck.sh -output info -d $LC_DBNAME -type core -species $LC_SPECIES ".
               "post_genebuild >>& $LC_workDIR/healthchecks/firstrun.out \n";
  print STDERR "\nyou can use 'backup' to backup your database, you should do this before making any changes.\n";
  
  # Fixes:
  #'select * from repeat_feature where repeat_start < 1'
  #'delete from repeat_feature where repeat_start < 1'
  #'select sf.* from simple_feature sf, seq_region sr where sf.seq_region_end > sr.length and sr.seq_region_id=sf.seq_region_id'
  #'delete simple_feature.* from simple_feature, seq_region where simple_feature.seq_region_end > seq_region.length and seq_region.seq_region_id=simple_feature.seq_region_id'
  
} elsif ($comm eq "backup") {
  print STDERR "\nYou have this amount of space on $LC_workDIR :\n";
  system("df -h $LC_workDIR | gawk '{print \$4}' | grep G ");
  print STDERR "Ready to make the backup? \nyes or no:";
  if (<STDIN> =~ /yes/i){
    if (system("mkdir $LC_workDIR"."/".$LC_DBNAME."_db_backup")){
      print STDERR "\ncannot create $LC_workDIR"."/".$LC_DBNAME."_db_backup\n";
      exit;
    }

    if (system("chmod a+w $LC_workDIR"."/".$LC_DBNAME."_db_backup")){
      print STDERR "\ncannot chmod $LC_workDIR"."/".$LC_DBNAME."_db_backup\n";
      exit;
    }
    $cmd = "mysqldump -h".$LC_DBHOST." -u".$LC_DBUSER." -p".$LC_DBPASS." -P".$LC_DBPORT." ".
           $LC_DBNAME." -T $LC_workDIR"."/".$LC_DBNAME."_db_backup";
    if (system($cmd)){
      print STDERR "Error backing up $LC_DBNAME in $LC_workDIR"."/".$LC_DBNAME."_db_backup\n";
      exit;
    }else{
      print STDERR "$LC_DBNAME is backed up\n".
      "you can tar and zip it to save space\n";
    }
  }else{
    print STDERR "you chose not to make the backup\n";
  }
   
  #Tar and zip afterwards to save lots of space!

  # If you need to recreate the database after mysqldump -T:
  #create a db
  #make the tables either with appropriate ensembl table.sql or by piping in all the .sql from the dump - go to directory of dump and do this:
  #cat *.sql |  mysql [db details]
  #Then import the data using mysqlimport.  --local seems to be required.
  #Like this:
  #mysqlimport  --local -h yourhost -P yourport -u ensadmin -p passwd yourdbname  *.txt
  #This will take quite a long time for large tables 



} elsif ($comm eq "clean"){

  #have kept all of the config files
  print STDERR "\n\n>> Cleaning up...will empty assembly and repeat_libraries and remove ".
               "seqdump directory from $LC_workDIR. \n** Is this ok? yes or no:";
  if (<STDIN> =~ /yes/i){
    #empty the assembly and repeat library directories - 
    if(!checkdir($LC_workDIR."/assembly", 1)){ warn "could not prepare directory! [".$LC_workDIR."/assembly"."]";}
    if(!checkdir($LC_workDIR."/seqdump", 0)){ warn "could not prepare directory! [".$LC_workDIR."/seqdump"."]";}
    if(!checkdir($LC_workDIR."/repeat_libraries", 1)){ warn "could not prepare directory! [".$LC_workDIR."/repeat_libraries"."]";}
    print STDERR "assembly, seqdump and repeat_libraries directories have been emptied/removed\n";
   
  }else{
    print STDERR "assembly, seqdump and repeat_libraries directories have not been emptied\n";
  }
  
  print STDERR "\n\nemptying $LC_scratchDIR/raw_computes\n"; #this takes a while - lots of files
  if(!checkdir($LC_scratchDIR."/raw_computes", 1)){ warn "could not prepare directory! [".$LC_scratchDIR."/raw_computes"."]";}
  
  if(-e $LC_scratchDIR."/".$LC_BUILD_VERSION."_repeatmasking.out"){
    $cmd = "rm " . $LC_scratchDIR."/".$LC_BUILD_VERSION."_repeatmasking.out";
    $status += system($cmd);
  }
  if(-e $LC_scratchDIR."/".$LC_BUILD_VERSION."_repeatmask_dependents.out"){
    $cmd = "rm " . $LC_scratchDIR."/".$LC_BUILD_VERSION."_repeatmask_dependents.out";
    $status += system($cmd);
  }

  if($status){ warn("Error deleting files.\n"); $status = 0; }

}
set_env();

print STDERR "\nDONE\n";
exit 1;

##############################

sub set_env{
  my $set = shift;
  #BlastdDB
  #setenv PATH /usr/local/ensembl/bin:/usr/local/ensembl/mysql/bin:/usr/apps/bin/:/usr/local/bin:/usr/local/pubseq/bin:/bin/:./:/usr/opt/java142/bin:${PATH} 

  if($set){
    #set some environment variables
	if (! ($ENV{"PERL5LIB"} =~ m|.*ensembl-config/$LC_SPECIES/$LC_BUILD_VERSION.*|) ){
      $perllib = $ENV{"PERL5LIB"};
	  $ENV{"PERL5LIB"} = "$LC_PERL5LIB:".$ENV{"PERL5LIB"};
    }
	
	$path = $ENV{"PATH"};
    if (! ($ENV{"PATH"} =~ m|.*/usr/opt/java142/bin.*|) ) {
      $ENV{"PATH"} = "/usr/opt/java142/bin:".$ENV{"PATH"}
    }
  }
  else{
    #re-set them to the initial values
	if($perllib){ $ENV{"PERL5LIB"} = $perllib; }    
    $ENV{"PATH"} = $path;
  }
}

sub createDB{
  my $status = 0;
  print STDERR "Creating new database $LC_DBNAME\n";
  system("mysql -h$LC_DBHOST -P$LC_DBPORT -u$LC_DBUSER -p$LC_DBPASS -e\"DROP database IF EXISTS $LC_DBNAME\"");
  $status += system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -e\"CREATE database $LC_DBNAME\"");
  $status += system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME < $LC_cvsDIR/ensembl-pipeline/sql/table.sql");
  $status += system("mysql -h$LC_DBHOST -u$LC_DBUSER -P$LC_DBPORT -p$LC_DBPASS -D$LC_DBNAME < $LC_cvsDIR/ensembl/sql/table.sql");
  if($status){ die("Error while building the database.") }
}

sub config_setup{

  my $ensconfDIR = $LC_cvsDIR."/ensembl-config/$LC_SPECIES/$LC_BUILD_VERSION";

   system("mkdir -p $ensconfDIR/Bio/EnsEMBL/Pipeline/Config/")
           && warn "\nProblem with making Pipeline in config setup\n";
	
   system("mkdir -p $ensconfDIR/Bio/EnsEMBL/Analysis/Config/")
           && warn "\nProblem with making Analysis in config setup\n";

   system ("mkdir -p $ensconfDIR/pipe_conf")
   	   && warn "\nProblem with making pipe_conf in config setup\n";

   system ("cp $LC_cvsDIR/ensembl-pipeline/scripts/LowCoverage/generic_config/Bio/EnsEMBL/Pipeline/Config/*.pm $ensconfDIR/Bio/EnsEMBL/Pipeline/Config/");
   
   system("cp $LC_cvsDIR/ensembl-pipeline/scripts/LowCoverage/generic_config/Bio/EnsEMBL/Analysis/Config/*.pm $ensconfDIR/Bio/EnsEMBL/Analysis/Config/");
   system ("cp $LC_cvsDIR/ensembl-pipeline/scripts/LowCoverage/generic_config/pipe_conf/*.conf $ensconfDIR/pipe_conf/");

	#have removed the warnings from "cp" because it warns you that it doesn't copy CVS - but that's a good thing 
	
	#need to change path to supplemental and ab-initio data for species in pipe_conf/analysis.conf
	open(IN, "<", "$ensconfDIR/pipe_conf/analysis.conf") or die ("Can't open $ensconfDIR/pipe_conf/analysis.conf");
	open (OUT, ">", "$ensconfDIR/pipe_conf/analysis.temp") or die ("Can't open $ensconfDIR/pipe_conf/analysis.temp");
	local $/ = '[';
	while(<IN>){
		my $entry = $_;
		if ($entry=~/Ab_initio_RepeatMask/){ #name of method
			$entry=~s/parameters=-lib [\/\w\.]+/parameters=-lib $LC_AB_INITIO_LIB/;
		}elsif ($entry=~/Supp_RepeatMask/){ #name of method
			$entry=~s/parameters=-lib [\/\w\.]+/parameters=-lib $LC_SUPP_LIB/;
		}
		print OUT $entry;
	}
	close IN;
	close OUT;
	#replace analysis.conf with analysis.temp
	system(`mv $ensconfDIR/pipe_conf/analysis.temp $ensconfDIR/pipe_conf/analysis.conf`);
	local $/ = "\n";

	
	print STDERR "The $ensconfDIR directories which have just been created have not been added to cvs\n";
}

sub checkdir{
  my $dirname = shift;
  my $option  = shift;
  #go through dirs recursively
  unless (opendir(DIR, $dirname)) {
    closedir(DIR);
    print STDERR "\ncan't open $dirname.\n";
    return 0;
  }
  my @files = grep (!/^\.\.?$/, readdir(DIR));
  foreach my $file (@files) {
    if( -f "$dirname/$file") {
      system("rm $dirname/$file");
    }
    if( -d "$dirname/$file") {
      checkdir("$dirname/$file");
    }
  }
  closedir(DIR);
  if(!$option){
    eval{ system("rmdir $dirname"); };
  }
  return 1;
}

sub find_matching_headers {
  my $num = shift;
  my @matches;
  open(FH, "$LC_workDIR/assembly/assembly.agp.fasta.60");
  my $re = qr/^\>(scaffold\_$num[^0-9]+)\s+.*$/;
  print STDERR "  Searching for $re in $LC_workDIR/assembly/assembly.agp.fasta.60\n";
  while ($_ = <FH>) {
    if ($_ =~ m/^\>(scaffold\_$num[^0-9]+\S+)/){
      push  @matches, $1;
      print STDERR "  Got $1\n";
      #last;
    }
  }
  close (FH);
  return \@matches;
}

__END__


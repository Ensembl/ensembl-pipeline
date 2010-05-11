#! /software/bin/perl -w


# Author: ck1@sanger.ac.uk
# Modified by ml6@sanger.ac.uk

use strict;
use Getopt::Long 'GetOptions';
use DBI;

#ensdb-1-11:5317:vega_homo_sapiens_ext_20080919_v49_NCBI36
#ensdb-1-11:5317:vega_homo_sapiens_ext_20080919_v49_for_ck.

# need to change a line (from throw() to warn() )in sub execute_sanity_check() of
# ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Utils/PipelineSanityChecks.pm
# to get rid of MSG: Some of your analyses don't have entries in the input_id_type_analysis table
# due to rulemanager (as there are other analysis in analysis table)


# typical parameters:
# (1a) preparing pipeline (1st step) : 
#	run_prot_pipeline.pl -dbuser ensadmin -dbhost xx -dbport xx -dbpass xx -dbname xx -prepare1
#       submit the translations dumping script to a LSF big memory machine and wait 
# (1b) preparing pipeline (2nd step) : 
#	run_prot_pipeline.pl -dbuser ensadmin -dbhost xx -dbport xx -dbpass xx -dbname xx -prepare2
#       run once big memory job sumitted in step one has completed
# (2) running pipeline   : 
#	run_prot_pipeline.pl -dbuser ensadmin -dbhost xx -dbport xx -dbpass xx -dbname xx [-start_from 200] [-analysis 202]
#       omitting both -start_from and -analysis options trigers the script to do all analyses
# (3) check duplicate prot. features :
#	run_prot_pipeline.pl -dbuser ensembl -dbhost xx -dbport xx -dbpass xx -dbname xx -finish


my ($dbhost, $dbuser, $dbpass, $dbport, $dbname, $prepare1, $prepare2, $analysis_id, $start_from, $help, $finishing);
my $species = 'human';

GetOptions('dbhost=s'     => \$dbhost,
	   'dbuser=s'     => \$dbuser,
	   'dbpass=s'     => \$dbpass,
	   'dbport=s'     => \$dbport,
	   'dbname=s'     => \$dbname,
       'prepare1'     => \$prepare1,
	   'prepare2'     => \$prepare2,
       'species=s'	  => \$species,
       'start_from=i' => \$start_from,  # dummy analysis_id
       'analysis=i'   => \$analysis_id,
       'finish'       => \$finishing
		);

exec('perldoc', $0) if !($dbhost && $dbuser && $dbpass && $dbport && $dbname);

my $dbh = DBI->connect("DBI:mysql:$dbname:$dbhost:$dbport", $dbuser, $dbpass, {RaiseError => 1, AutoCommit => 0})
        || die "cannot connect to $dbname, $DBI::errstr";

$dbh->debug();

my $perl        = "/software/bin/perl";
my $tmp_xref_db = "ml6_interpro_xref_for_vega";
my $REPOSITORY = "/nfs/anacode/protein_pipeline";
my $PIPE_DIR   = "/software/anacode/pipeline";

my $params1 = "-dbhost $dbhost -dbname $dbname -dbuser $dbuser -dbport $dbport -dbpass $dbpass";
my $params2 = "-user $dbuser -pass $dbpass -host $dbhost -port $dbport -dbname $tmp_xref_db";
my $sql_db  = "mysql --host=$dbhost --user=$dbuser --port=$dbport --pass=$dbpass";

my $ana_dir   = "$PIPE_DIR/ensembl-pipeline/scripts/Finished";
my $ens_dir   = "$PIPE_DIR/ensembl-pipeline/scripts/protein_pipeline";

my $chunk_dir = "$REPOSITORY/chunks";
my $bkp_dir   = "$REPOSITORY/backups";
my $protein   = "$REPOSITORY/proteins/dumped_proteins.fa";
my $log_file  = "$REPOSITORY/log/${dbname}_dump.$$";
my $ext_db_file  = "$PIPE_DIR/ensembl/misc-scripts/external_db/external_dbs.txt";
my $pipe_tables  = "$REPOSITORY/template/prot_pipeline_tables.sql";
my $drop_tables  = "$REPOSITORY/template/drop_pipeline_tables.sql";

my $xref_script        = "$PIPE_DIR/ensembl/misc-scripts/xref_mapping/xref_parser.pl";
my $dump_script        = "$perl $ens_dir/dump_translations.pl";
my $chunk_script       = "$perl $ens_dir/chunk_protein_file.pl";
my $fill_id_script     = "$perl $ana_dir/fillin_ids.pl";
my $rulemanager_script = "$perl $ana_dir/rulemanager.pl";
my $del_dup_script     = "$perl $ana_dir/delete_duplicates_from_protfeat_table.pl";


my $status;


if ( $prepare1 ){

  #-----------------------------------------------------
  # populate xref/interpro table for interpro proteins
  #-----------------------------------------------------

  # first prepare an intermediate interpro/xref database
  $status = system("$xref_script $params2 -source Interpro -species $species -create -download_path $REPOSITORY");
  check_status($status, "(1) Preparing interpro/xref database: $tmp_xref_db");

  #--------------------------------
  #   populates vega xref table
  #--------------------------------

  eval {
    $dbh->do(qq{
              INSERT INTO xref
              SELECT null, x.source_id, x.accession, x.label, x.version, x.description, NULL, NULL
              FROM $tmp_xref_db.xref x, $tmp_xref_db.source s
              WHERE x.source_id=s.source_id
            });

    my $getid = $dbh->prepare("SELECT source_id FROM $tmp_xref_db.source WHERE name = 'Interpro'");
    $getid->execute;
    my $ens_interpro_src_id = $getid->fetchrow;
    my $ext_interpro_src_id = `grep Interpro $ext_db_file | cut -f1`;

    #--------------------------------------------------------------------------------------------------
    #  Change Interpro source id which Ensembl xref system assigns to default in $interpro_src_id file
    #--------------------------------------------------------------------------------------------------
    $dbh->do("UPDATE xref SET external_db_id=$ext_interpro_src_id WHERE external_db_id=$ens_interpro_src_id");

    #--------------------------------
    #  populates vega interpro table
    #--------------------------------
    $dbh->do("INSERT IGNORE INTO interpro SELECT * FROM $tmp_xref_db.interpro");

  };

  $status = $@ ?  'failed' : 'successful'; # ?
  print "(2) Populating xref/interpro tables $status ...\n\n";
  die if $status eq "failed";

  #---------------------------------------
  # Create and populate 7 pipeline tables
  # (1) analysis,
  # (2) job,
  # (3) job_status,
  # (4) rule_goal,
  # (5) rule_conditions,
  # (6) input_id_analysis,
  # (7) input_id_type_analysis
  #---------------------------------------
  $status = system("$sql_db $dbname < $pipe_tables");
  check_status($status, "(3) Populating pipeline tables");

  #-----------------------------
  #  dump protein from database
  #-----------------------------

  # expected number of proteins: "select count(*) from translation;"
  $status = system("rm -f $protein");
  check_status($status, "(4.1) Remove all protein dump");
  
  $status = system("bsub -J dump_protein -e $log_file -sp 100 -q normal -M1000000 -R 'select[mem>1000] rusage[mem=1000]' $dump_script $params1 -db_id 1 -file $protein");
  check_status($status, "(4.2) Submit protein dumping script to the farm (type 'bjobs -J dump_protein' to check its status)");
}
if($prepare2){
  #---------------
  # chunk protein
  #---------------
  # current setup: 100 proteins / chunk
  $status = system("rm -f $chunk_dir/chunk*");
  check_status($status, "(5.1) Removing old chunk files");
  $status = system("$chunk_script $params1");
  check_status($status, "(5.2) Chunking protein");

  #---------------------------------
  # Fill in input_id_analysis table
  # TRANSLATONID => 100 (1 analysis)
  # FILENAME     => 200 (6 analyses)
  # PROTEOME     => 300 (1 analyses)
  #+-------------+-----------------------+
  #| analysis_id | logic_name            |
  #+-------------+-----------------------+
  #|         100 | SubmitTranscript      |
  #|         101 | Pfam                  |
  #|                                     |
  #|         200 | SubmitTranscriptChunk |
  #|         201 | Signalp               |
  #|         202 | Tmhmm                 |
  #|         203 | Prints                |
  #|         204 | Ncoils                |
  #|         205 | pfscan                |
  #|         206 | scanprosite*          | * Removed
  #|                                     |
  #|         300 | SubmitProteome        |
  #|         301 | Seg                   |
  #+-------------+-----------------------+

  # input_id_type TRANSLATIONID
  $status = system("$fill_id_script $params1 -id 100 -translation");
  check_status($status, "(7) TRANSLATIONID filled in");

  # input_id_type FILENAME
  $status = system("$fill_id_script $params1 -id 200 -file $chunk_dir");
  check_status($status, "(8) FILENAME filled in");

  # input_id_type PROTEOME
  my $ins = $dbh->prepare(qq{INSERT INTO input_id_analysis VALUES (?,?,?,now(),?,?,?)});;
  eval { $ins->execute('proteome', 'PROTEOME', 300, '', '', 0); };

  $status = $@ ? 'failed' : 'successful'; # ?
  print "(9) PROTEOME filled in $status ...\n\n";
  die if $status eq "failed";

  $dbh->disconnect;
}

#-----------------------
#  run rule_manager.pl
#-----------------------

if ( $start_from and $analysis_id ){
 $status = system("$rulemanager_script $params1 -verbose -start_from $start_from -analysis $analysis_id");
 check_status($status, "Rulemanager start_from $start_from analysis_id $analysis_id");
}
elsif ( $start_from ){
  $status = system("$rulemanager_script $params1 -verbose -start_from $start_from");
  check_status($status, "Rulemanager start_from $start_from");
}
elsif ( !$start_from and !$analysis_id && !$prepare1 && !$prepare2 && !$finishing ){
  # do all
  foreach my $id ( 100, 200, 300 ) {
    $status = system("$rulemanager_script $params1 -verbose -start_from $id");
    check_status($status, "Rulemanager start_from $id");
  }
}

# check duplicated protein_feature
elsif ( $finishing ){
  $status = system("$del_dup_script $params1 -bkp_dir $bkp_dir");
  check_status($status, "Deleting duplicated protein_features");
  $status = system("$sql_db $dbname < $drop_tables") unless $status != 0;
  check_status($status, "Dropping pipeline tables") unless $status != 0;
}

sub check_status {
  my ($status, $step) = @_;
  ($status == 0) ? (print "$step successful ...\n\n") : (warn "$step failed.\n");
}

__END__


#----------------------
#  running script
#----------------------
step 1a: /software/anacode/pipeline/run_prot_pipeline.pl -dbuser ensadmin -dbhost ensdb-1-11 -dbname vega_homo_sapiens_20080313_original -dbport 5317 -dbpass ensembl -prepare1

step 1b: /software/anacode/pipeline/run_prot_pipeline.pl -dbuser ensadmin -dbhost ensdb-1-11 -dbname vega_homo_sapiens_20080313_original -dbport 5317 -dbpass ensembl -prepare2

step 2: /software/anacode/pipeline/run_prot_pipeline.pl -dbuser ensadmin -dbhost ensdb-1-11 -dbname vega_homo_sapiens_20080313_original -dbport 5317 -dbpass ensembl

step 3: /software/anacode/pipeline/run_prot_pipeline.pl -dbuser ensadmin -dbhost ensdb-1-11 -dbname vega_homo_sapiens_20080313_original -dbport 5317 -dbpass ensembl -finish


#---------------------
# configs, templates
#---------------------
protein pipeline files:	/nfs/anacode/protein_pipeline/
table template: 	/nfs/anacode/protein_pipeline/template/prot_pipeline_tables.sql

# modified version of ProteinAnnotation.pm for Seq
/nfs/anacode/dna_pipeline/config/module/ProteinAnnotation.pm

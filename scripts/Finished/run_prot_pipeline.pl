#!/software/bin/perl -w


# Author: ck1@sanger.ac.uk

use strict;
use Getopt::Long 'GetOptions';
use DBI;

#ensdb-1-11:5317:vega_homo_sapiens_ext_20080919_v49_NCBI36 
#ensdb-1-11:5317:vega_homo_sapiens_ext_20080919_v49_for_ck. 

# need to change a line (from throw() to warn() )in sub execute_sanity_check() of
# /lustre/work1/sanger/ck1/NEW_PIPE/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Utils/PipelineSanityChecks.pm
# to get rid of MSG: Some of your analyses don't have entries in the input_id_type_analysis table
# due to rulemanager (as there are other analysis in analysis table)


# typical parameters:
# (a) preparing pipeline: run_prot_pipeline.pl -dbuser ensadmin -dbhost xx -dbport xx -dbpass xx -dbname xx -prepare
# (b) running pipeline  : run_prot_pipeline.pl -dbuser ensadmin -dbhost xx -dbport xx -dbpass xx -dbname xx [-start_from 200] [-analysis 202]
#                         omitting both -start_from and -analysis options trigers the script to do all analyses

# (c) check duplicate
#     prot. features    : run_prot_pipeline.pl -dbuser ensembl -dbhost xx -dbport xx -dbpass xx -dbname xx -finish


my ($dbhost, $dbuser, $dbpass, $dbport, $dbname, $prepare, $analysis_id, $start_from, $help, $finishing);

GetOptions('dbhost=s'     => \$dbhost,
		   'dbuser=s'     => \$dbuser,
		   'dbpass=s'     => \$dbpass,
		   'dbport=s'     => \$dbport,
		   'dbname=s'     => \$dbname,
           'prepare'      => \$prepare,
           'start_from=i' => \$start_from,  # dummy analysis_id
           'analysis=i'   => \$analysis_id,
           'finish'       => \$finishing
		  );
				
exec('perldoc', $0) if !($dbhost && $dbuser && $dbpass && $dbport && $dbname);

my $dbh = DBI->connect("DBI:mysql:$dbname:$dbhost:$dbport", $dbuser, $dbpass, {RaiseError => 1, AutoCommit => 0})
        || die "cannot connect to $dbname, $DBI::errstr";

$dbh->debug();

my $params       = "-dbhost $dbhost -dbname $dbname -dbuser $dbuser -dbport $dbport -dbpass $dbpass";
my $root         = "/lustre/work1/sanger/ck1";
my $PIPELINE     = "$root/NEW_PIPE";
my $protpipe     = "$PIPELINE/ensembl-pipeline/scripts/protein_pipeline";
my $chunk_dir    = "$root/PIPELINE/chunks";
my $protein      = "$root/PIPELINE/proteins/dumped_proteins.fa";
my $chunk_prot   = "perl $protpipe/chunk_protein_file.pl";
my $dump_prot    = "perl $protpipe/dump_translations.pl";
my $pipe_tables  = "$root/PIPELINE/template/prot_pipeline_tables.sql";
my $fill_in_id   = "$root/PIPELINE/scripts/fillin_ids.pl $params";
my $rulemanager3 = "perl $PIPELINE/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/RuleManager3.pl";
my $rulemanager  = "perl $PIPELINE/ensembl-pipeline/scripts/rulemanager.pl";
my $del_dup      = "perl $root/PIPELINE/scripts/delete_duplicates_from_protfeat_table.pl";
my $sql_db       = "mysql --host=$dbhost --user=$dbuser --port=$dbport --pass=$dbpass";
my $perl         = "/usr/local/ensembl/bin/perl";
my $xref_script  = "/lustre/work1/sanger/ck1/NEW_PIPE/ensembl-head/misc-scripts/xref_mapping/xref_parser.pl";
my $ext_db_file  = "$PIPELINE/ensembl-head/misc-scripts/external_db/external_dbs.txt";
my $status;


if ( $prepare ){

  #-------------------------------------------------------
  #  make sure all config files for pipeline are correct
  #-------------------------------------------------------

  my $ens_anal  = "$PIPELINE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Config";
  my $ens_pipe  = "$PIPELINE/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config";
  my $rdb       = "$PIPELINE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/RunnableDB";
  my $copy;

  # config files to update
  $copy  = system("cp $PIPELINE/configs/ensembl-analysis/Config/*.pm $ens_anal");
  $copy .= system("cp $PIPELINE/configs/ensembl-analysis/Config/GeneBuild/*.pm $ens_anal/GeneBuild/");
  $copy .= system("cp $PIPELINE/configs/ensembl-pipeline/Config/*pm $ens_pipe");
  $copy .= system("cp $PIPELINE/configs/ensembl-pipeline/Config/Protein_Annotation/General.pm $ens_pipe/Protein_Annotation/");
  $copy .= system("cp $PIPELINE/configs/ensembl-analysis/ProteinAnnotation.pm $rdb/"); # this is my modified version for Seq and scanprosite
  check_status($copy, "(1) Checking config files");

  #-----------------------------------------------------
  # populate xref/interpro table for interpro proteins
  #-----------------------------------------------------

  # first prepare an intermediate interpro/xref database
  $status = system("$perl $xref_script -user $dbuser -pass $dbpass -host $dbhost -port $dbport -dbname ck1_interpro_xref_for_vega -source Interpro -species human -create");
  check_status($status, "(2) Preparing interpro/xref database: ck1_interpro_xref_for_vega");

  #--------------------------------
  #   populates vega xref table
  #--------------------------------

  eval {
    $dbh->do(qq{
              INSERT INTO xref
              SELECT null, x.source_id, x.accession, x.label, x.version, x.description, NULL, NULL
              FROM ck1_interpro_xref_for_vega.xref x, ck1_interpro_xref_for_vega.source s
              WHERE x.source_id=s.source_id
            });

    my $getid = $dbh->prepare("SELECT source_id FROM ck1_interpro_xref_for_vega.source WHERE name = 'Interpro'");
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
    $dbh->do("INSERT IGNORE INTO interpro SELECT * FROM ck1_interpro_xref_for_vega.interpro");

  };

  $status = $@ ? 'failed' : 'successful';
  print "(3) Populating xref/interpro tables $status ...\n\n";
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
  check_status($status, "(4) Populating pipeline tables");

  #-----------------------------
  #  dump protein from database
  #-----------------------------

  # expected number of proteins: "select count(*) from translation;"
  $status = system("rm -f $protein");
  check_status($status, "(5.1) Remove all protein dump");
  $status = system("$dump_prot $params -db_id 1 -file $protein");
  check_status($status, "(5.2) Dump protein");

  #---------------
  # chunk protein
  #---------------

  # current setup: 100 proteins / chunk
  $status = system("rm -f $chunk_dir/chunk*");
  check_status($status, "(6.1) Removing old chunk files");
  $status = system("$chunk_prot $params");
  check_status($status, "(6.2) Chunking protein");

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
  #|         206 | scanprosite           |
  #|                                     |
  #|         300 | SubmitProteome        |
  #|         301 | Seg                   |
  #+-------------+-----------------------+

  # input_id_type TRANSLATIONID
  $status = system("$fill_in_id -id 100 -translation");
  check_status($status, "(7) TRANSLATIONID filled in");

  # input_id_type FILENAME
  $status = system("$fill_in_id -id 200 -file $chunk_dir");
  check_status($status, "(8) FILENAME filled in");

  # input_id_type PROTEOME
  my $ins = $dbh->prepare(qq{INSERT INTO input_id_analysis VALUES (?,?,?,now(),?,?,?)});;
  eval { $ins->execute('proteome', 'PROTEOME', 300, '', '', 0); };

  $status = $@ ? 'failed' : 'successful';
  print "(9) PROTEOME filled in $status ...\n\n";
  die if $status eq "failed";

  $dbh->disconnect;
}

#-----------------------
#  run rule_manager.pl
#-----------------------

if ( $start_from and $analysis_id ){
 $status = system("$rulemanager3 $params -verbose -start_from $start_from -analysis $analysis_id");
 check_status($status, "Rulemanager start_from $start_from analysis_id $analysis_id");
}
elsif ( $start_from ){
  $status = system("$rulemanager3 $params -verbose -start_from $start_from");
  check_status($status, "Rulemanager start_from $start_from");
}
elsif ( !$start_from and !$analysis_id && !$prepare && !$finishing ){
  # do all
  foreach my $id ( 100, 200, 300 ) {
    $status = system("$rulemanager3 $params -verbose -start_from $id");
    check_status($status, "Rulemanager start_from $id");
  }
}

# check duplicated protein_feature
elsif ( $finishing ){
  $status = system("$del_dup $params");
  check_status($status, "Deleting duplicated protein_features");
}

sub check_status {
  my ($status, $step) = @_;
  ($status == 0) ? (print "$step successful ...\n\n") : (warn "$step failed.\n");
}

__END__

#-------------------------
# Further documentation:
#-------------------------
(1) http://intwebdev.sanger.ac.uk/Teams/Team71/vertann/docs/running_prot_pipeline_on_new_schema.shtml
(2) /lustre/work1/sanger/ck1/NEW_PIPE/ensembl-doc/pipeline_docs/protein_annotation.txt

/lustre/work1/sanger/ck1/NEW_PIPE/configs/ensembl-analysis/Config/
/lustre/work1/sanger/ck1/NEW_PIPE/ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm

#----------------
# useful queries
#----------------

\. ~ck1/query/analysis             -------- show some columns of analysis table
\. ~/query/show_input_id_analysis  -------- show total jobs to do
\. ~ck1/query/protprog             -------- show finished jobs
\. ~ck1/query/protjob              -------- show comprehensive job status
\. ~ck1/query/failed               -------- show failed jobs
\. ~ck1/query/failed_c             -------- show failed jobs grouped by analysis
\. ~ck1/query/kill                 -------- empty job, job_status tables

# checks all external dbs in xref
select edb.db_name, count(*) from xref x, external_db edb where x.external_db_id = edb.external_db_id group by 
edb.db_name;


#-----------------
#  useful aliases
#-----------------
prot - goes to protein dump dir
chunk - goes to chunk file dir
lsfo  - goes to scratch dir

#----------------
#  setup ENV
#----------------
setup ENV: su ml6; source ~ck1/bin/set_pipecodes_csh newpipe


#----------------------
#  running script
#----------------------
step 1: /lustre/work1/sanger/ck1/PIPELINE/scripts/run_prot_pipeline.pl -dbuser ensadmin -dbhost ensdb-1-11 -dbname vega_homo_sapiens_20080313_original -dbport 5317 -dbpass ensembl -prepare

step 2: /lustre/work1/sanger/ck1/PIPELINE/scripts/run_prot_pipeline.pl -dbuser ensadmin -dbhost ensdb-1-11 -dbname vega_homo_sapiens_20080313_original -dbport 5317 -dbpass ensembl

step 3: /lustre/work1/sanger/ck1/PIPELINE/scripts/run_prot_pipeline.pl -dbuser ensadmin -dbhost ensdb-1-11 -dbname vega_homo_sapiens_20080313_original -dbport 5317 -dbpass ensembl -finish


#---------------------
# configs, templates
#---------------------
configfiles:    /lustre/work1/sanger/ck1/NEW_PIPE/configs/
table template: /lustre/work1/sanger/ck1/PIPELINE/template/prot_pipeline_tables.sql

# modified version of ProteinAnnotation.pm for Seq
/lustre/work1/sanger/ck1/NEW_PIPE/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/RunnableDB/ProteinAnnotation.pm

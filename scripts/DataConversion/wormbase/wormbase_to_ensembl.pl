#!/usr/bin/perl

=pod

=head1 NAME

wormbase_setup

set up script for the import of genome annotation for c.elegans 
from the wormbase database and the running of rawcomputes thereon.

=head1 SYNOPSIS / OPTIONS

Fill in the necessary variables in the script and in the config files.
Copy and adjust rule.conf & analysis.conf from ensembl-config/celegans.
Call the script with one of the following options, one after the other:
   setup    - get the data files, create the database,
              load sequence, parse and store genes & features,
              create input-id, dump translations.
   run      - use the rulemanager to start the rawcomputes
   sqlcheck - run some basic sql commands to check data sanity
   xrefs    - load xrefs for the worm
  (clean    - remove tmp file and clean-up the system.)

=head1 DESCRIPTION

The script sets up an ensembl database for c.elegans with
data from the wormbase in an automated fashion.

 1-fetch wormbase files from ftp sites
 2-create empty db
 3-insert sequences & initial analysis entries
 4-insert mtDNA
 5-insert all other features (genes, etc.)
 6-insert meta data & top-level
 7-prepare raw computes with pipeline (make input ids)
 8-dump translations from wormbase genes
 9-run the rulemanager
 10-check the data with some sql
 11-load the xrefs
 12-clean-up the system

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##  PREPARATION:  ##

#adjust configuration files:
#raw-computes
#ensembl-pipeline/scripts/DataConversion/wormbase/WormBaseConf.pm & WormBase.pm
#ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm, Blast.pm, General.pm
#ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config/Protein_Annotation/General.pm

#necessary fixes
# 1.
# sequence_store.pl--$clone_cs and $chromosome_cs objects: 
#			     need additional parameter RANK (1 for chromosome & 2 for clone)
# WormBase.pm--store_coord_system function: additional Hash key RANK, remove Hash key TOP_LEVEL
# 2.
# defined resource=linux for the protein annotations to avoidsystem-specific differences
# some machines dont have gawk installed: signalp fails!
# Fix: copy signalp-bin and set AWK=awk; point entry in analysis-table to this bin

#configuration variables:
my $release       = "CEL140";
my $workDIR       = "../wormbase140/";
my $cvsDIR        = "../cvs_checkout/";
my $analysisConf  = "./analysis.conf"; #to be created
my $ruleConf      = "./rule.conf";     #to be created

my $WB_DBHOST     = "";
my $WB_DBPORT     = "";
my $WB_DBNAME     = $ENV{'USER'}."_caenorhabditis_elegans_core_".$release;
my $WB_DBUSER     = "";
my $WB_DBPASS     = "";

##########################################

use lib '/nfs/acari/fsk/cvs_checkout/ensembl-pipeline/scripts/DataConversion/wormbase';

use warnings;
use strict;
use WormBaseConf;
use WormBase;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
$| = 1;

my $db;
my $perllib;
my $zoe;
my $blastbd;
my $MtDNA_length = 0;
my $cmd;
my $comm = $ARGV[0];
if(!$comm or ($comm ne "setup" or $comm ne "run" or $comm ne "xref" or $comm ne "sqlcheck")){
  exec('perldoc', $0);
  exit 1;
}

print "\nbuilding the worm...\n";

# 0-set environment
set_env(1);

if($comm eq "setup"){

  print "\npreparing for Wormbase data import\n";

  # 1-fetch wormbase files from ftp site
  getfiles();

  # 2-create empty db
  createDB();
  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $WB_DBHOST,
					   -user   => $WB_DBUSER,
					   -dbname => $WB_DBNAME,
					   -pass   => $WB_DBPASS,
					   -port   => $WB_DBPORT,
					  ) or die "cannot access the database!";

  # 3-create config-files
  # TODO

  # 4-insert  analysis entries
  if(insert_analysis($analysisConf, $ruleConf)){ die "error inserting analysis entries."; }

  # 5-insert mtDNA
  insert_mtDNA($db);

  # 6-store other chromosomes
  my $script_path = $cvsDIR."ensembl-pipeline/scripts/DataConversion/wormbase/sequence_store.pl";
  if( system("perl", $script_path) ){ die "could not store genomic sequence."; }

  # 7-insert meta data & toplevel
  insert_meta($db);

  # 8-insert other features
  # skipping Mt-DNA as no agp file??  $WB_TRNA_LOGIC_NAME all on Mt-DNA!
  my @steps = ($WB_LOGIC_NAME, $WB_OPERON_LOGIC_NAME, $WB_RNAI_LOGIC_NAME, $WB_EXPR_LOGIC_NAME, $WB_TRNA_LOGIC_NAME);
  foreach my $step (@steps){
    run_scripts($step);
  }

  # 9-create dummy entries
  make_input_ids("dna");
  # 9-create dummy entries
  make_input_ids("protein");

  # 10-dump translations for protein analysis'
  dump_translations();

  print "\n finnished first part.\n please check the setup.\n";

}
elsif($comm eq "run"){

  print "\nhanding over to rulemanager.".
        "\nPlease run 'perl set_wormbase.pl sqlcheck' after finnishing!\n";

  # 11-RUN RULEMANAGER
  run();

}
elsif($comm eq "sqlcheck"){

  print "\nexecuting some basic sql health-checks\n".
        "the standard health-checks need to be run independently later!\n";
  # 12-run some basic sql-checks
  sql_checks();

}
elsif($comm eq "xref"){

  # 13-CREATE & LOAD XREFS
  X_mapping();

}
elsif($comm eq "clean"){

  # 14-clean up the system
  print "\nstill has to be done manually! <g>\n\n";

}

# DO CHECKS

# 15-re-set environment
set_env(0);

print "\nDONE\n";
exit 1;

###########################################

sub set_env{
  my $set = shift;

  if($set){
    #set some environment variables
    if(! ($ENV{"PERL5LIB"} =~ m|.*ensembl-pipeline/scripts/DataConversion/wormbase.*|) ){
      $perllib = $ENV{"PERL5LIB"};
      $ENV{"PERL5LIB"} = $ENV{"PERL5LIB"}.":".$cvsDIR."ensembl-pipeline/scripts/DataConversion/wormbase";
    }
    $zoe = $ENV{"ZOE"};
    $ENV{"ZOE"} = "/usr/local/ensembl/Zoe";
    $blastbd = $ENV{BLASTDB};
    $ENV{BLASTDB} = "/data/blastdb";
  }
  else{
    #re-set them to the intital values
    if($perllib){ $ENV{"PERL5LIB"} = $perllib; }
    $ENV{"ZOE"} = $zoe;
    $ENV{BLASTDB} = $blastbd;
  }
}


sub getfiles{
    use LWP::Simple;

    #base FTP location
    my $baseurl = "ftp://ftp.sanger.ac.uk/pub/wormbase/$release/CHROMOSOMES/";
    #base file location
    my $basefile = $workDIR;
    my $file;

    my @cromosomes = qw(I II III IV V X);

    foreach my $chomosome (@cromosomes){
      #get agps
      $file = "CHROMOSOME_".$chomosome.".agp";
      print ">>getting file ".$file."\n";
      _fetchfile($baseurl, $basefile, $file);

      #get gffs
      $file = "CHROMOSOME_".$chomosome.".gff.gz";
      print ">>getting file ".$file."\n";
      _fetchfile($baseurl, $basefile, $file);
      _unzipfile($basefile, $file);
    }

    #get mt-dna
    $file = "CHROMOSOME_MtDNA.gff.gz";
    print ">>getting file ".$file."\n";
    _fetchfile($baseurl, $basefile, $file);
    _unzipfile($basefile, $file);
    #create agp file for Mt-DNA
    $file = "CHROMOSOME_MtDNA.agp";
    open(MTDNA, "<", $basefile.$file) or die "cant create file ".$basefile.$file;
    print MTDNA "MtDNA\t1\t$MtDNA_length\t1\tF\tMtDNA-Clone\t1\t$MtDNA_length\t+\n";
    close MTDNA;

    #get clones
    $baseurl = "ftp://ftp.sanger.ac.uk/pub/wormbase/sequences/ALL_COSMIDS/";
    $file = "allcmid.gz";
    print ">>getting file ".$file."\n";
    _fetchfile($baseurl, $basefile, $file);
    _unzipfile($basefile, $file);


    sub _fetchfile{
	my ($url, $path, $filename) = @_;
	my $status;
	$status = getstore( $url.$filename, $path.$filename );
	unless ( $status == 200) {die("Can not get the file: ".$url.$filename)}
    }

    sub _unzipfile{
	my  ($basefile, $filename) = @_;
	$filename = $basefile.$filename;
	system("gzip", "-df", "$filename");
    }
}


sub createDB{
  my $status = 0;
  print ">>creating new database $WB_DBNAME\n";
  system("mysql -h$WB_DBHOST -P$WB_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"DROP DATABASE IF EXISTS $WB_DBNAME;\"");
  $status += system("mysql -h$WB_DBHOST -P$WB_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"create database $WB_DBNAME;\"");
  $status += system("mysql -h$WB_DBHOST -P$WB_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_DBNAME < ".$cvsDIR."ensembl-pipeline/sql/table.sql");
  $status += system("mysql -h$WB_DBHOST -P$WB_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_DBNAME < ".$cvsDIR."ensembl/sql/table.sql");
  if($status){ die("Error while building the database.") }
}


sub insert_mtDNA{
  my $db = shift;
  my $length = 0;
  my $seqfile = $workDIR."CHROMOSOME_MtDNA.dna";
  print ">>inserting mitochondrial DNA\n";
  #fetching raw sequence string
  my $fh;
  my $newid = 1;
  open(IN, "<$seqfile") or die("could not open file $seqfile.");
  my $filecontent = <IN>;
  $filecontent = uc( join("", <IN>) );
  close(IN);
  $filecontent =~ s/\s//g;
  $MtDNA_length = length($filecontent);

  #get next valid id to use
  my $sql = "SELECT seq_region_id FROM seq_region ORDER BY seq_region_id DESC LIMIT 1;";
  my $sth = $db->dbc->prepare($sql);
  $sth->execute() or die("insert_mtDNA mySQL query problem");
  ($newid) = $sth->fetchrow_array();
  if(!$newid){ $newid = 1; }

  #fetch chr-length
  foreach my $chromosome_info(@{$WB_CHR_INFO}) {
    if( $chromosome_info->{'chr_name'} eq "MtDNA" ){
      $length = $chromosome_info->{'length'};
      if($length and $length != $MtDNA_length){
	print "\nMt-DNA length $length does not coorespond to seqence length $MtDNA_length!\n";
	#TODO: fix in file or create file dynamically in the first place...
      }
    }
  }
  if(!$length){ die("could not get MtDNA length") }

  #save as "Mt-chromosome"
  $sql = 'INSERT INTO seq_region SET seq_region_id='.$newid.
         ', name="MtDNA", coord_system_id=2, length='.$length.';';
  $db->dbc->do($sql) or die("insert_mtDNA mySQL query problem");
  #save as Mt-Sequence fragments
  $sql = 'INSERT INTO seq_region SET seq_region_id='.($newid+1).
         ', name="MtDNA-Clone", coord_system_id=1, length='.$length.';';
  $db->dbc->do($sql) or die("insert_mtDNA mySQL query problem");
  #save as Mt-Sequence fragments
  $sql = 'INSERT INTO assembly SET asm_seq_region_id='.$newid.
         ', cmp_seq_region_id='.($newid+1).', asm_start=1, asm_end='.
         $length.', cmp_start=1, cmp_end='.$length.', ori=1;';
  $db->dbc->do($sql) or die("insert_mtDNA mySQL query problem");
  #save sequence
  $sql = 'INSERT INTO dna SET seq_region_id='.($newid+1).', sequence=\''.$filecontent.'\'';
  $db->dbc->do($sql) or die("insert_mtDNA mySQL query problem");
  #set codon table
  $sql = 'INSERT INTO attrib_type SET attrib_type_id = 1, '.
         'code = "codon_table", '.
	 'name = "Codon Table", '.
	 'description = "Alternate codon table";';
  $db->dbc->do($sql) or die("insert_mtDNA mySQL query problem");
  $sql = 'INSERT INTO seq_region_attrib SET seq_region_id  = ?, '.
         'attrib_type_id = 1, '.
	 'value = "5";';
  $sth = $db->dbc->prepare($sql);
  $sth->execute($newid) or die("insert_mtDNA mySQL query problem");
  $sth->execute($newid+1) or die("insert_mtDNA mySQL query problem");
}


sub insert_meta{
  my $db = shift;
  my $cmd;
  my $dday;

  # adding assembly type to meta table
  my $meta_sql = "insert into meta(meta_key, meta_value) values(?,?)";
  my $meta_sth = $db->prepare($meta_sql);
  my @meta = (
	      'assembly.default', $release,
	      'species.taxonomy_id', '6239',
	      'species.common_name', 'C.elegans',
              'species.classification', 'elegans',
              'species.classification', 'Caenorhabditis',
	      'species.classification', 'Peloderinae',
	      'species.classification', 'Rhabditidae',
	      'species.classification', 'Rhabditoidea',
	      'species.classification', 'Rhabditida',
	      'species.classification', 'Chromadorea',
	      'species.classification', 'Nematoda',
	      'species.classification', 'Metazoa',
	      'species.classification', 'Eukaryota',
	      'genebuild.version', $dday.'Wormbase',
	      'assembly.mapping', 'chromosome:|clone',
	     );
  for(my $i=0; $i<scalar @meta-1; $i+=2){
    my $meta_key   = $meta[$i];
    my $meta_value = $meta[$i+1];
    $meta_sth->execute( $meta_key, $meta_value ) or
      die("insert_meta mySQL query problem");
  }

  #set toplevel
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/set_toplevel.pl ".
         "-dbhost $WB_DBHOST -dbuser $WB_DBUSER -dbpass $WB_DBPASS -dbname $WB_DBNAME";
  if(system($cmd)){ warn "error seting toplevel!"; }
}


sub run_scripts{
  my $analysis = shift;
  my $status = 0;
  print ">>running: $analysis\n";
  my $script_path = $cvsDIR."ensembl-pipeline/scripts/DataConversion/wormbase/feature_store.pl";
  #needs lib $cvsDIR.'ensembl-pipeline/scripts/DataConversion/wormbase' in PERL path.
  my $perl5lib = $ENV{"PERL5LIB"};
  if(! ($perl5lib =~ m|.*ensembl-pipeline/scripts/DataConversion/wormbase.*|) ){
    $ENV{"PERL5LIB"} = $ENV{"PERL5LIB"}.":".$cvsDIR."ensembl-pipeline/scripts/DataConversion/wormbase";
  }
  $status = system("perl", $script_path, $analysis);
  print "could not complete this step! " unless ! $status;
  if($@){ die "\nERROR: $@" }
}


sub check_translation{
  my $db = shift;
  print ">>running: translation check\n";
  $| = 1;
  foreach my $chromosome_info(@{$WB_CHR_INFO}) {
    print "CHR ".$chromosome_info->{chr_name};
    my $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', $chromosome_info->{chr_name});
    print "\nSLICE ".$slice->seq_region_name()."\n";
    my @genes = @{$slice->get_all_Genes};
    open(TRANSLATE, "+>>".$WB_NON_TRANSLATE) or die "couldn't open ".$WB_NON_TRANSLATE." $!";
    foreach my $gene(@genes){
      #print $gene->stable_id.", ";
      my $translation = &translation_check($gene, $db);
      if(!$translation){
	print TRANSLATE $gene->stable_id." from ".$chromosome_info->{'chr_name'}." doesn't translate\n";
	print "\n".$gene->stable_id." from ".$chromosome_info->{'chr_name'}." doesn't translate\n";
      }
    }
    close(TRANSLATE);
  }
}


sub insert_analysis{
  my ($analysis_file, $rule_file) = @_;
  my $cmd;
  my $status = 0;

  #insert analysis entries
  print ">>entering analysis & rule entries.\n";
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/analysis_setup.pl ".
         "-dbhost $WB_DBHOST -dbuser $WB_DBUSER -dbpass $WB_DBPASS -dbname $WB_DBNAME ".
	 "-dbport $WB_DBPORT -read -file ".$analysis_file;
  $status += system($cmd);
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/rule_setup.pl ".
         "-dbhost $WB_DBHOST -dbuser $WB_DBUSER -dbpass  $WB_DBPASS -dbname $WB_DBNAME ".
	 "-dbport $WB_DBPORT -read -file ".$rule_file;
  $status += system($cmd);
  if($status){ warn("Error entering analysis entries.") }
}


sub make_input_ids{
  my $set = shift;
  my $status = 0;

  if($set eq "dna"){
    #dna analysis
    print "\ncreating dummy entries for SubmitClone, SubmitChromosome.\n";
    $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/make_input_ids ".
        "-dbhost $WB_DBHOST -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS ".
	"-dbport $WB_DBPORT -slice -coord_system clone -logic SubmitClone";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/make_input_ids ".
        "-dbhost $WB_DBHOST -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass ".
	"$WB_DBPASS -dbport $WB_DBPORT -slice -coord_system chromosome -logic SubmitChromosome";
    $status += system($cmd);
  }
  elsif($set eq "protein"){
    #protein analysis
    print "\ncreating dummy entries for SubmitProteome, SubmitTranscript and SubmitTranscriptChunk.\n";
    $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/make_input_ids ".
        "-dbhost $WB_DBHOST -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass ".
	"$WB_DBPASS -dbport $WB_DBPORT -single -single_name proteome -logic SubmitProteome";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/make_input_ids ".
        "-dbhost $WB_DBHOST -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass ".
	"$WB_DBPASS -dbport $WB_DBPORT  -translation_id -logic SubmitTranscript";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/make_input_ids ".
        "-dbhost $WB_DBHOST -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass ".
	"$WB_DBPASS -dbport $WB_DBPORT -file -dir /ecs2/scratch1/fsk/wormbase140/peptide_chunks ".
	"-logic SubmitTranscriptChunk";
    $status += system($cmd);
  }
  if($status){ warn("Error entering dummy entries.") }
}


sub run{
  my $cmd= $cvsDIR."/ensembl-pipeline/scripts/rulemanager.pl -dbhost $WB_DBHOST -dbport ".
           "$WB_DBPORT -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS";
  exec($cmd);
}


sub dump_translations{
  my $cmd;
  my $status = 0;

  #write out all translations
  print ">>dumping translations.\n";
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/protein_pipeline/dump_translations.pl ".
         "-dbname $WB_DBNAME -dbport $WB_DBPORT -dbhost $WB_DBHOST -dbuser $WB_DBUSER -dbpass ".
	 "$WB_DBPASS -db_id 1 > ".$workDIR."peptide.fa";
  $status += system($cmd);
  #create small parts of the translation file
  print ">>chopping up translations.\n";
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/protein_pipeline/chunk_protein_file.pl";
  $status += system($cmd);
  if($status){ warn("Error dumping translations.") }
}


sub sql_checks{
  my $sql = "";
  my $res = 0;

  sub query{
    my $sql = shift;
    my $sth;
    $sth = $db->prepare($sql) or die "cant execute sql check.";
    return $sth->execute();
  }

  #1. all genes have transcript
  $sql = "select count(gene.gene_id) from gene left join transcript on gene.gene_id=transcript.gene_id ".
         "where transcript.gene_id is NULL;";
  $res = query($sql);
  if($res){
    print "\nError:$res genes have no transcript!\n";
  }

  #2. all transcripts have exons
  $sql = "select count(transcript.transcript_id) from transcript left join exon_transcript on ".
         "transcript.transcript_id=exon_transcript.transcript_id where exon_transcript.transcript_id is NULL;";
  $res = query($sql);
  if($res){
    print "\nError:$res transcripts have no exons!\n";
  }

  #3. all transcripts belong to a gene
  $sql = "select count(distinct t.transcript_id) from transcript t left join gene g on g.gene_id ".
         "= t.gene_id where g.gene_id is null;";
  $res = query($sql);
  if($res){
    print "\nError:$res transcripts have no gene!\n";
  }

  #4. all exons belong to a transcript
  $sql = "select count(distinct e.exon_id) from exon e left join exon_transcript et on e.exon_id ".
         "= et.exon_id where et.exon_id is null;";
  $res = query($sql);
  if($res){
    print "\nError:$res exons have no transcript!\n";
  }

  #5.check the number of single exon genes
  $sql = "select count(*) as c,transcript_id from exon_transcript group by transcript_id having c=1;";
  $res = query($sql);
  if($res){
    print "\nWarning: $res transcripts have only one exon!\n";
  }

  #6.check that start < end
  $sql = "select count(*) from translation where start_exon_id = end_exon_id and seq_start > seq_end;";
  $res = query($sql);
  if($res){
    print "\nError: $res translations with single exons have start > end!\n";
  }

  #7.starts and ends
  $sql = "SELECT count(*) FROM `simple_feature` where seq_region_start > seq_region_end;";
  $res = query($sql);
  if($res){
    print "\nError:$res simple-feature have start > end!\n";
  }
  $sql = "select count(*) FROM `dna_align_feature` where hit_start > hit_end;";
  $res = query($sql);
  if($res){
    print "\nError:$res dna_align_features have start > end!\n";
  }

  print "\n\n done with tests!\n"

}


sub X_mapping{
  my $XMAP_DIR    = $workDIR."xref_mapping";
  my $XREF_DIR    = $workDIR."xref_refs";
  my $xfile       = $workDIR."celegans.mapping";
  my $XREF_DBNAME = $ENV{'USER'}."_caenorhabditis_elegans_core_xrefs";

  #create mapping input file
  open(XDEF , ">$xfile") or die("cant create mappinf input file");
  print XDEF "xref\nhost=$WB_DBHOST\nport=$WB_DBPORT\ndbname=$XREF_DBNAME\n".
        "user=$WB_DBUSER\npassword=$WB_DBPASS\ndir=$XREF_DIR\n\n".
	"species=caenorhabditis_elegans\n".
	"host=$WB_DBHOST\nport=$WB_DBPORT\ndbname=$WB_DBNAME\n".
        "user=$WB_DBUSER\npassword=$WB_DBPASS\ndir=$XMAP_DIR\n\n";
  close XDEF;
  #do the parsing
  print "\nparsing xrefs. This will take a while.\n";
  $cmd = "perl ".$cvsDIR."/ensembl/misc-scripts/xref_mapping/xref_parser.pm ".
         "-host $WB_DBHOST -port $WB_DBPORT -user $WB_DBUSER -pass $WB_DBPASS ".
         "-dbname $XREF_DBNAME species caenorhabditis_elegans -create";
  if(system($cmd)){
    print "\nERROR PARSING XREFS.\n";
    exit 1;
  }
  #do the mapping, write directly into db
  print "\nmapping xrefs. This will take a while, too.\n";
  $cmd = "perl ".$cvsDIR."/ensembl/misc-scripts/xref_mapping/xref_mapper.pl ".
         "-file $xfile -upload -deleteexisting";
  if(system($cmd)){
    print "\nERROR MAPPING XREFS.\n";
    exit 1;
  }
  print "\ndone with xrefs. you can delete the database ".$ENV{'USER'}.
        "_caenorhabditis_elegans_core_xrefs if you like.\n";

}

__END__


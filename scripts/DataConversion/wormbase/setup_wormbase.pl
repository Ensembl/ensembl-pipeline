#!/usr/bin/perl

###########################################
# set up ensembl database for c.elegans with
# data from the wormbase.
#
# steps:
# 1-fetch wormbase files from ftp site
# 2-create empty db
# 3-insert mtDNA
# 4-insert sequences & other features
# 5-insert meta data & toplevel
# 6-create dummy entries
# 7- RUN RULEMANAGER: DNA-analysis'
# 8-dump translations for protein analysis'
# DO CHECKS
# 9-create dummy entries
# 10- RUN RULEMANAGER: Protein-Analysis
# DO CHECKS
# 11- CREATE & LOAD XREFS
# DO CHECKS
#
###########################################

##  PREPARATION:  ##

#adjust configuration files:
#raw-computes:
#ensembl-pipeline/scripts/DataConversion/wormbase/WormBaseConf.pm & WormBase.pm
#ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config/BatchQueue.pm, Blast.pm, General.pm
#ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config/Protein_Annotation/General.pm
#ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config/analysis.conf, rule.conf (see below)
#ensembl-pipeline/modules/Bio/EnsEMBL/Pipeline/Config/Protein_Annotation/General.pm
#xrefs:
#see below
#other:
#during gene-creation Selenocysteine is added. Check for errors during Protein analysis'.

#configuration variables:
my $workDIR       = "(data / work DIR)";
my $cvsDIR        = "(cvs-base-DIR)";
my $analysisConf  = "(path-to analysis.conf)";
my $ruleConf      = "(path-to rule.conf)";

##########################################

use warnings;
use strict;
use WormBaseConf;
use WormBase;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $db;
my $perllib;
my $zoe;
my $blastbd;

print "\nbuilding the worm...\n";

# 0-set environment
set_env(1);

# 1-fetch wormbase files from ftp site
#getfiles();

# 2-create empty db
#createDB();
$db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $WB_DBHOST,
					 -user   => $WB_DBUSER,
					 -dbname => $WB_DBNAME,
					 -pass   => $WB_DBPASS,
					 -port   => $WB_DBPORT,
					);

# 3-insert mtDNA
#insert_mtDNA($db);

# 3-insert sequences and other features
my @steps = qw(sequence_store.pl gene_store.pl operon_store.pl expr_store.pl rnai_store.pl pseudo_store.pl tRNA_store.pl);
foreach my $step (@steps){
  #run_scripts($step);
}

# 4-insert meta data & toplevel
#insert_meta($db);

# 5-insert analyis and rule definitions
#insert_analysis($analysisConf, $ruleConf);

# 6-create dummy entries
#make_input_ids("dna");

# 7- RUN RULEMANAGER: DNA-Analysis'

# 8-dump translations for protein analysis'
#dump_translations();

# DO CHECKS

# 9-create dummy entries
#make_input_ids("protein");

# 10- RUN RULEMANAGER: Protein-Analysis

# DO CHECKS

# 11-re-set environment
set_env(0);

# 12- CREATE & LOAD XREFS

# DO CHECKS

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
    my $baseurl = "ftp://ftp.sanger.ac.uk/pub/wormbase/live_release/CHROMOSOMES/";
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
	unless ( $status == 200) {die("Can not get the file: ".$url.$path)}
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
  eval{
    system("mysql -h$WB_DBHOST -P$WB_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"drop database $WB_DBNAME;\"");
  };
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

  #get next valid id to use
  my $sql = "SELECT seq_region_id FROM seq_region ORDER BY seq_region_id DESC LIMIT 1;";
  my $sth = $db->prepare($sql);
  $sth->execute() or die("insert_mtDNA mySQL query problem");
  ($newid) = $sth->fetchrow_array();
  if(!$newid){ $newid = 1; }

  #fetch chr-length
  foreach my $chromosome_info(@{$WB_CHR_INFO}) {
    if( $chromosome_info->{'chr_name'} eq "MtDNA" ){
      $length = $chromosome_info->{'length'};
    }
  }
  if(!$length){ die("could not get MtDNA length") }

  #save as "Mt-chromosome"
  $sql = 'INSERT INTO seq_region SET seq_region_id='.$newid.
    ', name="MtDNA", coord_system_id=2, length='.$length.';';
  $db->do($sql) or die("insert_mtDNA mySQL query problem");
  #save as Mt-Sequence fragments
  $sql = 'INSERT INTO seq_region SET seq_region_id='.($newid+1).
    ', name="MtDNA-Clone", coord_system_id=1, length='.$length.';';
  $db->do($sql) or die("insert_mtDNA mySQL query problem");
  #save as Mt-Sequence fragments
  $sql = 'INSERT INTO assembly SET asm_seq_region_id='.$newid.
    ', cmp_seq_region_id='.($newid+1).', asm_start=1, asm_end='.
    $length.', cmp_start=1, cmp_end='.$length.', ori=1;';
  $db->do($sql) or die("insert_mtDNA mySQL query problem");
  #save sequence
  $sql = 'INSERT INTO dna SET seq_region_id='.($newid+1).', sequence=\''.$filecontent.'\'';
  $db->do($sql) or die("insert_mtDNA mySQL query problem");

  #set codon table
  $sql = 'INSERT INTO attrib_type SET attrib_type_id = 1, 
                            code = "codon_table", 
                            name = "Codon Table", 
                            description = "Alternate codon table";';
  $db->do($sql) or die("insert_mtDNA mySQL query problem");

  $sql = 'INSERT INTO seq_region_attrib SET seq_region_id  = ?, 
                                  attrib_type_id = 1, 
                                  value = "5";';
  $sth = $db->prepare($sql);
  $sth->execute($newid) or die("insert_mtDNA mySQL query problem");
  $sth->execute($newid+1) or die("insert_mtDNA mySQL query problem");
}


sub insert_meta{
  my $db = shift;
  my $cmd;
  # adding assembly type to meta table
  my $meta_sql = "insert into meta(meta_key, meta_value) values(?,?)";
  my $meta_sth = $db->prepare($meta_sql);
  my %meta = (
	      '00.assembly.mapping'       => 'Chromosome:|Clone',
	      '01.species.taxonomy_id'    => '6239',
	      '02.species.common_name'    => 'C.elegans',
              '03.species.classification' => 'elegans',
              '04.species.classification' => 'Caenorhabditis',
	      '05.species.classification' => 'Peloderinae',
	      '06.species.classification' => 'Rhabditidae',
	      '07.species.classification' => 'Rhabditoidea',
	      '08.species.classification' => 'Rhabditida',
	      '09.species.classification' => 'Chromadorea',
	      '10.species.classification' => 'Nematoda',
	      '11.species.classification' => 'Metazoa',
	      '12.species.classification' => 'Eukaryota',
	      '13.assembly.default'       => 'CEL130',
	      '14.genebuild.version'      => '0409Wormbase',
	     );
  foreach my $meta_key (keys %meta){
    $meta_sth->execute( substr($meta_key,3,length($meta_key)-3), $meta{$meta_key} ) or
      die("insert_meta mySQL query problem");
  }
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/set_toplevel.pl ".
         "-dbhost $WB_DBHOST -dbuser $WB_DBUSER -dbpass $WB_DBPASS -dbname $WB_DBNAME";
  if(system($cmd)){ warn "error seting toplevel!"; }
}


sub run_scripts{
  my $analysis = shift;
  my $status = 0;
  print ">>running: $analysis\n";
  my $basepath = $cvsDIR."ensembl-pipeline/scripts/DataConversion/wormbase/";
  $analysis = $basepath.$analysis;
  #needs lib $cvsDIR.'ensembl-pipeline/scripts/DataConversion/wormbase' in PERL path.
  my $perl5lib = $ENV{"PERL5LIB"};
  if(! ($perl5lib =~ m|.*ensembl-pipeline/scripts/DataConversion/wormbase.*|) ){
    $ENV{"PERL5LIB"} = $ENV{"PERL5LIB"}.":".$cvsDIR."ensembl-pipeline/scripts/DataConversion/wormbase";
  }
  eval{ $status = system("perl", "$analysis"); };
  warn "could not complete this step! " unless ! $status;
  if($@){ print "\nERROR: $@" }
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
  my $cmd;
  my $status = 0;

  if($set eq "dna"){
    #dna analysis
    print "\ncreating dummy entries for SubmitClone, SubmitChromosome.\n";
    $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/make_input_ids ".
        "-dbhost $WB_DBHOST -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS ".
	"-dbport $WB_DBPORT -slice -coord_system Clone -logic SubmitClone";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/make_input_ids ".
        "-dbhost $WB_DBHOST -dbname $WB_DBNAME -dbuser $WB_DBUSER -dbpass ".
	"$WB_DBPASS -dbport $WB_DBPORT -slice -coord_system Chromosome -logic SubmitChromosome";
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
	"$WB_DBPASS -dbport $WB_DBPORT -file -dir /ecs2/work3/fsk/wormbase/chunks -logic SubmitTranscriptChunk";
    $status += system($cmd);
  }
  if($status){ warn("Error entering dummy entries.") }
}


sub dump_translations{
  my $cmd;
  my $status = 0;

  #write out all translations
  print ">>dumping translations.\n";
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/protein_pipeline/dump_translations.pl ".
         "-dbname $WB_DBNAME -dbhost $WB_DBHOST -dbuser $WB_DBUSER -dbpass $WB_DBPASS -db_id 1 > ".
	 $workDIR."peptide.fa";
  $status += system($cmd);
  #create small parts of the translation file
  print ">>chopping up translations.\n";
  $cmd = "perl ".$cvsDIR."ensembl-pipeline/scripts/protein_pipeline/chunk_protein_file.pl";
  $status += system($cmd);
  if($status){ warn("Error dumping translations.") }
}


sub protein_X_mapping{

  print "\nDO BE DONE MANUALLY\nGET FILES & SET VARIABLES FIRST.";
  #load external DB definitions
  #-->load data infile $cvsDIR.'ensembl/misc-scripts/external_db/external_dbs.txt' into table external_db;
  #get current proteome / xref files
  #-->ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprotfiles_files/proteomes/6239.SPC.gz
  #-->ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/fasta_files/proteomes/6239.FASTAC.gz
  #-->ftp://ftp.sanger.ac.uk/pub/databases/wormpep/wormpep.table
  #-->ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/{I-X}/*.faa
  #-->ftp://ftp.ncbi.nih.gov/refseq/release/invertebrate/{1 - 12}.protein.gpff

  #set all paths in mappping_conf.pl!!!
  #also see: protein_mapping.txt
  #needs files in Bio::EnsEMBL::Pipeline::Config::GeneBuild::{Targetted, GeneBuild, Similarity, GeneBuilder, Scripts}

  #combine wormpep and RefSeq data
  #-->new_prepare_proteome.pl
  #prepare SwissProt data
  #-->perl get_Xmapping.pl
  #run exonerate to apply this data to proteome
  # might need specific version (0.7.1)
  #-->perl ensembl/misc-scripts/protein_match/exonerate_wrapper.pl -d
  #combine Refs and write to db
  #-->perl ensembl/misc-scripts/protein_match/maps2db.pl
  #set display_ids in gene & transcript tables
  #-->perl ensembl/misc-scripts/protein_match/load_transcript_display_id.pl

}


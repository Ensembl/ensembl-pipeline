#!/usr/bin/perl

=pod

=head1 NAME

cDNA_setup

Procedure for getting the latest cDNAs aligned to the human build
in an automated fashion.

=head1 SYNOPSIS / OPTIONS

  A. Fill in config variables, start script with argument 'run':
     "perl cDNA_setup.pl run"
  B. Start script again with argument 'clean' after finnishing the pipeline
     run to clean up: "perl cDNA_setup.pl clean", removing tmp files

=head1 DESCRIPTION

This is a set-up script for generating new cDNA alignments as an isolated step of the
build process with the pipeline. We re using the current build as a basis to add the
latest cDNA / EST information as an additional track. This is done using exonerate in
the pipeline, fasta files and repeat-masked chromosome files.
The configuration variables at the beginning of the script must ALL be filled in.
The whole process usually finishes in ~24h if there are no complications.
The results are dna_align_features in an ESTGENE database.
Use the clean-up option when finnished to leave the system in the original state.

The steps the script performes:
  1. config_setup: check config valiables & files
  2. DB_setup: partly copy current db, insert analysis etc.,
     create TARGET database, synchronise
  3. fastafiles: get & read input files
  4. run_analysis: run exonerate using the pipeline

  5. cleanup: post-process result DB, restore config files, remove tmp files and dbs

What YOU need to do:
  1. Fill in the config variables in this script (just below this).
  2. Run it, check the set-up and re-run if there are errors.
  3. Check the results
  4. Clean up any mess by running 'clean'.
  5. Hand over target-database.

If there is an error and the script dies, the original config files are restored
without removing the data files and databases, allowing the re-run of the script.

The setup of scripts and databases runs for ~ 10 min, the exonerate pipeline needs 
between 5 and 24 h, depending on farm usage.

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# configuration variables, adjust to your needs:

# personal base DIR for ensembl perl libs
# expects to find directories 'ensembl' & 'ensembl-analysis' here
$cvsDIR             = "/nfs/acari/fsk/cvs_checkout";

# personal data dir (for temporaty & result/error files)
$dataDIR            = "/ecs2/work3/fsk/cDNA_update";

# sequence data files, which are used for the update
# if in doubt, ask Hans
$vertrna            = "embl_vertrna-1";
$vertrna_update     = "emnew_vertrna-1";
$refseq             = "hs.fna";
$sourceHost         = "cbi1";
$sourceDIR          = "/data/blastdb";
$masked_genome      = "/data/blastdb/Ensembl/Human/NCBI35/softmasked_dusted";

# db parameters
#admin rights required
$WB_DBUSER          = "";
$WB_DBPASS          = "";
# reference db (current build)
$WB_REF_DBNAME      = "homo_sapiens_core_29_35b";
$WB_REF_DBHOST      = "ecs2";
$WB_REF_DBPORT      = "3364";
# new source db (PIPELINE)
$WB_PIPE_DBNAME     = $ENV{'USER'}."_cDNA_pipe";
$WB_PIPE_DBHOST     = "ecs1a";
$WB_PIPE_DBPORT     = "3306";
# new target db (ESTGENE)
$WB_TARGET_DBNAME   = $ENV{'USER'}."_cDNA_update";
$WB_TARGET_DBHOST   = "ecs2";
$WB_TARGET_DBPORT   = "3362";
# new EST db (might not be really necessary)
$WB_EST_DBNAME      = $ENV{'USER'}."_EST";
$WB_EST_DBHOST      = "ecs2";
$WB_EST_DBPORT      = "3366";


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;

my %saved_files;
my $cmd;
my $status;
$config_file       = "config_files.txt";
$newfile           = "cdna_update";
$configDIR         = $dataDIR."/configbackup";
$chunkDIR          = $dataDIR."/chunks";
$outDIR            = $dataDIR."/output";
my @configvars     = qw(cvsDIR dataDIR chunkDIR outDIR vertrna vertrna_update refseq 
		     configDIR sourceDIR newfile config_file masked_genome WB_DBUSER WB_DBPASS 
		     WB_REF_DBNAME WB_REF_DBHOST WB_REF_DBPORT WB_PIPE_DBNAME WB_PIPE_DBHOST 
		     WB_PIPE_DBPORT WB_TARGET_DBNAME WB_TARGET_DBHOST WB_TARGET_DBPORT 
		     WB_EST_DBNAME WB_EST_DBHOST WB_EST_DBPORT);


my $option = $ARGV[0];
if(!$option or ($option ne "run" and $option ne "clean")){
   exec('perldoc', $0);
   exit 1;
}
if($option eq "run"){
  print "\nstarting cDNA-update for current human build.\n";

  config_setup();

  if(! fastafiles() ){ unclean_exit(); }

  if(! DB_setup()   ){ unclean_exit(); }

  run_analysis();

}
elsif($option eq "clean"){
  print "\ncleanign up after cDNA-update.\n";

  clean_up(0);

}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Write required config files for analysis.
# Stores content and location of original files for later restoration.
# The required config values are written into placeholders in the config skeleton file (config_files.txt)
#  using the variables defined above.

sub config_setup{
  $status = 0;
  my $filecount = 0;
  my ($header, $filename, $path, $tempdir);
  #import function, to be included in all config files
  my $import_sub = '
  sub import {
    my ($callpack) = caller(0);
    my $pack = shift;
    my @vars = @_ ? @_ : keys(%Config);
    return unless @vars;
    eval "package $callpack; use vars qw(".
      join(\' \', map { \'$\'.$_ } @vars) . ")";
    die $@ if $@;
    foreach (@vars) {
      if (defined $Config{ $_ }) {
        no strict \'refs\';
	*{"${callpack}::$_"} = \$Config{ $_ };
      }else {
	die "Error: Config: $_ not known\n";
      }
    }
  }
  1;
  ';

  #check files & directories
  foreach my $configvar (@configvars){
    $configvarref = "$$configvar";
    if(!$configvarref){ die "please define all configuration variables! [$configvarref]\n"; }
    if($configvar =~ m/.+DIR.*/){
      if(!-e $configvarref){
	if(system("mkdir $configvarref")){ die "could not create directory! [$configvarref]\n"; }
      }
      if(!-r $configvarref){ die "directory not writeable! [$configvarref]";}
    }
  }

  #check existence of source databases
  if(!connect_db($WB_REF_DBHOST, $WB_REF_DBPORT, $WB_REF_DBNAME, $WB_DBUSER, $WB_DBPASS)){
    die "could not find $WB_REF_DBNAME.";
  }
  if(!connect_db($WB_PREF_DBHOST, $WB_PREF_DBPORT, $WB_PREF_DBNAME, $WB_DBUSER, $WB_DBPASS)){
    die "could not find $WB_PREF_DNA_DBNAME.";
  }
  if(!connect_db($WB_PREF_DNA_DBHOST, $WB_PREF_DNA_DBPORT, $WB_PREF_DNA_DBNAME, $WB_DBUSER, $WB_DBPASS)){
    die "could not find $WB_PREF_DNA_DBNAME.";
  }
  #go thru config info to create defined files,
  #back-up the original, write version with filled-in variables
  open(RP, "<", "$config_file") or die("can't open config file definitions $config_file");
  local $/ = '>>';
  <RP>;
  while(my $content = <RP>){
    #get specific config-file name
    $content  =~ s/([\w\/_\-\.]+)\n//;
    $header   = $1;
    $header   =~ m/(.+\/)([\w\._\-]*)$/g;
    $path     = $1;
    $filename = $2;
    $content  =~ s/>>//;
    #replace variables in config file
    foreach my $configvar (@configvars){
      $configvarref = "$$configvar";
      $content =~ s/\<$configvar\>/$configvarref/g;
    }
    #backup file if exists
    if(-e $cvsDIR."/".$header){
      $cmd = "mv ".$cvsDIR."/".$header." ".$configDIR."/".$filecount;
      if(system($cmd)){
	die "could not backup config file $header.\n";
      }
    }
    else{
      $path = 0;
    }
    #store file location
    $saved_files{$filecount} = $header;
    #write modified file
    $filename = $cvsDIR."/".$path.$filename;
    open(WP, ">", $filename) or die("can't create new config file.\n");
    print WP $content."\n".$import_sub;
    close WP;
    $filecount++;
  }
  close(RP);

  #save data dump with config paths
  $Data::Dumper::Purity = 1;
  open(WP, "> config_paths.perldata") or die "\ncan't create file for data dumping!\n";
  print WP Data::Dumper->Dump([\%saved_files], ['*saved_files']);
  close(WP);
  $/ = "\n";
  print "created backup of current config files, new config files written.\n";
}


# delete old content if any for a given directory
# (recuresively, as some tmp-dir get pretty crowded...)

sub checkdir{
  my $dirname = shift;
  my $option  = shift;
  #go trough dirs recursively
  unless (opendir(DIR, $dirname)) {
    closedir(DIR);
    print "NOPE";
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


# fetch fasta files, combine them, chop them up

sub fastafiles{
  use Net::SSH qw(sshopen2);

  $status = 0;
  my %EMBL_ids;
  my $header;
  my $vertrna_ver     = 1;
  my $vertrna_upd_ver = 1;
  my $refseq_ver      = 1;
  my $update          = 0;
  my @filestamp;

  eval{
    #check file versions, copy only if changed
    $cmd  = "cd $sourceDIR; ls -n ".$vertrna." ".$vertrna_update." ".$refseq;
    sshopen2("$user\@$host", *READER, *WRITER, "$cmd") || die "ssh: $!";
    while(<READER>){
      @filestamp = split(" ", $_);
      my $stampA = join("-", @filestamp[5..7]);
      $cmd  = "cd ".$dataDIR."; "."ls -n ".$filestamp[8];
      @filestamp = split(" ", `$cmd`);
      my $stampB = join("-", @filestamp[5..7]);
      if($stampA eq $stampB){
	#no changes...
	if($filestamp[8] eq $vertrna){ $vertrna_ver = 0; }
	elsif($filestamp[8] eq $vertrna_update){ $vertrna_upd_ver = 0; }
	elsif($filestamp[8] eq $refseq){ $refseq_ver = 0; }
      }
    }
    close(READER);
    close(WRITER);

    #copy files
    if($vertrna_ver){
      $cmd = "scp -p " . $sourceHost.":".$sourceDIR."/".$vertrna . " " . $dataDIR."/".$vertrna;
      $status += system($cmd);
    }
    if($vertrna_upd_ver){
      $cmd = "scp -p " . $sourceHost.":".$sourceDIR."/".$vertrna_update . " " . $dataDIR."/".$vertrna_update;
      $status += system($cmd);
    }
    if($refseq_ver){
      $cmd = "scp -p " . $sourceHost.":".$sourceDIR."/".$refseq . " " . $dataDIR."/".$refseq;
      $status += system($cmd);
    }
    if($status){ die("Error while copying files.\n") }
    print "copied necessary files.\n";

    if($vertrna_upd_ver or $vertrna_ver or $refseq_ver){
      $update = 1;
      #get human entries, combine base file & update file
      #read update file
      local $/ = "\n>";
      open(RP, "<", $dataDIR."/".$vertrna_update) or die("can t read $vertrna_update\n");
      open(WP, ">", $dataDIR."/".$newfile) or die("can t create $newfile\n");
      <RP>;
      while (my $entry = <RP>){
	if($entry =~ m/Homo sapiens/){
	  #extract & save id
	  $entry =~ s/^[\w\d]+\s([\w\.\d]+)\s.+\n{1}?/$1\n/;
	  if(!$1){ die "\n$vertrna_update: unmatched id pattern:\n$entry\n"; }
	  $EMBL_ids{$1} = 1;
	  #re-write fasta entry
	  $entry =~ s/\>//g;
	  print WP '>'.$entry;
	}
      }
      close(RP);
      print "read update EMBL file.\n";

      #read base file
      open(RP, "<", $dataDIR."/".$vertrna) or die("can t read $vertrna\n");
      <RP>;
      while (my $entry = <RP>){
	if($entry =~ m/Homo sapiens/){
	  #extract & save id
	  $entry =~ s/^[\w\d]+\s([\w\.\d]+).+\n{1}?/$1\n/;
	  if(!$1){ die "\n$vertrna: unmatched id pattern:\n$entry\n"; }
	  if( !defined($EMBL_ids{$1}) ){
	    #add fasta entry for unchanged id
	    $entry =~ s/\>//g;
	    print WP '>'.$entry;
	  }
	}
      }
      close(RP);
      print "read base EMBL file.\n";

      #read RefSeq file
      open(RP, "<", $dataDIR."/".$refseq) or die("can t read $refseq.\n");
      <RP>;
      while (my $entry = <RP>){
	#we're not using 'predicted' XM entries for now
	if($entry =~ m/^gi.+ref\|(NM_.+)\| Homo sapiens.*/){
	  $header = $1;
	}
	elsif($entry =~ m/^gi.+ref\|(NR_.+)\| Homo sapiens.*/){
	  $header = $1;
	}
	else{
	  next;
	}
	$entry =~ s/\>//g;
	if($header){
	  #reduce header to accesion number
	  $entry =~ s/^gi.+\n{1}?/$header\n/g;
	  print WP '>'.$entry;
	}
      }
      print "read RefSeq file.\n";
      close(RP);
      close(WP);
      local $/ = "\n";
    }

    if($update){
      #split fasta files, store into CHUNKDIR
      print "splitting fasta file.\n";
      my $newfasta = $dataDIR."/".$newfile;
      my $outdir   = $chunkDIR;
      my $chunknum = 1000;   #(<300 sequences / file)
      if(system("/nfs/acari/searle/progs/fastasplit/fastasplit $newfasta $chunknum $outdir")){
	die("couldn t split file.$@\n");
      }
      print "chopped up file.\n";
    }
  };
  if($@){
    print STDERR "\nERROR: $@";
    return 0;
  }
  return 1;
}


# prepare required databases: EST_DB, result DB,
# fill required tables with data

sub DB_setup{
  #use strict;
  $status = 0;

  eval{
    #create dbs, deleting if existing
    $status  = system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"DROP DATABASE IF EXISTS $WB_PIPE_DBNAME;\"");
    if($status && $status != 256){ die("couldnt drop old database $WB_PIPE_DBNAME!\n"); }
    $status  = system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"CREATE DATABASE $WB_PIPE_DBNAME;\"");
    $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < ".$cvsDIR."/ensembl/sql/table.sql");
    $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < ".$cvsDIR."/ensembl-pipeline/sql/table.sql");

    $status  = system("mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"DROP DATABASE IF EXISTS $WB_TARGET_DBNAME;\"");
    if($status && $status != 256){ die("couldnt drop old database $WB_TARGET_DBNAME!\n"); }
    $status += system("mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"CREATE DATABASE $WB_TARGET_DBNAME;\"");
    $status += system("mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_TARGET_DBNAME < ".$cvsDIR."/ensembl/sql/table.sql");

    $status  = system("mysql -h$WB_EST_DBHOST -P$WB_EST_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"DROP DATABASE IF EXISTS $WB_EST_DBNAME;\"");
    if($status && $status != 256){ die("couldnt drop old database $WB_EST_DBNAME!\n"); }
    $status += system("mysql -h$WB_EST_DBHOST -P$WB_EST_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"create database $WB_EST_DBNAME;\"");
    $status += system("mysql -h$WB_EST_DBHOST -P$WB_EST_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_EST_DBNAME < ".$cvsDIR."/ensembl/sql/table.sql");

    #copy defined db tables from current build
    $cmd = "mysqldump -u$WB_DBUSER -p$WB_DBPASS -h$WB_REF_DBHOST -P$WB_REF_DBPORT".
      " --add-drop-table $WB_REF_DBNAME".
      " analysis assembly attrib_type coord_system exon exon_stable_id exon_transcript gene gene_stable_id meta meta_coord".
      " seq_region seq_region_attrib transcript transcript_stable_id translation translation_stable_id".
      " > ".$dataDIR."/import_tables.sql";
    $status += system($cmd);
    $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < ".$dataDIR."/import_tables.sql");
    #copy dna table from current build
    $cmd = "mysqldump -u$WB_DBUSER -p$WB_DBPASS -h$WB_REF_DBHOST -P$WB_REF_DBPORT".
           " --add-drop-table $WB_REF_DBNAME dna"." > ".$dataDIR."/import_tables2.sql";
    $status += system($cmd);
    $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < ".$dataDIR."/import_tables2.sql");
    if($status){ die("couldnt create databases!\n"); }
    print "created databases.\n";

    #insert analysis entries
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/add_Analysis ".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
           " -logic_name Exonerate_cDNA_update -program exonerate -program_version 0.8.3".
	   " -program_file /usr/local/ensembl/bin/exonerate-0.8.3 -module Exonerate2Genes".
	   " module_version 1 -gff_source Exonerate -gff_feature similarity -input_id_type FILENAME";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/add_Analysis ".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
           " -logic_name SubmitcDNAChunk -module dummy -input_id_type FILENAME";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/RuleHandler.pl".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
	   " -insert -goal Exonerate_cDNA_update -condition SubmitcDNAChunk";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/make_input_ids".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
	   " -file -dir $chunkDIR -logic_name SubmitcDNAChunk";
    $status += system($cmd);
    if($status){ die("Error while setting up the database.\n") }
    print "database set up.\n";

    #copy analysis entries (and others, just to make sure)
    $cmd = "mysqldump -u$WB_DBUSER -p$WB_DBPASS -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT".
           " --add-drop-table $WB_PIPE_DBNAME".
           " analysis assembly attrib_type coord_system meta meta_coord seq_region seq_region_attrib ".
           "> ".$dataDIR."/import_tables3.sql";
    $status += system($cmd);
    $status += system("mysql -h$WB_EST_DBHOST -P$WB_EST_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_EST_DBNAME < ".$dataDIR."/import_tables3.sql");
    $status += system("mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_TARGET_DBNAME < ".$dataDIR."/import_tables3.sql");
    if($status){ die("Error while synchronising databases.\n") }
    print "databases in sync.\n";
  };
  if($@){
    print STDERR "\nERROR: $@";
    return 0;
  }
  return 1;
}


# call rulemanager to start the exonerate run, leaving the set-up script

sub run_analysis{
  $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/rulemanager.pl ".
         "-dbhost $WB_PIPE_DBHOST -dbport $WB_PIPE_DBPORT -dbuser $WB_DBUSER -dbpass $WB_DBPASS -dbname $WB_PIPE_DBNAME";
  print "\nSTARTING PIPELINE.\nusing the command:\n".$cmd."\n\nPlease monitor results/errors of the pipeline.\n\n";
  exec($cmd);
}


# remove files and database leftovers after analysis,
# restore original config files

sub clean_up{
  my $option = shift;
  $status = 0;
  #read data dump
  open(RP, "< config_paths.perldata") or $status = 1;
  if(!$status){
    undef $/;
    eval <RP>;
    die "\ncant recreate data dump.\n$@\n" if $@;
    close(RP);
    $/ = "\n";
  }
  else{
    warn "\ncan't open data dumping file! Already cleaned?\n";
    $status = 0;
  }

  if(!$option){
    #remove files (fasta, chunks, sql)
    if(-e $dataDIR."/".$vertrna){
      $cmd = "rm " . $dataDIR."/".$vertrna;
      $status += system($cmd);
    }
    if(-e $dataDIR."/".$vertrna_update){
      $cmd = "rm " . $dataDIR."/".$vertrna_update;
      $status += system($cmd);
    }
    if(-e $dataDIR."/".$refseq){
      $cmd = "rm " . $dataDIR."/".$refseq;
      $status += system($cmd);
    }
    if(-e $dataDIR."/".$newfile){
      $cmd = "rm " . $dataDIR."/".$newfile;
      $status += system($cmd);
    }
    if(-e $dataDIR."/import_tables.sql"){
      $cmd = "rm " . $dataDIR."/import_tables.sql";
      $status += system($cmd);
    }
    if(-e $dataDIR."/import_tables2.sql"){
      $cmd = "rm " . $dataDIR."/import_tables2.sql";
      $status += system($cmd);
    }
    if(-e $dataDIR."/import_tables3.sql"){
      $cmd = "rm " . $dataDIR."/import_tables3.sql";
      $status += system($cmd);
    }
    #clean output directories
    if(!checkdir($chunkDIR, 1)){ warn "could not prepare directory! [".$chunkDIR."]";}
    if(!checkdir($outDIR, 1))  { warn "could not prepare directory! [".$outDIR."]";}

    if($status){ warn("Error deleting files.\n"); $status = 0; }

    #remove dbs
    print "\n\nshould we remove the pipeline database? (y/n)   ";
    my $ant = <>;
    if($ant eq "y" or $ant eq "Y" or $ant eq "yes"){
      $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"drop database IF EXISTS $WB_PIPE_DBNAME;\"");
    }
    $status += system("mysql -h$WB_EST_DBHOST -P$WB_EST_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"drop database IF EXISTS $WB_EST_DBNAME;\"");
    if($status){ warn("Error deleting databases.\n"); $status = 0; }

    print "cleaned out databases, removed temporary files.\n";
  }

  if(%saved_files){
    #restore original config files
    foreach my $config_file (keys %saved_files){
      if($saved_files{$config_file}){
	$cmd = "mv ".$dataDIR."/configbackup/".$config_file." ".$cvsDIR."/".$saved_files{$config_file};
      }
      else{
	$cmd = "rm ".$configDIR."/".$config_file;
      }
      $status += system($cmd);
    }
  }
  if($status){ warn("Error restoring config files.\n") }
  print "restored original config files.\n\n";
  if((-e config_paths.perldata) and (system("rm config_paths.perldata"))){
    warn "\ncould not remove perldata file.\n";
  }
}


#partly clean-up after error

sub unclean_exit{
  clean_up(1);
  print "\n Restored original config files.\n Check for errors and restart script.\n";
  exit 1;
}


#connect to a given database, optional with attached DNA db

sub connect_db{

  my $host       = shift;
  my $port       = shift;
  my $dbname     = shift;
  my $user       = shift;
  my $pass       = shift;
  my $dnadb      = shift;
  my $dbObj;
print "$host $port $dbname $user $pass\n";
  if($dnadb){
    $dbObj      = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host   => $host,
						     -port   => $port,
						     -user   => $user,
						     -pass   => $pass,
						     -dbname => $dbname,
						     -dnadb  => $dnadb
						    );
  }
  else{
    $dbObj      = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						     -host    => $host,
						     -port   => $port,
						     -user    => $user,
						     -pass   => $pass,
						     -dbname  => $dbname
						    );
  }
  if(!$dbObj){
    return 0;
  }
  return $dbObj;
}


1;


__END__


The following errors are reported, but can be ignored:

"OSTYPE: Undefined variable"

"You should also add a rule that has SubmitcDNAChunk as its goal,
or this rule will never have its conditions fulfilled."


"-------------------- WARNING ----------------------
MSG: Some of your analyses don t have entries in the input_id_type_analysis table
FILE: Pipeline/Utils/PipelineSanityChecks.pm LINE: 99
CALLED BY: Pipeline/Utils/PipelineSanityChecks.pm  LINE: 73
---------------------------------------------------"

"MSG: Could not find fasta file for 'MT in '/data/blastdb/Ensembl/Human/NCBI35/softmasked_dusted'"

_________________________________________________________


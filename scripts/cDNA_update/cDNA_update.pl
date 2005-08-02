#!/usr/local/ensembl/bin/perl

=pod

=head1 NAME

cDNA_setup

Procedure for getting the latest cDNAs aligned to the human build
in an automated fashion.

=head1 SYNOPSIS / OPTIONS

  A. Fill in config variables, start script with argument 'prepare':
     "perl cDNA_setup.pl prepare".
     After the preperation, the modified genome files will have to be pushed across
     the farm. The previous files should be removed before this! Please mail systems
     about these two things.
  B. Start script again with argument 'run' to start the pipeline.
  B. Check the results by comparing them to the previous alignment
     by calling "perl cDNA_setup.pl compare"
  C. Start script again with argument 'clean' after finnishing the pipeline
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
  4. chop off the polyA tails and chunk the fasta files into smaller pieces
  5. copy & modify the softmasked genome files: thte coordinates of the DR fragments
     need to be corrected according to the assembly table..
  6. run_analysis: run exonerate using the pipeline

  7. comparison: check some numbers by comparing the results to previsous alignments
  8. cleanup: post-process result DB, restore config files, remove tmp files and dbs

What YOU need to do:
  1. Fill in the config variables in this script (just below this).
  2. Check for the existance of two additional programs needed:
	fastasplit, splitting a fasta file into a number of chunks
 	polyA_clipping, removing poly-A tails from sequences.
  3. Ask systems to push the genome files across the farm.
  4. Run it, check the set-up and re-run if there are errors.
  5. Check the results directly and by running 'compare'.
  6. Clean up any mess by running 'clean'.
  7. Hand over target-database (patch to new version if neccessary).

If there is an error and the script dies, the original config files are restored
without removing the data files and databases, allowing the re-run of the script.

The setup of scripts and databases runs for ~ 10 min, the exonerate pipeline needs
around 12 h, depending on farm usage.
Set resource => 'select[mem>2500] rusage[mem=2500]' in BatchQueue.pm and re-run
the pipeline commandif jobs fail or take too long.


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# configuration variables, adjust to your needs:
# all directories without trailing '/'

# personal base DIR for ensembl perl libs
# expects to find directories 'ensembl' & 'ensembl-analysis' here
$cvsDIR               = "";

# personal data dir (for temporaty & result/error files)
$dataDIR              = "";

# sequence data files, which are used for the update
# if in doubt, ask Hans
$vertrna              = "embl_vertrna-1";
$vertrna_update       = "emnew_vertrna-1";
$refseq               = "hs.fna";
$sourceHost           = "cbi1";
$sourceDIR            = "/data/blastdb";
$assembly_version     = "NCBI35";
$org_masked_genome    = ""; #"/data/blastdb/Ensembl/Human/".$assembly_version."/softmasked_dusted";
$target_masked_genome = "/data/blastdb/Ensembl/Human/".$assembly_version."/modified/softmasked_dusted";

# external programs needed (absolute paths):
$fastasplit           = "/nfs/acari/searle/progs/fastasplit/fastasplit";
$polyA_clipping       = "/nfs/acari/fsk/projects/cDNA_update/steve_clip_ployA.pl";

# db parameters
#admin rights required
$WB_DBUSER            = "ensadmin";
$WB_DBPASS            = "ensembl";
# reference db (current build)
$WB_REF_DBNAME        = "homo_sapiens_core_33_35f";
$WB_REF_DBHOST        = "ecs2";
$WB_REF_DBPORT        = "3364";
# new source db (PIPELINE)
$WB_PIPE_DBNAME       = $ENV{'USER'}."_cDNA_pipe";
$WB_PIPE_DBHOST       = "ecs1a";
$WB_PIPE_DBPORT       = "3306";
# new target db (ESTGENE)
$WB_TARGET_DBNAME     = $ENV{'USER'}."_cDNA_update";
$WB_TARGET_DBHOST     = "ia64g";
$WB_TARGET_DBPORT     = "3306";

#use & adjust assembly exception sequences (DR52 & DR53)
$adjust_assembly      = 0;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#no changes should be nesessary below this

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;
use lib '~fsk/perls/';
use Net::SSH qw(sshopen2);

my %saved_files;
my $cmd;
my $status;
#temp. dirs & files:
$config_file       = $cvsDIR."/ensembl-pipeline/scripts/cDNA_update/config_files.txt";
$newfile           = "cdna_update";
$configDIR         = $dataDIR."/configbackup";
$chunkDIR          = $dataDIR."/chunks";
$outDIR            = $dataDIR."/output";
$masked_genome     = $target_masked_genome;
my $oldFeatureName = "Exonerate_cDNA";
my $newFeatureName = "Exonerate_cDNA_update"; #also used as analysis name
my $submitName     = "SubmitcDNAChunk";
my @configvars     = qw(cvsDIR dataDIR chunkDIR outDIR vertrna vertrna_update refseq 
		     configDIR sourceDIR newfile config_file masked_genome fastasplit
                     polyA_clipping WB_DBUSER WB_DBPASS WB_REF_DBNAME WB_REF_DBHOST 
                     WB_REF_DBPORT WB_PIPE_DBNAME WB_PIPE_DBHOST WB_PIPE_DBPORT 
                     WB_TARGET_DBNAME WB_TARGET_DBHOST WB_TARGET_DBPORT);
#fasta chunk specifications:
my $chunknum       = 1000;   #(<300 sequences / file)
my $maxseqlenght   = 20000;
$tmp_masked_genome = "/ecs2/scratch1/fsk/cDNA_update/data/genome";
#program specifications:
my $program_name    = "exonerate";
my $program_version = "0.9.0";
my $program_file    = "/usr/local/ensembl/bin/exonerate-0.9.0";
my $module_name     = "Exonerate2Genes";


my $option = $ARGV[0];
if(!$option or ($option ne "prepare" and $option ne "run" and $option ne "clean" and $option ne "compare")){
   exec('perldoc', $0);
   exit 1;
}
if($option eq "prepare"){
  print "\nstarting cDNA-update for current human build.\n";

  config_setup();

  if(! fastafiles() ){ unclean_exit(); }

  if($adjust_assembly){
    adjust_assembly();
  }

  if(! DB_setup()   ){ unclean_exit(); }

  print "\n\nFinnished setting up the analysis.\n";
  if($adjust_assembly){
    print "The genome files' directory will have to be distributed across the farm!\n".
	  "SOURCE PATH: ".$tmp_masked_genome."\nTARGET PATH: ".$target_masked_genome."\n\n";
  }
}
elsif($option eq "run"){

  run_analysis();

}
elsif($option eq "clean"){
  print "\ncleanign up after cDNA-update.\n";

  clean_up(0);
}
elsif($option eq "compare"){
  print "\nrunning checks after cDNA-update.\n".
        "checking through alignments & genes.\n";

  compare();
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
  #set env var to avoid warnings
  if(!defined $ENV{"OSTYPE"} ){
    $ENV{"OSTYPE"} = "";
  }
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
    #store file location
    $saved_files{$filecount} = $header;
    #write modified file
    $filename = $cvsDIR."/".$path.$filename;
    open(WP, ">", $filename) or die("can't create new config file $filename.\n");
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
    print "\ncan't open $dirname.\n";
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
  $status = 0;
  my %EMBL_ids;
  my $header;
  my $vertrna_ver     = 1;
  my $vertrna_upd_ver = 1;
  my $refseq_ver      = 1;
  my $update          = 1;
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
      #clip ployA tails
      my $newfile2 = $dataDIR."/".$newfile.".clipped";
      $cmd = "$polyA_clipping -mRNA ".
	     $dataDIR."/".$newfile." -out ".$newfile2." -clip";
      if(system($cmd)){
	die("couldn t clip file.$@\n");
      }
      #split fasta files, store into CHUNKDIR
      print "splitting fasta file.\n";
      my $outdir   = $chunkDIR;
      $cmd = "$fastasplit $newfile2 $chunknum $outdir";
      if(system($cmd)){
	die("couldn t split file.$@\n");
      }

      #isolate biggest sequences
      check_chunksizes();

      print "\nchopped up file.\n";
    }
  };
  if($@){
    print STDERR "\nERROR: $@";
    return 0;
  }
  return 1;
}


#adjust sequence files for assembly-exceptions
#currently looks for DR*.fa files

sub adjust_assembly{
  my $filename;
  #move original genome files to defined temporary location
  $cmd = 'ln -s '.$org_masked_genome.'/* '.$tmp_masked_genome.'/';
  if(system($cmd)){
    die("couldn t copy masked genome files.$@\n");
  }
  #get the correct location of DR-assembly pieces
  my $db  = db_connector($WB_PIPE_DBHOST, $WB_PIPE_DBPORT, $WB_PIPE_DBNAME, 'ensro');
  my $sql = 'select s.name, s.seq_region_id, ae.seq_region_start, ae.seq_region_end '.
            'from seq_region s, assembly_exception ae where s.seq_region_id=ae.seq_region_id '.
	    'and s.name like "DR%";';
  my $sth = $db->prepare($sql) or die "sql error!";
  $sth->execute();
  while( my ($name, $seq_region_id, $seq_region_start, $seq_region_end) = $sth->fetchrow_array ){
    #read original file (link)
    $filename = $tmp_masked_genome."/".$name.".fa";
    open(FASTAFILE, "<$filename") or die("cant open fasta file $filename.");
    my $headerline = <FASTAFILE>;
    $headerline =~ s/^\>(\w+)\:([\w\d]*)\:([\w\d]*)\:(\d+)\:(\d+)\:(.+)//;
    $headerline = ">".$1.":".$2.":".$3.":".$seq_region_start.":".$seq_region_end.":".$6;
    local $/ = '';
    my $seq = <FASTAFILE>;
    local $/ = '\n';
    close(FASTAFILE);
    #remove file (link)
    $cmd = 'rm '.$filename;
    if($cmd){ die 'can t remove link '.$filename.'!'; }
    #write modified file
    open(FASTAFILE, ">$filename") or die("cant open fasta file $filename for writing.");
    print FASTAFILE $headerline."\n";
    print FASTAFILE $seq;
    close(FASTAFILE);
  }
  #create README for new genome directory
  my $date = "";
  ($day, $month, $year) = (localtime)[3,4,5];
  my $datestring = printf("%04d %02d %02d", $year+1900, $month+1, $day);
  my $readme = $masked_genome."/README";
  open(README, ">$readme") or die "can t create README file.";
  print README "Directory ".$target_masked_genome."\n\n".
        "These are the softmasked dusted genome files from human ".$assembly_version.
	" assembly with two small modifications:\nThe coordinates of the DR52 & DR53 ".
	"contigs were adjusted according to the assembly table of the database ".
	$WB_REF_DBNAME.
        ".\nThey are used for the cDNA-update procedure to produce an up-to-date cDNA track ".
	"every month.\nCreated by ".$ENV{USER}." on ".$datestring.".\n\n";
  close(README);
}


#find the really big sequences & put them into seperate chunks

sub check_chunksizes{
  local $/ = '>';
  my $allseqs;
  my $file;
  my $toolongs;
  my $seqname;
  my $newfile;

  unless ( opendir( DIR, $chunkDIR ) ) {
    die "can t read $chunkDIR";
  }
  foreach( readdir(DIR) ){
    if(($_ =~ /^\.+$/) or ($_ =~ /^newchunk.+$/)){ next; }
    $file = $chunkDIR."/".$_;
    $toolongs = 0;
    $allseqs = "";

    open(CHUNKFILE, "<$file") or die("can t open file $file.");
    <CHUNKFILE>;
    while(my $seq = <CHUNKFILE>){
      $seq =~ s/\>//;
      $seq =~ m/(.+)\n/;
      $seqname = $1;
      if(length($seq) > $maxseqlenght){
	print "\nTOO LONG: $seqname";
	if(!$toolongs){
	  $toolongs = 1;
	}
	$newfile = $chunkDIR."/newchunk_".$seqname;
	open(NEWFILE, ">$newfile") or die "cant create new fasta file $newfile!";
	print NEWFILE ">".$seq;
	close(NEWFILE);
      }
      else{
	$allseqs .= ">".$seq;
      }
    }
    close(CHUNKFILE);

    if($toolongs){
      open(CHUNKFILE, ">$file") or die("can t open file $file.");
      print CHUNKFILE $allseqs;
      close(CHUNKFILE);
    }
  }
  closedir(DIR);
  local $/ = "\n";
}


# prepare required databases: EST_DB, result DB,
# fill required tables with data

sub DB_setup{
  $status = 0;
  my $program_name    = "exonerate";
  my $program_version = "0.9.0";
  my $program_file    = "/usr/local/ensembl/bin/exonerate-0.9.0";
  my $module_name     = "Exonerate2Genes";

  eval{
    #create dbs, deleting if existing
    $status  = system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"DROP DATABASE IF EXISTS $WB_PIPE_DBNAME;\"");
    if($status && $status != 256){ die("couldnt drop old database $WB_PIPE_DBNAME!\n"); }
    $status  = system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"CREATE DATABASE $WB_PIPE_DBNAME;\"");
    $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < "."/nfs/acari/fsk/cvs_checkout/ensembl/sql/table.sql");
    $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < "."/nfs/acari/fsk/cvs_checkout/ensembl-pipeline/sql/table.sql");
    print ".";
    $status  = system("mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"DROP DATABASE IF EXISTS $WB_TARGET_DBNAME;\"");
    if($status && $status != 256){ die("couldnt drop old database $WB_TARGET_DBNAME!\n"); }
    $status += system("mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e\"CREATE DATABASE $WB_TARGET_DBNAME;\"");
    $status += system("mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_TARGET_DBNAME < "."/nfs/acari/fsk/cvs_checkout/ensembl/sql/table.sql");
    print ".";
    #copy defined db tables from current build
    $cmd = "mysqldump -u$WB_DBUSER -p$WB_DBPASS -h$WB_REF_DBHOST -P$WB_REF_DBPORT -t $WB_REF_DBNAME".
      " analysis assembly attrib_type coord_system exon exon_stable_id exon_transcript gene gene_stable_id meta meta_coord".
      " assembly_exception seq_region seq_region_attrib transcript transcript_stable_id translation translation_stable_id".
      " > ".$dataDIR."/import_tables.sql";
    $status += system($cmd);
    print ".";
    $cmd = "mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e  '".
	   "DELETE FROM analysis; DELETE FROM assembly; DELETE FROM attrib_type; DELETE FROM coord_system;" .
           "DELETE FROM exon; DELETE FROM exon_stable_id; DELETE FROM exon_transcript; DELETE FROM gene; ".
	   "DELETE FROM gene_stable_id; DELETE FROM meta; DELETE FROM meta_coord;  DELETE FROM assembly; ".
           "DELETE FROM assembly_exception; DELETE FROM seq_region_attrib; DELETE FROM transcript; DELETE FROM transcript_stable_id; ".
           "DELETE FROM translation; DELETE FROM translation_stable_id; DELETE FROM assembly_exception;' $WB_PIPE_DBNAME ";
    $status += system($cmd);
    print ".";
    $status += system("mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < ".$dataDIR."/import_tables.sql");
    print ".";
    #copy dna table from current build
    $cmd = "mysqldump -u$WB_DBUSER -p$WB_DBPASS -h$WB_REF_DBHOST -P$WB_REF_DBPORT".
           " -t $WB_REF_DBNAME dna"." > ".$dataDIR."/import_tables2.sql";
    $status += system($cmd);
    print ".";
    $cmd = "mysql -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -u$WB_DBUSER -p$WB_DBPASS $WB_PIPE_DBNAME < ".
            $dataDIR."/import_tables2.sql";
    $status += system($cmd);
    if($status){ die("couldnt create databases!\n"); }
    print "created databases.\n";
    #insert analysis entries
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/add_Analysis ".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
           " -logic_name $newFeatureName -program $program_name -program_version $program_version".
	   " -program_file $program_file -module $module_name".
	   " module_version 1 -gff_source Exonerate -gff_feature similarity -input_id_type FILENAME";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/add_Analysis ".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
           " -logic_name SubmitcDNAChunk -module dummy -input_id_type FILENAME";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/RuleHandler.pl".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
	   " -insert -goal $newFeatureName -condition SubmitcDNAChunk";
    $status += system($cmd);
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/make_input_ids".
           " -dbhost $WB_PIPE_DBHOST -dbname $WB_PIPE_DBNAME -dbuser $WB_DBUSER -dbpass $WB_DBPASS".
           " -file -dir $chunkDIR -logic_name SubmitcDNAChunk";
    $status += system($cmd);
    if($status){ die("Error while setting up the database.\n") }
    print "database set up.\n";
    #copy analysis entries (and others, just to make sure)
    $cmd = "mysqldump -u$WB_DBUSER -p$WB_DBPASS -h$WB_PIPE_DBHOST -P$WB_PIPE_DBPORT -t $WB_PIPE_DBNAME".
           " analysis assembly assembly_exception attrib_type coord_system meta meta_coord seq_region seq_region_attrib ".
           "> ".$dataDIR."/import_tables3.sql";
    $status += system($cmd);
    $cmd = "mysql -h$WB_TARGET_DBHOST -P$WB_TARGET_DBPORT -u$WB_DBUSER -p$WB_DBPASS -e '".
           "DELETE FROM analysis; DELETE FROM assembly; DELETE FROM attrib_type; DELETE FROM coord_system; ".
           "DELETE FROM meta; DELETE FROM meta_coord; DELETE FROM seq_region; DELETE FROM seq_region_attrib; ".
           "DELETE FROM assembly_exception;' $WB_TARGET_DBNAME";
    $status += system($cmd);

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
  #running a test first
  print "\nRunning the test-RunnablDB first.\nPlease monitor the output.\nShould we start? (y/n)";
  my $ant = "";
  chomp($ant = <STDIN>);
  if($ant eq "y" or $ant eq "Y"){
    #get one input id for testing
    my $db = db_connector($WB_PIPE_DBHOST, $WB_PIPE_DBPORT, $WB_PIPE_DBNAME, "ensro");
    my $sql = 'SELECT input_id FROM input_id_analysis i, analysis a WHERE i.analysis_id=a.analysis_id '.
              'AND a.logic_name="'.$submitName.'" LIMIT 1;';
    my $sth = $db->prepare($sql) or die "sql error getting an input-id!";
    $sth->execute();
    my ($input_id) = $sth->fetchrow_array;
    if(!$input_id){
      die "\nCould not get an input id from database!\nQuery used: $sql\n\n";
    }
    $cmd = "perl ".$cvsDIR."/ensembl-analysis/scripts/test_RunnableDB ".
           "-dbhost $WB_PIPE_DBHOST -dbport $WB_PIPE_DBPORT -dbuser $WB_DBUSER -dbpass $WB_DBPASS -dbname $WB_PIPE_DBNAME ".
	   "-input_id $input_id -logic_name $newFeatureName -verbose -nowrite";
    print $cmd."\n";
    system($cmd);
  }
  #start the real process
  print "\n\nShould we start the actual analysis? (y/n)";
  chomp($ant = <STDIN>);
  if($ant eq "y" or $ant eq "Y"){
    $cmd = "perl ".$cvsDIR."/ensembl-pipeline/scripts/rulemanager.pl ".
           "-dbhost $WB_PIPE_DBHOST -dbport $WB_PIPE_DBPORT -dbuser $WB_DBUSER -dbpass $WB_DBPASS -dbname $WB_PIPE_DBNAME";
    print "\nSTARTING PIPELINE.\nusing the command:\n".$cmd."\n\nPlease monitor results/errors of the pipeline.\n\n";
    exec($cmd);
  }
  else{
    print "\nProcess interrupted. Not running pipeline.\n\n";
  }
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
    if($@){
      $/ = "\n";
      die "\ncant recreate data dump.\n$@\n";
    }
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
    print "\n\nshould we remove the clipped fasta file? (y/n)   ";
    my $ant = <>;
    if($ant eq "y" or $ant eq "Y" or $ant eq "yes"){
      if(-e $dataDIR."/".$newfile){
	$cmd = "rm " . $dataDIR."/".$newfile;
	$status += system($cmd);
      }
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


#compare results to previous data as a health-check
#bsubs a function call for every chromosome

sub compare{
  my (%chromosomes_1, %chromosomes_2);
  my ($sql, $sql2, $sth1, $sth2);
  my $exclude_NT = 1;
  my (%hits_per_chrom_1, %hits_per_chrom_2);
  my $hitcount1 = 0;
  my $hitcount2 = 0;
  my (%chromosomes_1, %chromosomes_2);
  my ($sql, $sql2, $sth1, $sth2);
  my $exclude_NT = 1;
  my (%hits_per_chrom_1, %hits_per_chrom_2);
  my $hitcount1 = 0;
  my $hitcount2 = 0;

  #get db connectors
  #old alignments
  my $db1     = db_connector("ecs2", 3365, "homo_sapiens_cdna_30_35c", "ensro");
  #new alignments
  my $db2     = db_connector($WB_TARGET_DBHOST,$WB_TARGET_DBPORT, $WB_TARGET_DBNAME , "ensro");

  #get chromsome names / ids
  $sql = 'select coord_system_id from coord_system where name="chromosome"';
  $sth1 = $db1->prepare($sql) or die "sql error!";
  $sth1->execute();
  my ($coord_system_id) = $sth1->fetchrow_array;
  $sql = 'select seq_region_id, name from seq_region where coord_system_id = '.$coord_system_id;
  if ($exclude_NT) {
    $sql .= ' and name not like "%NT%"';
  }
  $sth1 = $db1->prepare($sql) or die "sql error!";
  $sth1->execute();
  while (my ($seq_region_id, $name) = $sth1->fetchrow_array) {
    $chromosomes_1{$name} = $seq_region_id;
  }
  $sth2 = $db2->prepare($sql) or die "sql error!";
  $sth2->execute();
  while (my ($seq_region_id, $name) = $sth2->fetchrow_array) {
    $chromosomes_2{$name} = $seq_region_id;
  }

#  #create LSF jobs for in-depth analysis
#  foreach my $chomosome (keys %chromosomes_1){
#    $cmd = "bsub -q normal -o ".$dataDIR."/".$chomosome.".out perl ".$cvsDIR.
#           "/ensembl-pipeline/scripts/cDNA_update/comparison.pl ".
#           $chomosome." ".$oldFeatureName." ".$newFeatureName." ".$dataDIR;
#    print $cmd."\n";
#    `$cmd`;
#  }

  print "\nGetting hits per chromosome\n";
  #check hits per chromosome
  $sql = "select count(*) from  dna_align_feature dnaa, analysis a where a.logic_name='human_cDNA_update'".
    "and a.analysis_id=dnaa.analysis_id and dnaa.seq_region_id=?";
  #$sth1 = $db1->prepare($sql) or die "sql error!";
  $sth2 = $db2->prepare($sql) or die "sql error!";
  foreach my $chromosome (keys %chromosomes_1) {
    #$sth1->execute($chromosomes_1{$chromosome});
    $sth2->execute($chromosomes_2{$chromosome});
    #$hits_per_chrom_1{$chromosome} = $sth1->fetchrow_array;
    $hits_per_chrom_2{$chromosome} = $sth2->fetchrow_array;
  }

  my @sorted_chromosomes = sort bychrnum keys %chromosomes_1;
  foreach my $chromosome (@sorted_chromosomes) {
    print "\n$chromosome:".
      "\t".$hits_per_chrom_1{$chromosome}.
      "\t".$hits_per_chrom_2{$chromosome};
    $hitcount1 += $hits_per_chrom_1{$chromosome};
    $hitcount2 += $hits_per_chrom_2{$chromosome};
  }
  print "\ntotal sum:".
      "\t".$hitcount1.
      "\t".$hitcount2;
}


#sort chroms by name.

sub bychrnum {
  my @awords = split /_/,$a;
  my @bwords = split /_/,$b;

  my $anum = $awords[0];
  my $bnum = $bwords[0];

  $anum =~ s/chr//;
  $bnum =~ s/chr//;

  if ($anum !~ /^[0-9]*$/) {
    if ($bnum !~ /^[0-9]*$/) {
      return $anum cmp $bnum;
    } else {
      return 1;
    }
  }
  if ($bnum !~ /^[0-9]*$/) {
    return -1;
  }

  if ($anum <=> $bnum) {
    return $anum <=> $bnum;
  } else {
    if ($#awords == 0) {
      return -1;
    } elsif ($#bwords == 0) {
      return 1;
    } else {
      return $awords[1] cmp $bwords[1];
    }
  }
}


#get a db connection

sub db_connector{
  my $host       = shift;
  my $port       = shift;
  my $dbname     = shift;
  my $user       = shift;
  my $dbCon      = new Bio::EnsEMBL::DBSQL::DBConnection(
							-host    => $host,
							-port    => $port,
							-user    => $user,
							-dbname  => $dbname
						       );
  if(!$dbCon){ die "\ncould not connect to \"$dbname\".\n"; }
  return $dbCon;
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
						     -port    => $port,
						     -user    => $user,
						     -pass    => $pass,
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


---FIX DISPLAY PROBLEM---

[features are not diplayed with exon/iontron structure.]

mysql $READARGS -D$DBNAME -e'select daf.dna_align_feature_id from dna_align_feature daf, transcript_supporting_feature tsf, supporting_feature sf where daf.dna_align_feature_id=tsf.feature_id and daf.dna_align_feature_id!=sf.feature_id;' > cdnaupdate

cat cdnaupdate | awk '{  print "delete from dna_align_feature where dna_align_feature_id="$1";" }' > cdnaupdate.fix

cvs co -r branch-ensembl-32 ensembl
cvs co -r branch-ensembl-32 ensembl-pipeline
cvs co -r branch-ensembl-32 ensembl-analysis

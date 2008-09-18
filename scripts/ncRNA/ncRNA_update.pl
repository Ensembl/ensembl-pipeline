#!/usr/local/ensembl/bin/perl


use ncRNA_update_config; 
use strict;
use vars qw(%Config);
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;
use Bio::SeqIO;

my $usage = "perl ncRNA_update.pl 
-pass *  
-verbose 
-dbsetup (create the dbs) 
-refresh (refresh/create the RFAM/miRBasefiles) 
-species (list of species to run on)
-norfam  (Only run miRNA annotation)
writes the rulemanager command and path to a shell script species.csh 
* = required\n";
my $pass;
my $verbose;
my $config_file = "config_files.txt";
my $dbsetup;
my $refresh;
my $norfam;
my @species_list;
$| = 1;
&GetOptions(
	    'pass=s'       => \$pass,
	    'verbose!'     => \$verbose,
	    'config=s'     => \$config_file,
	    'dbsetup!'     => \$dbsetup,
	    'refresh!'     => \$refresh,
	    'norfam!'      => \$norfam,
            'species=s' => \@species_list,
	   );
die "$usage\n" unless ($config_file && $pass);

if ($refresh){
  print "Updating RFAM and miRNA to latest version\n";
  prepare_RFAM();
} 

print "Looking at the config file\n";
my @speciess;
if (scalar(@species_list)) {
  @species_list = split(/,/,join(',',@species_list));
  foreach my $species (@species_list){
    if ($CONFIG->{$species}){
      push @speciess, $species;
    } else {
    print "Skipping species $species\n";
    }
  }
} else {
  @speciess = (keys %$CONFIG);
}

my @globalvars = ($DATADIR, $CVSDIR, $WRITEUSER, $CONFIG, $BLASTDIR);
check_vars(\@globalvars);
foreach my $species (@speciess){
###########################################################
# gonna need seperate config and checkouts of the code    #
# for each species?....                                   #
# actually the code dont matter but the config will       #
# databases.pm will need to change as will the perl5 path #
# and also the batchqueue                                 #
###########################################################

  system ("mkdir $DATADIR") unless -e $DATADIR;
  system ("mkdir $DATADIR/$species") unless -e "$DATADIR/$species";
  system ("mkdir $DATADIR/$species/Bio") unless -e "$DATADIR/$species/Bio";
  system ("mkdir $DATADIR/$species/Bio/EnsEMBL") unless -e "$DATADIR/$species/Bio/EnsEMBL" ;
  system ("mkdir $DATADIR/$species/Bio/EnsEMBL/Analysis") unless -e "$DATADIR/$species/Bio/EnsEMBL/Analysis";
  system ("mkdir $DATADIR/$species/Bio/EnsEMBL/Pipeline") unless -e "$DATADIR/$species/Bio/EnsEMBL/Pipeline";
  system ("mkdir $DATADIR/$species/Bio/EnsEMBL/Analysis/Config") unless -e "$DATADIR/$species/Bio/EnsEMBL/Analysis/Config";
  system ("mkdir $DATADIR/$species/Bio/EnsEMBL/Analysis/Config/GeneBuild") unless -e "$DATADIR/$species/Bio/EnsEMBL/Analysis/Config/GeneBuild";
  system ("mkdir $DATADIR/$species/Bio/EnsEMBL/Pipeline/Config") unless -e "$DATADIR/$species/Bio/EnsEMBL/Pipeline/Config";
  # set up batchqueue output dirs
  # start the db set up if required
  DB_setup($species,
	   $CONFIG->{$species}->{"WRITEHOST"},
	   $CONFIG->{$species}->{"WRITEPORT"},
	   $WRITEUSER,
	   $pass,
	   $CONFIG->{$species}->{"DBHOST"},
	   $CONFIG->{$species}->{"DBPORT"},
	   $CONFIG->{$species}->{"DBNAME"},
	   $CONFIG->{$species}->{"WRITENAME"}) if $dbsetup;

  print "Checking config\n";
  my @localconfigvars =qw(WRITEHOST WRITEPORT DBNAME DBPORT DBHOST 
			  REFINS WRITEINS WRITELOAD REFLOAD WRITENAME );
  config_setup(\@localconfigvars,$species);
}
# once thats all done sucessfully start the rule managers
my $perlpath = $ENV{"PERL5LIB"};
foreach my $species (@speciess){
  open (SPEC,">$species.csh") or die ("Cannot open file $species.csh");
  print SPEC "#!/bin/csh\n\n";
  $ENV{"PERL5LIB"} = "$DATADIR/$species:$CVSDIR/ensembl-analysis/modules:$CVSDIR/ensembl-analysis/scripts:".
     "$CVSDIR/ensembl-pipeline/scripts:$CVSDIR/ensembl-pipeline/modules:".
      "$CVSDIR/ensembl/scripts:$CVSDIR/ensembl/modules:".
        "$BIOPERL_LIVE_PATH:$BIOPERL_RUN_PATH";
  print SPEC "setenv PERL5LIB ".$ENV{"PERL5LIB"}."\n";
  
  system ("perl $CVSDIR/ensembl-pipeline/scripts/setup_batchqueue_outputdir.pl"); 
  # if all be well, run the rulemanager
  my $cmd_rulemanager = "bsub -o $species.out -q normal perl $CVSDIR/ensembl-pipeline/scripts/rulemanager.pl ".
    "-dbname  $CONFIG->{$species}->{\"WRITENAME\"} ".
      "-dbport $CONFIG->{$species}->{\"WRITEPORT\"} ".
	"-dbhost $CONFIG->{$species}->{\"WRITEHOST\"} ".
	  "-dbuser ensadmin -dbpass $pass -once\n";

  print SPEC "$cmd_rulemanager\n";
  print "Monitor:\n";
  my $cmd = "perl $CVSDIR/ensembl-pipeline/scripts/monitor ".
    "-dbname  $CONFIG->{$species}->{\"WRITENAME\"} ".
      "-dbport $CONFIG->{$species}->{\"WRITEPORT\"} ".
	"-dbhost $CONFIG->{$species}->{\"WRITEHOST\"} ".
	  "-dbuser ensadmin -dbpass $pass -current";
  
  print "$cmd\n";
  close SPEC;
};

# set it back to previous
 $ENV{"PERL5LIB"}= $perlpath;
# next need to make sure paths output directories and config files are set up correctly.
# felixes script does this beautifully, I shall pinch methods from that.
exit;

# Write required config files for analysis.
# Stores content and location of original files for later restoration.
# The required config values are written into placeholders in the config skeleton file (config_files.txt)
#  using the variables defined above.

sub config_setup{
  my ($vars,$species) =@_;
  print "SPECIES = $species\n";
  my $status = 0;
  my $filecount = 0;
  my ($header, $filename, $path, $tempdir);
  my $db;
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


  check_vars($vars,$species);
  #check existence of source databases
  print "$species\n";
  print $CONFIG->{$species}->{"DBHOST"}," ".
          $CONFIG->{$species}->{"DBPORT"}," ".
	    $CONFIG->{$species}->{"DBNAME"}," "."\n";
			
  unless( connect_db(
		 $CONFIG->{$species}->{"DBHOST"},
		 $CONFIG->{$species}->{"DBPORT"},
		 $CONFIG->{$species}->{"DBNAME"},
		 'ensro')){
      die "could not find ".$CONFIG->{$species}->{"DBNAME"}."\n";
  }
  $db = connect_db(
		 $CONFIG->{$species}->{"WRITEHOST"},
		 $CONFIG->{$species}->{"WRITEPORT"},
		 $CONFIG->{$species}->{"WRITENAME"},
		 'ensro');

  unless  ( $db ) {
      die "could not find ".$CONFIG->{$species}->{"WRITENAME"}."\n";
  }
# add global variables
  push @$vars,"DATADIR";
  push @$vars,"species";
  # Numbers of slices vary wildly esp with 2x genomes so I am counting slice number for each db
  # So I can sumbit 50 Blast RFAM and 25 BLast miRNA jobs per species
  my $slices = sql("SELECT COUNT(*) FROM input_id_analysis WHERE analysis_id = 7",$db)->[0];
  my $miRNA_batch = int( $slices / 25 );
  my $RFAM_batch  = int( $slices / 50 );
  # put a hard limit on it - batches can be *too* big
  $RFAM_batch  = 5000 if  $RFAM_batch  > 5000;
  $miRNA_batch = 5000 if  $miRNA_batch > 5000;  
  # don't allow a batch size of 0
  $miRNA_batch = 1 unless $miRNA_batch;
  $RFAM_batch  = 1 unless $RFAM_batch;
  print "Have $slices slices, using a batch size of $RFAM_batch for RFAM and " . 
     " $miRNA_batch for miRNA\n";
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
    foreach my $configvar (@$vars){
      my $configvarref;
      if ($configvar eq "DATADIR"){     
      	$configvarref = $DATADIR;
      } 
      elsif ($configvar eq "species") {
       	$configvarref = $species;
      }
       else {
        $configvarref = $CONFIG->{$species}->{$configvar};
      }
      $content =~ s/\<PASS\>/$pass/g;
      $content =~ s/\<$configvar\>/$configvarref/g;
      $content =~ s/\<MIRNABATCH\>/$miRNA_batch/;
      $content =~ s/\<RFAMBATCH\>/$RFAM_batch/;
    }
    #write modified file
    $filename = "$DATADIR/$species/$path/$filename";
    open(WP, ">", $filename) or die("can't create new config file $filename.\n");
    print WP $content."\n".$import_sub;
    close WP;
    $filecount++;
  }
  close(RP);
  print "new config files written.\n";
}


#check files & directories, create if necessary

sub check_vars{
  my ($vars,$species) = @_;
  foreach my $configvar (@$vars){
    print "CONFIG:\t \"$configvar\"\t" if $verbose;
    my $configvarref;
    $configvarref = $configvar unless $species;
    $configvarref = $CONFIG->{$species}->{$configvar} if $species;
    print "$configvarref\n" if $verbose;
    if(!$configvarref){ die "please define all configuration variables! [$configvar]\n"; }
    if($configvar =~ /\//){
      if(!-e $configvarref){
	print "making directory $configvar\n" if $verbose;
	if(system("mkdir $configvarref")){ die "could not create directory! [$configvarref]\n"; }
      }
      if(!-r $configvarref){ die "directory not accessible! [$configvarref]";}
    }
  }
}

sub DB_setup{
  my ($species,$WRITEHOST,$WRITEPORT,$WRITEUSER,$pass,$REFDBHOST,$REFDBPORT,$REFDBNAME,$WRITENAME) = @_;
  my $status = 0;
  print "Creating and loading database $WRITENAME\@$WRITEHOST:$WRITEPORT\n" if $verbose;
  eval{
    #create dbs, deleting if existing
    # need to check if the db exists first
    if(connect_db(
		 $WRITEHOST,
		 $WRITEPORT,
		 $WRITENAME,
		 'ensro'))
    {
      print "Database $WRITENAME @ $WRITEHOST : $WRITEPORT already exists do you want to delete it? (y/n)\n";
      my $reply = <>;
      chomp;
      if ($reply =~ /^Y/ or $reply =~ /^y/){    
          $status  = system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass -e\"DROP DATABASE IF EXISTS $WRITENAME;\"");
          if($status && $status != 256){ die("couldnt drop old database $WRITENAME!\n"); 
	  }
      } else {
      die "Exiting\n";
	}
    }
    print "Creating database $WRITENAME @ $WRITEHOST : $WRITEPORT \n" if $verbose;
    $status  =  system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass -e\"CREATE DATABASE $WRITENAME;\"");
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < "."$CVSDIR/ensembl/sql/table.sql");
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < "."$CVSDIR/ensembl-pipeline/sql/table.sql");
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < "."$CVSDIR/ensembl-pipeline/sql/flag.sql");
    print ".";;
    #copy defined db tables from current build
    my $cmd = "mysqldump -u$WRITEUSER -p$pass -h$REFDBHOST -P$REFDBPORT -t $REFDBNAME".
      " assembly attrib_type coord_system meta meta_coord".
      " assembly_exception seq_region seq_region_attrib ".
      " > ".$DATADIR."/$species/import_tables.sql";
    $status += system($cmd);
    print ".";
    $cmd = "mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass -e  '".
	   "DELETE FROM analysis; DELETE FROM assembly; DELETE FROM attrib_type; DELETE FROM coord_system;" .
	   "DELETE FROM meta; DELETE FROM meta_coord;  DELETE FROM assembly; ".
           "DELETE FROM assembly_exception; DELETE FROM seq_region_attrib;; ".
           "DELETE FROM assembly_exception;' $WRITENAME ";
    $status += system($cmd);
    print ".";
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME -e\"load data local infile \'$CVSDIR/ensembl/misc-scripts/unmapped_reason/unmapped_reason.txt\' into table unmapped_reason;\"");
    print ".";
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME -e\"load data local infile \'$CVSDIR/ensembl/misc-scripts/external_db/external_dbs.txt\' into table external_db;\"");
    print ".";
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < ".$DATADIR."/$species/import_tables.sql");
    print ".";
     #insert analysis entries
    # dont automatically run rfam jobs uless specifically requested
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT -dbpass $pass".
	" -logic_name ncRNA -program cmsearch -program_file /software/ensembl/bin/cmsearch_0.72 -database Rfam -database_version $RFAMVERSION -database_file $BLASTDIR ".
	  " -module Bio::EnsEMBL::Analysis::RunnableDB::Infernal".
	    " module_version 1 -gff_source ensembl -gff_feature gene -input_id_type FLAG";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -logic_name RfamBlast -program wublastn -program_file wublastn -database Rfam ".
	  " -database_file $BLASTDIR/high_copy.fasta -parameters \'lowcopy => $BLASTDIR/low_copy.fasta\' ".
	    " -module Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam".
	      " module_version 1 -gff_source ensembl -gff_feature gene -input_id_type SLICE";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -logic_name BlastmiRNA -program wublastn -program_file wublastn -database  all_mirnas.fa ".
	  "-database_file $BLASTDIR/all_mirnas.fa ".
	    " -module Bio::EnsEMBL::Analysis::RunnableDB::BlastmiRNA".
	      " module_version 1 -gff_source ensembl -gff_feature gene -input_id_type SLICE";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -logic_name DummyFlag -module Dummy -input_id_type  FLAG";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -logic_name miRNA -database miRBase -database_file $BLASTDIR/all_mirnas.embl  ".
	  " -database_version $MIRBASEVERSION -module Bio::EnsEMBL::Analysis::RunnableDB::miRNA -input_id_type GENOME";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER   -dbport $WRITEPORT -dbpass $pass".
	" -logic_name SubmitmiRNA -module Dummy -input_id_type GENOME";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER   -dbport $WRITEPORT -dbpass $pass".
	" -logic_name DummySlice -module Dummy -input_id_type SLICE";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER   -dbport $WRITEPORT -dbpass $pass".
	" -logic_name DummyFlag -module Dummy -input_id_type FLAG";
    # rules 
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/RuleHandler.pl ".
      "-dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER   -dbport $WRITEPORT -dbpass $pass".
	" -insert -goal BlastmiRNA -condition DummySlice";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/RuleHandler.pl ".
      "-dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER   -dbport $WRITEPORT -dbpass $pass".
	" -insert -goal RfamBlast -condition DummySlice";
    $status += system($cmd)    unless ($norfam);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/RuleHandler.pl ".
      "-dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER   -dbport $WRITEPORT -dbpass $pass".
	" -insert -goal ncRNA -condition DummyFlag";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/make_input_ids ".
      "-dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -logic_name DummySlice -slice -slice_size 200000 -coord_system toplevel";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/RuleHandler.pl ".
      "-dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -insert -goal miRNA -condition SubmitmiRNA ";
    $status += system($cmd);
    print "database set up.\n";
  };
  if($@){
    print STDERR "\nERROR: $@";
    return 0;
  }
  return 1;
}
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

  eval {
  	$dbObj->dbc->connect;
       };
  if ($@){
	  return undef;
	 }
  return $dbObj;
}

sub prepare_RFAM{
  # high copy domains
  my %high_copy = 
    (RF00001 => 1,
     RF00003 => 1,
     RF00004 => 1,
     RF00015 => 1,
     RF00017 => 1,
     RF00019 => 1,
     RF00020 => 1,
     RF00023 => 1,
     RF00024 => 1,
     RF00026 => 1,
     RF00045 => 1,
     RF00066 => 1,
     RF00100 => 1,
     RF00108 => 1,
     RF00177 => 1,
    );
  my %families;
  
  # create Rfam.seed file
  my $exit;
  print "Updating RFAM descriptions file ...\n";
  system ("mkdir $BLASTDIR") unless -e "$BLASTDIR";
  $exit =  system ("wget ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz  -O $BLASTDIR/Rfam.tar.gz");
  die ("Error with obtaining Rfam covariance model  file from ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.tar.gz\n") if $exit > 0;
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.tar.gz");
  die ("Error decompressing Rfam.tar.gz\n") if $exit > 0;
  $exit =  system ("tar -xf $BLASTDIR/Rfam.tar -C $BLASTDIR");
  die ("Error extracting Rfam covariance models  file from $BLASTDIR/Rfam.tar\n") if $exit > 0;
  $exit =  system ("wget ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz  -O $BLASTDIR/Rfam.seed.gz");
  die ("Error with obtaining Rfam.seed file from ftp://ftp.sanger.ac.uk/pub/databases/Rfam/Rfam.seed.gz\n") if $exit > 0;
  $exit =   system ("wget ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.fasta.gz -O $BLASTDIR/Rfam.fasta.gz");
  die ("Error with obtaining Rfam.fasta file from ftp://ftp.sanger.ac.uk/pub/databases/Rfam/Rfam.fasta.gz\n") if $exit > 0;
  $exit =   system ("wget ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.thr.gz -O $BLASTDIR/Rfam.thr.gz");
  die ("Error with obtaining Rfam.thr file from ftp://ftp.sanger.ac.uk/pub/databases/Rfam/Rfam.thr.gz\n") if $exit > 0;
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.seed.gz");
  die ("Error decompressing Rfam.seed.gz\n") if $exit > 0;
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.fasta.gz");
  die ("Error decompressing Rfam.fasta.gz\n") if $exit > 0;   
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.thr.gz");
  die ("Error decompressing Rfam.thr.gz\n") if $exit > 0; 
  print "done\nUpdating miRNA file...";
  $exit =   system ("wget ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/miRNA.dat.gz -O $BLASTDIR/all_mirnas.embl.gz");
  $exit =   system ("gzip -d  $BLASTDIR/all_mirnas.embl.gz");
  die ("Error with obtaining miRNA.dat  file from ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/miRNA.dat\n") if $exit > 0;
  $exit =   system ("wget ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/hairpin.fa.gz  -O $BLASTDIR/all_mirnas.fa.gz");
  $exit =   system ("gzip -d  $BLASTDIR/all_mirnas.fa.gz");
  die ("Error with obtaining hairpin.fa  file from ftp://ftp.sanger.ac.uk/pub/mirbase/sequences/CURRENT/hairpin.fa\n") if $exit > 0;
  print "done\n";
  # use bioperl to parse it
    my $miRNA = Bio::SeqIO-> new
    (
     -file => "$BLASTDIR/all_mirnas.embl",
     -format => "embl",
    );
  my $miRNA_fasta = Bio::SeqIO-> new
    (
     -file => ">$BLASTDIR/all_mirnas.fa",
     -format => "Fasta",
    );
  # write to fasta file 
  while (my $seq = $miRNA->next_seq){
    $seq->desc($seq->accession." ".$seq->desc);
    $miRNA_fasta->write_seq($seq);
  }

  print "Fetching all Rfam fasta sequences...";
  my $RFAM = Bio::SeqIO-> new
    (
     -file => "$BLASTDIR/Rfam.fasta",
     -format => "Fasta",
    );
  
  my $high_copy = Bio::SeqIO-> new
    (
     -file => ">$BLASTDIR/high_copy.fasta",
     -format => "Fasta",
    );
  my $low_copy = Bio::SeqIO-> new
    (
     -file => ">$BLASTDIR/low_copy.fasta",
     -format => "Fasta",
    );
  
  die ("Cannot open  /data/blastdb/Rfam/Rfam.fasta\n") unless $RFAM;
  die ("Cannot open  $BLASTDIR/high_copy.fasta\n") unless $high_copy;
  die ("Cannot open  $BLASTDIR/low_copy.fasta\n") unless $low_copy;
  print "Done\n";
  open (DESC,"$BLASTDIR/Rfam.seed") or die "Cannot open description file $BLASTDIR/Rfam.seed\n";
  
  # determine if the ncRNA is a gene or cis acting etc
  # if it is a gene we want it, if not we dont!
  my $domain;
  while (<DESC>){
    chomp;
    $domain = $1 if ($_ =~ /^\#=GF AC   (RF.+)/);
    if ($_ =~ /^\#=GF TP   Gene;/){
      next if $_ =~ /miRNA/;
      next if $_ =~ /tRNA/;
      next if $_ =~ /antisense/;
      $families{$domain} = 1;
    }
  }
  close DESC;
  
  
  while (my $seq = $RFAM->next_seq){
    my $domain = $1 if $seq->desc =~ /^(RF\d+);.+/;
    if ($families{$domain}){
      if ($high_copy{$domain}){
	$high_copy->write_seq($seq);
      }
      else {
	$low_copy->write_seq($seq);
      }
    }
  }
  print "Done\n Formatting blast databases....\n";
  
  system ("xdformat -n $BLASTDIR/all_mirnas.fa");
  system ("xdformat -n $BLASTDIR/high_copy.fasta");
  system ("formatdb -i $BLASTDIR/low_copy.fasta -p f ");

return 0;
}

sub sql {
  my ($query,$db) = @_;
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  my @array = $sth->fetchrow_array;
  return \@array;
}
__END__

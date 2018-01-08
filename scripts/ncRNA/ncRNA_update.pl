#!/usr/bin/env perl


# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


use warnings ;
use ncRNA_update_config; 
use strict;
use vars qw(%Config);
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;
use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Tools::ConfigWriter;

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
GetOptions(
	    'pass|dbpass|p=s'       => \$pass,
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

my @globalvars = ($DATADIR, $CVSDIR, $WRITEUSER, $ROUSER, $CONFIG, $BLASTDIR);
check_vars(\@globalvars);
foreach my $species (@speciess) {
    print "SPECIES = $species\n";
###########################################################
# gonna need seperate config and checkouts of the code    #
# for each species?....                                   #
# actually the code dont matter but the config will       #
# databases.pm will need to change as will the perl5 path #
# and also the batchqueue                                 #
###########################################################

  config_setup($species);
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

# once thats all done sucessfully start the rule managers
  open (SPEC,">$species.csh") or die ("Cannot open file $species.csh");
  print SPEC "#!/bin/csh\n\n";
  my $perl5lib = "$DATADIR/$species:$CVSDIR/ensembl-analysis/modules".
     ":$CVSDIR/ensembl-pipeline/modules:$CVSDIR/ensembl/modules".
     ":$BIOPERL_PATH";
  print SPEC "setenv BLASTMAT /software/ensembl/genebuild/usrlocalensembldata/blastmat\n";
  print SPEC 'setenv PERL5LIB '.$perl5lib."\n";
  
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
	  "-dbuser ensadmin -dbpass $pass -current -finishedpercent";
  
  print "$cmd\n";
  close SPEC;
}

# next need to make sure paths output directories and config files are set up correctly.
# felixes script does this beautifully, I shall pinch methods from that.
exit;

# Write required config files for analysis.

sub config_setup {
    my $species = shift;

    my @analysis_list = ('RfamBlast', 'BlastmiRNA', 'miRNA', 'BlastWait', 'ncRNA');
    my @repeat_masking = ('repeatmask', 'dust');
    my $output_dir = $DATADIR.'/'.$species;
    #set env var to avoid warnings
    if(!defined $ENV{"OSTYPE"} ){
      $ENV{"OSTYPE"} = "";
    }

    # We can simply do a cp of the example file here
    my $pipe_general_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new( -modulename => 'Bio::EnsEMBL::Pipeline::Config::General',
        -is_example => 1);
    $pipe_general_cfg->moduledir($output_dir);
    $pipe_general_cfg->create_path;
    $pipe_general_cfg->write_config(1);
    my $bq_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new( -modulename => 'Bio::EnsEMBL::Pipeline::Config::BatchQueue',
        -is_example => 1,);
    $bq_cfg->moduledir($output_dir);
    $bq_cfg->empty_queue_config(\@analysis_list);
    my $token;
    $token = $CONFIG->{$species}->{DBTOKEN}.'=20:duration=10' if ($CONFIG->{$species}->{DBTOKEN});
    $token = $CONFIG->{$species}->{WRITETOKEN}.'=20' if ($CONFIG->{$species}->{WRITETOKEN});
    if ($token) {
        foreach my $analysis (@analysis_list) {
            my $value = $bq_cfg->analysis_from_batchqueue($analysis);
            if (exists $value->{resource} and $value->{resource} =~ /rusage/) {
                $value->{resource} =~ s/(rusage\[)/$1$token, /;
            }
            else {
                $value->{resource} .= " rusage[$token]";
            }
            $bq_cfg->analysis_from_batchqueue($analysis, $value);
        }
    }
    $bq_cfg->root_value('DEFAULT_OUTPUT_DIR', $output_dir);
    $bq_cfg->write_config(1);
    # We can simply do a cp of the example file here
    my $analysis_general_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new( -modulename => 'Bio::EnsEMBL::Analysis::Config::General',
        -is_example => 1);
    # RepeatMask and Dust were in the config file template so I want to be sure that we still have them
    $analysis_general_cfg->root_value('ANALYSIS_REPEAT_MASKING', \@repeat_masking);
    $analysis_general_cfg->moduledir($output_dir);
    $analysis_general_cfg->create_path;
    $analysis_general_cfg->write_config(1);
    my $database_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new( -modulename => 'Bio::EnsEMBL::Analysis::Config::Databases',
        -is_example => 1);
    $database_cfg->moduledir($output_dir);
    $database_cfg->delete_databases(['REFERENCE_DB', 'GENEBUILD_DB']);
    $database_cfg->key_by_parent('REFERENCE_DB', '-dbname', $CONFIG->{$species}->{"DBNAME"});
    $database_cfg->key_by_parent('REFERENCE_DB', '-host', $CONFIG->{$species}->{"DBHOST"});
    $database_cfg->key_by_parent('REFERENCE_DB', '-user', $ROUSER);
    $database_cfg->delete_key_by_parent('REFERENCE_DB', '-pass');
    $database_cfg->key_by_parent('REFERENCE_DB', '-port', $CONFIG->{$species}->{"DBPORT"});
    $database_cfg->key_by_parent('GENEBUILD_DB', '-dbname', $CONFIG->{$species}->{"WRITENAME"});
    $database_cfg->key_by_parent('GENEBUILD_DB', '-host', $CONFIG->{$species}->{"WRITEHOST"});
    $database_cfg->key_by_parent('GENEBUILD_DB', '-user', $WRITEUSER);
    $database_cfg->key_by_parent('GENEBUILD_DB', '-pass', $pass);
    $database_cfg->key_by_parent('GENEBUILD_DB', '-port', $CONFIG->{$species}->{"WRITEPORT"});
    $database_cfg->write_config(1);
    # We can simply do a cp of the example file here
    my $blast_general_cfg = Bio::EnsEMBL::Analysis::Tools::ConfigWriter->new( -modulename => 'Bio::EnsEMBL::Analysis::Config::Blast',
        -is_example => 1);
    $blast_general_cfg->moduledir($output_dir);
    $blast_general_cfg->write_config(1);
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
      chomp($reply);
      if ($reply =~ /^Y/ or $reply =~ /^y/){    
          $status  = system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass -e\"DROP DATABASE IF EXISTS $WRITENAME;\"");
          if($status && $status != 256){ die("couldnt drop old database $WRITENAME!\n"); 
	  }
      } else {
      die "Exiting\n";
	}
    }
    print "Creating database $WRITENAME @ $WRITEHOST : $WRITEPORT \n" if $verbose;
    $status  = system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass -e\"CREATE DATABASE $WRITENAME;\"");
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < "."$CVSDIR/ensembl/sql/table.sql");
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < "."$CVSDIR/ensembl-pipeline/sql/table.sql");
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < "."$CVSDIR/ensembl-pipeline/sql/flag.sql");
    print ".";;
    #copy defined db tables from current build
    my $cmd = "mysqldump --add-drop-table -u$WRITEUSER -p$pass -h$REFDBHOST -P$REFDBPORT $REFDBNAME ".
      " assembly  coord_system meta meta_coord".
      " assembly_exception seq_region seq_region_attrib ".
      " unmapped_reason  " .
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
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < ".$DATADIR."/$species/import_tables.sql");
    print ".";

    $cmd = "mysqldump -u ensro -hens-staging1 ensembl_production master_attrib_type master_external_db > " .$DATADIR."/$species/master.sql";
    $status += system($cmd);
    print ".";
    $status += system("mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass $WRITENAME < ".$DATADIR."/$species/master.sql");
    print ".";        
    $cmd = "mysql -h$WRITEHOST -P$WRITEPORT -u$WRITEUSER -p$pass -e  '".
    "drop table attrib_type; drop table external_db; rename table master_attrib_type to attrib_type; rename table master_external_db to external_db;' $WRITENAME ";
    $status += system($cmd);
    print ".";

     #insert analysis entries
    # dont automatically run rfam jobs uless specifically requested
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT -dbpass $pass".
	" -logic_name ncRNA -program cmsearch -program_file /software/ensembl/bin/cmsearch_1.0 -database Rfam -database_version $RFAMVERSION -database_file $BLASTDIR ".
	  " -module Bio::EnsEMBL::Analysis::RunnableDB::Infernal".
	    " module_version 1 -gff_source ensembl -gff_feature gene -input_id_type FLAG";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -logic_name RfamBlast -program wublastn -program_file wublastn -database Rfam ".
	  " -database_file $BLASTDIR/filtered.fasta  -parameters  \'W=12 B=10000 V=10000 -hspmax 0 -gspmax 0 -kap -cpus=1\'" .
	    " -module Bio::EnsEMBL::Analysis::RunnableDB::BlastRfam -database_version $RFAMVERSION ".
	      " module_version 1 -gff_source ensembl -gff_feature gene -input_id_type SLICE";
    $status += system($cmd);
    $cmd = "perl ".$CVSDIR."/ensembl-pipeline/scripts/add_Analysis ".
      " -dbhost $WRITEHOST -dbname $WRITENAME -dbuser $WRITEUSER  -dbport $WRITEPORT  -dbpass $pass".
	" -logic_name BlastmiRNA -program wublastn -program_file wublastn -database  all_mirnas.fa ".
	  "-database_file $BLASTDIR/all_mirnas.fa -database_version $MIRBASEVERSION ".
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
  # big families
  my %avoid = ( "RF00177" => 1, #ssu
		"RF00028" => 1, #GroupI Intron
		"RF00029" => 1, #Group II Intron
		"RF00005" => 1
	      );
  my %families;
  
  # create Rfam.seed file
  my $exit;
  print "Updating RFAM descriptions file ... using Sanger FTP site (consider using EBI mirror)\n";
  system ("mkdir $BLASTDIR") unless -e "$BLASTDIR";
  $exit =  system ("wget ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz  -O $BLASTDIR/Rfam.cm.gz");
  die ("Error with obtaining Rfam covariance model  file from
  ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz\n") if $exit > 0;
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.cm.gz");
  die ("Error decompressing Rfam.tar.gz\n") if $exit > 0;
  $exit =  system ("wget ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz  -O $BLASTDIR/Rfam.seed.gz");
  die ("Error with obtaining Rfam.seed file from ftp://ftp.sanger.ac.uk/pub/databases/Rfam/Rfam.seed.gz\n") if $exit > 0;
  $exit =   system ("wget ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.fasta.gz -O $BLASTDIR/Rfam.fasta.gz");
  die ("Error with obtaining Rfam.fasta file from ftp://ftp.sanger.ac.uk/pub/databases/Rfam/Rfam.fasta.gz\n") if $exit > 0;
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.seed.gz");
  die ("Error decompressing Rfam.seed.gz\n") if $exit > 0;
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.fasta.gz");
  die ("Error decompressing Rfam.fasta.gz\n") if $exit > 0;   
  $exit =   system ("gzip -d  $BLASTDIR/Rfam.thr.gz");
  print "done\nUpdating miRNA file...";
  $exit =   system ("wget ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz -O $BLASTDIR/all_mirnas.embl.gz");
  $exit =   system ("gzip -d  $BLASTDIR/all_mirnas.embl.gz");
  die ("Error with obtaining miRNA.dat  file from ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat\n") if $exit > 0;
  $exit =   system ("wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz  -O $BLASTDIR/all_mirnas.fa.gz");
  $exit =   system ("gzip -d  $BLASTDIR/all_mirnas.fa.gz");
  die ("Error with obtaining hairpin.fa  file from ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa\n") if $exit > 0;
  print "done\n";
  # use bioperl to parse it
#    my $miRNA = Bio::SeqIO-> new
#    (
#     -file => "$BLASTDIR/all_mirnas.embl",
#     -format => "embl",
#    );
#  my $miRNA_fasta = Bio::SeqIO-> new
#    (
#     -file => ">>$BLASTDIR/all_mirnas.fa",
#     -format => "Fasta",
#    );
#  # write to fasta file 
#  while (my $seq = $miRNA->next_seq){
#    $seq->desc($seq->accession." ".$seq->desc);
#    $miRNA_fasta->write_seq($seq);
#  }

  print "Fetching all Rfam fasta sequences...";
  my $RFAM = Bio::SeqIO-> new
    (
     -file => "$BLASTDIR/Rfam.fasta",
     -format => "Fasta",
    );
  
  my $filtered_sequences = Bio::SeqIO-> new
    (
     -file => ">$BLASTDIR/filtered.fasta",
     -format => "Fasta",
    );
  
  die ("Cannot open  /data/blastdb/Rfam/Rfam.fasta\n") unless $RFAM;
  die ("Cannot open  $BLASTDIR/filtered.fasta\n") unless $filtered_sequences;
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
    my $domain = $1 if $seq->display_id =~ /^(RF\d+);.+/;
    if ($families{$domain}){
      unless ($avoid{$domain}){
	$filtered_sequences->write_seq($seq);
      }
    }
  }
  print "Done\n Formatting blast databases....\n";
  
  system ("xdformat -n $BLASTDIR/all_mirnas.fa");
  system ("xdformat -n $BLASTDIR/filtered.fasta");
  
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

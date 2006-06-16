#!/usr/local/ensembl/bin/perl

=head1 NAME

extra_database_setup.pl 

=head1 SYNOPSIS

This will dump specified tables from database

perl extra_database_setup.pl -dbhost host -dbuser user -dbpass **** -dbport 3306 -dbname yourdatabase -genebuild_database -output_dir $MYSQL_OUT/mouse_test_dump/ -verbose -dump

this will load specified tables into database

perl extra_database_setup.pl -dbuser user -dbpass ** -target_info database@host:3306 -genebuild_database -output_dir $MYSQL_OUT/mouse_test_dump/ -verbose -load

=head1 DESCRIPTION

This script automates the process of dumping tables from and existing 
database and loading tables to a new database

=head1 OPTIONS
   
  -dbhost    host name for database (gets put as host= in locator)
  
  -dbport    For RDBs, what port to connect to (port= in locator)
  
  -dbname    For RDBs, what name to connect to (dbname= in locator)
  
  -dbuser    For RDBs, what username to connect as (user= in locator)
  
  -dbpass    For RDBs, what password to use (pass= in locator)
  
  -help      Displays script documentation with PERLDOC

  -target_info This should be information about a database to load data 
               into. The format of this string should be 
               dbname@host:port:user:pass. Port, user and pass are 
               optional. Port if not specified will be assumed to be 3306 
               and pass and user will be taken from -dbuser and -dbpass. 

  -dump flag to dump the specified tables from the source database

  -load flag to load the specified tables into the target database/s

  -genebuild_database use a predefined list of tables as needed by a
                      non sequence genebuild database
 
  -genes uses a predefined list of tables for genes

  -sequence uses a predefined list of tables for sequence

  -table_name the name of a table to dump or load

  -output_dir the directory to dump into/load from this must be read and 
  writable by the mysqluser

  -verbose print out commands used

  -local This is needed when you are using mysql 4.1 and it is not setup to
         to allow imports from filesystems. You will see this error is this
         
         mysqlimport: Error: The used command is not allowed with this 
         MySQL  version, when using table: meta


=head1 EXAMPLES

This will dump specified tables from database

perl extra_database_setup.pl -dbhost host -dbuser user -dbpass **** -dbport 3306 -dbname yourdatabase -table_name seq_region -table_name seq_region_attrib -table_name attrib_type -output_dir $MYSQL_OUT/mouse_test_dump/ -verbose -dump

this will load specified tables into database

perl extra_database_setup.pl -target_info database@host:3306:user:pass
-genebuild_database -output_dir $MYSQL_OUT/mouse_test_dump/ -verbose -load


This will load the specified tables into all three databases

perl extra_database_setup.pl -load -output_dir /path/to/output_dir 
-verbose -table_name assembly -table_name coord_system -table_name 
seq_region -table_name analysis -target_info database1@host:3306:user:pass
 -target_info database2@host:3306:user:pass
 -target_info database3@host:3306:user:pass

=head1 NOTES

There are a couple of important notes to add to this. First it is a good
idea to both run the script and have the output directory mounted on the
same filesystem the mysqlinstance you are loading into is on to avoid
nfs issues 

Secondly this script runs mysqlimport with the -i option. This means if it
finds any duplicated lines between your dumped files and the database you
are loading into it will not load the duplicate lines but instead skip them

=cut
use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my @target_info;
my $dump;
my $load;
my @table_names;
my $genebuild_database;
my $help;
my $output_dir;
my @command_args = @ARGV;
my $verbose;
my $local;
my $genes;
my $sequence;

&GetOptions(
            'dbhost=s'            => \$dbhost,
            'dbname=s'            => \$dbname,
            'dbuser=s'            => \$dbuser,
            'dbpass=s'            => \$dbpass,
            'dbport=s'            => \$dbport,
            'help!'               => \$help,
            'target_info=s@'      => \@target_info,
            'dump!'               => \$dump,
            'load!'               => \$load,
            'genebuild_database!' => \$genebuild_database,
            'table_name=s@'       => \@table_names,
            'output_dir=s'        => \$output_dir,
            'verbose!'            => \$verbose,
            'local!'              => \$local,
            'genes!'              => \$genes,
            'sequence!'           => \$sequence,
           ) or perldocs("Couldnt get your options");

my $error_msg;
print join("\t", @command_args)."\n";
if($dump && (!$dbhost || !$dbuser || !$dbname)){
  $help = 1;
  $error_msg = "If you are dumping tables you must provide database ".
    "information from where you want to dump from\n";
}

if($load && !(@target_info >= 1)){
  $help = 1;
  $error_msg .= "If you want to load data into databases you must ".
    "provide target database information in the format ".
      'name @host:port:user:pass\n';
}

if($genebuild_database){
  push(@table_names, 'assembly', 'assembly_exception', 'analysis', 
       'attrib_type', 'coord_system', 'meta', 'meta_coord', 'seq_region',
       'seq_region_attrib');
}
if($genes){
  push(@table_names, 'gene', 'transcript', 'transcript_supporting_feature',
       'translation', 'exon_transcript', 'exon', 'supporting_feature',
       'protein_align_feature', 'dna_align_feature', 'gene_stable_id',
       'transcript_stable_id', 'translation_stable_id', 'exon_stable_id',
       'analysis', 'meta_coord');
}
if($sequence){
  push(@table_names, 'dna', 'assembly', 'assembly_exception', 'analysis', 
       'attrib_type', 'coord_system', 'meta', 'meta_coord', 'seq_region',
       'seq_region_attrib');
}
my %tables;
foreach my $table(@table_names){
  $tables{$table} = 1;
}
@table_names = keys(%tables);
if(!(@table_names  >= 1)){
  $help = 1;
  $error_msg .= "If you want to load or dump tables you must provide a ".
    "list of tables to load or dump\n";
}

if(!$output_dir || ! -e $output_dir || ! -d $output_dir){
  $help = 1;
  $error_msg .= "Your outputdir must be defined and exist and ".
    "be a directory\n";
}

if($help){
  perldocs($error_msg);
}

if($dump){
  my $source_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host   => $dbhost,
     -dbname => $dbname,
     -user   => $dbuser,
     -pass   => $dbpass,
     -port   => $dbport,
    );
  foreach my $table(@table_names){
    my $filename = $output_dir."/".$table.".".$dbname;
    my $dump_sql = "select * from $table into outfile '$filename'";
    print $dump_sql."\n" if($verbose);
    my $dump_sth = $source_db->dbc->prepare($dump_sql);
    $dump_sth->execute;
  }
}

if($load){
  foreach my $info(@target_info){
    my ($target_dbname, $host, $port, $user, $pass) = 
      db_info_parse($info, $dbuser, $dbpass);
    my $import_command = "mysqlimport -i -h$host -u$user -p$pass -P$port";
    $import_command .= "--local " if($local);
    $import_command .= " $target_dbname ";
    foreach my $table(@table_names){
      my $filename = $output_dir."/".$table.".*";
      my $cmd = $import_command." ".$filename;
      print $cmd."\n" if($verbose);
      system($cmd) == 0 or throw("Failed to run ".$cmd);
    }
  }
}


sub db_info_parse{
  my ($info, $source_user, $source_pass) = @_;
  my ($dbname, $other_info) = split /\@/, $info;
  my ($host, $port, $user, $pass) = split /\:/, $other_info;
  $user = $source_user if(!$user);
  $pass = $source_pass if(!$pass);
  $port = 3306 if(!$port);

  return ($dbname, $host, $port, $user, $pass);
}



sub perldocs{
  my ($msg) = @_;
  print $msg."\n" if($msg);
  exec('perldoc', $0);
  exit(0);
}

package TestDB;

use strict;
use warnings;
use lib '../test_system';

use vars qw(@ISA);
use Bio::EnsEMBL::Root;
use DBI;
use File::Path;

@ISA = qw(Bio::EnsEMBL::Root);

use Bio::EnsEMBL::Utils::Exception qw(throw warning verbose);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;


#homo sapiens is used if no species is specified
my $DEFAULT_SPECIES  = 'homo_sapiens';
my $CONF_FILE = 'TestDB.conf';


=head2 new

  Arg [1]   : string, TestDB name of package
  Arg [2]   : string, species name
  Arg [3]   : int, boolean flag for verbosity
  Arg [4]   : string, path to conf file
  Function  : create a TestDB object
  Returntype: TestDB
  Exceptions: none
  Example   : my $testdb = TestDB->new;

=cut



sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($species, $verbose, $conf_file, $local) = 
    rearrange (['SPECIES', 'VERBOSE', 'CONF_FILE', 'LOCAL'], @args);
  my $curr_dir = $ENV{'PWD'};
  $self->curr_dir($curr_dir);
  $self->species($species);  
  $self->verbosity($verbose);
  $self->conf_file($conf_file);
  $self->islocal($local);
  $self->load_databases;
  return $self;
}



=head2 Container methods

  Arg [1]   : TestDB
  Arg [2]   : string to be contained
  Function  : holder for a string, the exception to this is db which holds
  a DBI object
  Returntype: string
  Exceptions: none
  Example   : 

=cut



#containers
sub verbosity{
  my $self = shift;
  $self->{'verbose'} = shift if(@_);
  return $self->{'verbose'};
}

sub curr_dir{
  my $self = shift;
  $self->{'curr_dir'} = shift if(@_);
  return $self->{'curr_dir'};
}

sub islocal{
  my $self = shift;
  $self->{'islocal'} = shift if(@_);
  return $self->{'islocal'};
}

sub species{
  my $self = shift;
  $self->{'species'} = shift if(@_);
  return $self->{'species'} || $DEFAULT_SPECIES;
}
sub conf_file{
  my $self = shift;
  $self->{'conf_file'} = shift if(@_);
  return $self->{'conf_file'} || $self->curr_dir."/".$CONF_FILE;
}
sub data_dir{
  my $self = shift;
  $self->{'data_dir'} = shift if(@_);
  return $self->{'data_dir'};
}

sub dbi_connection{
  my $self = shift;
  $self->{'dbi'} = shift if(@_);
  return $self->{'dbi'};
}




=head2 conf_hash

  Arg [1]   : TestDB
  Arg [2]   : hashref for conf file
  Function  : holds the hashref of the test configuration
  Returntype: hashref
  Exceptions: throws if not passed a hashref
  Example   : my $dbname = $self->conf_hash->{'dbname'};

=cut


sub conf_hash{
  my ($self, $conf_hash) = @_;
  if($conf_hash){
    throw($conf_hash." must be a hash ref ") 
      unless(ref($conf_hash) eq 'HASH');
    $self->{'conf_hash'} = $conf_hash;
  }
  return $self->{'conf_hash'};
}



=head2 load_databases

  Arg [1]   : TestDB
  Function  : Creates a new database with a random name unless
  one is specifies in config and then loads the sql specified in config
  Returntype: none
  Exceptions: throws if config file doesn't exist or if can't create the
  database specifed
  Example   : 

=cut


sub load_databases {
  my ($self) = shift;

  print "\nTrying to load [$self->{'species'}] databases\n" 
    if($self->verbosity);

  #create database from conf and from zip files 
  
  if(!-e $self->conf_file){
    throw("Can't load databases from ".$self->conf_file." if doesn't ".
          "exist");
  }
  my $db_conf = do $self->conf_file;

  my $port = $db_conf->{'port'};
  my $driver = $db_conf->{'driver'};
  my $host = $db_conf->{'host'};
  my $pass = $db_conf->{'pass'};
  my $user = $db_conf->{'user'};
  my $dbname = $db_conf->{'dbname'};
  my $sql_files = $db_conf->{'sql_files'};
  my $preloaded = $db_conf->{'preloaded_tables'};
  $self->data_dir($db_conf->{'data_dir'});

  $self->conf_hash($db_conf);

  #connect to the database
  my $locator = "DBI:".$driver.":host=".$host.";port=".$port;
  print "Connecting to ".$locator."\n" if($self->verbosity);
  my $db = DBI->connect($locator, $user, $pass, {RaiseError => 1});
  unless($db) {
    throw("Can't connect to database $locator");
    return;
  }
  $self->dbi_connection($db);
  if(!$db_conf->{'dont_unzip'}){
    $self->unzip_data_files();
  }
  if(!$dbname){
    if($preloaded){
      $self->exception("How can you database have preloaded tables if ".
                       "you didn't specify the database name in ".
                       $self->conf_file);
    }
    $dbname = $self->_create_db_name;
    $self->conf_hash->{'dbname'} = $dbname;
  }
  if(!$preloaded){
    my $sql = "CREATE DATABASE $dbname";
    print "Creating database ".$sql."\n" if($self->verbosity);
    unless($self->dbi_connection->do("CREATE DATABASE $dbname")) {
      $self->exception("Couldn't not create database");
    }
  }
  $db->do("use $dbname");
  if(!$preloaded){
    $self->do_sql_files($sql_files);
  }
}



=head2 do_sql_files

  Arg [1]   : TestDB
  Arg [2]   : arrayref, list of sql files to load
  Function  : load sql table definitions into the established DBI 
  connection
  Returntype: int, the number of sql statements executed
  Exceptions: throws if can't open sql file'
  Example   : 

=cut



sub do_sql_files {
    my( $self, $files ) = @_;
    local *SQL;
    my $i = 0;
    my $dbh = $self->dbi_connection;

    my $comment_strip_warned=0;

    foreach my $file (@$files){
      my $sql = '';
      open SQL, $file or $self->exception("Can't read SQL file ".
                                          "'$file' : $!");
      while (<SQL>) {
        # careful with stripping out comments; quoted text
        # (e.g. aligments) may contain them. Just warn (once) and ignore
        if (    /'[^']*#[^']*'/ 
          || /'[^']*--[^']*'/ ) {
            if ( $comment_strip_warned++ ) { 
              # already warned
            } else {
              warning("#################################\n".
                      "# found comment strings inside quoted string; ".
                      "not stripping, too complicated: $_\n".
                      "# (continuing, assuming all these they are simply ".
                      "valid quoted strings)\n".
                      "#################################\n");
            }
          } else {
            s/(#|--).*//;       # Remove comments
             }
            next unless /\S/;   # Skip lines which are all space
            $sql .= $_;
            $sql .= ' ';
          }
      close SQL;
      
      #Modified split statement, only semicolumns before end of line,
      #so we can have them inside a string in the statement
      #\s*\n, takes in account the case when there is space before the 
      #new line
      foreach my $s (grep /\S/, split /;[ \t]*\n/, $sql) {
        $self->validate_sql($s);
        $dbh->do($s);
        $i++
      }
    }
    return $i;
} 


=head2 db

  Arg [1]   : TestDB
  Arg [2]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor (optional)
  Function  : holds a pipeline dbadaptor and creates one from config
  if one is requested but not defined
  Returntype: Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Exceptions: throws if not passed a 
  Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Example   : 

=cut



sub db{
  my ($self, $db) = @_;
  if($db){
    if(!$db->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor')){
      $self->exception("Can't run the RuleManager with $db you need a ".
                       "Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor");
    }
    $self->{'pipeline_adaptor'} = $db;
  }
  if(!$self->{'pipeline_adaptor'}){
    $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new
      (-dbname => $self->conf_hash->{'dbname'},
       -user   => $self->conf_hash->{'user'},
       -pass   => $self->conf_hash->{'pass'},
       -port   => $self->conf_hash->{'port'},
       -host   => $self->conf_hash->{'host'},
       -driver => $self->conf_hash->{'driver'}
      );
    $self->{'pipeline_adaptor'} = $db;
  }
  return  $self->{'pipeline_adaptor'};
}


=head2 validate_sql

  Arg [1]   : TestDB
  Arg [2]   : string, sql statement
  Function  : checks sql statement is an insert with explicit column 
  definitions if it is an insert
  Returntype: none
  Exceptions: throws if is an insert without explicit column definitions
  Example   : 

=cut



sub validate_sql {
    my ($self, $statement) = @_;
    if ($statement =~ /insert/i)
    {
        $statement =~ s/\n/ /g; #remove newlines
        $self->exception("INSERT should use explicit column names ".
                         "(-c switch in mysqldump)\n$statement\n")
          unless ($statement =~ /insert.+into.*\(.+\).+values.*\(.+\)/i);
    }
}


=head2 load_data

  Arg [1]   : TestDB
  Arg [2]   : string file name must be in the format table_name.ext
  Function  : constructs and calls a mysqlimport statement to import
  data in database specified by config
  Returntype: none
  Exceptions: throws if the mysqlimport doesn't return a 0 exit status'
  Example   : $self->load_data('seq_region.sql');

=cut


sub load_data{
  my ($self, $file) = @_;
  if($self->{'preloaded_data'}){
    return 1;
  }
  my $db_conf = $self->conf_hash;
  my $port = $db_conf->{'port'};
  my $host = $db_conf->{'host'};
  my $pass = $db_conf->{'pass'};
  my $user = $db_conf->{'user'};
  my $dbname = $db_conf->{'dbname'};
  my $cmd = "mysqlimport -i -h$host -u$user -p$pass -P$port ".
    ( $self->islocal ? "--local " : "") . " $dbname $file";
  print $cmd."\n" if($self->verbosity);
  system($cmd) == 0 or $self->exception("Failed to run ".$cmd);
  return 1;
}


=head2 load_*table_set*_tables

  Arg [1]   : TestDB
  Arg [2]   : string path to directory where files are located
  Function  : A set of methods which loads a predefined set of tables
  into the database
  Returntype: none
  Exceptions: throws if data directory doesnt exist
  Example   : 

=cut



sub load_core_tables{
  my ($self, $data_dir) = @_;
  my @tables = ('meta', 'meta_coord', 'coord_system', 'analysis', 
                'attrib_type', 'seq_region', 'assembly', 'dna', 
                'seq_region_attrib', 'marker');
  if(!$data_dir){
    $data_dir = $self->curr_dir."/".$self->species;
  }
  if(! -d $data_dir){
    $self->exception("Can't import core tables if ".$data_dir.
                     " doesnt exist");
  }
  $self->load_tables(\@tables, $data_dir);
}

sub load_pipeline_tables{
  my ($self, $data_dir) = @_;
  my @tables = ('rule_goal', 'rule_conditions', 'input_id_type_analysis',
                'input_id_analysis');
  if(!$data_dir){
    $data_dir = $self->curr_dir."/".$self->species;
  }
  if(! -d $data_dir){
    $self->exception("Can't import core tables if ".$data_dir." doesnt exist");
  }
  $self->load_tables(\@tables, $data_dir);
}

sub load_gene_tables{
  my ($self, $data_dir) = @_;
  my @tables = ('gene', 'exon', 'transcript', 'translation',
                'exon_transcript', 'supporting_feature',
                'protein_align_feature', 'dna_align_feature',
                'gene_stable_id', 'exon_stable_id',
                'translation_stable_id', 'transcript_stable_id');
  if(!$data_dir){
    $data_dir = $self->curr_dir."/".$self->species;
  }
  if(! -d $data_dir){
    $self->exception("Can't import core tables if ".$data_dir." doesnt exist");
  }
  $self->load_tables(\@tables, $data_dir);
}

sub load_raw_compute_tables{
  my ($self, $data_dir) = @_;
  my @tables = ('repeat_consensus', 'repeat_feature',
                'prediction_exon', 'prediction_transcript',
                'dna_align_feature', 'protein_align_feature',
                'simple_feature');
  if(!$data_dir){
    $data_dir = $self->curr_dir."/".$self->species;
  }
  if(! -d $data_dir){
    $self->exception("Can't import core tables if ".$data_dir." doesnt exist");
  }
  $self->load_tables(\@tables, $data_dir);
}


=head2 load_tables

  Arg [1]   : TestDB
  Arg [2]   : arrayref list of table names
  Arg [3]   : string path to data directory
  Function  : constructs a full path for each table_name and calls to
  load_data to run the mysqlimport
  Returntype: none
  Exceptions: throws if filepath constructed doesnt exist
  Example   : 

=cut


sub load_tables{
  my ($self, $tables, $data_dir) = @_;
  foreach my $table(@$tables){
    my $file = $data_dir."/".$table;
    #print STDERR "Trying to load ".$file."\n";
    if(! -e $file){
      $self->exception("Can't load ".$data_dir."/".$table." it doesn't exist");
    }else{
      $self->load_data($file);
    }
  }
}


=head2 unzip_data_files

  Arg [1]   : TestDB
  Arg [2]   : string, path to directory containing file to unzip
  Function  : unzips a file containing all the sql data. The filename
  is based on species_name as is the directory where the data is unzipped 
  to
  Returntype: none
  Exceptions: throws if data dir doesnt exist, throws if zip file doesnt
  exist, also throws if the destination dir already exists
  Example   : 

=cut



sub unzip_data_files{
  my ($self, $data_dir) = @_;
  if(!$data_dir){
    $data_dir = $self->data_dir;
  }
  if(! -d $data_dir){
    $self->exception("Can't unzip data if directory data lives in: $data_dir ".
          "doesn't exist");
  }
  my $zip_file = $self->species.".zip";
  my $path = $data_dir."/".$zip_file;
  my $dest_dir = $self->curr_dir."/".$self->species;
  if(!-e $path){
    $self->exception("Can't unzip ".$zip_file." if ".$path." doesn't exist");
  }
  if(! -d $dest_dir){
    mkdir($dest_dir);
  }elsif(-d $dest_dir){
    $self->exception("Directory ".$dest_dir." already exists the data might ".
          "not be correct");
  }
  my $cmd = "unzip -q $path -d $dest_dir ";
  print $cmd."\n" if($self->verbosity);
  system($cmd) == 0 or $self->exception("Error running ".$cmd);
}



=head2 _create_db_name

  Arg [1]   : TestDB
  Function  : creates a random name for the database based on the username
  date and time
  Returntype: string, dbname
  Exceptions: none
  Example   : my $dbname = $self->_create_db_name;

=cut



sub _create_db_name {
    my( $self) = @_;

    my @t_info = localtime;

    my $date = join ( "_", $t_info[3],$t_info[4]+1);  
    my $time = join ( "", $t_info[2],$t_info[1],$t_info[0]);  

    my $species = $self->species;

    # create a unique name using host and date / time info
    my $db_name = "_test_db_${species}_pipeline_".$ENV{'USER'}."_".$date.
      "_".$time;
    return $db_name;
}


=head2 cleanup

  Arg [1]   : TestDB
  Function  : deletes the directory the data was unzipped into and
  drops the database
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub cleanup{
  my ($self) = @_;
  my $data_dir = $self->curr_dir."/".$self->species;
  print "Deleting data from ".$data_dir."\n" if($self->verbosity);
  rmtree($data_dir) unless($self->conf_hash->{'dont_unzip'});
  print "Dropping database ".$self->conf_hash->{'dbname'}."\n" 
    if($self->verbosity);
  $self->dbi_connection->do("Drop database ".$self->conf_hash->{'dbname'})
    unless($self->conf_hash->{'preloaded_tables'});
  $self->dbi_connection->disconnect;
}

=head2 database_args

  Arg [1]   : TestDB
  Function  : create a string of the standard database args for
  ensembl-pipeline scripts for the script RunTest runs
  Returntype: string
  Exceptions: none
  Example   : 

=cut


sub database_args{
  my ($self) = @_;
  my $db_conf = $self->conf_hash;
  my $dbport = $db_conf->{'port'};
  my $dbhost = $db_conf->{'host'};
  my $dbpass = $db_conf->{'pass'};
  my $dbuser = $db_conf->{'user'};
  my $dbname = $db_conf->{'dbname'};
  my $db_args = " -dbhost ".$dbhost." -dbuser ".$dbuser;
  $db_args .= " -dbpass ".$dbpass if($dbpass);
  $db_args .= " -dbport ".$dbport if($dbport);
  $db_args .= " -dbname ".$dbname." ";
  return $db_args;
}

=head2 cleanup_command

  Arg [1]   : TestDB
  Function  : creates a command and prints it to screen to cleanup after
  a test run. This is to allow easy removal of test databases and output
  directories if you wish to keep the output for a little while for 
  investigation
  Returntype: none
  Exceptions: 
  Example   : 

=cut
sub cleanup_command{
  my ($self) = @_;
  my $db_args;
  if($self->conf_hash->{'dbname'}){
    $db_args = $self->database_args;
  }
  my $data_dir = $self->curr_dir."/".$self->species;
  my $cleanup_command = "cleanup_output.pl ";
  $cleanup_command .= $db_args." " if($db_args);
  $cleanup_command .= " -sql_data_dir ".$data_dir;
  print "You have specifed -dont_cleanup when running your test \n".
    "If you want to delete your output you can run this script ".
      "ensembl-pipeline/test_system/cleanup_output.pl\n".
        "this is the command you should use \n".$cleanup_command."\n".
          "If you don't want any of the data sets deleted remove either".
              " -dbname, -sql_data_dir or -output_dir options from the ".
                "commandline\n";
}



sub exception{
  my ($self, $msg) = @_;
  $self->cleanup_command;
  throw($msg);
}

1;

#
# Config.pm - Object which reads and stores configuration information
#
# 
# You may distribute this module under the same terms as perl itself
#

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Config - Reads and stores configuration information

=head1 SYNOPSIS

  use Bio::EnsEMBL::Pipeline::Pipeline::Config;

  #
  # Read config from files and store it in the db (as specified in files)
  #
  my $config = Bio::EnsEMBL::Pipeline::Config->new(-files => \@_);

  #print out the configuration 
  foreach my $header ($config->get_headers) {
    print "[$header]\n";
    foreach my $key ($config->get_keys) {
      print "$key = ", $config->get_parameter($key), "\n";
    }
  }
  
  # Read config which has been stored in the database
  $config = Bio:EnsEMBL::Pipeline::Config->new
    (-dbname => 'pipedb',
     -host   => 'ecs1g',
     -user   => 'ensadmin',
     -pass   => 'ensembl');

=head1 DESCRIPTION

This module reads a configuration when instantiated, either from a set of 
files or from a database.  

If read from files the configuration is then stored in the database specified 
within the configuration file.  A required header in the configuration file is
therefore [PIPELINE_DATABASE].  If the config table in the database pointed to 
by the files is already populated with data the constructor will throw an 
exception and will not load the new configuration or alter the existing 
configuration in the database.

Configuration information can be retrieved from a configuration object via the 
accessor methods get_headers, get_keys, get_parameter.

Configuration files have the following structure:

  [DEFAULT_DATABSE]
  user = ensadmin
  host=ecs1g
  user=ensadmin
  pass=ensembl
  port=3306

  [PIPELINE_DATABASE]
  dbname = mcvicker_test_briggsae

  #comment
  [DEFAULT_HEADER]
  key1 = value
  key2 = value

  [MY_HEADER]
  key1 = another value
  key3 = test #not a comment

The headers which begin with DEFAULT_ are special headers and their keys and
values are propogated to all subsequent headers with the same postfix.  Default
keys may be overridden, however.  

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Pipeline::Config;

use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

@ISA = ('Bio::EnsEMBL::Root');



=head2 new

  Arg [-FILES]: listref $filenames
  Arg [-DB]   : a dbadaptor for the pipeline
  Example    : $conf = Bio::EnsEMBL::Pipeline::Config->new(-FILE => $filename);
  Description: Creates a new pipeline configuration object.  Configuration
               may be read from files or from the database.
  Returntype : Bio::EnsEMBL::Pipeline::Config
  Exceptions :
  Caller     :

=cut

sub new {

  my $class = shift;

  my $self = bless {}, $class;

 
  $self->{'config'} = {};  #hash of hashes - will hold config as we build it up

  my ($files, $db) = $self->_rearrange([qw( FILES DB )], @_);

  # Files and DB should not be defined
  if ($files && @$files && $db) {
    $self->throw("Cannot read config from files and db at the same time.");
  }

  if($files && @$files) {
    $self->_parse_files(@$files);

  } elsif($db) {
    # store the DB
    $self->{'_dbobj'} = $db;

    my $config_rows = $self->_read_db($db);
  }

  $self->_update_all_defaults();

  # if using a file, write the contents to the database
  # note this needs to be done *after* the call to _update_all_defaults
   if ($files ) {
     # Store contents of $self->{'config'} in database
    $self->_write_config_to_db();
   }

  return $self;

}

sub _read_db {

  my ($self, $dbobj) = @_;

  my $stmt = $dbobj->prepare("SELECT * FROM config");
  my $result = $stmt->execute;

  my $rows = 0;

  while (my $ref = $stmt->fetchrow_hashref) {
    my $header  = lc($ref->{'header'});
    my $keyname = lc($ref->{'key_name'});
    my $value   = $ref->{'value'};
    $rows++;
    #print "|$header"  . "\t|$keyname" . "\t|$value" . "\t|\n"
    $self->{'config'}->{$header}->{$keyname} = $value;
  }

  $stmt->finish();

  return $rows;
}



sub _write_config_to_db {

  my $self = shift;

  # get the config values that will be required to set up the DBAdaptor
  my @required_params = ( "host", "user", "pass", "dbname", "port" );
  my %params;
  foreach my $param (@required_params) {
    $params{$param} =  $self->get_parameter('pipeline_database', $param);
    #print "Got pipeline database parameter $param, value = $params{$param}\n";
  }

  # create the DBAdaptor
  my $dbobj = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
    ( '-host'   => $params{'host'},
      '-port'   => $params{'port'},
      '-dbname' => $params{'dbname'},
      '-user'   => $params{'user'},
      '-pass'   => $params{'pass'});

  # check if there is any existing config and throw if there is
  my $stmt = $dbobj->prepare("SELECT count(*) FROM config");
  my $result = $stmt->execute;
  my ($count) = $stmt->fetchrow_array;
  $stmt->finish();
  if ($count != 0) {
    $self->throw("config table in " . $params{'dbname'} . " on " . 
        $params{'host'} . " already contains $count rows; it should be empty");
  }

  # do the insert
  $stmt = $dbobj->prepare("INSERT INTO config VALUES (?,?,?)");

  foreach my $header ($self->get_headers()) {

    foreach my $key ($self->get_keys($header)) {

      my $value = $self->get_parameter($header, $key);
      my $result = $stmt->execute($header, $key, $value);

    }

  }

  $stmt->finish();

  $self->{'_dbobj'} = $dbobj;

}


sub _parse_files {


  my $self = shift;
  my @files = shift;

  my %headers;     # will store names of headers and number of keys for each

  # read each file

  foreach my $file (@files) {

    if (! -e $file) {
      $self->throw("Configuration file $file not found\n");
    }
    my $header = "";

    open (FILE, $file);
    while (<FILE>) {
      chomp();

      # Comment or blank line
      next if (/^\s$/ || /^\#/);

      # [HEADER]
      if (/^\[(.*)\]$/) {         # $1 will be the header name, without the []
	$header = lc($1);
	$headers{$header} = 0;
	#print "Reading stanza $header\n";
      }

      # key=value
      if (/^(\S+)\s*=\s*(\S+)/) {   # $1 = key, $2 = value

	my $key = lc($1);           # keys stored as all lowercase, values have case preserved
	my $value = $2;

	if (length($header) == 0) {
	  $self->throw("Found key/value pair $key/$value outside stanza");
	}

	#print "Key: $key Value: $value\n";

	# Check if this header/key is already defined
	if (exists($self->{'config'}->{$header}->{$key})) {
	  $self->throw("$key is already defined for [$header]; cannot be redefined");
	} else {
	  # store them in the config hash
	  $self->{'config'}->{$header}->{$key} = $value;
	  #print "$header:$key=$value\n";
	  $headers{$header}++;  # will be used to check for headers with no keys
	}

      }

    } # while <FILE>

    close FILE;
  }


  # add a blank key/value for any headers that have no keys
  foreach my $h (keys (%headers)) {
    if ($headers{$h} == 0) {
      #print "$h has no keys; adding blank key/value\n";
      $self->{'config'}->{$h}->{""} = "";
    }
  }

}

# modify the config hash so that default values are stored where no value is specified
sub _update_all_defaults {

  my $self = shift;

  my @headers = $self->get_headers();
  foreach my $header (@headers) {

    next if lc($header) =~ "^default";

    if ($header =~ /^\w+_database$/i) { $self->_update_single_header_defaults($header, "DEFAULT_DATABASE"); }

    elsif ($header =~ /^\w+_task$/i) { $self->_update_single_header_defaults($header, "DEFAULT_TASK"); }

    else { $self->_update_single_header_defaults($header, "DEFAULT"); }
	
  }

}

# put in defaults for all the keys for a specific header
sub _update_single_header_defaults {

  my $self = shift;
  my $header = lc(shift);
  my $default_header = lc(shift);

  if (exists($self->{'config'}->{$default_header})) {  # don't break if no default header is defined

    my @default_keys = $self->get_keys($default_header);

    my @current_keys = $self->get_keys($header);

    foreach my $key (@default_keys) {

      # if the key is not defined in current_keys, add the default
      if (!grep /^$key$/, @current_keys) {
	my $value = $self->get_parameter($default_header, $key);
	$self->{'config'}->{$header}->{$key} = $value;
      }
    }

  }

}

=head2 get_parameter

  Arg [1]    : string $header
  Arg [2]    : string $key
  Example    : 
  Description: Gives a value, provided with a config header and key
  Returntype : string
  Exceptions : 
  Caller     : 

=cut

sub get_parameter {

  my $self      = shift;
  my $header = lc(shift);
  my $key    = lc(shift);

  if (!exists($self->{'config'}->{$header}->{$key})) {
    my ($p, $f, $l) = caller;
    $self->warn("Config key [$key] is not defined in header [$header] $f:$l");
    return undef;
  }
  return $self->{'config'}->{$header}->{$key};

}



=head2 get_keys

  Arg [1]    : string $header
  Example    : 
  Description: Gives a list of keys for a given header.  Keys for defaults
               provided by a DEFAULT class are also provided.
  Returntype : list of strings
  Exceptions : 
  Caller     : 

=cut

sub get_keys {

  my $self      = shift;
  my $header = lc(shift);

  if (!exists($self->{'config'}->{$header})) {
    $self->throw ("No such header: $header");
  }
  my $keys = $self->{'config'}->{$header};
  return keys(%$keys);
}



=head2 get_headers

  Arg [1]    : None
  Example    : 
  Description: Returns a list of headers present in the configuration,
               including default headers
  Returntype : list of strings
  Exceptions : none
  Caller     : general

=cut

sub get_headers {

  my $self = shift;

  my $headers = $self->{'config'};
  return keys(%$headers);

}



=head2 get_DBAdaptor

  Arg [1]    : (optional)
  Example    : $db = $conf->db();
  Description: Getter/Setter for the dbadaptor to the pipeline database, 
               which is either specified directly in the constructor or
               read from a configuration file.
  Returntype : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub get_DBAdaptor {

  my $self = shift;

  return $self->{'_dbobj'};

}

1;

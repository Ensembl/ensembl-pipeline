use strict;
use warnings;


package Bio::EnsEMBL::Pipeline::Config;

use vars qw(@ISA);

use Bio::EnsEMBL::Root;

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

  my ($files) = $self->_rearrange([qw( FILES )], @_);

  my @files = split(/ /, $files);

  $self->{'config'} = {};


  # read each file

  foreach my $file (@files) {

    print "Parsing $file \n";

    my $header = "";

    open (FILE, $file);
    while (<FILE>) {
      chomp();

      # Comment or blank line
      next if (/^\s$/ || /^\#/);

      # [HEADER]
      if (/^\[(.*)\]$/) {         # $1 will be the header name, without the []
	$header = lc($1);
	print "Reading stanza $header\n";
      }

      # key=value
      if (/^(\w+)\s*=\s*(\w+)/) {   # $1 = key, $2 = value

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
	  # store them
	  $self->{'config'}->{$header}->{$key} = $value;
	}

      }

    }
    close FILE;
  }


 # TODO DB Stuff

  return $self;

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
    $self->throw("Key $key is not defined in header $header");
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

  Arg [1]    : string $header
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



=head2 db

  Arg [1]    : (optional)
  Example    : $db = $conf->db();
  Description: Getter/Setter for the dbadaptor to the pipeline database, 
               which is either specified directly in the constructor or
               read from a configuration file.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

sub db {

}



1;

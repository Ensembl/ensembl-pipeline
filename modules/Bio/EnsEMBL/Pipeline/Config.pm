use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::Config;



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
	my $self = shift;
	my $header = shift;
	my $key = shift;
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
	my $self = shift;
	my $header = shift;

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

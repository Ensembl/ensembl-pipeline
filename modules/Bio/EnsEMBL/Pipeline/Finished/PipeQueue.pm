package Bio::EnsEMBL::Pipeline::Finished::PipeQueue;

=head1 NAME

Bio::EnsEMBL::Pipeline::Finished::PipeQueue - access to pipe_queue.queue

=head1 DESCRIPTION

This module wraps up access to the otterlive.pipe_queue database.

It is still not in the Ensembl canonical "DBAdaptor" style, but it is
better than having multiple copies of the C<< DBI->connect >> code.

=head1 METHODS

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Pipeline::Config::BatchQueue qw( QUEUE_HOST QUEUE_NAME );

use DBI;


=head2 configure($queue_host, $queue_name)

Replace the default configuration, which is taken from
L<Bio::EnsEMBL::Pipeline::Config::BatchQueue>

=cut

{
    # State used by the class.  Overridden by dequeuer.pl
    my ($queue_host, $queue_name) = ($QUEUE_HOST, $QUEUE_NAME);

    sub configure {
        my $called = shift;
        ($queue_host, $queue_name) = @_;
    }

    sub _get_config {
        return ($queue_host, $queue_name);
    }
}


=head2 qdbh()

Class method.  Returns a new L<DBI> pointing at the currently
configured C<pipe_queue>.

Password comes from L<Net::Netrc> per team practice.

=cut

sub qdbh {
    my ($called) = @_;
    my ( $dbhost, $dbname ) = $called->_get_config;

    my ( $dbuser, $dbpass, $dbport ) = &get_db_param( $dbhost );
    my $dsn = "DBI:mysql:host=$dbhost;dbname=$dbname;port=$dbport";

    return DBI->connect( $dsn, $dbuser, $dbpass );
}

# XXX:DRY there are already many. fix later.
use Net::Netrc;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
sub get_db_param {
	my ( $dbhost ) = @_;
	my ( $dbuser, $dbpass, $dbport );

	my $ref = Net::Netrc->lookup($dbhost);
	throw("$dbhost entry is missing from ~/.netrc") unless ($ref);
	$dbuser = $ref->login;
	$dbpass = $ref->password;
	$dbport = $ref->account;
	throw(
		"Missing parameter in the ~/.netrc file:\n
			machine " .  ( $dbhost || 'missing' ) . "\n
			login " .    ( $dbuser || 'missing' ) . "\n
			password " . ( $dbpass || 'missing' ) . "\n
			account "
		  . ( $dbport || 'missing' )
		  . " (should be used to set the port number)"
	  )
	  unless ( $dbuser && $dbpass && $dbport );

	return ( $dbuser, $dbpass, $dbport );
}


1;

=head1 AUTHOR

Matthew Astley mca@sanger.ac.uk

=cut

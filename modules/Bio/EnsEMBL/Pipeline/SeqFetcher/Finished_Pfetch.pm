#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new(
    								-pfetch_server => '',
    								-pfetch_port => '',
    								-options => ''
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using pfetch.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;

use strict;
use Bio::Root::RootI;
use Bio::DB::RandomAccessI;
use Bio::Seq;
use IO::Socket;

use Bio::EnsEMBL::Pipeline::Tools::Embl;
use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);

BEGIN{
	print "Finished_Pfetch Loaded\n";
}

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;

    my ( $server, $port, $options ) = $self->_rearrange( [ 'PFETCH_SERVER', 'PFETCH_PORT', 'OPTIONS' ], @args );

    if ( !defined $server ) {
        $server = 'pubseq.internal.sanger.ac.uk';
    }
    $self->server($server);

    if ( !defined $port ) {
        $port = '22100';
    }
    $self->port($port);

    if ( defined $options ) {
        $self->options($options);
    }

    return $self;    # success - we hope!
}

=head2 server

  Title   : server
  Usage   : $self->server('address.of.server');
  Function: Get/set for the path to the server being used by the module. 
  Returns : string
  Args    : string

=cut

sub server {
    my ( $self, $server ) = @_;
    if ($server) {
        $self->{'_server'} = $server;
    }
    return $self->{'_server'};
}

=head2 port

  Title   : port
  Usage   : $self->port('port');
  Function: Get/set for the port to the pfetch server. 
  Returns : string
  Args    : string

=cut

sub port {
    my ( $self, $port ) = @_;
    if ($port) {
        $self->{'_port'} = $port;
    }
    return $self->{'_port'};
}

=head2 options

  Title   : options
  Usage   : $self->options('tc');
  Function: Get/set for options to pfetch
  Returns : string
  Args    : string

=cut

sub options {

    my ( $self, $options ) = @_;
    if ($options) {
        $self->{'_options'} = $options;
    }
    return $self->{'_options'};

}

sub get_server {
    my ($self) = @_;

    local $^W = 0;

    my $host = $self->server;
    my $port = $self->port;

    my $server = IO::Socket::INET->new(
        PeerAddr => $host,
        PeerPort => $port,
        Proto    => 'tcp',
        Type     => SOCK_STREAM,
        Timeout  => 10,
    );
    if ($server) {
        $server->autoflush(1);
        return $server;
    }

}

=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_eq_by_acc($accession);
  Function: Does the sequence retrieval via pfetch
  Returns : Bio::Seq
  Args    : 

=cut

sub get_Seq_by_acc {
    my ( $self, @id_list ) = @_;

    #confess "No names provided" unless @id_list;
    unless (@id_list) {
        $self->throw("No accession input");
    }

    my $server = $self->get_server();
    print $server "-q @id_list\n";
    my (@seq_list);
    for ( my $i = 0 ; $i < @id_list ; $i++ ) {
        chomp( my $seq_string = <$server> );
        eval {
            if ( defined $seq_string && $seq_string ne 'no match' )
            {

                my $seq = new Bio::Seq(
                    '-seq'              => $seq_string,
                    '-accession_number' => $id_list[$i],
                    '-display_id'       => $id_list[$i]
                );

                $self->throw("Could not pfetch sequence for $id_list[[$i]]\n") unless defined $seq;
                $seq_list[$i] = $seq;

            }

        };
        if ($@) {
            print STDERR "$@\n";
        }

    }

    if (wantarray) {

        # an array was passed - return array of Bio:Seq objects
        return @seq_list;

    }
    else {

        # one acc was passed then return the first(and only) element of the array    
        return $seq_list[0];
    }
}
sub write_descriptions {
    my ( $self, $dbobj, @ids ) = @_;
    my $sth = $dbobj->prepare( qq{ 
                                REPLACE DELAYED INTO 
                                hit_description (hit_name, hit_description, hit_length, hit_taxon, hit_db)
                                VALUES (?,?,?,?,?)}
    );
	my $embl_parser = Bio::EnsEMBL::Pipeline::Tools::Embl->new();
	my $server = $self->get_server();
	print  $server "-F " . join(" ", @ids) . "\n";
	local $/ = "//\n";
	while(<$server>){
		#print "**************************************\n";
		next unless $_ !~ /no match/;
		$embl_parser->parse($_);
		my @hit_row = (
		               $embl_parser->sequence_version || shift(@{$embl_parser->accession}),
		               $embl_parser->description,
		               $embl_parser->seq_length,
		               $embl_parser->taxon,
		               $embl_parser->which_database
		);
		#print "Name                 :" . $hit_row[0] . "\n";
		#print "Description      : " . $hit_row[1] . "\n";
		#print "Sequence length: " . $hit_row[2] . "\n";
		#print "TaxonomyID         : " .$hit_row[3] . "\n";
		#print "Database           : " . $hit_row[4] . "\n";
		#print "**************************************\n";
		$sth->execute(@hit_row);
	}
	
}
1;
__END__

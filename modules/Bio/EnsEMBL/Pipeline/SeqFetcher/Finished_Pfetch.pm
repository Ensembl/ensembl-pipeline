
package Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;

use strict;
use Bio::Root::RootI;
use Bio::DB::RandomAccessI;
use Bio::Seq;
use IO::Socket;

use Bio::EnsEMBL::Pipeline::Tools::Embl;
use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;

    my ( $server, $port, $options, $archive_port ) = $self->_rearrange(
        [ 'PFETCH_SERVER', 'PFETCH_PORT', 'OPTIONS', 'ARCHIVE_PORT' ], @args );

    $self->server($server || 'cbi2.internal.sanger.ac.uk');
    $self->port($port || 22100);
    $self->archive_port($archive_port || 23100);
    $self->options($options);

    return $self;
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

sub archive_port {
    my( $self, $archive_port ) = @_;
    
    if ($archive_port) {
        $self->{'_archive_port'} = $archive_port;
    }
    return $self->{'_archive_port'};
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

sub get_archive_server {
    my ($self) = @_;

    my $host = $self->server;
    my $port = $self->archive_port;

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
            if (defined $seq_string && $seq_string ne 'no match') {

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

    my $sth = $dbobj->prepare(qq{ 
        REPLACE INTO hit_description (hit_name
              , hit_description
              , hit_length
              , hit_taxon
              , hit_db)
        VALUES (?,TRIM(?),?,?,?)
        });

    my $count = 100;
    my( @failed_to_fetch );
    while (my @hundred_ids = splice(@ids, @ids > $count ? -$count : 0)) {
    
        printf STDERR "Pfetching %d EMBL sequences\n", scalar(@hundred_ids);
    
	my $embl_parser = Bio::EnsEMBL::Pipeline::Tools::Embl->new();
	my $server = $self->get_server();
	print  $server join(' ', '-F', @hundred_ids) . "\n";
        my %ids_to_fetch = map {$_, 1} @hundred_ids;
        {
	    local $/ = "//\n";  #"
	    while (<$server>) {
                s/^no match\n//mg;
                next unless $_;
                $embl_parser->parse($_);
                my $name = $embl_parser->sequence_version || $embl_parser->accession->[0];
                delete($ids_to_fetch{$name})
                    or die "Got identifier '$name' from EMBL entry, but it was not in list of sequences to fetch:\n@hundred_ids";
                $sth->execute(
                    $name,
                    $embl_parser->description,
                    $embl_parser->seq_length,
                    $embl_parser->taxon,
                    $embl_parser->which_database
                    );
	    }
        }
        
        if (my @missing = keys %ids_to_fetch) {
            # Get lengths from the archvie server
            my $server = $self->get_archive_server;
            print $server join(" ", '-l', @missing), "\n";
            my( %acc_sv_lengths );
            my $i = 0;
            while (<$server>) {
                chomp;
                unless ($_ eq 'no match') {
                    my $acc_sv = $missing[$i];
                    #warn "Found $_ for $acc_sv\n";
                    $acc_sv_lengths{$acc_sv} = $_;
                }
                $i++;
            }

            # Make a hash of accessions without the .SV for
            # everything we got a length for.
            my( %acc_missing );
            foreach my $acc_sv (keys %acc_sv_lengths) {
                my ($acc) = $acc_sv =~ /^(.+)\.\d+$/;
                unless ($acc) {
                    warn "Can't make accession from '$acc_sv'";
                    next;
                }
                $acc_missing{$acc} = $acc_sv;
            }
            
            # Now get the rest of the info from the latest version
            # the sequence from the live server.
            {
                my $server = $self->get_server;
                print $server join(" ", '-F', keys %acc_missing), "\n";
                local $/ = "//\n";  #"
                while (<$server>) {
                    s/^no match\n//mg;
                    next unless $_;
                    $embl_parser->parse($_);
                    my $acc = $embl_parser->accession->[0];
                    my $acc_sv = $acc_missing{$acc}
                        or die "Can't see ACCESSION.SV for ACCESSION '$acc'";
                    my $length = $acc_sv_lengths{$acc_sv};
                    delete($ids_to_fetch{$acc_sv})
                        or die "Got identifier '$acc_sv' from EMBL entry, but it was not in list of sequences to fetch:\n@hundred_ids";
                    $sth->execute(
                        $acc_sv,
                        $embl_parser->description,
                        $embl_parser->seq_length,
                        $embl_parser->taxon,
                        $embl_parser->which_database
                        );
	        }
            }
        }
        push(@failed_to_fetch, keys %ids_to_fetch);
    }
    return @failed_to_fetch;
}

1;

__END__

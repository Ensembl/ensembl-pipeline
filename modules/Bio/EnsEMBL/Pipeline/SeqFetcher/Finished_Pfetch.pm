
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
    my( $self, $dbobj, $id_list, $chunk_size ) = @_;

    $chunk_size ||= 200;

	my ($descriptions,$failed) = fetch_descriptions($self, $id_list, $chunk_size);

    my $sth = $self->prepare_hit_desc_sth($dbobj);
	for my $accession (keys %$descriptions) {
		$sth->execute(
			$accession,
			map { $descriptions->{$accession}{$_}; } ('description','length','taxonid','database')
		);
	}
}

sub min {
	my ($a,$b) = @_;
	return ($a<$b) ? $a : $b;
}

sub fetch_descriptions {
    my( $self, $id_list, $chunk_size ) = @_;
    
	my $descriptions = {};	# results are stored here

		# first pass fails when the sequence version has changed:
    my $failed_first_pass = [];
	my $start = 0;
	do {
		my $next = min($chunk_size, scalar(@$id_list));
        my $chunk = [@$id_list[$start..$next-1]];

		push @$failed_first_pass, @{ $self->fetch_descriptions_by_accession($chunk, $descriptions) };
		$start = $next;
	} while($start<@$id_list);

	my $failed_second_pass = [];
	$start = 0;
	do {
		my $next = min($chunk_size, scalar(@$failed_first_pass));
        my $chunk = [@$failed_first_pass[$start..$next-1]];

		my ($succeeded_arch_len, $failed_arch_len)
			= $self->fetch_lengths_from_archive($failed_first_pass, $descriptions);

		my $failed_desc = $self->fetch_descriptions_by_accession($succeeded_arch_len, $descriptions);

		push @$failed_second_pass, @$failed_arch_len, @$failed_desc;
		$start = $next;
	} while($start<@$failed_first_pass);

    return ($descriptions, $failed_second_pass);
}

sub fetch_descriptions_by_accession {
    my( $self, $id_list, $descriptions ) = @_;
    
    my $srv = $self->get_server;
    print $srv join(' ', '-F', @$id_list), "\n";

    # Set the input record separator to split on EMBL
    # entries, which end with "//\n" on a line.
    local ${/} = "//\n";  # ${/} is the same as $/, but doesn't mess up
                          # syntax highlighting in my editor!

    my $embl_parser = Bio::EnsEMBL::Pipeline::Tools::Embl->new;

	my $failed = [];

    for (my $i = 0; $i < @$id_list; $i++) {
        my $entry = <$srv>;

        # Each entry may have one or more "no match" lines at the start
        $i += $entry =~ s/^no match\n//mg;

        # The end of the loop
        last unless $entry;
        
        my $full_name = $id_list->[$i] or die "No id at '$i' in list:\n@$id_list";
        $embl_parser->parse($entry);

		my $name_without_version = (split('.',$full_name))[0];

		my $found = 0;
		NAMES: for my $one_of_names (@{ $embl_parser->accession }) {
			if ($name_without_version eq $one_of_names) {
				$found = 1;
				warn "Found '$name_without_version'";

					# NB: do not redefine any of these if already defined
				$descriptions->{$full_name}{description}||= $embl_parser->description;
				$descriptions->{$full_name}{length}		||= $embl_parser->seq_length;
				$descriptions->{$full_name}{taxonid}	||= $embl_parser->taxon;
				$descriptions->{$full_name}{database}	||= $embl_parser->which_database;

				last NAMES;
			}
		}
		if(!$found) {
			my $acc_sv = $embl_parser->sequence_version;
            warn "Expecting '$full_name' but got '$acc_sv'";
			push @$failed, $full_name;
        }
    }
    return $failed;
}

sub fetch_lengths_from_archive {
    my( $self, $id_list, $descriptions ) = @_;

	my $server = $self->get_archive_server;
	print $server join(" ", '-l', @$id_list), "\n";

	my $succeeded = [];
	my $failed    = [];
	my $i = 0;
	while (<$server>) {
		chomp;
		my $full_name = $id_list->[$i];
		if($_ ne 'no match') {
			$descriptions->{$full_name}{length}	= $_;
			push @$succeeded, $full_name;
		} else {
			warn "${full_name}'s length is not in the archive";
			push @$failed, $full_name;
		}
		$i++;
	}
	return ($succeeded,$failed);
}

sub prepare_hit_desc_sth {
    my( $self, $dbobj ) = @_;
    
    my $sth = $dbobj->prepare(qq{ 
        REPLACE INTO hit_description (hit_name
              , hit_description
              , hit_length
              , hit_taxon
              , hit_db)
        VALUES (?,TRIM(?),?,?,?)
        });
    return $self->{'_hit_desc_sth'} = $sth; # James, do we really need to save it?
}

sub get_hit_desc_sth {
    my( $self ) = @_;
    
    reuturn $self->{'_hit_desc_sth'};	# ------,,----- and restore it, if it's only used in 1 place?
}

sub OLD_write_descriptions {
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
    
        printf STDERR "Pfetching %d sequences\n", scalar(@hundred_ids);
    
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
                    or die "Got identifier '$name' from entry, but it was not in list of sequences to fetch:\n@hundred_ids";
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
                # May be missing more than one SV from the same accession!
                my $list = $acc_missing{$acc} ||= [];
                push(@$list, $acc_sv);
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
                    my $list = $acc_missing{$acc}
                        or die "Can't see ACCESSION.SV for ACCESSION '$acc'";
                    foreach my $acc_sv (@$list) {
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
        }
        push(@failed_to_fetch, keys %ids_to_fetch);
    }
    return @failed_to_fetch;
}

1;

__END__

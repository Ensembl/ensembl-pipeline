
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

	my ($descriptions,$failed) = $self->fetch_descriptions($id_list, $chunk_size);

    my $sth = $self->prepare_hit_desc_sth($dbobj);
	for my $accession (keys %$descriptions) {
		$sth->execute(
			$accession,
			map { $descriptions->{$accession}{$_}; } ('description','length','taxonid','database')
		);
	}
    return $failed;
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
    for (my $i = 0; $i < @$id_list; $i += $chunk_size) {
        my $j = $i + $chunk_size - 1;
        # Set second index to last element if we're off the end of the array
        $j = $#$id_list if $#$id_list < $j;
        # Take a slice from the array
        my $chunk = [@$id_list[$i..$j]];
        push @$failed_first_pass, @{ $self->fetch_descriptions_by_accession($chunk, $descriptions) };
    }
    
    my $failed_second_pass = [];
    for (my $i = 0; $i < @$failed_first_pass; $i += $chunk_size) {
        my $j = $i + $chunk_size - 1;
        $j = $#$failed_first_pass if $#$failed_first_pass < $j;

        my $chunk = [@$failed_first_pass[$i..$j]];
        my ($succeeded_arch_len, $failed_arch_len)
			= $self->fetch_lengths_from_archive($chunk, $descriptions);
        my $failed_desc = $self->fetch_descriptions_by_accession($succeeded_arch_len, $descriptions, 1);
        push @$failed_second_pass, @$failed_arch_len, @$failed_desc;
    }
    
    return ($descriptions, $failed_second_pass);
}

sub fetch_descriptions_by_accession {
    my( $self, $id_list, $descriptions, $trim_sv_flag ) = @_;
    
    my $srv = $self->get_server;
    warn "Fetching full entries for ", scalar(@$id_list), " identifiers\n";
    if ($trim_sv_flag) {
        my @req = @$id_list;
        foreach (@req) {
            s/\.\d+$//;
        }
        print $srv join(' ', '-F', @req), "\n";
    } else {
        print $srv join(' ', '-F', @$id_list), "\n";
    }

    # Set the input record separator to split on EMBL
    # entries, which end with "//\n" on a line.
    local ${/} = "//\n";  # ${/} is the same as $/, but doesn't mess up
                          # syntax highlighting in my editor!

    my $embl_parser = Bio::EnsEMBL::Pipeline::Tools::Embl->new;

	my $failed = [];
    for (my $i = 0; $i < @$id_list; $i++) {
        my $entry = <$srv>;

        # Each entry may have one or more "no match" lines at the start
        while ($entry =~ s/^no match\n//m) {
            push(@$failed, $id_list->[$i]);
            $i++;
        }

        # The end of the loop
        last unless $entry;

        my $full_name = $id_list->[$i] or die "No id at '$i' in list:\n@$id_list";
        $embl_parser->parse($entry);

        my $name_without_version = $full_name;
        $name_without_version =~ s/\.\d+$//;

	    my $found = 0;
	    NAMES: for my $one_of_names (@{ $embl_parser->accession }) {
		    if ($name_without_version eq $one_of_names) {
			    $found = 1;
			    # warn "Found '$name_without_version'";

			    # NB: do not redefine any of these if already defined
			    $descriptions->{$full_name}{description} ||= $embl_parser->description;
			    $descriptions->{$full_name}{length}      ||= $embl_parser->seq_length;
			    $descriptions->{$full_name}{taxonid}     ||= $embl_parser->taxon;
			    $descriptions->{$full_name}{database}    ||= $embl_parser->which_database;

			    last NAMES;
		    }
	    }
	    unless ($found) {
            my $all_accs = $embl_parser->accession;
            if (my $sv = $embl_parser->sequence_version) {
                unshift @$all_accs, $sv;
            }
            warn "Expecting '$full_name' but got (@$all_accs)";
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

1;

__END__

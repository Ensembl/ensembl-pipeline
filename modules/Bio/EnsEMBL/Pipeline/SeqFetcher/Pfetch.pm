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
							      -executable => $exe
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
package Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;

use strict;
use Bio::Root::RootI;
use Bio::DB::RandomAccessI;
use Bio::Seq;
use IO::Socket;

use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);

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

    if ( wantarray ) {

        # an array was passed - return array of Bio:Seq objects
        return @seq_list; 
        
    }
    else {

        # one acc was passed then return the first(and only) element of the array    
        return $seq_list[0];
    }
}

sub get_descriptions {
    my ( $self, @id_list ) = @_;
    
    unless (@id_list) {
        $self->throw("No accession input");
    }

    my $server = $self->get_server();
    print $server "-D @id_list\n";
    my (@desc_list);
    for ( my $i = 0 ; $i < @id_list ; $i++ ) {
        chomp( my $desc = <$server> );
        if ( $desc eq 'no match' ) {
            $desc_list[$i] = undef;
        }
        else {
            $desc_list[$i] = $desc;
        }
    }

   if ( wantarray ) {

        # an array was passed - return array of Bio:Seq objects
        return @desc_list; 
        
    }
    else {

        # one acc was passed then return the first(and only) element of the array    
        return $desc_list[0];
    }

}

sub get_lengths {
    my ( $self, @id_list ) = @_;

    unless (@id_list) {
        $self->throw("No accession input");
    }
    my $server = $self->get_server();
    print $server "-l @id_list\n";
    my (@length_list);
    for ( my $i = 0 ; $i < @id_list ; $i++ ) {
        chomp( my $length = <$server> );
        if ( $length eq 'no match' ) {
            $length_list[$i] = undef;
        }
        else {
            $length_list[$i] = $length;
        }
    }

   if ( wantarray ) {

        # an array was passed - return array of Bio:Seq objects
        return @length_list; 
        
    }
    else {

        # one acc was passed then return the first(and only) element of the array    
        return $length_list[0];
    }
}

sub write_descriptions {
    my ( $self, $dbobj, @ids,  ) = @_;

    my ( $hid, $desc );
    
    my @desc_line = $self->get_descriptions(@ids);
    my @lengths   = $self->get_lengths(@ids);

    if ( scalar(@lengths) != scalar(@ids) ) {
        die qq{scalar(@ids) elements in id_list\t scalar(@lengths) elements in Lengths};
    }
    if ( scalar(@desc_line) != scalar(@ids) ) {
        die qq{scalar(@ids) elements in id_list\t scalar(@desc_line) elements in Descriptions};
    }

    for ( my $i = 0 ; $i < @ids ; $i++ ) {

        # parse description from dbEST
        if ( $desc_line[$i] =~ /^[^\|]+\|[^\|]+\|[^\|]+\|(\w+\.\d+)\|\w+\s+(.*)/ ) {
            ( $hid, $desc ) = ( $1, $2 );
            if ( $hid ne $ids[$i] ) {
                warn qq{Hid : $hid    does not match    Query_hid : $i};
                next;
            }
            $desc =~ s/\n//g;
            $desc =~ s/"/\\"/g;
            $desc =~ s/'/\\'/g;
            my $hit_db_prefix = 'Em:';
            my $l             = $lengths[$i];

            my $sth = $dbobj->prepare( qq{ 
                                INSERT DELAYED INTO 
                                hit_descriptions_test.hit_description
                                VALUES ('$hid','$hit_db_prefix','$l','$desc')}
            );

            $sth->execute;
            next;

        }

        # parse description from vertrna, dbSTS and dbGSS
        elsif ( $desc_line[$i] =~ /^\S+\s+(\S+\.\d+)(.*)/ ) {
            ( $hid, $desc ) = ( $1, $2 );
            if ( $hid ne $ids[$i] ) {
                warn qq{Hid : $hid    does not match    Query_hid : $i};
                next;
            }

            $desc =~ s/\n//g;
            $desc =~ s/"/\\"/g;
            $desc =~ s/'/\\'/g;
            my $hit_db_prefix = 'Em:';
            
            my $l             = $lengths[$i];

            my $sth = $dbobj->prepare( qq{ 
                                INSERT DELAYED INTO 
                                hit_descriptions_test.hit_description
                                VALUES ('$hid','$hit_db_prefix','$l','$desc')}
            );

            $sth->execute;
            next;

        }

        # parse description from Swall
        elsif ( $desc_line[$i] =~ /Desc:\s(.*)/ ) {

            $desc = $1;

            # is the id present anywhere in the swall description line? 
            if ( grep /$ids[$i]/, $desc_line[$i] ) {
                $hid = $ids[$i];
            }
            else {
                next;
            }

            # next decide whether this swall hit is Swissprot or Trembl.
            # logic in operation here is that only swall ids will 
            # ever have a '_'. So grep fro '_' in ids leading up to 
            # the descripton line.
            my ($pre_desc) = $desc_line[$i] =~ /(^.*)\sDesc:/;
            my $hit_db_prefix;
            if ( grep /_/, $pre_desc ) {

                # swall hit 
                $hit_db_prefix = 'Sw:';
            }
            else {

                # trembl hit 
                $hit_db_prefix = 'Tr:';
            }

            $desc =~ s/\n//g;
            $desc =~ s/"/\\"/g;
            $desc =~ s/'/\\'/g;

            my $l = $lengths[$i];

            my $sth =
              $dbobj->prepare(qq{ 
                            INSERT DELAYED INTO 
                            hit_descriptions_test.hit_description
                            VALUES ('$hid','$hit_db_prefix','$l','$desc')}
            );

            $sth->execute;
            next;
        }

    }

}

1;

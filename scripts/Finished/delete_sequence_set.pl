#!/usr/local/bin/perl -w

=head1 NAME

delete_sequence_set.pl

=head1 SYNOPSIS

delete_sequence_set.pl -host otterpipe2 -port 3303 -dbname pipe_zebrafish -user otta -pass **** -show

delete_sequence_set.pl -host otterpipe2 -port 3303 -dbname pipe_zebrafish -user otta -pass **** -delete -set chr1

=head1 DESCRIPTION

This script is used to either delete sequence sets or show the list of loaded sequence sets from
a new schema database. The Connection parameters (login, password and port) are
retrieved from the ~/.netrc file if it exists or provided as command line arguments.

=head1 OPTIONS

    -help         - display this pod

    -host       - (default: otterpipe1) host name for database (gets put as host= in locator)
    -port       - (check the ~/.netrc file) For RDBs, what port to connect to (port= in locator)
    -dbname       - For RDBs, what name to connect to (p_name= in locator)
    -user       - (check the ~/.netrc file) For RDBs, what username to connect as (user= in locator)
    -pass       - (check the ~/.netrc file) For RDBs, what password to use (pass= in locator)

    -cs_version	- the coordinate system version of the chromosome to use for deletion (default Otter)
    -delete       - delete the set
     -set|name     - set the set name to delete (could be used several times)
     or
    -show         - show the list of sequence sets present in the pipeline database

=head1 AUTHOR

Mustapha Larbaoui, B<email> ml6@sanger.ac.uk

=cut

use strict;
use Getopt::Long;
use Net::Netrc;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host = 'otterpipe1';
my $user;
my $pass;
my $port;
my $dbname;
my $delete;
my @sets;
my $cs_version = 'Otter';
my $show;
my $help;

my $table_header = '+' . ( '-' x 4 ) . '+' . ( '-' x 7 ) . '+' . ( '-' x 20 ) . '+' . ( '-' x 90 ) . '+' .( '-' x 12 ) . '+' .( '-' x 6 ) . '+' . "\n";
my $line_format  = "|%-4s|%-7s|%-20s|%-90s|%-12s|%-6s|\n";

&GetOptions(
	'host:s'   => \$host,
	'port:n'   => \$port,
	'user:s'   => \$user,
	'pass:s'   => \$pass,
	'dbname:s'   => \$dbname,
	'delete!'    => \$delete,
	'show!'      => \$show,
	'set|name:s' => \@sets,
	'cs_version:s' => \$cs_version,
	'h|help'     => \$help
);

if ($help) {
	exec( 'perldoc', $0 );
}

if ( !$user || !$pass || !$port ) {
	my $ref = Net::Netrc->lookup($host);
	throw(
		"~/.netrc file unavailable;
			need to provide missing parameter:
			user [$user]; password [$pass]; port [$port]"
	  )
	  unless ($ref);
	$user = $ref->login    unless $user;
	$pass = $ref->password unless $pass;
	$port = $ref->account  unless $port;
	throw(
		"Missing parameter in the ~/.netrc file:\n
			machine " .  ( $host || 'missing' ) . "\n
			login " .    ( $user || 'missing' ) . "\n
			password " . ( $pass || 'missing' ) . "\n
			account "
		  . ( $port || 'missing' ) . " (should be used to set the port number)"
	  )
	  unless ( $user && $pass && $port );
}

if ( !$dbname ) {
	throw("You must specify a database name (-dbname)");
}

my $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(
	-host   => $host,
	-user   => $user,
	-pass   => $pass,
	-port   => $port,
	-dbname => $dbname
);

if ($show) {
	show_sequence_sets($db);
}

if ( $delete && !@sets ) {
	throw("You must specify at least one sequence set to delete");
}
elsif ( !$delete && @sets ) {
	throw("What am I supposed to do with @sets ? You must specify delete");
}
elsif ( $delete && @sets ) {
	delete_sequence_sets( $db, \@sets );
}

if ( !$delete && !$show ) {
	throw("You must specify a command, either show or delete");
}

#
# Methods
#

sub show_sequence_sets {
	my ( $db, $s_sets ) = @_;
	print STDOUT $table_header;
	print STDOUT
	  sprintf( $line_format, 'chr', 'version', 'sequence_set', 'description','write_access','hidden' );
	print STDOUT $table_header;
	my $sets = $s_sets || get_sequence_sets($db);
	foreach my $sr_id ( sort keys %$sets ) {
		my $set  = $sets->{$sr_id}->[0];
		my $version  = $sets->{$sr_id}->[1];
		my $chr  = $sets->{$sr_id}->[2];
		my $desc = $sets->{$sr_id}->[3];
		my $wa = defined($sets->{$sr_id}->[4]) ? $sets->{$sr_id}->[4] : 'null';
		my $hidden= defined($sets->{$sr_id}->[5]) ? $sets->{$sr_id}->[5] : 'null';;
		print STDOUT sprintf( $line_format, $chr, $version, $set, $desc, $wa , $hidden );
	}
	print STDOUT $table_header;
}

sub delete_sequence_sets {
	my ( $db, $del_sets ) = @_;

	my $sql = qq{
		DELETE a,s,at
		FROM assembly a, seq_region s,  coord_system c, seq_region_attrib at
		WHERE s.name = ?
		AND ( a.asm_seq_region_id = s.seq_region_id
		OR a.cmp_seq_region_id = s.seq_region_id )
		AND at.seq_region_id = s.seq_region_id
		AND c.version = '$cs_version'
		AND c.name = 'chromosome'
		AND c.coord_system_id = s.coord_system_id
	};
	my $sth     = $db->prepare($sql);

	foreach my $set (@$del_sets) {
		my $nb = $sth->execute($set);
		if($nb >= 1 ) {
			print STDOUT "Sequence set $set on chromosome:$cs_version deleted from $dbname\n";
		} else {
			print STDOUT "There is no sequence set $set on chromosome:$cs_version in $dbname\n";
		}
	}

	show_sequence_sets($db);
}

sub get_sequence_sets {
	my ($db) = @_;
	my $hash;
	my $query = qq{
				SELECT s.seq_region_id,c.version,
				            a1.value AS chromosome,
				            s.name AS seq_set,
				            a2.value AS description,
				            a3.value AS write_access,
				            a4.value AS hidden
				FROM attrib_type t1, attrib_type t2,attrib_type t3, attrib_type t4,coord_system c, seq_region s
				LEFT JOIN seq_region_attrib a1 ON (a1.seq_region_id =  s.seq_region_id AND t1.attrib_type_id = a1.attrib_type_id)
				LEFT JOIN seq_region_attrib a2 ON (a2.seq_region_id =  s.seq_region_id AND t2.attrib_type_id = a2.attrib_type_id)
				LEFT JOIN seq_region_attrib a3 ON (a3.seq_region_id =  s.seq_region_id AND t3.attrib_type_id = a3.attrib_type_id)
				LEFT JOIN seq_region_attrib a4 ON (a4.seq_region_id =  s.seq_region_id AND t4.attrib_type_id = a4.attrib_type_id)
				WHERE c.name = 'chromosome'
				AND c.coord_system_id = s.coord_system_id
				AND t1.code = 'chr'
				AND t2.code = 'description'
				AND t3.code = 'write_access'
				AND t4.code = 'hidden'
				ORDER BY  s.seq_region_id DESC
					};
	my $sth = $db->prepare($query);
	$sth->execute();
	while ( my ( $sr_id, $version, $chr, $set, $desc, $wa, $hidden ) = $sth->fetchrow ) {
		$hash->{$sr_id} = [ $set, $version, $chr, $desc, $wa, $hidden ];
	}

	return $hash;
}

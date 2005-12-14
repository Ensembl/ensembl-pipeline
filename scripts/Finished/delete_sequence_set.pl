#!/usr/local/bin/perl -w

=head1 NAME 

delete_sequence_set.pl

=head1 SYNOPSIS

delete_sequence_set.pl -p_host otterpipe2 -p_port 3303 -p_name pipe_zebrafish -p_user otta -p_pass **** -show

delete_sequence_set.pl -p_host otterpipe2 -p_port 3303 -p_name pipe_zebrafish -p_user otta -p_pass **** -delete -set chr1

=head1 DESCRIPTION

This script is used to either B<delete> sequence sets or B<show> the list of loaded sequence sets from 
a new B<pipeline> database. The Connection parameters (login, password and port) are 
retrieved from the ~/.netrc file if it exists or provided as command line arguments.

=head1 OPTIONS

    -help         - display this pod
    
    -p_host       - (default: otterpipe1) host name for database (gets put as host= in locator)
    -p_port       - (check the ~/.netrc file) For RDBs, what port to connect to (port= in locator) 
    -p_name       - For RDBs, what name to connect to (p_name= in locator)
    -p_user       - (check the ~/.netrc file) For RDBs, what username to connect as (user= in locator)
    -p_pass       - (check the ~/.netrc file) For RDBs, what password to use (pass= in locator)
    
    -delete       - delete the B<set>
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
my $p_name;
my $delete;
my @sets;
my $show;
my $help;

&GetOptions(
	'p_host:s'   => \$host,
	'p_port:n'   => \$port,
	'p_user:s'   => \$user,
	'p_pass:s'   => \$pass,
	'p_name:s'   => \$p_name,
	'delete!'    => \$delete,
	'show!'      => \$show,
	'set|name:s' => \@sets,
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

if ( !$p_name ) {
	throw("You must specify a database name (-p_name)");
}

my $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(
	-host   => $host,
	-user   => $user,
	-pass   => $pass,
	-port   => $port,
	-dbname => $p_name
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
	my $line =
	  '+' . ( '-' x 9 ) . '+' . ( '-' x 29 ) . '+' . ( '-' x 38 ) . '+' . "\n";
	print STDOUT $line;
	print STDOUT
	  sprintf( "|%-9s|%-29s|%-38s|\n", 'chr', 'sequence_set', 'description' );
	print STDOUT $line;
	my $sets = $s_sets || get_sequence_sets($db);
	foreach my $set ( keys %$sets ) {
		my $chr  = $sets->{$set}->[0];
		my $desc = $sets->{$set}->[1];
		print STDOUT sprintf( "|%-9s|%-29s|%-38s|\n", $chr, $set, $desc );
	}
	print STDOUT $line;
}

sub delete_sequence_sets {
	my ( $db, $del_sets ) = @_;
	my $stored_sets = get_sequence_sets($db);
	foreach my $set (@$del_sets) {
		if ( !$stored_sets->{$set} ) {
			show_sequence_sets( $db, $stored_sets );
			throw("Sequence set [$set] is not in pipeline database");
		}
	}
	my $del_attrib = "DELETE FROM seq_region_attrib
					    WHERE seq_region_id IN (
					  		SELECT seq_region_id 
							FROM seq_region
							WHERE name = ? )";
	my $att_sth      = $db->prepare($del_attrib);
	my $del_assembly = "DELETE FROM assembly
						WHERE asm_seq_region_id IN (
							SELECT seq_region_id 
							FROM seq_region
							WHERE name = ? )";
	my $asm_sth     = $db->prepare($del_assembly);
	my $del_seq_reg = "DELETE FROM seq_region
				    	WHERE name = ?";
	my $seq_sth = $db->prepare($del_seq_reg);

	foreach my $set (@$del_sets) {
		$att_sth->execute($set);
		$asm_sth->execute($set);
		$seq_sth->execute($set);
		print STDOUT "Set $set deleted from pipeline\n";
		show_sequence_sets($db);
	}
}

sub get_sequence_sets {
	my ($db) = @_;
	my $hash;
	my $query =
	  "SELECT a1.value AS chromosome, s.name AS seq_set, a2.value AS description
	 FROM seq_region_attrib a1, seq_region_attrib a2, 
		  attrib_type t1, attrib_type t2, seq_region s, coord_system c
	 WHERE c.name = 'chromosome'
	 AND c.version = 'Otter'
	 AND c.coord_system_id = s.coord_system_id
	 AND t1.code = 'chr'
	 AND t2.code = 'desc'
	 AND t1.attrib_type_id = a1.attrib_type_id 
	 AND t2.attrib_type_id = a2.attrib_type_id 
	 AND a1.seq_region_id =  s.seq_region_id
	 AND a2.seq_region_id =  s.seq_region_id
	 ORDER BY chromosome";
	my $sth = $db->prepare($query);
	$sth->execute();
	while ( my ( $chr, $set, $desc ) = $sth->fetchrow ) {
		$hash->{$set} = [ $chr, $desc ];
	}

	return $hash;
}

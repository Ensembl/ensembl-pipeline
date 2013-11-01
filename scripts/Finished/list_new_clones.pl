#!/usr/bin/env perl

=pod

=head1 NAME
	list_new_clones.pl

=head1 DESCRIPTION
	compare current (visible & writable) and previous chromosome assemblies and return the list of new clones

=head1 OPTION

     Database options

    -host (default:otterlive)   host name of the target database
    -dbname (no default)  For RDBs, what database to connect to
    -user (check the ~/.netrc file)  For DBAs, what username to connect as
    -pass (check the ~/.netrc file)  For DBAs, what password to use
    -port (check the ~/.netrc file)   For DBAs, what port to use

     Other options

    -help prints out the perl docs

=head1 EXAMPLES
	list_new_clones.pl -dbname loutre_human

=head1 AUTHOR

	ml6@sanger.ac.uk

=cut


use warnings ;
use strict;
use Net::Netrc;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $host = 'otterlive';
my $user;
my $pass;
my $port;
my $name;
my $useage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
            'host=s'      => \$host,
            'dbname=s'      => \$name,
            'user=s'      => \$user,
            'pass=s'      => \$pass,
            'port=s'      => \$port,
            'h|help!'       => $useage,
	    ) or $useage->();

if ( !$user || !$pass || !$port ) {
	my @param = &get_db_param($host);
	$user = $param[0] unless $user;
	$pass = $param[1] unless $pass;
	$port = $param[2] unless $port;
}

my $loutre_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-user   => $user,
	-dbname => $name,
	-host   => $host,
	-port   => $port,
	-pass   => $pass
);

my $dbc = $loutre_dba->dbc;

my $current_chr_list = $dbc->prepare(
qq{
    SELECT s.name
    FROM (attrib_type t1
          , attrib_type t2
          , attrib_type t3
          , attrib_type t4
          , attrib_type t5
          , coord_system c
          , seq_region s)
    LEFT JOIN seq_region_attrib a1
      ON (a1.seq_region_id = s.seq_region_id
          AND t1.attrib_type_id = a1.attrib_type_id)
    LEFT JOIN seq_region_attrib a2
      ON (a2.seq_region_id = s.seq_region_id
          AND t2.attrib_type_id = a2.attrib_type_id)
    LEFT JOIN seq_region_attrib a3
      ON (a3.seq_region_id = s.seq_region_id
          AND t3.attrib_type_id = a3.attrib_type_id)
    LEFT JOIN seq_region_attrib a4
      ON (a4.seq_region_id = s.seq_region_id
          AND t4.attrib_type_id = a4.attrib_type_id)
    LEFT JOIN seq_region_attrib a5
      ON (a5.seq_region_id = s.seq_region_id
          AND t5.attrib_type_id = a5.attrib_type_id)
    WHERE c.name IN ('chromosome','subregion')
      AND c.coord_system_id = s.coord_system_id
      AND t1.code = 'chr'
      AND t2.code = 'description'
      AND t3.code = 'write_access'
      AND t4.code = 'hidden'
      AND t5.code = 'equiv_asm'
      AND a4.value = 0
      AND a3.value = 1
      AND s.name LIKE 'chr%'
    ORDER BY s.name
} );

my $chr_clone_list = $dbc->prepare(
qq{
    SELECT DISTINCT(sc.name)
    FROM seq_region sa
      , seq_region sc
      , assembly a
      , coord_system ca
      , coord_system cc
    WHERE sa.seq_region_id = a.asm_seq_region_id
      AND sc.seq_region_id = a.cmp_seq_region_id
      AND ca.name = 'chromosome'
      AND ca.coord_system_id = sa.coord_system_id
      AND cc.name = 'contig'
      AND cc.coord_system_id = sc.coord_system_id
      AND sa.name = ?
} );

$current_chr_list->execute;

while(my $row = $current_chr_list->fetchrow_arrayref) {
	my ($current_chr) = @$row;
	my ($chr,$version) = $current_chr =~ /chr(\w+)-(\d+)/;
	$version--;
	my $previous_chr = sprintf("chr%s-%.2d",$chr,$version);
	print STDOUT "compare $current_chr to $previous_chr\n";

	my $current_hash = &get_clone_hash($current_chr);
	my $previous_hash = &get_clone_hash($previous_chr);

	foreach my $new_clone (keys %$current_hash) {
		print "\t$new_clone\n" unless $previous_hash->{$new_clone};
	}
}

sub get_clone_hash {
	my ($seq_name) = @_;
	my $hash;
	$chr_clone_list->execute($seq_name);
	while(my $r = $chr_clone_list->fetch) {
		my ($clone) = @$r;
		$hash->{$clone} = 1;
	}

	return $hash;
}


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


#!/software/bin/perl -w

=head1 NAME

make_subregion_agp.pl

=head1 SYNOPSIS

make_subregion_agp.pl

=head1 DESCRIPTION

This script is used to produce an agp for a chromosome subregion defined with
the chromosome name, the first accession and last accession in the subergion.
If the login, password and port parameters of the loutre connexion are not provided, they will be
recovered from the ~/.netrc file. See the Net::Netrc module for more details.


=head1 OPTIONS

    -host   (default:otterlive)        host name for the loutre database (gets put as phost= in locator)
    -dbname (no default)               for RDBs, what name to connect to (pname= in locator)
    -user   (check the ~/.netrc file)  for RDBs, what username to connect as (puser= in locator)
    -pass   (check the ~/.netrc file)  for RDBs, what password to use (ppass= in locator)
    -port   (check the ~/.netrc file)  for RDBs, what port to use (pport= in locator)

    -cs_name (default:chromosome) the name of the coordinate system being stored
    -cs_version (default:Otter) the version of the chromosome coordinate system being stored
    -set	 the sequence set name
    -first	first accession in the subregion
    -last	last accession in the subregion
    -chr	(optional) the chromosome name to be put in the agp
    -help|h      displays this documentation with PERLDOC


=head1 DEPENDENCE

regions_to_agp needs to be in the PATH

=head1 EXAMPLES

Here is a command line example:

    ./make_subregion_agp.pl \
	-dbname loutre_zebrafish \
	-set	chr1_20090930 \
	-first	FP016183 \
	-last	FP015984 \
	-chr	1

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use strict;

use Getopt::Long;
use Net::Netrc;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# loutre connexion parameters, default values.
my $host = 'otterlive';
my $port = '';
my $name = '';
my $user = '';
my $pass = '';

my $cs_name = 'chromosome';
my $cs_version = 'Otter';
my $chromosome;
my $set;
my $first;
my $last;

my $usage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
	'host:s'                => \$host,
	'port:n'                => \$port,
	'dbname:s'              => \$name,
	'user:s'                => \$user,
	'pass:s'                => \$pass,
	'cs_name:s' 			=> \$cs_name,
	'cs_version:s'			=> \$cs_version,
	'set=s'                 => \$set,
	'first=s'				=> \$first,
	'last=s'				=> \$last,
	'chr=s'			=> \$chromosome,
	'h|help!'               => $usage
  )
  or $usage->();

throw("No sequence set name given") unless ($set);

throw("Must provide first and last accessions of the subregion") unless $first && $last;

throw("Can't connect to a database without a dbname") unless $name;

if ( !$user || !$pass || !$port ) {
	my @param = &get_db_param($host);
	$user = $param[0] unless $user;
	$pass = $param[1] unless $pass;
	$port = $param[2] unless $port;
}

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
		-user   => $user,
		-dbname => $name,
		-host   => $host,
		-port   => $port,
		-pass   => $pass);

my $agp_sql = qq{
	SELECT	chr_ctg.asm_start,
			chr_ctg.asm_end,
			substring_index(ctg.name,".",2),
			chr_ctg.cmp_start,
			chr_ctg.cmp_end,
			chr_ctg.ori
	FROM seq_region chr, seq_region ctg, assembly chr_ctg, coord_system cs_chr, coord_system cs_ctg
	WHERE chr.seq_region_id = chr_ctg.asm_seq_region_id
	AND ctg.seq_region_id = chr_ctg.cmp_seq_region_id
	AND cs_chr.name = 'chromosome'
	AND cs_chr.coord_system_id = chr.coord_system_id
	AND cs_ctg.name = 'contig'
	AND cs_ctg.coord_system_id = ctg.coord_system_id
	AND chr.name = ?
	ORDER BY chr_ctg.asm_start ASC
};

my $sth = $dba->dbc->prepare($agp_sql);

$sth->execute($set);

my ($region,$s_chr_end,$f,$l);

while(my ($chr_start, $chr_end, $clone, $clone_start, $clone_end, $ori) = $sth->fetchrow) {
	next unless $clone =~ /$first/ || $f;$f =1;

	if($s_chr_end && (my $gap = $chr_start-$s_chr_end-1)) {
		$region .= "N\t$gap\n";
	}
	$region .= "F\t$clone\t$clone_start\t$clone_end\t$ori\n";

	if($clone =~ /$last/) { $l = 1 ; last }
}

throw("First accession $first not found in chromosome $set") unless ($f);
throw("Last accession $last not found in chromosome $set") unless ($l);

open my $agp, "| regions_to_agp ".($chromosome ? "-chr $chromosome " : "") or die "Can't open pipe to regions_to_agp; $!";
print $agp $region;
close $agp or die "Error running regions_to_agp; exit($?)";

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



#!/usr/bin/env perl

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
    -first	first accession in the subregion (regexp)
    -last	last accession in the subregion (regexp)
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


To find out whether the contigs from clones in the assembly are sorted
or jumbled, try

  make_subregion_agp.pl -host otterlive -port 3324 -user ottadmin -pass mumble \
   -dbname loutre_pig -set chrX-05 -chr pigX5 -first FP236733 -last CU928280 \
    |cut -f6 | cut -d. -f1 | uniq | sort | uniq -c | sort -rn
  # Take the accession names and uniq them twice.
  # When each clone is used in one place it has prefix "1" in output.

To count the contigs per clone just remove the first uniq:

  make_subregion_agp.pl -host otterlive -port 3324 -user ottadmin -pass mumble \
   -dbname loutre_pig -set chrX-05 -chr pigX5 -first FP236733 -last CU928280 \
    |cut -f6 | cut -d. -f1 | sort | uniq -c | sort -rn


=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use warnings ;
use strict;

use Getopt::Long;
use Net::Netrc;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


die <<B0RK;
Instead of this script, use
  ensembl-otter/scripts/lace/show_sequence_set -agp -dataset pig -set chrX-05


Problems with $0,

  It cannot output gaps.

  The "-last" clone will be truncated to one contig.

  It assumes the contig seq_region.name is accession.version .

 The feature this script has which show_sequence_set lacks is the
 -first / -last limit.  If this is needed, to correctly output the
 whole of the last clone it will be necessary either to peek the
 assembly for (max(cmp_end) where clonename = ?) OR cut the tail after
 reading the whole assembly.

B0RK

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

# If -first or -last are not provided, we should read the chromosome
# anyway and tell what we saw
foreach ($first, $last) { $_ = "^notprovided" unless defined $_ }

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
my ($input_first, $input_last); # currently just for diagnostics

while(my ($chr_start, $chr_end, $clone, $clone_start, $clone_end, $ori) = $sth->fetchrow) {

	$input_first = $clone unless defined $input_first;
	$input_last = $clone;

	next unless $clone =~ /$first/ || $f;$f =1;

	if($s_chr_end && (my $gap = $chr_start-$s_chr_end-1)) {
		$region .= "N\t$gap\n";
	}
	$region .= "F\t$clone\t$clone_start\t$clone_end\t$ori\n";

	if($clone =~ /$last/) { $l = 1 ; last }
}

throw("First accession $first not found in chromosome $set\n (input contained first=$input_first, last=$input_last)") unless ($f);
throw("Last accession $last not found in chromosome $set\n (input contained first=$input_first, last=$input_last)") unless ($l);

print $region; die;
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



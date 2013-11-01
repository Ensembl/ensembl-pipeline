#!/usr/bin/env perl

=head1 NAME

create_agp_from_scaffold_xtrop.pl

=head1 SYNOPSIS

create_agp_from_scaffold_xtrop.pl

=head1 DESCRIPTION

This script is used to refactor xtrop scaffold into clone sized pieces (150Kb).
It will generate an agp of the new assembly into the chromosome coordinates,
the clones sequence files (fasta) and GFF files for the internal gaps.

here is an example commandline

./create_agp_from_scaffold.pl \
-seq_region_name chr1.47 \
-host ens-livemirror \
-port 3306 \
-name bos_taurus_core_43_3 \
-user ensadmin \
-pass *****	\
-outdir ./tmp/

=head1 OPTIONS

    -host (default:ens-livemirror)   host name for the database (gets put as phost= in locator)
    -name (no default)  For RDBs, what name to connect to (pname= in locator)
    -user (check the ~/.netrc file)  For RDBs, what username to connect as (puser= in locator)
    -pass (check the ~/.netrc file)  For RDBs, what password to use (ppass= in locator)
    -port (check the ~/.netrc file)   For RDBs, what port to use (pport= in locator)

    -cs_name	(default:scaffold) coordinate system of the seq_region
    -seq_region_name	seq_region name
    -start	The start of the slice on the sequence region
    -end	The end of the slice on the sequence region
    -outdir		output directory
    -verbose 		display the chromosome -> scaffold -> contig coordinates mapping
    -help|h		displays this documentation with PERLDOC

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use warnings ;
use strict;

use Getopt::Long;
use Net::Netrc;
use GFF;
use GFF::GeneFeature;
use GFF::GeneFeatureSet;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# Ensembl database connection parameters
my $phost = 'ens-livemirror';
my $pport = '3306';
my $pname = '';
my $puser = '';
my $ppass = '';

my $coord_sytem = 'scaffold';
my $seq_region_name;
my $start;
my $end;
my $output = '';
my $verbose;

my $usage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
	'host:s'                => \$phost,
	'port:n'                => \$pport,
	'name=s'                => \$pname,
	'user:s'                => \$puser,
	'pass:s'                => \$ppass,
	'cs_name:s'				=> \$coord_sytem,
	'seq_region_name=s'     => \$seq_region_name,
	'start:s'				=> \$start,
	'end:s'					=> \$end,
	'outdir:s'				=> \$output,
	'verbose!'				=> \$verbose,
	'h|help!'               => $usage
  )
  or $usage->();

throw("No seq_region name given") unless ($seq_region_name);

$output .= '/' if($output && $output !~ /\/$/);

if ( !$puser || !$ppass || !$pport ) {
	my @param = &get_db_param($phost);
	$puser = $param[0] unless $puser;
	$ppass = $param[1] unless $ppass;
	$pport = $param[2] unless $pport;
}

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-user   => $puser,
	-dbname => $pname,
	-host   => $phost,
	-port   => $pport,
	-pass   => $ppass
);

my $cs_a                 = $dba->get_CoordSystemAdaptor();
my $slice_a              = $dba->get_SliceAdaptor();


my $seq_region_cs = $cs_a->fetch_by_name($coord_sytem);
my $slice = $slice_a->fetch_by_region($seq_region_cs->name,$seq_region_name,$start,$end);
my $slice_length = $slice->length;

my $loop = 1;
my $i = 1;
$start ||= 1;
$end   ||= $slice_length;
my $e;

my $agp_out;
open($agp_out, ">${output}${seq_region_name}.agp") or die "Can't create file $seq_region_name : $!";

while($loop)
{
	last unless($start < $end);
	$e = $start + 149999;
	if($e > $end) {
		$e = $end;
		$loop = 0;
	}
	print2file($slice->sub_Slice($start,$e),$i) if(($e-$start) > 0);
	$start += 150000;
	$i++;
}

sub print2file {
	my ($s,$i) = @_;
	my $seq_out;
	# print AGP row here
	printf $agp_out "chrU\t%d\t%d\t%d\tF\t%s\t1\t%d\t+\n",$s->start,$s->end,$i,"${seq_region_name}_${i}.1",$s->length;
	open($seq_out, ">${output}${seq_region_name}_${i}.1.seq") or die "Can't create file ${seq_region_name}_${i}.1.seq : $!";
	print $seq_out ">${seq_region_name}_${i}.1\n".$s->seq;

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
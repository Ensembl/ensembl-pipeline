#!/software/bin/perl -w

=head1 NAME

create_agp_from_scaffold.pl

=head1 SYNOPSIS

create_agp_from_scaffold.pl

=head1 DESCRIPTION

This script is used to refactor a scaffold which is a set of short sequences into clone sized pieces (150Kb).
( It can also be used to refactor any coord_system type ) It will generate an agp of the new assembly
into the chromosome coordinates, the clones sequence files (fasta) and GFF files for the internal gaps.

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

my $chromosome_cs = $cs_a->fetch_by_name( "chromosome" );;
my $contig_cs = $cs_a->fetch_by_name("contig");
my $seq_region_cs = $cs_a->fetch_by_name($coord_sytem);


my $slice = $slice_a->fetch_by_region($seq_region_cs->name,$seq_region_name,$start,$end);
my $chromosome;
my $contig_projection = $slice->project($contig_cs->name);

my $contigs = {};

foreach my $segment (@$contig_projection) {
      my $contig_slice = $segment->to_Slice();
      my $chr_projection = shift @{ $contig_slice->project($chromosome_cs->name) };
      my $chr = $chr_projection->to_Slice();
      $chromosome ||= $chr->seq_region_name();
      my $name = $contig_slice->seq_region_name();
      $contigs->{$name} = [$chr->start(),$chr->end()];
      print 'chr',$chr->seq_region_name(), ':', $chr->start(), '-',$chr->end(), ':', $chr->strand(),"\t(",$chr->length,")\t|\t",
      		$slice->seq_region_name(), ':', $segment->from_start(), '-',$segment->from_end(), "\t(",$segment->from_end()-$segment->from_start()+1,")\t|\t",
            $contig_slice->seq_region_name(), ':', $contig_slice->start(), '-',$contig_slice->end(), ':', $contig_slice->strand(), "\t(",$contig_slice->length,")\n"
      if $verbose;
}

my $clones = &merge_contigs_into_clone_1($contigs);

my $l;
my $agp_out;
open($agp_out, ">${output}${seq_region_name}.agp") or die "Can't create file $seq_region_name : $!";

foreach my $clone (sort {$clones->{$a}->{coord}->[0] <=> $clones->{$b}->{coord}->[0]} keys %$clones) {
	my $asm_start = $clones->{$clone}->{coord}->[0];
	my $asm_end = $clones->{$clone}->{coord}->[1];
	my $cmp_start =	1;
	my $cmp_end = $clones->{$clone}->{coord}->[2];
	$l++;

	my $seq_out;
	my $gap_out;

	open($seq_out, ">${output}${clone}.1.seq") or die "Can't create file ${clone}.1.seq : $!";
	my $clone_slice = $slice_a->fetch_by_region('chromosome',$chromosome,$asm_start,$asm_end);
	print $seq_out ">${clone}.1\n".$clone_slice->seq;

	open($gap_out, ">${output}${clone}.1.gff") or die "Can't create file ${clone}.1.gap : $!";
	my $gaps = GFF::GeneFeatureSet->new();

	print $agp_out "chr$chromosome\t$asm_start\t$asm_end\t$l\tF\t${clone}.1\t$cmp_start\t$cmp_end\t+\n";

	while(my $gap = shift @{$clones->{$clone}->{gaps}}) {
		my $gff = GFF::GeneFeature->new();
		$gff->seqname("${clone}.1");
		$gff->feature('Gap');
		$gff->source($pname);
		$gff->start($gap->[0]);
		$gff->end($gap->[1]);
		$gff->strand('+');
		$gaps->addGeneFeature($gff);
	}

	$gaps->dump($gap_out);
}

sub merge_contigs_into_clone_1 {
	my ($contigs) = @_;
	my $clones = {};
	my $contigs_count = scalar (keys %$contigs);
	my $count;
	my $clone_length;
	my $clone_start;
	my $clone_end;
	my $gaps;
	my $i;
	foreach my $contig (sort { $a->[0] <=> $b->[0] } values %$contigs) {
		$count++;
		my $contig_start   = $contig->[0];
		my $contig_end     = $contig->[1];
		my $contig_length  = $contig_end-$contig_start+1;
		my $gap_start;
		my $gap_end;
		my $gap_length;

		# save clone internal gaps
		if($clone_end && ($clone_end+1) != $contig_start) {
			$gap_length = $contig_start - $clone_end -1;
			$gap_start = $clone_length+1;
			$gap_end   = $gap_start + $gap_length -1;

			push @$gaps, [$gap_start,$gap_end];
		}

		$clone_start  ||= $contig_start;
		$clone_end      = $contig_end;
		$clone_length  += $contig_length;
		$clone_length  += $gap_length if $gap_length;

		# create clone if clone length >= 140Kb or last contig of the loop
		if($clone_length && ($clone_length>=140000 || $count == $contigs_count) ) {
			$i++;
			my $clone_name = $seq_region_name.'_'.$i;
			$clones->{$clone_name}->{coord} = [$clone_start,$clone_end,$clone_length];
			$clones->{$clone_name}->{gaps} = $gaps;
			$gaps = [];
			$clone_length = 0;
			$clone_start  = 0;
			$clone_end    = 0;
		}
	}

	return $clones;
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
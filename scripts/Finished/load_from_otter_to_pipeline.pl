#!/usr/local/ensembl/bin/perl -w

=head1 NAME

load_from_otter_to_pipeline.pl

=head1 SYNOPSIS

load_from_otter_to_pipeline.pl

=head1 DESCRIPTION

This script is used to load either one or several sequence sets from an organisme specific Otter database into the seq_region, seq_region_attrib, attrib_type, assembly and dna tables of the Pipeline database (new schema). 
Appropriate option like the coord_system version for the specific chromosome has to be given. The sequence that is loaded into the dna table is Pfetched.
If no values are provided for the database connexion (login, password and port), the ~/.netrc file will be checked. See Net::Netrc module for more details.

here is an example commandline

./load_from_otter_to_pipeline.pl -set chr11 -o_host ecs4 -o_port 3352 -o_name otter_human -o_user pipuser -o_pass ***** -p_host otterpipe2 -p_port 3352 -p_name pipe_human -p_user pipuser -p_pass *****


=head1 OPTIONS

    -o_host (default:humsrv1)   host name for the dataset specific Otter database (gets put as ohost= in locator)
    -o_name (default:otter_human)   For RDBs, what name to connect to (oname= in locator)
    -o_user (check the ~/.netrc file) For RDBs, what username to connect as (ouser= in locator)
    -o_pass (check the ~/.netrc file) For RDBs, what password to use (opass= in locator)
    -o_port (check the ~/.netrc file)   For RDBs, what port to use (oport= in locator)
    -p_host (default:otterpipe1)   host name for the Pipeline database (gets put as phost= in locator)
    -p_name (no default)  For RDBs, what name to connect to (pname= in locator)
    -p_user (check the ~/.netrc file)  For RDBs, what username to connect as (puser= in locator)
    -p_pass (check the ~/.netrc file)  For RDBs, what password to use (ppass= in locator)
    -p_port (check the ~/.netrc file)   For RDBs, what port to use (pport= in locator)

    -chromosome_cs_version (default:Otter) the version of the coordinate system being stored
    -set|chr the sequence set to load 
    -help|h      displays this documentation with PERLDOC

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use strict;

use Getopt::Long;
use Net::Netrc;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::DBSQL::SequenceAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# Otter and Pipeline connexion parameters, default values.
my $phost = 'otterpipe1';
my $pport = '';
my $pname = '';
my $puser = '';
my $ppass = '';
my $ohost = 'humsrv1';
my $oport = '';
my $oname = 'otter_human';
my $ouser = '';
my $opass = '';

my $chromosome_cs_version = 'Otter';
my @seq_sets;

my $usage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
	'p_host:s'                 => \$phost,
	'p_port:n'                 => \$pport,
	'p_name:s'                 => \$pname,
	'p_user:s'                 => \$puser,
	'p_pass:s'                 => \$ppass,
	'o_host:s'                 => \$ohost,
	'o_port:n'                 => \$oport,
	'o_name:s'                 => \$oname,
	'o_user:s'                 => \$ouser,
	'o_pass:s'                 => \$opass,
	'chromosome_cs_version:s' => \$chromosome_cs_version,
	'chr|set=s'               => \@seq_sets,
	'h|help!'                 => $usage
  )
  or $usage->();

if ( !$puser || !$ppass || !$pport ) {
	my $ref = Net::Netrc->lookup($phost);
	throw(
		"~/.netrc file unavailable; 
			need to provide missing parameter: 
			user [$puser]; password [$ppass]; port [$pport]"
	  )
	  unless ($ref);
	$puser = $ref->login    unless $puser;
	$ppass = $ref->password unless $ppass;
	$pport = $ref->account  unless $pport;
	throw(
		"Missing parameter in the ~/.netrc file:\n
			machine " .  ( $phost || 'missing' ) . "\n
			login " .    ( $puser || 'missing' ) . "\n
			password " . ( $ppass || 'missing' ) . "\n
			account "
		  . ( $pport || 'missing' ) . " (should be used to set the port number)"
	  )
	  unless ( $puser && $ppass && $pport );
}

if ( !$ouser || !$opass || !$oport ) {
	my $ref = Net::Netrc->lookup($ohost);
	throw(
		"~/.netrc file unavailable; 
			need to provide missing parameter: 
			user [$ouser]; password [$opass]; port [$oport]"
	  )
	  unless ($ref);
	$ouser = $ref->login    unless $ouser;
	$opass = $ref->password unless $opass;
	$oport = $ref->account  unless $oport;
	throw(
		"Missing parameter in the ~/.netrc file:\n
			machine " .  ( $ohost || 'missing' ) . "\n
			login " .    ( $ouser || 'missing' ) . "\n
			password " . ( $opass || 'missing' ) . "\n
			account "
		  . ( $oport || 'missing' ) . " (should be used to set the port number)"
	  )
	  unless ( $ouser && $opass && $oport );
}

if ( !$pname ) {
	print STDERR
	  "Can't load sequence set without a target pipeline database name\n";
	print STDERR "-p_host $phost -p_user $puser -p_pass $ppass\n";
	$usage->();
}

if ( !scalar(@seq_sets) ) {
	print STDERR "Need chr|set to be able to run\n";
	$usage->();
}

my $otter_dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
	-user   => $ouser,
	-dbname => $oname,
	-host   => $ohost,
	-port   => $oport,
	-pass   => $opass
);
my $pipe_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-user   => $puser,
	-dbname => $pname,
	-host   => $phost,
	-port   => $pport,
	-pass   => $ppass
);
my $pipe_dbc = $pipe_dba->dbc();

my %end_value;
my $contigs_hashref = {};
my $seqset_info     = {};

#
# Get the sequence data set from Otter database and store it in a Hashtable.
#
{
	print STDOUT
"Getting Sequence set @seq_sets from Otter database: $oname ($ohost:$oport)\n";
	my $contigs_sth = $otter_dbc->prepare(
		q{
        SELECT chr.name, a.chr_start, a.chr_end
          , a.contig_start, a.contig_end, a.contig_ori
          , c.embl_acc, c.embl_version, a.superctg_name
        FROM chromosome chr
          , assembly a
          , contig g
          , clone c
        WHERE chr.chromosome_id = a.chromosome_id
          AND a.contig_id = g.contig_id
          AND g.clone_id = c.clone_id
          AND a.type = ?
        ORDER BY a.chr_start
        }
	);
	my $description_sth = $otter_dbc->prepare(
		q{
		SELECT description
		FROM sequence_set
		WHERE assembly_type = ?
		}
	);
	SET:foreach my $sequence_set (@seq_sets) {
		my $contig_number = 0;
		my $chr_href;

		# fetch all the contigs for this sequence set
		$contigs_sth->execute($sequence_set);
		CONTIG:while (
			my (
				$chr_name,     $chr_start,  $chr_end,
				$contig_start, $contig_end, $strand,
				$acc,          $sv,         $superctg_name
			)
			= $contigs_sth->fetchrow
		  )
		{
			if ( !$end_value{$sequence_set} ) {
				$end_value{$sequence_set} = $chr_end;
			}
			else {
				if ( $chr_end > $end_value{$sequence_set} ) {
					$end_value{$sequence_set} = $chr_end;
				}
			}
			$contigs_hashref->{ $acc . $sv.$sequence_set } = [
				$chr_name,     $chr_start,  $chr_end,
				$contig_start, $contig_end, $strand,
				$acc,          $sv,         $sequence_set
			];
			$chr_href->{$chr_name} = 1;
			$contig_number++;

		}

		# fetch the sequence set description from the sequence_set table
		$description_sth->execute($sequence_set);
		my ($desc) = $description_sth->fetchrow;
		$seqset_info->{$sequence_set} = [ $desc, keys %$chr_href ];

		throw("No data for '$sequence_set' in $oname")
		  unless ($contig_number);
		print STDOUT $contig_number
		  . " contigs retrieved for $sequence_set sequence set\n";
	}
}

#
# Load the sequence data set into the Pipeline database
#
{
	print STDOUT
	  "Writing data into Pipeline database: $pname ($phost:$pport)\n";

	my %asm_seq_reg_id;

	my $attr_a  = $pipe_dba->get_AttributeAdaptor();
	my $cs_a    = $pipe_dba->get_CoordSystemAdaptor();
	my $slice_a = $pipe_dba->get_SliceAdaptor();
	my $seq_a   = $pipe_dba->get_SequenceAdaptor();

	my $chromosome_cs;
	my $clone_cs;
	my $contig_cs;

	eval {
		$chromosome_cs =
		  $cs_a->fetch_by_name( "chromosome", $chromosome_cs_version );
		$clone_cs  = $cs_a->fetch_by_name("clone");
		$contig_cs = $cs_a->fetch_by_name("contig");

	};
	if ($@) {
		throw(
"A coord_system matching the arguments does not exsist in the cord_system table, please ensure you have the right coord_system entry in the database [$@]\n"
		);
	}
	my $slice;

# insert chromosome(s) in seq_region table and his attributes in seq_region_attrib & attrib_type tables
	foreach my $name ( keys(%end_value) ) {
		my $endv = $end_value{$name};
		eval {
			$slice =
			  $slice_a->fetch_by_region( 'chromosome', $name, undef, undef,
				undef, $chromosome_cs_version );
		};
		if ($slice) {
			warn(
"Sequence set <$name> is already in pipeline database <$pname>\n"
			);
			throw(  "There is a difference in size for $name: stored ["
				  . $slice->length
				  . "] =! new ["
				  . $endv
				  . "]" )
			  unless ( $slice->length eq $endv );
			$asm_seq_reg_id{$name} = $slice->get_seq_region_id;
		}
		else {
			$slice = &make_slice( $name, 1, $endv, $endv, 1, $chromosome_cs );
			$asm_seq_reg_id{$name} = $slice_a->store($slice);
			$attr_a->store_on_Slice( $slice,
				&make_seq_set_attribute( $seqset_info->{$name} ) );
		}
	}

# insert clone & contig in seq_region, seq_region_attrib, dna and assembly tables
	my $insert_query = qq {
				INSERT INTO assembly 
				(asm_seq_region_id, cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori) 
				values 
				(?,?,?,?,?,?,?)};
	my $sth = $pipe_dbc->prepare($insert_query);
	while ( my ( $k, $v ) = each %$contigs_hashref ) {
		my @values       = @$v;
		my $chr_name     = $values[0];
		my $sequence_set = $values[8];
		my $chr_start    = $values[1];
		my $chr_end      = $values[2];
		my $ctg_start    = $values[3];
		my $ctg_end      = $values[4];
		my $ctg_ori      = $values[5];
		my $acc          = $values[6];
		my $ver          = $values[7];
		my $acc_ver      = $acc . "." . $ver;

		#if($acc_ver ne 'AC004775.1') { next ; }
		#if(!($acc_ver eq 'BX571961.11' || $acc_ver eq 'BX649427.7')) { next ; }

		my $clone;
		my $clone_seq_reg_id;
		my $contig;
		my $contig_name = $acc_ver . "." . "1" . ".";
		my $ctg_seq_reg_id;
		my $seqlen;
		eval {
			$clone  = $slice_a->fetch_by_region( 'clone', $acc_ver );
			$seqlen = $clone->length;
			$contig =
			  $slice_a->fetch_by_region( 'contig', $contig_name . $seqlen );
		};
		if ( $clone && $contig ) {
			warn "clone and contig <"
			  . $clone->name . " ; "
			  . $contig->name
			  . "> are already in the pipeline database\n";
			$clone_seq_reg_id = $clone->get_seq_region_id;
			$ctg_seq_reg_id   = $contig->get_seq_region_id;
		}
		elsif ( $clone && !$contig ) {
			### code to be added to create contigs related to the clone
			throw(  "Missing contig entry in the database for clone <"
				  . $clone->name
				  . ">\n" );
		}
		else {
			##fetch the dna sequence from pfetch server with acc_ver id
			my $seqobj ||= &pfetch_acc_sv($acc_ver);
			my $seq = $seqobj->seq;
			$seqlen = $seqobj->length;
			$contig_name .= $seqlen;

			##make clone and insert clone to seq_region table
			$slice = &make_slice( $acc_ver, 1, $seqlen, $seqlen, 1, $clone_cs );
			$clone_seq_reg_id = $slice_a->store($slice);
			if ( !defined($clone_seq_reg_id) ) {
				print
"clone seq_region_id has not been returned for the accession $acc_ver";
				exit;
			}
			##make attribute for clone and insert attribute to seq_region_attrib & attrib_type tables
			$attr_a->store_on_Slice( $slice,
				&make_clone_attribute( $acc, $ver ) );

			##make contig and insert contig, and associated dna sequence to seq_region & dna table
			$slice =
			  &make_slice( $contig_name, 1, $seqlen, $seqlen, 1, $contig_cs );
			$ctg_seq_reg_id = $slice_a->store( $slice, \$seq );
			if ( !defined($ctg_seq_reg_id) ) {
				print
"contig seq_region_id has not been returned for the contig $contig_name";
				exit;
			}
		}
		##insert chromosome to contig assembly data into assembly table
		$sth->execute( $asm_seq_reg_id{$sequence_set},
			$ctg_seq_reg_id, $chr_start, $chr_end, $ctg_start, $ctg_end,
			$ctg_ori );
		##insert clone to contig assembly data into assembly table
		$sth->execute( $clone_seq_reg_id, $ctg_seq_reg_id, 1, $seqlen, 1,
			$seqlen, $ctg_ori );

	}

}

#
# Methods part
#
sub make_clone_attribute {
	my ( $acc, $ver ) = @_;
	my @attrib;
	my $attrib = &make_attribute(
		'htgs_phase',                              'HTGS Phase',
		'High Throughput Genome Sequencing Phase', '3'
	);
	push @attrib, $attrib;
	push @attrib,
	  &make_attribute( 'intl_clone_name', 'International Clone Name',
		'', $acc . "." . $ver );
	push @attrib,
	  &make_attribute( 'embl_accession', 'EMBL Accession', '', $acc );
	push @attrib, &make_attribute( 'embl_version', 'EMBL Version', '', $ver );
	return \@attrib;
}

sub make_seq_set_attribute {
	my ($arr_ref) = @_;
	my $desc = shift @$arr_ref;
	my @attrib;
	push @attrib,
	  &make_attribute(
		'desc',
		'Assembly Description',
		'Assembly Description for  this Sequence Set', $desc
	  );
	foreach my $ch (@$arr_ref) {
		push @attrib,
		  &make_attribute(
			'chr',
			'Chromosome Name',
			'Chromosome Name Contained in the Assembly', $ch
		  );
	}

	return \@attrib;
}

sub make_attribute {
	my ( $code, $name, $description, $value ) = @_;
	my $attrib = Bio::EnsEMBL::Attribute->new(
		-CODE        => $code,
		-NAME        => $name,
		-DESCRIPTION => $description,
		-VALUE       => $value
	);
	return $attrib;
}

sub make_slice {

	my ( $name, $start, $end, $length, $strand, $coordinate_system ) = @_;
	my $slice = Bio::EnsEMBL::Slice->new(
		-seq_region_name   => $name,
		-start             => $start,
		-end               => $end,
		-seq_region_length => $length,
		-strand            => $strand,
		-coord_system      => $coordinate_system,
	);
	return $slice;
}

#
# Pfetch the sequences
#-------------------------------------------------------------------------------
#  If the sequence isn't available from the default pfetch
#  the archive pfetch server is tried.
#
{
	my ( $pfetch, $pfetch_archive );

	sub pfetch_acc_sv {
		my ($acc_ver) = @_;
		print "Fetching '$acc_ver'\n";
		$pfetch ||= Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
		$pfetch_archive ||=
		  Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new(
			-PFETCH_PORT => 23100, );
		my $seq = $pfetch->get_Seq_by_acc($acc_ver);
		unless ($seq) {
			$seq = $pfetch_archive->get_Seq_by_acc($acc_ver);
		}
		unless ($seq) {
			my $seq_file = "$acc_ver.seq";
			warn
			  "Attempting to read fasta file <$acc_ver.seq> in current dir.\n";
			my $in = Bio::SeqIO->new(
				-file   => $seq_file,
				-format => 'FASTA',
			);
			$seq = $in->next_seq;
			my $name = $seq->display_id;
			unless ( $name eq $acc_ver ) {
				die "Sequence in '$seq_file' is called '$name' not '$acc_ver'";
			}
		}
		return $seq;
	}
}

1;

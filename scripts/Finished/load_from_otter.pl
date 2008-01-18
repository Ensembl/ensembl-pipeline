#!/software/bin/perl -w

=head1 NAME

load_from_otter.pl

=head1 SYNOPSIS

load_from_otter.pl

=head1 DESCRIPTION

This script is used to load one or several sequence sets from an organisme
specific otter database into either a pipeline or a loutre database, both based
on ensembl schema version 20+. The loaded sequences are either Pfetched or fetched
from fasta files in the current directory. If the login, password and port parameters
are not provided, they will be recovered from the ~/.netrc file.
See the Net::Netrc module for more details.

here is an example commandline

./load_from_otter.pl
-set chr11-02
-o_host ecs4
-o_port 3352
-o_name otter_human
-o_user pipuser
-o_pass *****
-host otterpipe2
-port 3352
-name pipe_human
-user pipuser
-pass *****


=head1 OPTIONS

    -o_host (default:otterlive)   host name for the dataset specific otter database (gets put as ohost= in locator)
    -o_name (default:otter_human)   For RDBs, what name to connect to (oname= in locator)
    -o_user (check the ~/.netrc file) For RDBs, what username to connect as (ouser= in locator)
    -o_pass (check the ~/.netrc file) For RDBs, what password to use (opass= in locator)
    -o_port (check the ~/.netrc file)   For RDBs, what port to use (oport= in locator)
    -host (default:otterpipe1)   host name for the target database (gets put as phost= in locator)
    -name (no default)  For RDBs, what name to connect to (pname= in locator)
    -user (check the ~/.netrc file)  For RDBs, what username to connect as (puser= in locator)
    -pass (check the ~/.netrc file)  For RDBs, what password to use (ppass= in locator)
    -port (check the ~/.netrc file)   For RDBs, what port to use (pport= in locator)

    -cs_version (default:Otter) the version of the chromosome coordinate system being stored
    -set|chr	the sequence set to load
    -nosubmit	don't prime the pipeline with the SubmitContig analysis in case of pipeline
    		database
    -help|h	displays this documentation with PERLDOC

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use strict;

use Getopt::Long;
use Net::Netrc;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# otter and destination database connexion parameters, default values.
my $host = 'otterpipe1';
my $port = '';
my $name = '';
my $user = '';
my $pass = '';

my $ohost = 'otterlive';
my $oport = '';
my $oname = 'otter_human';
my $ouser = '';
my $opass = '';

my $cs_version = 'Otter';
my @seq_sets;
my $do_submit = 1;	  # Set if we want to prime the pipeline with the SubmitContig analysis

my $usage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
	'host:s'                 => \$host,
	'port:n'                 => \$port,
	'name=s'                 => \$name,
	'user:s'                 => \$user,
	'pass:s'                 => \$pass,
	'o_host:s'                 => \$ohost,
	'o_port:n'                 => \$oport,
	'o_name:s'                 => \$oname,
	'o_user:s'                 => \$ouser,
	'o_pass:s'                 => \$opass,
	'cs_version:s' 				=> \$cs_version,
	'chr|set=s'               => \@seq_sets,
	'no_submit!'			  => sub{ $do_submit = 0 },
	'submit!'			  	  => \$do_submit,
	'h|help!'                 => $usage
  )
  or $usage->();

if ( !$user || !$pass || !$port ) {
	my @param = &get_db_param($host);
	$user = $param[0] unless $user;
	$pass = $param[1] unless $pass;
	$port = $param[2] unless $port;
}

if ( !$ouser || !$opass || !$oport ) {
	my @param = &get_db_param($ohost);
	$ouser = $param[0] unless $ouser;
	$opass = $param[1] unless $opass;
	$oport = $param[2] unless $oport;
}

if ( !$name ) {
	print STDERR
	  "Can't load sequence set without a target loutre/pipeline database name\n";
	print STDERR "-host $host -user $user -pass $pass\n";
}

if ( !scalar(@seq_sets) ) {
	print STDERR "Need chr|set to be able to run\n";
}

my $otter_dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
	-user   => $ouser,
	-dbname => $oname,
	-host   => $ohost,
	-port   => $oport,
	-pass   => $opass
);
my $module = $name =~ /pipe_/ ? 'Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor' :
								'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $target_dba = $module->new(
	-user   => $user,
	-dbname => $name,
	-host   => $host,
	-port   => $port,
	-pass   => $pass
);
my $target_dbc = $target_dba->dbc();

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
          , c.embl_acc, c.embl_version
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
		SELECT description, hide
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
				$acc,          $sv
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
		my ($desc, $hide) = $description_sth->fetchrow;
		$hide = $hide eq 'N' ? 0 : 1;
		my @chrs = keys %$chr_href;
		throw("Set $sequence_set should contains (only) one chromosome: [@chrs]")
		  unless (scalar(@chrs) == 1);

		$seqset_info->{$sequence_set} = [ shift @chrs,$desc, $hide ];

		throw("No data for '$sequence_set' in $oname")
		  unless ($contig_number);

		print STDOUT $contig_number." contigs retrieved for $sequence_set\n";
	}
}

#
# Load the sequence data set into the loutre/pipeline database
#
{
	print STDOUT
	  "Writing data into database: $name ($host:$port)\n";

	my $dbh = $target_dbc->db_handle;
	$dbh->begin_work;

	my %asm_seq_reg_id;

	my $attr_a		= $target_dba->get_AttributeAdaptor();
	my $cs_a		= $target_dba->get_CoordSystemAdaptor();
	my $slice_a		= $target_dba->get_SliceAdaptor();
	my $seq_a		= $target_dba->get_SequenceAdaptor();
	my $analysis_a	= $target_dba->get_AnalysisAdaptor;
	my $state_info_container = $target_dba->get_StateInfoContainer if($name =~ /pipe_/);

	my $chromosome_cs;
	my $clone_cs;
	my $contig_cs;

	eval {
		$chromosome_cs =
		  $cs_a->fetch_by_name( "chromosome", $cs_version );
		$clone_cs  = $cs_a->fetch_by_name("clone");
		$contig_cs = $cs_a->fetch_by_name("contig");

	};
	if ($@) {
		throw(
			qq{
				A coord_system matching the arguments does not exsist in the cord_system table,
				please ensure you have the right coord_system entry in the database [$@]
			}
		);
	}
	my $slice;
	my $ana = $analysis_a->fetch_by_logic_name('SubmitContig');

	# insert the sequence sets as chromosomes in seq_region
	# table and their attributes in seq_region_attrib & attrib_type tables
	foreach my $set_name ( keys(%end_value) ) {
		my $endv = $end_value{$set_name};
		eval {
			$slice =
			  $slice_a->fetch_by_region( 'chromosome', $set_name, undef, undef,
				undef, $cs_version );
		};
		if ($slice) {
			print STDOUT "Sequence set <$set_name> is already in database <$name>\n";
			throw(  "There is a difference in size for $set_name: stored ["
				  . $slice->length. "] =! new [". $endv. "]" )
			unless ( $slice->length eq $endv );
			$asm_seq_reg_id{$set_name} = $slice->get_seq_region_id;
		}
		else {
			$slice = &make_slice( $set_name, 1, $endv, $endv, 1, $chromosome_cs );
			$asm_seq_reg_id{$set_name} = $slice_a->store($slice);
			$attr_a->store_on_Slice( $slice,
				&make_seq_set_attribute( $seqset_info->{$set_name} ) );
		}
	}

	# insert clone & contig in seq_region, seq_region_attrib,
	# dna and assembly tables
	my $insert_query = qq{
				INSERT IGNORE INTO assembly
				(asm_seq_region_id, cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)
				values
				(?,?,?,?,?,?,?)};
	my $insert_sth = $target_dbc->prepare($insert_query);
	for ( values %$contigs_hashref ) {
		my $chr_name     = $_->[0];
		my $sequence_set = $_->[8];
		my $chr_start    = $_->[1];
		my $chr_end      = $_->[2];
		my $ctg_start    = $_->[3];
		my $ctg_end      = $_->[4];
		my $ctg_ori      = $_->[5];
		my $acc          = $_->[6];
		my $ver          = $_->[7];
		my $acc_ver      = $acc . "." . $ver;

		my $clone;
		my $clone_seq_reg_id;
		my $contig;
		my $contig_name = $acc_ver . "." . "1" . ".";
		my $ctg_seq_reg_id;
		my $seqlen;
		eval {
			$clone  = $slice_a->fetch_by_region( 'clone', $acc_ver );
			$seqlen = $clone->length;
			$contig_name .= $seqlen;
			$contig = $slice_a->fetch_by_region( 'contig', $contig_name );
		};
		if ( $clone && $contig ) {
			print STDOUT "\tclone and contig < ${acc_ver} ; ${contig_name} > are already in the database\n";
			$clone_seq_reg_id = $clone->get_seq_region_id;
			$ctg_seq_reg_id   = $contig->get_seq_region_id;
		}
		elsif ( $clone && !$contig ) {
			$dbh->rollback;
			### code to be added to create contigs related to the clone
			throw(
				"Missing contig entry in the database for clone ${acc_ver}"
			);
		}
		else {
			##fetch the dna sequence from pfetch server with acc_ver id
			my $seqobj;
			eval {
				$seqobj = &pfetch_acc_sv($acc_ver);
				#$seqobj = &pfetch_acc_sv($acc);
			};
			if($@) {
				$dbh->rollback;
				throw($@);
				#print STDERR "$sequence_set => $acc_ver\n";
				#next;
			}
			my $seq = $seqobj->seq;
			$seqlen = $seqobj->length;
			$contig_name .= $seqlen;

			##make clone and insert clone to seq_region table
			$clone = &make_slice( $acc_ver, 1, $seqlen, $seqlen, 1, $clone_cs );
			$clone_seq_reg_id = $slice_a->store($clone);
			if(!$clone_seq_reg_id) {
				$dbh->rollback;
				throw("clone seq_region_id has not been returned for the accession $acc_ver");
			}

			##make attribute for clone and insert attribute to seq_region_attrib & attrib_type tables
			$attr_a->store_on_Slice( $clone,
				&make_clone_attribute( $acc, $ver ) );

			##make contig and insert contig, and associated dna sequence to seq_region & dna table
			$contig =
			  &make_slice( $contig_name, 1, $seqlen, $seqlen, 1, $contig_cs );
			$ctg_seq_reg_id = $slice_a->store( $contig, \$seq );
			if(!$ctg_seq_reg_id){
				$dbh->rollback;
				throw("contig seq_region_id has not been returned for the contig $contig_name");
			}
		}
		##insert chromosome to contig assembly data into assembly table
		$insert_sth->execute( $asm_seq_reg_id{$sequence_set}, $ctg_seq_reg_id,
			$chr_start, $chr_end, $ctg_start, $ctg_end, $ctg_ori );
		##insert clone to contig assembly data into assembly table
		$insert_sth->execute( $clone_seq_reg_id, $ctg_seq_reg_id, 1, $seqlen, 1,
			$seqlen, 1 );

		##prime the input_id_analysis table
		$state_info_container->store_input_id_analysis( $contig->name(), $ana,'' )
			if($do_submit && $name =~ /pipe_/);


	}
	#$dbh->rollback;
	$dbh->commit;

}

#
# Methods
#
sub make_clone_attribute {
	my ( $acc, $ver ) = @_;
	my @attrib;
	my $attrib = &make_attribute(
		'htg',                              'htg',
		'High Throughput phase attribute', '3'
	);
	push @attrib, $attrib;
#	push @attrib,
#	  &make_attribute( 'intl_clone_name', 'International Clone Name',
#		'', '' );
	push @attrib,
	  &make_attribute( 'embl_acc', 'EMBL accession', '', $acc );
	push @attrib, &make_attribute( 'embl_version', 'EMBL Version', '', $ver );
	return \@attrib;
}

sub make_seq_set_attribute {
	my ($arr_ref) = @_;
	my ($chr,$desc,$hide,$write) =  @$arr_ref;
	my @attrib;

	$hide = defined($hide) ? $hide : 1;
	$write = defined($write) ? $write : 1;

	push @attrib,
	  &make_attribute(
		'description',
		'Description',
		'A general descriptive text attribute', $desc
	  );
	push @attrib,
		&make_attribute(
			'chr',
			'Chromosome Name',
			'Chromosome Name Contained in the Assembly', $chr
	);
	push @attrib,
		&make_attribute(
			'write_access',
			'Write access for Sequence Set',
			'1 for writable , 0 for read-only', $write
	);
	push @attrib,
		&make_attribute(
			'hidden',
			'Hidden Sequence Set', '',$hide
	);

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

#
# Pfetch the sequences
#-------------------------------------------------------------------------------
#  If the sequence isn't available from the default pfetch
#  the archive pfetch server is tried.
#
{
	my ( $pfetch );

	sub pfetch_acc_sv {
		my ($acc_ver) = @_;
		print "Fetching '$acc_ver'\n";
		$pfetch ||= Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
		my $seq = $pfetch->get_Seq_by_acc($acc_ver);
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

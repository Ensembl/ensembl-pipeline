#!/software/bin/perl -w

=head1 NAME

load_from_loutre_or_pipe.pl

=head1 SYNOPSIS

load_from_loutre_or_pipe.pl

=head1 DESCRIPTION

This script is used to copy one or several sequence sets between two new schema databases (loutre/pipeline).
If the login, password and port parameters are not provided, they will be recovered from the ~/.netrc file.
See the Net::Netrc module for more details.

here is an example commandline

./load_from_loutre_or_pipe.pl
-set chr11-02
-shost ecs4
-sport 3352
-sname loutre_human
-suser pipuser
-spass *****
-thost otterpipe2
-tport 3352
-tname pipe_human
-tuser pipuser
-tpass *****


=head1 OPTIONS

    -shost (default:otterlive)   source database host name
    -sname (default:loutre_human)   For RDBs, what name to connect to
    -suser (check the ~/.netrc file) For RDBs, what username to connect as
    -spass (check the ~/.netrc file) For RDBs, what password to use
    -sport (check the ~/.netrc file)   For RDBs, what port to use
    -thost (default:otterpipe1)   target database host name
    -tname (no default)  For RDBs, what name to connect to
    -tuser (check the ~/.netrc file)  For RDBs, what username to connect as
    -tpass (check the ~/.netrc file)  For RDBs, what password to use
    -tport (check the ~/.netrc file)   For RDBs, what port to use

    -chromosome_cs_version (default:Otter) the version of the coordinate system being stored
    -set|chr	the sequence set(s) to load
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
my $shost = 'otterlive';
my $sport = '';
my $sname = 'loutre_human';
my $suser = '';
my $spass = '';

my $thost = 'otterpipe1';
my $tport = '';
my $tname = '';
my $tuser = '';
my $tpass = '';

my $chromosome_cs_version = 'Otter';
my @seq_sets;
my $do_submit = 1;	  # Set if we want to prime the pipeline with the SubmitContig analysis

my $usage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
	'shost:s'                 => \$shost,
	'sport:n'                 => \$sport,
	'sname:s'                 => \$sname,
	'suser:s'                 => \$suser,
	'spass:s'                 => \$spass,
	'thost:s'                 => \$thost,
	'tport:n'                 => \$tport,
	'tname=s'                 => \$tname,
	'tuser:s'                 => \$tuser,
	'tpass:s'                 => \$tpass,
	'chromosome_cs_version:s' => \$chromosome_cs_version,
	'chr|set=s'               => \@seq_sets,
	'no_submit!'			  => sub{ $do_submit = 0 },
	'submit!'			  	  => \$do_submit,
	'h|help!'                 => $usage
  )
  or $usage->();

if ( !$tuser || !$tpass || !$tport ) {
	my @param = &get_db_param($thost);
	$tuser = $param[0] unless $tuser;
	$tpass = $param[1] unless $tpass;
	$tport = $param[2] unless $tport;
}

if ( !$suser || !$spass || !$sport ) {
	my @param = &get_db_param($shost);
	$suser = $param[0] unless $suser;
	$spass = $param[1] unless $spass;
	$sport = $param[2] unless $sport;
}

if ( !$tname ) {
	print STDERR
	  "Can't load sequence set without a target database name\n";
	print STDERR "-thost $thost -tuser $tuser -tpass $tpass\n";
}

if ( !scalar(@seq_sets) ) {
	print STDERR "Need chr|set to be able to run\n";
}

my $source_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-user   => $suser,
	-dbname => $sname,
	-host   => $shost,
	-port   => $sport,
	-pass   => $spass
);
my $source_dbc = $source_dba->dbc();

my $module = $tname =~ /pipe_/ ? 'Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor' :
								'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $target_dba = $module->new(
	-user   => $tuser,
	-dbname => $tname,
	-host   => $thost,
	-port   => $tport,
	-pass   => $tpass
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
"Getting Sequence set @seq_sets from source database: $sname ($shost:$sport)\n";
	my $contigs_sth = $source_dbc->prepare(
		q{
        SELECT t.value, a.asm_start, a.asm_end, a.cmp_start, a.cmp_end, a.ori,
            substring_index(s2.name,'.',1) ,
            substring_index(substring_index(s2.name,'.',-3),'.',1) ,
            t.value
		FROM assembly a, seq_region_attrib t, seq_region s1, seq_region s2
		WHERE s1.name = ?
		AND s2.seq_region_id = a.cmp_seq_region_id
		AND s1.seq_region_id = a.asm_seq_region_id
		AND t.seq_region_id = s1.seq_region_id
		AND t.attrib_type_id = 97
		ORDER BY  a.asm_start
		}
	);
	my $description_sth = $source_dbc->prepare(
		q{
		SELECT t.value
		FROM seq_region_attrib t, seq_region s
		WHERE s.name = ?
		AND t.seq_region_id = s.seq_region_id
		AND t.attrib_type_id = 49
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
		my @chrs = keys %$chr_href;
		$seqset_info->{$sequence_set} = [ $desc, @chrs ];

		throw("No data for '$sequence_set' in $sname")
		  unless ($contig_number);
		throw("Set $sequence_set should contains (only) one chromosome: [@chrs]")
		  unless (scalar(@chrs) == 1);
		print STDOUT $contig_number." contigs retrieved for $sequence_set\n";
	}
}

#
# Load the sequence data set into the target database
#
{
	print STDOUT
	  "Writing data into target database: $tname ($thost:$tport)\n";

	my $dbh = $target_dbc->db_handle;
	$dbh->begin_work;

	my %asm_seq_reg_id;

	my $attr_a		= $target_dba->get_AttributeAdaptor();
	my $cs_a		= $target_dba->get_CoordSystemAdaptor();
	my $slice_a		= $target_dba->get_SliceAdaptor();
	my $seq_a		= $target_dba->get_SequenceAdaptor();
	my $analysis_a	= $target_dba->get_AnalysisAdaptor;
	my $state_info_container = $target_dba->get_StateInfoContainer if($tname =~ /pipe_/);

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
				undef, $chromosome_cs_version );
		};
		if ($slice) {
			print STDOUT "Sequence set <$set_name> is already in target database <$tname>\n";
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
	my $insert_query = qq {
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
			print STDOUT "\tclone and contig < ${acc_ver} ; ${contig_name} > are already in the pipeline database\n";
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
			if($do_submit && $tname =~ /pipe_/);


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
	my ($desc,@chr) =  @$arr_ref;
	my @attrib;
	push @attrib,
	  &make_attribute(
		'description',
		'Description',
		'A general descriptive text attribute', $desc
	  );
	foreach my $ch (@chr) {
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

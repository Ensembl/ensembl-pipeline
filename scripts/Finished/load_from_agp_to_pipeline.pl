#!/usr/local/ensembl/bin/perl -w

=head1 NAME

load_from_agp_to_pipeline.pl

=head1 SYNOPSIS

load_from_agp_to_pipeline.pl

=head1 DESCRIPTION

This script is used to load one sequence set from an agp file into the 
following tables of a pipeline database (new schema):
seq_region, seq_region_attrib, attrib_type, assembly and dna. 
The sequence that is loaded into the dna table is either Pfetched or fetched from a raw file.
If the login, password and port parameters are not provided, they will be 
recovered from the ~/.netrc file. See the Net::Netrc module for more details.

here is an example commandline

./load_from_agp_to_pipeline.pl \
-set chr11-02 \
-description 'chromosome 11' \
-p_host otterpipe2 \
-p_port 3352 \
-p_name pipe_human \
-p_user pipuser \
-p_pass *****

=head1 OPTIONS

    -p_host (default:otterpipe1)   host name for the Pipeline database (gets put as phost= in locator)
    -p_name (no default)  For RDBs, what name to connect to (pname= in locator)
    -p_user (check the ~/.netrc file)  For RDBs, what username to connect as (puser= in locator)
    -p_pass (check the ~/.netrc file)  For RDBs, what password to use (ppass= in locator)
    -p_port (check the ~/.netrc file)   For RDBs, what port to use (pport= in locator)

    -chromosome_cs_version (default:Otter) the version of the coordinate system being stored
    -set	the sequence set name
    -description the sequence set description
    -nosubmit	Used to avoid the pipeline priming with the SubmitContig analysis 
    -help|h		displays this documentation with PERLDOC

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use strict;

use Getopt::Long;
use Net::Netrc;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::SequenceAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
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

my $chromosome_cs_version = 'Otter';
my $set;
my $description;
my $do_submit = 1; # Set if we want to prime the pipeline with the SubmitContig analysis

my $usage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
	'p_host:s'                => \$phost,
	'p_port:n'                => \$pport,
	'p_name:s'                => \$pname,
	'p_user:s'                => \$puser,
	'p_pass:s'                => \$ppass,
	'chromosome_cs_version:s' => \$chromosome_cs_version,
	'set=s'                   => \$set,
	'description=s'           => \$description,
	'submit!'                 => \$do_submit,
	'h|help!'                 => $usage
  )
  or $usage->();

my $agp_file =
  $ARGV[0];    # takes the remaining argument as the filename to be read
throw("cannot load assembly details, as there is no agp file")
  unless ( defined $agp_file );

throw("no description given") unless ( defined $description );

throw("No sequence set name given") unless ($set);

if ( !$puser || !$ppass || !$pport ) {
	my @param = &get_db_param($phost);
	$puser = $param[0] unless $puser;
	$ppass = $param[1] unless $ppass;
	$pport = $param[2] unless $pport;
}

if ( !$pname ) {
	print STDERR
	  "Can't load sequence set without a target pipeline database name\n";
	print STDERR "-p_host $phost -p_user $puser -p_pass $ppass\n";
}

my $pipe_dba = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(
	-user   => $puser,
	-dbname => $pname,
	-host   => $phost,
	-port   => $pport,
	-pass   => $ppass
);

my %end_value;
my $contigs_hashref = {};
my $seqset_info     = {};

#
# Get the sequence data set from agp file and store it in a Hashtable.
#
{
	print STDOUT "Getting data from agp file $agp_file\n";

	open( my $fh, "$agp_file" ) or die "Can't read '$agp_file' : $!";

	my $agp_chr_name;
	my $contig_number = 0;
	my $chr_href;

	while (<$fh>) {
		chomp;
		next if $_ =~ /^\#/;
		my (
			$input_type, $agp_chr_name, $chr_start, $chr_end,
			$n,          $type,         $acc_ver,   $ctg_start,
			$ctg_end,    $ctg_ori
		  )
		  = check_line($_);
		if ( $input_type ne 'AGP' ) {
			#print STDOUT "Skipping $input_type type line\n";
			next;
		}

		$agp_chr_name =~ s/^chr//i;

		# translate orientation to integer
		if ( $ctg_ori eq '-' ) {
			$ctg_ori = -1;
		}
		elsif ( $ctg_ori eq '+' ) {
			$ctg_ori = 1;
		}
		else {
			throw("Invalid orientation '$ctg_ori'\n");
		}

		#split into accesion number and version number
		my ( $acc, $sv ) = $acc_ver =~ /^(.+)\.(\d+)$/;

		if ( !$end_value{$set} ) {
			$end_value{$set} = $chr_end;
		}
		else {
			if ( $chr_end > $end_value{$set} ) {
				$end_value{$set} = $chr_end;
			}
		}
		$contigs_hashref->{ $acc . $sv . $set } = [
			$agp_chr_name, $chr_start, $chr_end,
			$ctg_start,    $ctg_end,   $ctg_ori,
			$acc,          $sv,        $set
		];
		$chr_href->{$agp_chr_name} = 1;
		$contig_number++;
	}
	my @chrs = keys %$chr_href;
	$seqset_info->{$set} = [ $description, @chrs ];

	close $fh;
	
	throw("AGP should contains (only) one chromosome: [@chrs]")
		  unless (scalar(@chrs) == 1);
	print STDOUT $contig_number." contigs retrieved from $agp_file\n";
}

#
# Load the sequence data set into the Pipeline database
#
{
	print STDOUT
	  "Writing data into Pipeline database: $pname ($phost:$pport)\n";

	my %asm_seq_reg_id;

	my $attr_a               = $pipe_dba->get_AttributeAdaptor();
	my $cs_a                 = $pipe_dba->get_CoordSystemAdaptor();
	my $slice_a              = $pipe_dba->get_SliceAdaptor();
	my $seq_a                = $pipe_dba->get_SequenceAdaptor();
	my $analysis_a           = $pipe_dba->get_AnalysisAdaptor;
	my $state_info_container = $pipe_dba->get_StateInfoContainer;

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

	# insert the sequence set as a chromosome in seq_region
	# table and its attributes in seq_region_attrib & attrib_type tables
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
				  . $slice->length. "] =! new [". $endv. "]" )
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

	# insert clone & contig in seq_region, seq_region_attrib, 
	# dna and assembly tables
	my $insert_query = qq {
			INSERT IGNORE INTO assembly 
			(asm_seq_region_id, cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori) 
			values 
			(?,?,?,?,?,?,?)};
	my $insert_sth = $pipe_dba->prepare($insert_query);
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
			$clone = $slice_a->fetch_by_region( 'clone', $acc_ver );
			$seqlen = $clone->length;
			$contig_name .= $seqlen;
			$contig = $slice_a->fetch_by_region( 'contig', $contig_name );
		};
		if ( $clone && $contig ) {
			warn(
				qq{
					clone and contig < ${acc_ver} ; ${contig_name} >
					are already in the pipeline database
				}
			);
			$clone_seq_reg_id = $clone->get_seq_region_id;
			$ctg_seq_reg_id   = $contig->get_seq_region_id;
		}
		elsif ( $clone && !$contig ) {
			### code to be added to create contigs related to the clone
			throw(  
				"Missing contig entry in the database for clone ${acc_ver}" 
			);
		}
		else {
			##fetch the dna sequence from pfetch server with acc_ver id
			my $seqobj ||= &pfetch_acc_sv($acc_ver);
			my $seq = $seqobj->seq;
			$seqlen = $seqobj->length;
			$contig_name .= $seqlen;

			##make clone and insert clone to seq_region table
			$clone = &make_slice( $acc_ver, 1, $seqlen, $seqlen, 1, $clone_cs );
			$clone_seq_reg_id = $slice_a->store($clone);
			throw(
				"clone seq_region_id has not been returned for the accession $acc_ver"
			) unless $clone_seq_reg_id;

			##make attribute for clone and insert attribute to seq_region_attrib & attrib_type tables
			$attr_a->store_on_Slice( $clone,
				&make_clone_attribute( $acc, $ver ) );

			##make contig and insert contig, and associated dna sequence to seq_region & dna table
			$contig =
			  &make_slice( $contig_name, 1, $seqlen, $seqlen, 1, $contig_cs );
			$ctg_seq_reg_id = $slice_a->store( $contig, \$seq );
			throw(
				"contig seq_region_id has not been returned for the contig $contig_name"
			) unless $ctg_seq_reg_id;
		}
		##insert chromosome to contig assembly data into assembly table
		$insert_sth->execute( $asm_seq_reg_id{$sequence_set}, $ctg_seq_reg_id, 
			$chr_start, $chr_end, $ctg_start, $ctg_end, $ctg_ori );
		##insert clone to contig assembly data into assembly table
		$insert_sth->execute( $clone_seq_reg_id, $ctg_seq_reg_id, 1, $seqlen, 1,
			$seqlen, 1 );

		##prime the input_id_analysis table
		$state_info_container->store_input_id_analysis( $contig->name(), $ana,'' ) 
			if($do_submit);

	}

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
	my $desc = shift @$arr_ref;
	my @attrib;
	push @attrib,
	  &make_attribute(
		'description',
		'Description',
		'A general descriptive text attribute', $desc
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
# check line
#

sub check_line {
	my ($line) = @_;
	my $input_type;

# 0       1       2       3       4       5               6       7       8     9
# chr_20  2808333 2934911 29      F       AL121905.0      101     126679  +     Optional comment
# splits each line into its component parts - puts line in a temporary array (splits the line on whitespace)
	my @line_in_array = split /\s+/, $line, 10;

	if ( $line_in_array[4] && $line_in_array[4] =~ /^[AUF]$/ ) {
		$input_type = "AGP";
	}
	else {
		return ('SKIP');
	}
	return $input_type, @line_in_array;
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


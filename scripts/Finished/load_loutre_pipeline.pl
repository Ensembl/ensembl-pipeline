#!/software/bin/perl -w

=head1 NAME

load_loutre_pipeline.pl

=head1 SYNOPSIS

load_loutre_pipeline.pl

=head1 DESCRIPTION

This script is used to load an agp file into a loutre/ensembl database and its corresponding pipeline database if requested.
The loutre key 'pipeline_db_rw_head' in the meta table is used to retrieve the pipeline connexion parameters.
The sequence loaded into the dna table is either Pfetched or fetched from a raw file.
If the login, password and port parameters of the loutre connexion are not provided, they will be
recovered from the ~/.netrc file. See the Net::Netrc module for more details.


=head1 OPTIONS

    -host   (default:otterlive)        host name for the loutre database (gets put as phost= in locator)
    -dbname (no default)               for RDBs, what name to connect to (pname= in locator)
    -user   (check the ~/.netrc file)  for RDBs, what username to connect as (puser= in locator)
    -pass   (check the ~/.netrc file)  for RDBs, what password to use (ppass= in locator)
    -port   (check the ~/.netrc file)  for RDBs, what port to use (pport= in locator)

    -skip_type comma separated list of clone status to ignore (e.g. A,D,F,G,O,P,W)
    -assembly attach a equiv_asm attribute (e.g. NCBI36, GRCh37)
    -cs_name (default:chromosome) the name of the coordinate system being stored
    -cs_version (default:Otter) the version of the chromosome coordinate system being stored
    -set	 the sequence set name
    -description the sequence set description
    -do_pipe	(default: true) populate the pipeline sattelite db
    -nosubmit	 (default: true) don't prime the analysis pipeline
		i.e. won't add submitcontig rows in the input_id_anlysis table
    -help|h      displays this documentation with PERLDOC

=head1 EXAMPLES

Here are some command line examples:

=head2 1. Populate a loutre db and its pipeline sattelite with a GRCh38 chromosome assembly

    ./load_loutre_pipeline.pl \
	-set chr11-02 \
	-description 'chromosome 11 GRCh38' \
	-assembly GRCh38
	-dbname loutre_human \
	[AGP_FILE]

=head2 2. Load a NCBIM36 mouse chromosome in loutre only (needed for the creation of the mapping)

    ./load_loutre_pipeline.pl \
	-set 11 \
	-description 'chromosome 11 NCBIM36' \
	-cs_version NCBIM36
	-dbname loutre_mouse \
	[AGP_FILE]

=head2 3. Populate one database only

    ./load_loutre_pipeline.pl \
	-set chr11-02 \
	-description 'chromosome 11' \
	-dbname loutre_human \
	-do_pipe 0 \
	[AGP_FILE]

=head2 4. Load a subregion AGP

    ./load_loutre_pipeline.pl \
	-set subregion-01 \
	-cs_name subregion \
	-cs_version '' \
	-description 'This subregion description will be replaced anyway' \
	-dbname loutre_human \
	-do_pipe 0 \
	[AGP_FILE]


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

# loutre connexion parameters, default values.
my $host = 'otterlive';
my $port = '';
my $name = '';
my $user = '';
my $pass = '';

my $cs_name = 'chromosome';
my $cs_version = 'Otter';
my $set;
my $description;
my $assembly;
my $do_pipe = 1; # Set to load loutre and pipeline dbs
my $do_submit = 1; # Set if we want to prime the pipeline with the SubmitContig analysis
my @skip;

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
	'description=s'         => \$description,
	'skip_type=s'           => \@skip,
	'assembly=s'			=> \$assembly,
	'do_pipe!'				=> \$do_pipe,
	'submit!'               => \$do_submit,
	'h|help!'               => $usage
  )
  or $usage->();

@skip = split(/,/,join(',',@skip));

my $agp_file =
  $ARGV[0];    # takes the remaining argument as the filename to be read
throw("cannot load assembly details, as there is no agp file")
  unless ( defined $agp_file );

throw("No description given") unless ( defined $description );

throw("No sequence set name given") unless ($set);

throw("Can't load a sequence set without a target datbase name") unless $name;

if ( !$user || !$pass || !$port ) {
	my @param = &get_db_param($host);
	$user = $param[0] unless $user;
	$pass = $param[1] unless $pass;
	$port = $param[2] unless $port;
}

my $dbas;

my $dba = $name =~ /pipe_/ ?
	Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(
		-user   => $user,
		-dbname => $name,
		-host   => $host,
		-port   => $port,
		-pass   => $pass) :
	Bio::EnsEMBL::DBSQL::DBAdaptor->new(
		-user   => $user,
		-dbname => $name,
		-host   => $host,
		-port   => $port,
		-pass   => $pass);

push @$dbas, $dba;

my $pipe_dba;
if($do_pipe) {
	my $meta_container = $dba->get_MetaContainer();
	my ($pipe_param) = @{$meta_container->list_value_by_key('pipeline_db_rw_head')};
	if($pipe_param) {
		$pipe_dba = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(eval $pipe_param);
		push @$dbas, $pipe_dba;
	} else {
		throw("You need to add the meta key 'pipeline_db_rw_head' in ${host}/${name} or unset the do_pipe option to populate ${host}/${name} only\n");
	}
}

my %end_value;
my $contigs_hashref = {};
my $seqset_info     = {};

#
# Get the sequence set data from the agp and store it in a Hashtable.
#
{
	print STDOUT "Getting data from the agp file $agp_file\n";

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
		my $contig_key = $acc.$sv.$ctg_start.$ctg_end.$set;
		$contigs_hashref->{ $contig_key } = [
			$agp_chr_name, $chr_start, $chr_end,
			$ctg_start,    $ctg_end,   $ctg_ori,
			$acc,          $sv,        $set
		];
		$chr_href->{$agp_chr_name} = 1;
		$contig_number++;
	}
	my @chrs = keys %$chr_href;
	throw("AGP should contains (only) one chromosome: [@chrs]")
		unless (scalar(@chrs) == 1);

	$seqset_info->{$set} = [ shift @chrs ,$description, undef, undef, $assembly  ];

	close $fh;


	print STDOUT $contig_number." contigs retrieved from $agp_file\n";
}

#
# Load the sequence data set into the databases
#
{
	my %objects;

	foreach (@$dbas) { $_->dbc->db_handle->begin_work};

	foreach my $dba (@$dbas) {
		print STDOUT "Writing data into database: ".$dba->dbc->dbname." (".$dba->dbc->host.":".$dba->dbc->port.")\n";
		my %asm_seq_reg_id;

		my $attr_a               = $dba->get_AttributeAdaptor();
		my $cs_a                 = $dba->get_CoordSystemAdaptor();
		my $slice_a              = $dba->get_SliceAdaptor();
		my $seq_a                = $dba->get_SequenceAdaptor();
		my $analysis_a           = $dba->get_AnalysisAdaptor;
		my $state_info_container = $dba->get_StateInfoContainer if ($dba->dbc->dbname =~ /pipe_/);

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
			foreach (@$dbas) { $_->dbc->db_handle->rollback};
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
		my $stored_sets = [];
		foreach my $name ( keys(%end_value) ) {
			my $endv = $end_value{$name};
			eval {
				$slice =
				  $slice_a->fetch_by_region( 'chromosome', $name, undef, undef,
					undef, $cs_version );
			};
			if ($slice) {
				push @$stored_sets, $name;
			}
			else {
				$slice = &make_slice( $name, 1, $endv, $endv, 1, $chromosome_cs );
				$asm_seq_reg_id{$name} = $slice_a->store($slice);
				$attr_a->store_on_Slice( $slice,
					&make_seq_set_attribute( $seqset_info->{$name} ) );
			}
		}
		if(@$stored_sets) {
			foreach (@$dbas) { $_->dbc->db_handle->rollback};
			throw("Sequence set ".join(" ",@$stored_sets)." is/are already in database <".$dba->dbc->dbname.">");
		}

		# insert clone & contig in seq_region, seq_region_attrib,
		# dna and assembly tables
		my $insert_query = qq{
				INSERT IGNORE INTO assembly
				(asm_seq_region_id, cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)
				values
				(?,?,?,?,?,?,?) };
		my $insert_sth = $dba->dbc->prepare($insert_query);
		my $error = '';
		CLONE: for ( keys %$contigs_hashref ) {

			my $chr_name     = $contigs_hashref->{$_}->[0];
			my $sequence_set = $contigs_hashref->{$_}->[8];
			my $chr_start    = $contigs_hashref->{$_}->[1];
			my $chr_end      = $contigs_hashref->{$_}->[2];
			my $contig_start = $contigs_hashref->{$_}->[3];
			my $contig_end   = $contigs_hashref->{$_}->[4];
			my $contig_ori   = $contigs_hashref->{$_}->[5];
			my $acc          = $contigs_hashref->{$_}->[6];
			my $ver          = $contigs_hashref->{$_}->[7];
			my $acc_ver      = $acc . "." . $ver;

			my $clone;
			my $clone_length;
			my $clone_seq_reg_id;
			my $contig;
			my $contig_name = $acc_ver . ".1.";
			my $ctg_seq_reg_id;


			eval {
				$clone = $slice_a->fetch_by_region( 'clone', $acc_ver );
				$clone_length = $clone->length;
				$contig_name .= $clone_length;
				$contig = $slice_a->fetch_by_region( 'contig', $contig_name );
			};
			if ( $clone && $contig ) {
				print STDOUT "\tclone and contig < ${acc_ver} ; ${contig_name} > are already in database\n";
				$clone_seq_reg_id = $clone->get_seq_region_id;
				$ctg_seq_reg_id   = $contig->get_seq_region_id;
			} else {
				##fetch the dna sequence from pfetch server or a local fasta file with acc_ver id
				my $seqobj;

				if($objects{$acc_ver}){
					$seqobj = $objects{$acc_ver};
				}else{
					eval {
						$seqobj = &pfetch_acc_sv($acc_ver);
						$objects{$acc_ver} = $seqobj;
					};
					if($@) {
						$error .= $@;
						next CLONE;
					}
				}

				$clone_length = $seqobj->length;
				$contig_name .= $clone_length unless $clone;
				my $contig_seq = $seqobj->seq;


				if(! $clone){
					##make clone and insert clone to seq_region table
					$clone = &make_slice( $acc_ver, 1, $clone_length, $clone_length, 1, $clone_cs );
					$clone_seq_reg_id = $slice_a->store($clone);
					if(!$clone_seq_reg_id) {
						foreach (@$dbas) { $_->dbc->db_handle->rollback};
						throw("clone seq_region_id has not been returned for the accession $acc_ver");
					}
					##make attribute for clone and insert attribute to seq_region_attrib & attrib_type tables
					$attr_a->store_on_Slice( $clone,
						&make_clone_attribute( $acc, $ver ) );
				} else {
					print STDOUT "\tclone < ${acc_ver} > is already in database\n";
					$clone_seq_reg_id = $clone->get_seq_region_id;
				}

				if(!$contig) {
					##make contig and insert contig, and associated dna sequence to seq_region & dna table
					$contig =
					  &make_slice( $contig_name, 1, $clone_length,$clone_length, 1, $contig_cs );
					$ctg_seq_reg_id = $slice_a->store( $contig, \$contig_seq );
					if(!$ctg_seq_reg_id){
						foreach (@$dbas) { $_->dbc->db_handle->rollback};
						throw("contig seq_region_id has not been returned for the contig $contig_name");
					}
				}
			}

			##insert chromosome to contig assembly data into assembly table
			$insert_sth->execute( $asm_seq_reg_id{$sequence_set}, $ctg_seq_reg_id,
				$chr_start, $chr_end, $contig_start, $contig_end, $contig_ori );
			##insert clone to contig assembly data into assembly table
			$insert_sth->execute( $clone_seq_reg_id, $ctg_seq_reg_id, 1, $clone_length, 1,
				$clone_length, 1 );

			##prime the input_id_analysis table
			$state_info_container->store_input_id_analysis( $contig->name(), $ana,'' )
				if($do_submit and $dba->dbc->dbname =~ /pipe_/);

		}
		if($error) {
			foreach (@$dbas) { $_->dbc->db_handle->rollback};
			throw($error);
		}
	}
	foreach (@$dbas) { $_->dbc->db_handle->commit};
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
	my ( $chr, $desc, $hide, $write, $assembly) =  @$arr_ref;
	my @attrib;

	$hide = defined($hide) ? $hide : 1;
	$write = defined($write) ? $write : 0;

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

	if($assembly) {
		push @attrib,
			&make_attribute(
				'equiv_asm',
				'Equivalent EnsEMBL assembly',
				'For full chromosomes made from NCBI AGPs',
				$assembly
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

## AGP file specification 1.1
# The sequencing status of the component.
# These typically correspond to keywords in
# the International Sequence Database
# (GenBank/EMBL/DDBJ) submission.
#  Current acceptable values are:
#    A=Active Finishing
#    D=Draft HTG (often phase1 and phase2
#    are called Draft, whether or not they have
#    the draft keyword).
#    F=Finished HTG (phase 3)
#    G=Whole Genome Finishing
#    N=gap with specified size
#    O=Other sequence (typically means no
#    HTG keyword)
#    P=Pre Draft
#    U= gap of unknown size, typically
#    defaulting to predefined values.
#    W=WGS contig
##

	my $known_type = "^[ADFGNOPUW]\$";

	my $accepted_type = "ADFGOPW";
	if(@skip) {
		my $skip_pattern = "[".join('',@skip)."]";
		$accepted_type =~ s/$skip_pattern//g;
	}

	my $accepted_type_p = "^[$accepted_type]\$";

# 0       1       2       3       4       5               6       7       8     9
# chr_20  2808333 2934911 29      F       AL121905.0      101     126679  +     Optional comment
# splits each line into its component parts - puts line in a temporary array (splits the line on whitespace)
	my @line_in_array = split /\s+/, $line, 10;

	# throw exception if new type
	throw("Unknown sequence status in AGP row\n$line")
		if ($line_in_array[4] && $line_in_array[4] !~ /$known_type/);

	if ( $line_in_array[4] && $line_in_array[4] =~ /$accepted_type_p/ ) {
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


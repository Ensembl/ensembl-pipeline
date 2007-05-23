#!/software/bin/perl -w

=head1 NAME

load_clone_features.pl

=head1 SYNOPSIS

load_clone_features.pl

=head1 DESCRIPTION

This script is used to load clone features from a GFF file into a pipeline database.
If the login, password and port parameters are not provided, they will be
recovered from the ~/.netrc file. See the Net::Netrc module for more details.

here is an example commandline

./load_clone_features.pl -p_host otterpipe2 -p_port 3352 -p_name pipe_human -p_user pipuser -p_pass ***** -type SimpleFeature	features.gff


=head1 OPTIONS

	-p_host (default:otterpipe1)   host name for the Pipeline database (gets put as phost= in locator)
	-p_name (no default)  For RDBs, what name to connect to (pname= in locator)
	-p_user (check the ~/.netrc file)  For RDBs, what username to connect as (puser= in locator)
	-p_pass (check the ~/.netrc file)  For RDBs, what password to use (ppass= in locator)
	-p_port (check the ~/.netrc file)   For RDBs, what port to use (pport= in locator)

	-type		the features are of this type (default: SimpleFeature)
	-verbose
	-help|h		displays this documentation with PERLDOC

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use strict;

use Getopt::Long;
use Net::Netrc;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use GFF;
use GFF::GeneFeature;


# pipeline connexion parameters, default values.
my $phost = 'otterpipe1';
my $pport = '';
my $pname = '';
my $puser = '';
my $ppass = '';

my $type = 'SimpleFeature';
my $verbose;
my $usage = sub { exec( 'perldoc', $0 ); };

&GetOptions(
	'p_host:s'                 => \$phost,
	'p_port:n'                 => \$pport,
	'p_name=s'                 => \$pname,
	'p_user:s'                 => \$puser,
	'p_pass:s'                 => \$ppass,
	'type:s'                   => \$type,
	'verbose!'				   => \$verbose,
	'h|help!'                  => $usage
  )
  or $usage->();

my @gff_files = @ARGV;
throw("cannot load features as there is no gff files")
  unless ( @gff_files );

throw("need to specify a feature type")
  unless ( $type );

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
my $pipe_dbc = $pipe_dba->dbc();

# read GFF files
my $gff_arr;
foreach my $file (@gff_files) {
	print STDOUT "parse $file\n" if $verbose;
	open(my $fh,"$file")  or die "Can't read '$file' : $!";
	while(my $line = <$fh>) {
		print STDOUT $line;
		push @$gff_arr, GFF::GeneFeature->new_from_line($line);
	}
}

# write features into database
my $dbh = $pipe_dbc->db_handle;
$dbh->begin_work;

my $slice_a = $pipe_dba->get_SliceAdaptor();
my $analysis_a = $pipe_dba->get_AnalysisAdaptor();
my $method = "get_${type}Adaptor";
my $feature_a = $pipe_dba->${method};

my $feature_module = "Bio::EnsEMBL::${type}";
my $file = "$feature_module.pm";
$file =~ s{::}{/}g;
require "$file";

my %analysis;
my %slices;
my $features;

foreach my $gff (@$gff_arr) {
	my ($seq_name,$analysis,$start,$end,$score,$strand) = map { $gff->$_ }
		('seqname','feature','start','end','score','strand');
	$strand = $strand eq '+' ? 1 : -1;
	$score  = 0 if $score eq '.';
	# get the analysis object
	my $ana = $analysis{$analysis};
	if(!$ana) {
		$ana = $analysis_a->fetch_by_logic_name($analysis);
		$analysis{$analysis} = $ana;
	}
	# get the slice object
	my $slice = $slices{$seq_name};
	if(!$slice) {
		my $slice_clone = $slice_a->fetch_by_region('',$seq_name);
		$slice = $slice_clone->project('contig')->[0]->to_Slice();
		$slices{$seq_name} = $slice;
	}

	if(!$slice || !$ana) {
		$dbh->rollback;
		throw("Cannot find Slice $seq_name [$slice] or Analysis $analysis [$ana]\n");
	}

	# save features
	my $feature = $feature_module->new
                        (-start    => $start,
                         -end      => $end,
                         -strand   => $strand,
                         -slice    => $slice,
                         -analysis => $ana,
                         -display_label => 'Gap',
                         -score => $score);
    push @$features,$feature;
}

eval {
	$feature_a->store(@$features);
};
if($@){
	$dbh->rollback;
	throw($@);
}

# make the changes permanent
$dbh->commit;

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



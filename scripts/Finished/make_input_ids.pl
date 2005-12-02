#!/usr/local/ensembl/bin/perl

=head1 NAME

make_input_ids.pl  

=head1 SYNOPSIS

make_input_ids.pl -p_host ecs1a -p_user user -p_pass **** -p_name pipe_mouse -seq_region chr12-03

=head1 DESCRIPTION

This script allows to prime the pipeline database with the input_ids of the seq_region name(s) provided.
If no values are provided for the database connexion (login and password), the ~/.netrc file will be checked. 
See Net::Netrc module for more details.

=head1 OPTIONS

    -p_host    (default: otterpipe1) host name for database (gets put as host= in locator)
    -p_port    (default: 3302) For RDBs, what port to connect to (port= in locator) 
    -p_name    For RDBs, what name to connect to (p_name= in locator)
    -p_user    (check the ~/.netrc file) For RDBs, what username to connect as (user= in locator)
    -p_pass    (check the ~/.netrc file) For RDBs, what password to use (pass= in locator)

    -logic_name (default: SubmitContig) the logic_name of the analysis object which needs to be 
                associated with these entries 

	-seqreg_name the seq_region name you want to prime (could be used several times)
	-cs (default: chromosome) the coordinate system associated with the seq_region name 
    -cs_version (default: Otter) the version of the coord system you want 
    -target_cs (default: contig) the target coordinate system you want slices in 
    -target_cs_version the version of the target coord system you want

    -verbose if you want more information about what the script is doing
    -help      Displays script documentation with PERLDOC

=head1 CONTACT

Mustapha Larbaoui B<email> ml6@sanger.ac.uk

=cut

use strict;
use Getopt::Long;
use Net::Netrc;
use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host = 'otterpipe1';
my $user;
my $pass;
my $port = 3302;
my $p_name;
my $logic_name = 'SubmitContig';
my $cs         = 'chromosome';
my $cs_version = 'Otter';
my $target_cs  = 'contig';
my $target_cs_version;
my @seqreg_name;
my $verbose;
my $help = 0;
&GetOptions(
	'p_host:s'            => \$host,
	'p_port:n'            => \$port,
	'p_user:s'            => \$user,
	'p_pass:s'            => \$pass,
	'p_name:s'            => \$p_name,
	'seqreg_name:s'       => \@seqreg_name,
	'cs:s'                => \$cs,
	'cs_version:s'        => \$cs_version,
	'target_cs:s'         => \$target_cs,
	'target_cs_version:s' => \$target_cs_version,
	'logic_name:s'        => \$logic_name,
	'verbose!'            => \$verbose,
	'h|help'              => \$help
);

if ($help) {
	exec( 'perldoc', $0 );
}

if ( !$user || !$pass ) {
	my $ref = Net::Netrc->lookup($host);
	throw("Need to provide pipeline login ( -p_user ) and password ( -p_pass )")
	  unless ($ref);
	$user = $ref->login;
	$pass = $ref->password;
}

if ( !$p_name ) {
	throw("You must specify a database name");
}

if ( !@seqreg_name ) {
	throw("You must at least specify one seq_region name");
}

my $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(
	-host   => $host,
	-user   => $user,
	-pass   => $pass,
	-port   => $port,
	-dbname => $p_name
);

my $slice_a              = $db->get_SliceAdaptor;
my $analysis_a           = $db->get_AnalysisAdaptor;
my $state_info_container = $db->get_StateInfoContainer;

my $ana = $analysis_a->fetch_by_logic_name($logic_name);

foreach my $name (@seqreg_name) {
	print STDOUT "Storing seq_region name $name in pipeline [$p_name]\n"
	  if ($verbose);
	my $slice =
	  $slice_a->fetch_by_region( $cs, $name, undef, undef, undef, $cs_version );
	if ( !$slice ) {
		warn(
"No seq_region [$name] found in database [$p_name] for coord_system [$cs] and cs_version [$cs_version]"
		);
		next;
	}
	my $target_projection = $slice->project($target_cs);
	foreach my $ct (@$target_projection) {
		my $target_sclice = $ct->to_Slice();
		my $target        =
		  $slice_a->fetch_by_region( $target_cs,
			$target_sclice->seq_region_name,
			undef, undef, undef, $target_cs_version );

		$state_info_container->store_input_id_analysis( $target->name(), $ana,
			'' );
		print STDOUT $target->name() . "\t"
		  . $ana->input_id_type . "\t"
		  . $ana->logic_name
		  . "\tstored\n"
		  if ($verbose);
	}
}

#!/software/bin/perl

=head1 NAME

list_input_ids.pl - generate region restriction info for rulemanager.pl

=head1 SYNOPSIS

list_input_ids.pl -dataset human -set chr14-04

=head1 DESCRIPTION

This script lists the input_ids of the seq_region name(s) provided.

=head1 OPTIONS

    -dataset            the species or dataset to connect to
    -set|name           the seq_region name you want to prime (could be used several times)
    -cs                 the coordinate system associated with the seq_region name (default: chromosome)
    -cs_version         the version of the coord system you want                  (default: Otter)
    -target_cs          the target coordinate system you want slices in           (default: contig)
    -target_cs_version  the version of the target coord system you want           (optional)

Switches taking no argument,
    -add-target-cs      include in the output an extra column showing the target_cs
    -verbose            (unused so far)
    -help      Displays script documentation with PERLDOC

=head1 CONTACT

Michael Gray B<email> mg13@sanger.ac.uk

=cut

use strict;
use warnings;

use Bio::Otter::Lace::Defaults;
use Bio::Otter::Lace::PipelineDB;

use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $dataset_name;
my $cs         = 'chromosome';
my $cs_version = 'Otter';
my $target_cs  = 'contig';
my $target_cs_version;
my $add_target_cs;
my $seqreg_name;
my $verbose;
my $help = 0;
Bio::Otter::Lace::Defaults::do_getopt(
        'dataset=s'           => \$dataset_name,
        'set|name:s'          => \$seqreg_name,
        'cs:s'                => \$cs,
        'cs_version:s'        => \$cs_version,
        'target_cs:s'         => \$target_cs,
        'target_cs_version:s' => \$target_cs_version,
        'add-target-cs!'      => \$add_target_cs,
        'verbose!'            => \$verbose,
        'h|help'              => \$help
);

if ($help) {
        exec( 'perldoc', $0 );
}

if ( !$dataset_name ) {
        throw("You must specify a dataset name (-dataset option");
}

if ( !$seqreg_name ) {
        throw("You must specify a seq_region name (-set option)");
}

# Client communicates with otter HTTP server
my $cl = Bio::Otter::Lace::Defaults::make_Client();

# DataSet interacts directly with an otter database
my $ds = $cl->get_DataSet_by_name($dataset_name);

my $otter_dba = $ds->get_cached_DBAdaptor;
my $pipe_dba = Bio::Otter::Lace::PipelineDB::get_pipeline_DBAdaptor($otter_dba);

#my $db = new Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor(
#        -host   => $host,
#        -user   => $user,
#        -pass   => $pass,
#        -port   => $port,
#        -dbname => $p_name
#);

my $slice_a              = $pipe_dba->get_SliceAdaptor;

my $slice = $slice_a->fetch_by_region( $cs, $seqreg_name, undef, undef, undef, $cs_version );
if ( !$slice ) {
    warn(
        "No seq_region [$seqreg_name] found in dataset [$dataset_name] for coord_system [$cs] and cs_version [$cs_version]"
    );
}
my $target_projection = $slice->project($target_cs);
foreach my $ct (@$target_projection) {
    my $target_slice = $ct->to_Slice();
    my $target        =
        $slice_a->fetch_by_region( $target_cs,
                                   $target_slice->seq_region_name,
                                   undef, undef, undef, $target_cs_version );

    if ($add_target_cs) {
        print STDOUT $target->name(), "\t", uc $target_cs, "\n";
    } else {
        print STDOUT $target->name(), "\n";
    }
}

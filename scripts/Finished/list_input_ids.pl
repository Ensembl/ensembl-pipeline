#!/usr/bin/env perl

=head1 NAME

list_input_ids.pl - generate region restriction info for rulemanager.pl

=head1 SYNOPSIS

 list_input_ids.pl -dataset human [ -purpose <label> ] -set chr14-04,chr15-06
 list_input_ids.pl -dataset human [ -purpose <label> ]      chr14-04 chr15-06

=head1 DESCRIPTION

This script lists the input_ids of the seq_region name provided.

Multiple seq_regions may be listed with commas or separate words, or
by calling the script for each one piped through C<sort -u>.

=head2 Duplicates

Contigs are listed at most once per input seq_region.

=head1 OPTIONS

    -dataset            the species or dataset to connect to
    -set|name           the seq_region name(s) you want to prime
    -cs                 the coordinate system associated with the seq_region name (default: chromosome)
    -cs_version         the version of the coord system you want                  (default: Otter)
    -target_cs          the target coordinate system you want slices in           (default: contig)
    -target_cs_version  the version of the target coord system you want           (optional)

    -purpose            Write input_ids to the input_id_purpose table.
                        If used, you must supply the purpose field.

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

use Bio::EnsEMBL::Pipeline::Analysis;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Pipeline::Utils::InputIDFactory;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub main {
    my $dataset_name;
    my $cs         = 'chromosome';
    my $cs_version = 'Otter';
    my $target_cs  = 'contig';
    my $target_cs_version;
    my $add_target_cs;
    my $seqreg_name;
    my $verbose;
    my $purpose;
    my $help = 0;
    Bio::Otter::Lace::Defaults::do_getopt(
            'dataset=s'           => \$dataset_name,
            'set|name:s'          => \$seqreg_name,
            'cs:s'                => \$cs,
            'cs_version:s'        => \$cs_version,
            'target_cs:s'         => \$target_cs,
            'target_cs_version:s' => \$target_cs_version,
            'purpose=s'           => \$purpose,
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

    my @seqregion = @ARGV;
    push @seqregion, split ',', $seqreg_name if defined $seqreg_name;
    if (!@seqregion) {
        throw("You must specify a seq_region name (-set option or trailing arguments)");
    }

    # Client communicates with otter HTTP server
    my $cl = Bio::Otter::Lace::Defaults::make_Client();

    # DataSet interacts directly with an otter database
    my $ds = $cl->get_DataSet_by_name($dataset_name);
    my $pipe_dba = $ds->get_pipeline_DBAdaptor(1);

    my $exit = 0; # set bits on fail, per slice
    foreach my $seqreg_name (@seqregion) {
        $exit |= do_slice
          ($dataset_name, $pipe_dba,
           $cs, $cs_version, $seqreg_name,
           $purpose,
           $add_target_cs, $target_cs, $target_cs_version);
    }

    return $exit;
}


sub do_slice {
    my ($dataset_name, $pipe_dba,
        $cs, $cs_version, $seqreg_name,
        $purpose,
        $add_target_cs, $target_cs, $target_cs_version) = @_;

    my $slice_a              = $pipe_dba->get_SliceAdaptor;

    # This table is experimental, but proved itself useful to me.
#  CREATE TABLE `input_id_purpose` (
#    `iip_id` int(10) unsigned NOT NULL AUTO_INCREMENT COMMENT 'useful for taking range slices of some work',
#    `input_id` varchar(100) NOT NULL,
#    `purpose` varchar(20) NOT NULL COMMENT 'some programmer-given label',
#    PRIMARY KEY (`iip_id`),
#    KEY `purp_iid` (`purpose`,`input_id`)
#  ) COMMENT='bring input_id list in, to ease mauling of jobs by SQL'
    my $purph = $pipe_dba->prepare
      (q{INSERT INTO input_id_purpose (input_id, purpose) VALUES (?,?)});

    my $slice = $slice_a->fetch_by_region( $cs, $seqreg_name, undef, undef, undef, $cs_version );
    if ( !$slice ) {
        warn "No seq_region [$seqreg_name] found in dataset [$dataset_name] ".
          "for coord_system [$cs] and cs_version [$cs_version] - skipping\n";
        return 4;
    }

    my %seen; # key = slice name
    my $target_projection = $slice->project($target_cs);
    foreach my $ct (@$target_projection) {
        my $target_slice = $ct->to_Slice();
        my $target =
          $slice_a->fetch_by_region( $target_cs,
                                     $target_slice->seq_region_name,
                                     undef, undef, undef, $target_cs_version );

        # once per chromosome
        next if $seen{ $target->name } ++;

        if (defined $purpose) {
            $purph->execute($target->name, $purpose);
        }

        if ($add_target_cs) {
            print STDOUT $target->name(), "\t", uc $target_cs, "\n";
        } else {
            print STDOUT $target->name(), "\n";
        }
    }
    return 0;
}


exit main();

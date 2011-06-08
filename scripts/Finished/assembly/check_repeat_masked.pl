#!/software/bin/perl

=head1 NAME

check_repeat_masked.pl  

=head1 SYNOPSIS

check_repeat_masked.pl [arguments]

=head1 DESCRIPTION

This script checks that the repeat mask analyses have been run on the selected sets.
This should be done before running align_non_ident.pl.

=head1 OPTIONS

[see align_by_component_identity for now]

=head1 CONTACT

Michael Gray B<email> mg13@sanger.ac.uk

=cut

use strict;
use warnings;

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../../..";
    unshift(@INC, "$Bin");
}

use Pod::Usage;
use Readonly;

use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;

my %logic_names = map { $_ => 1 } @$ANALYSIS_REPEAT_MASKING;
my $n_logic_names = scalar  @$ANALYSIS_REPEAT_MASKING;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

$support->param( 'verbose',     1 );    # throw away all that garbage
$support->param( 'interactive', 0 );    # stop that garbage from coming up

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'assembly=s',
    'altassembly=s',
    'altdbname=s',
    'chromosomes|chr=s@',
    'altchromosomes|altchr=s@',
    'ref_start=i', 
    'ref_end=i',
    'alt_start=i',
    'alt_end=i',
    );
$support->allowed_params( $support->get_common_params,
                          'assembly',
                          'altassembly',
                          'altdbname',
                          'chromosomes',
                          'altchromosomes',
    );

if ( $support->param('help') or $support->error ) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list( 'chromosomes', 'altchromosomes' );

#FIXME: params if appropriate
Readonly my $TARGET_CS  => 'contig';
Readonly my $TARGET_CS_VERSION;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params( 'assembly', 'altassembly' );

# first set connection parameters for alternative db
# both databases have to be on the same host, so we don't need to configure
# them separately
for my $prm (qw(host port user pass dbname)) {
    $support->param( "alt$prm", $support->param($prm) )
        unless ( $support->param("alt$prm") );
}

# FIXME dup with align_by_component_identity.pl

# reference database
my $R_dba = $support->get_database( 'ensembl', '' );
my $R_dbc = $R_dba->dbc;
my $R_sa  = $R_dba->get_SliceAdaptor;
my $R_asm = $support->param('assembly');
my $Ref_start = $support->param('ref_start') || undef;
my $Ref_end = $support->param('ref_end') || undef;

# database containing the alternative assembly
my $A_dba = $support->get_database( 'ensembl', 'alt' );
my $A_sa  = $A_dba->get_SliceAdaptor;
my $A_asm = $support->param('altassembly');
my $Alt_start = $support->param('alt_start') || undef;
my $Alt_end = $support->param('alt_end') || undef;

$support->log_stamped("Looping over chromosomes...\n");

my @R_chr_list = $support->param('chromosomes');
if ( !scalar(@R_chr_list) ) {
    @R_chr_list = $support->sort_chromosomes;

    if ( scalar( $support->param('altchromosomes') ) ) {
        die "AltChromosomes list is defined while Chromosomes list is not!";
    }
}

my @A_chr_list = $support->param('altchromosomes');
if ( !scalar(@A_chr_list) ) {
    @A_chr_list = @R_chr_list;
}
elsif ( scalar(@R_chr_list) != scalar(@A_chr_list) ) {
    die "Chromosome lists do not match by length";
}

my $ok = 1;

for my $i ( 0 .. scalar(@R_chr_list) - 1 ) {
    my $R_chr = $R_chr_list[$i];
    my $A_chr = $A_chr_list[$i];

    $support->log_stamped( "Chromosome $R_chr/$A_chr ...\n", 1 );

    # fetch chromosome slices
    my $R_slice =
        $R_sa->fetch_by_region( 'chromosome', $R_chr, $Ref_start, $Ref_end, undef, $R_asm );
    $support->log($R_slice->seq_region_name." ".$R_slice->start." -> ".$R_slice->end."\n", 2);

    my $r_ok = check_repeat_analyses($R_slice);
    $ok &&= $r_ok;

    my $A_slice =
        $A_sa->fetch_by_region( 'chromosome', $A_chr, $Alt_start, $Alt_end, undef, $A_asm );
    $support->log($A_slice->seq_region_name." ".$A_slice->start." -> ".$A_slice->end."\n", 2);

    my $a_ok = check_repeat_analyses($A_slice);
    $ok &&= $a_ok;
}

$support->log_warning("\n*** Detected missing repeat masking ***\n") if not $ok;
$support->finish_log;

exit ($ok ? 0 : 1);

sub check_repeat_analyses {
    my ($slice, $analyses) = @_;

    my $input_ids = list_slice_input_ids($slice, $TARGET_CS, $TARGET_CS_VERSION);

    my $missing = 0;
    my $total   = 0;

    # FIXME - duplication with load_loutre_pipeline.pl
    my $meta_container = $slice->adaptor->db->get_MetaContainer();
    my ($pipe_param) = @{$meta_container->list_value_by_key('pipeline_db_rw_head')};
    my $state_info_container;
    if($pipe_param) {
        my $pipe_dba = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(eval $pipe_param);
        $state_info_container = $pipe_dba->get_StateInfoContainer;
    } else {
        throw("Missing meta key 'pipeline_db_rw_head'\n");
    }

    foreach my $id (@$input_ids) {
        ++$total;
        my $anals = $state_info_container->fetch_analysis_by_input_id($id);
        my @matches = map { $_->logic_name } grep { $logic_names{$_->logic_name} } @$anals;
        my $l_names = join(',', @matches) || '';
        if (scalar @matches == $n_logic_names) {
            $support->log_verbose("$id\tOK [$l_names]\n", 3);
        } else {
            $support->log_warning("$id\tNOT FOUND [$l_names]\n", 3);
            ++$missing;
        }
    }
    $support->log(sprintf("%s: %d / %d clones missing repeat analyses\n",
                          $slice->seq_region_name, $missing, $total), 2);
    return ($missing == 0);
}

sub list_slice_input_ids {
    my ($slice, $target_cs, $target_cs_version) = @_;

    my %input_ids = ();

    my $target_projection = $slice->project($target_cs);
    foreach my $ct (@$target_projection) {
        my $target_slice = $ct->to_Slice();
        my $target        =
            $slice->adaptor()->fetch_by_region( $target_cs,
                                                $target_slice->seq_region_name,
                                                undef, undef, undef, $target_cs_version );
        $input_ids{$target->name()} = 1;
    }

    return [ keys %input_ids ];
}

# EOF


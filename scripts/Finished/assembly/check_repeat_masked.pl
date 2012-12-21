#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/Finished/assembly/check_repeat_masked.pl,v $
# $Revision: 1.6 $

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
use lib "$Bin";

use AssemblyMapper::Support;

use Pod::Usage;
use Readonly;

use Bio::EnsEMBL::Utils::Exception qw(throw); # supply via AssemblyMapper::Support?
use Bio::EnsEMBL::Analysis::Config::General;

my %logic_names = map { $_ => 1 } @$ANALYSIS_REPEAT_MASKING;
my $n_logic_names = scalar  @$ANALYSIS_REPEAT_MASKING;

my $support = AssemblyMapper::Support->new(
    extra_options => [
        'ref_start=i',
        'ref_end=i',
        'alt_start=i',
        'alt_end=i',
    ],
    );

unless ($support->parse_arguments(@_)) {
    warn $support->error if $support->error;
    pod2usage(1);
}

#FIXME: params if appropriate
Readonly my $TARGET_CS  => 'contig';
Readonly my $TARGET_CS_VERSION;

$support->connect_dbs;

my $ok = $support->iterate_chromosomes(
    prev_stage =>     '10-align_identical',
    this_stage =>     '20-check_repeat_masked',
    worker =>         \&do_check_repeat_masked,
    do_all => 1, # find out about all input_ids at once
    );

$support->log_warning("\n*** Detected missing repeat masking ***\n") if not $ok;
if ($support->has_output_info) {
    my $oi = $support->output_info;
    # This is a nudge in the right direction for when we're mauling it
    # by hand.  WIBNI more integrated.
    $oi->{input_id_file} = join '',
      map {qq{$_\tCONTIG\n}} @{$oi->{input_id} };
    my $ana_id = $oi->{submitcontig};
    my @vals = map {qq{('$_', 'CONTIG', $ana_id, now(), '')}}
      @{ $oi->{input_id} };
    $oi->{insert} = q{/* Use this to insert SubmitContig jobs */
INSERT INTO input_id_analysis (input_id, input_id_type, analysis_id, created, runhost)  values
}.(join ",\n", @vals).";\n/* and tell the rulemanager */";
    $support->log_warning( $support->output_info_as_yaml );
}
$support->finish_log;

exit ($ok ? 0 : 1);

sub do_check_repeat_masked {
    my ($asp, $data) = @_;

    my $r_ok = check_repeat_analyses($asp->ref_slice, $asp->align_support);
    my $a_ok = check_repeat_analyses($asp->alt_slice, $asp->align_support);

    return ($r_ok and $a_ok);
}

sub check_repeat_analyses {
    my ($slice, $support) = @_;

    my $input_ids = list_slice_input_ids($slice, $TARGET_CS, $TARGET_CS_VERSION);

    my $missing = 0;
    my $total   = 0;

    my $dba = $support->get_pipe_db($slice->adaptor->db);
    $dba ||= $slice->adaptor->db;

    my $state_info_container = $dba->get_StateInfoContainer;

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
            push @{ $support->output_info->{input_id} }, $id;
        }
    }

    if ($missing) {
        my $anaA = $dba->get_AnalysisAdaptor;
        # stash it for the INSERT generator
        $support->output_info->{submitcontig} =
          $anaA->fetch_by_logic_name('SubmitContig')->dbID;
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


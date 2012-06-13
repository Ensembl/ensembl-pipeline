#! /software/bin/perl-5.12.2

use strict;
use warnings;

use YAML 'Dump';

=head1 NAME

postqc_logfilter.pl - post-processing of Vega QC log

=head1 DESCRIPTION

This was supposed to be a quick hack to rearrange the QC log output.
It may not be worth maintaining.

The aim is to collect together all problems of one type, so the
classes and size of problems can be evaluated more quickly.

=head1 AUTHOR

Matthew Astley mca@sanger.ac.uk

=cut

sub main {
    my $MEMSTAMP_RE = qr{\[\d{4}-\d{2}-\d{2} [0-9:]{8}, mem +\d+\]};
    my $W = qr{^WARNING:};
    my $CWHO = qr{\[ ?chromosome = (\S+?)\. (?:gene|transcript) author = (\S+?)\] ?};
    my $DS = qr{\[(\d{2}/\d{2}/\d{4})\]};
    my $GG = qr{(\S+ \(\S+\))};

    my @re =
      (# patterns keyed with _* are administrative, internal state or report layout
       _head_script => qr{^(Script|Date|User|Parameters):( .+)?$},
       _head_kvp_head => qr{^ {4}PARAMETER {12}VALUE\s*$|^ {4}-{71}$},
       _head_kvp => qr{^ {3}\S+( *\t+\S+)?$},
       _blank => qr{^\s*$},
       _foot => qr{^All done\. \d+ warnings\. Runtime \d+h \d+min \d+sec $MEMSTAMP_RE$},

       # state changes are tracked
       _state_chr_species => qr{^Chromosome: (\S+) Species: (\S+) $MEMSTAMP_RE$},
       _state_chr => qr{^Checking chromosome (\S+) $MEMSTAMP_RE$},
       _state_species => qr{^={30}Tests for: (\S+)={30}$},
       _state_test => qr{^-{30}Test: (\S+)-{30}$},

       # showing test internal state?
       _inprog => qr{^Ignoring gene (OTT\S+) since 'in progress' \(author = .*\)$},
       _readthrough => qr{^ Ignoring some tests for transcript $GG since it has a remark of readthrough$},
       _chr_hid => qr{^chromosome chromosome:Otter:\S+ is hidden$},
       _biotype_new => qr{$W Found new gene biotype \S+},
       _species_hashref => qr{$W  Can't match species name from database \S+ with 'regions' hashref. Exiting$},


       # VQC._* ignored as biological issues
       VQCT_pot_selC => qr{^VQCT_pot_selC POTENTIAL SELENO },
       VQCT_CDS_remark => qr{$W VQCT_CDS_remark (\S+) has (.*?) 'CDS (start|end) not found' remark (.*?) $CWHO$},
       VQCT_no_stop => qr{$W VQCT_no_stop NO STOP: (\S+) pig gene (\S+) has a (\S+) transcript $GG that has a translation without a stop codon\. $CWHO$},
       VQCT_internal_stop => qr{$W VQCT_internal_stop INTERNAL STOPS HAVANA: Transcript (\S+?) .* from gene (\S+) has non '\w+' stop codons.*?  $DS $CWHO\.\)$},
       VQCG_NPC_CDS => qr{$W VQCG_NPC_CDS MISCLASSIFIED GENE: ig_segment gene $GG has translations\. $CWHO$},
       VQCT_biotype_inconst_CDS => qr{$W VQCT_biotype_inconst_CDS MISCLASSIFIED: (\S+) transcript $GG doesn't have a translation\. $CWHO$},
       VQCG_biotype_CDS_mismatch => qr{$W VQCG_biotype_CDS_mismatch MISCLASSIFIED TRANS: (\S+) transcript $GG has a translation\. $CWHO$},

       VQCG_PC_NO_CDS => qr{$W VQCG_PC_NO_CDS MISCLASSIFIED GENE: protein coding gene $GG has no translations\. $CWHO$},
       VQCG_transcript_gaps => qr{$W VQCG_transcript_gaps GAP between transcripts found in gene $GG\. $CWHO$},
       VQCG_inapp_NMD => qr{$W VQCG_inapp_NMD INCORRECT NMD TRANSCRIPT: (\S+) gene $GG has an NMD transcript \(not readthrough\)\. $CWHO$},

       VQCT_zero_length_CDS => qr{$W VQCT_zero_length_CDS POTENTIALLY SHORT - $GG has a translation reported as being of (\d+) aa length\. $CWHO$},
       VQCT_short_CDS => qr{$W VQCT_short_CDS POTENTIALLY SHORT - $GG has a translation that is (\d+) and has no CDS start/end not found remark\. $CWHO$},
       VQCT_short_exon => qr{$W VQCT_short_exon SHORT exon - $GG has an exon (\S+) with a length of (\d+ bp)\. $CWHO$},
       VQCT_short_intron => qr{$W VQCT_short_intron SHORT intron - $GG has an intron between exon (\S+) and (\S+) with a length of (\d bp)\. $CWHO$},
       VQCx_long => qr{$W VQC[GT]_long_(gene|transcript) LONG (gene|transcript) - $GG is (\d+)\. $CWHO$},
       VQCG_strand_mismatch => qr{$W VQCG_strand_mismatch STRAND mismatch between gene $GG and its transcripts\. $CWHO$},

       VQCG_no_description => qr{$W VQCG_no_description gene - $GG has no description\. $CWHO$},
       VQCG_duplicated_name => qr{$W VQCG_duplicated_name DUPLICATED - genes (\S+) and (\S+) share a name \((\S+)\) and are prob duplicates$},

       # bio?
       name_VQCG_duplicated_name_region => qr{$W VQCG_duplicated_name DUPLICATED - havana genes (\S+) and (\S+) share coordinates \((\d+:\d+)\)$},

       # issues to fix
       coords_VQCT_CDS_STARTgtEND => qr{$W VQCT_CDS_STARTgtEND NEGATIVE: $GG starts \((\d+)\) after it ends \((\d+)\)\. $CWHO$},
       VQCG_duplicated_root_name => qr{$W VQCG_duplicated_root_name SPLIT: havana transcript base name (\S+) is found in genes $GG and $GG$},
       VQCG_multi_root_name => qr{$W VQCG_multi_root_name DUPLICATED: havana gene $GG has transcripts with more than one name root $CWHO$},

       name_VQCG_wrong_name => qr{$W VQCG_wrong_name UNEXPECTED extension for havana GENE (\S+)$},
       name_VQCT_wrong_name => qr{$W VQCT_wrong_name UNEXPECTED name for havana transcript $GG $CWHO$},
       name_VQCT_wrong_name___new => qr{$W VQCT_wrong_name UNEXPECTED new format extension for havana (\S+) $CWHO$},
       name_VQCT_wrong_name___old => qr{$W VQCT_wrong_name UNEXPECTED old format extension for havana (\S+) $CWHO$},

       # bio_* are ignored
       bio_biotype_comb => qr{$W Combination of gene biotype '(\S+)' and transcript biotype '(\S+)' is not allowed \((\d+) in total\)\. Cases are (.*)$},
      );

    my %fnum =
      (# keyed regexp's capture indices, to their function
       VQCG_PC_NO_CDS => { g_gg => 0, chr => 1 },
       VQCT_biotype_inconst_CDS => { t_gg => 1, chr => 2 },

       # chr from state should be right for these (seen from the layout of the raw log)
       name_VQCG_wrong_name => { g_name => 0 },
       VQCG_duplicated_root_name => {},
       name_VQCT_wrong_name => { t_gg => 0, chr => 1 },
       name_VQCT_wrong_name___new => { t_name => 0, chr => 1 },
       name_VQCT_wrong_name___old => { t_name => 0, chr => 1 },
      );

    my %re;
    while (my ($k, $v) = splice @re, 0, 2) {
        die "No regexp for $k" unless $v;
        die "Dup regexp for $k" if $re{$k};
        $re{$k} = $v;
    }

    my %count;  # keys from %re, values count regexp hits
    my %state;  # keys from %re, values written
    my %unused; # keys from %re; a set
    my %hit; # keys from %re, values are match captures

  LINE: while (<>) {
        my $M = match_keyed($_, \%re);

        warn "NO MATCH: $_" unless keys %$M;

        my (@unused, $skip, $push);
        while (my ($mk, $mv) = each %$M) {
            $count{$mk} ++;
            $skip = 1 if $mk =~ m{^_}; # process noise
            $skip = 1 if $mk =~ m{^(VQC[TGx]|bio)_}; # biological issue

            if ($mk =~ m{^_state_(.*)$}) {
                $state{$1} = $mv;
            } elsif ($mk =~ m{^(coords|name)_} || exists $fnum{$mk}) {
                $push = $mk;
                $skip = 0;
            } else {
                $unused{$mk} ++ unless $skip;
            }
        }
        die "match breakage" if $push && $skip;
        next LINE if $skip;

        chomp;
        if ($push) {
            push @{ $hit{$push} },
              { TEXT => $_, LINE => $.,
                %$M,
                state => { %state },
              };
        }
    }

    print qq{-*- org -*-\n
This file format is most usefully viewed with Emacs.

Move to a line and press TAB to show headings / show all / close\n};

    my $find_chr = sub {
        my ($k, $v) = @_; # key,value from %hit
        my $info = $fnum{$k};
        die "Source of 'chr' unknown for regex-key $k" unless defined $info;
        my $chr_txt = $info->{chr} ? $v->{$k}->[ $info->{chr} ] : undef;
        my $chr_state = $v->{state}->{chr}->[0];
## seen some, but they are OK (due to not clearing state when test changes)
#  ...they need checking per type
#            if (defined $chr_txt && $chr_txt ne $chr_state) {
#                # sanity check this, because not every line tells its chr
#                die Dump({ fail => 'chr_txt and chr_state mismatch',
#                           k=> $k, v => $v,
#                           chr_txt => $chr_txt, chr_state => $chr_state });
#            }
        return $chr_txt || $chr_state;
    };

    print "\n* Summary\nHit types starting with '_' are internal stuff.\n";
    my %chr_count; # key = regexp, value = { chr => hitcount }
    while (my ($k, $info) = each %fnum) {
        $chr_count{$k} = {};
        foreach my $v (@{ $hit{$k} }) {
            $chr_count{$k}{ $find_chr->($k, $v) } ++;
            $chr_count{$k}{total} ++;
        }
    }
    print Dump({ line_hit_counts => \%count, per_chr => \%chr_count });

    print "\n* Issues\n";
    my %dump; # key = **heading, value = \@text
    foreach my $k (keys %hit) {
        foreach my $v (@{ $hit{$k} }) {
            my $section =
              ($fnum{$k}
#               ? sprintf('** %-30s on %10s', $k, $find_chr->($k, $v))
               ? sprintf('** %10s, %30s', $find_chr->($k, $v), $k)
               : sprintf('**  [all] %s', $k));
            push @{ $dump{$section} },
              sprintf qq{*** %s\n%s\n}, (join " ", @{ $v->{$k} }), $v->{TEXT};
        }
    }
    nested_dump(\%dump);
}


sub nested_dump {
    my ($dumph) = @_;
    foreach my $H (sort keys %$dumph) {
        printf "%-25s (%4d)\n", $H, scalar @{ $dumph->{$H} };
        print @{ $dumph->{$H} };
    }
    %$dumph = ();
    return ();
}

sub match_keyed {
    my ($txt, $re) = @_;
    my %hit;
    while (my ($k, $v) = each %$re) {
        my @match = ($txt =~ $v);
        $hit{$k} = \@match if @match;
    }
    return \%hit;
}

main();

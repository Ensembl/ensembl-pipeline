#!/software/bin/perl

use strict;
use warnings;
use Getopt::Long;

# This is an expansion of a one-liner, which may be useful on
# subsequent runs of transfer_annotation.pl

sub main {
    my %opt = (verbose => 0, quiet => 0);
    die "Syntax: $0 [--quiet] [--verbose] [ <logfile>... or STDIN]\n" unless
      GetOptions(\%opt, 'verbose|v!', 'quiet|q!');

    my $date_re = qr{[-0-9]{10} [:0-9]{8}};
    my $biotype_re = qr{(?:[35]')?\w+};

    my @skip =
      (qr{^\*Changes observed$},
       qr{^\S+\ now\ has\ changed\ main\ key\ (
              seq_region_name-seq_region_start-seq_region_end-seq_region_strand-phase-end_phase|
              seq_region_name-seq_region_start-seq_region_end-seq_region_strand-biotype-status-exon_count-description-evidence_count-attrib_string|
              start_exon_vega_hashkey-end_exon_vega_hashkey-tl_start-tl_end)$}x,
       qr{^ (Before| After)-key: \w+[-0-9]+(-\w+)?(-\w+[-0-9]+)?(-\w+=[^=]+)*$},
#      1  Before-key: chr1-03-58468844-58480910--1-nonsense_mediated_decay-UNKNOWN-11--1-cds_end_NF=0-cds_start_NF=0-hidden_remark=frame-shifted-mRNA_end_NF=0-mRNA_start_NF=0-name=AC121091.1-005
       qr{^$},

       # Multiline warning. MSG is not excluded here.
       qr{^(-{20} WARNING -{22}|-{51})$},
       qr{^(FILE|CALLED BY): Bio/EnsEMBL/Transcript\.pm\s+LINE: \d+$},
       qr{^Ensembl API version = \d+$},

       # Already
       qr{^SKIP (TATA-box|polyA_signal|polyA_site|pseudo_polyA|heptamer|nonamer|\d+ exons? phase \d from \S+( \(?O?T?T?\w*\)?)?|EUCOMM exon\(s\)|OTT\w{3}E\d+_\S+) \d+ \d+ -?1, already saved on chr\S+$},

       # Boilerplate top&tail
       qr{^Script: \S+:/\S+$},
       qr{^Date: $date_re$},
       qr{^User: \w+$},
       qr{^Parameters:$},
       qr{^ {4}(PARAMETER {12}VALUE *|-{71})$},
       qr{^ {3}[a-z_]{4,14} ?\t{1,2}\S+$},
       qr{^\t{3}\S+$}, # linewrap of key:value
       qr{^( {3}chromosomes\t|\t\t| {3}altchromosomes)\t(chr\w+-\d+)(, chr\w+-\d+)*,?$}, # special case key:value, and its linewrap
       qr{^ {8}(Ref|Alt): \S+, seq_region: \d+$},
       qr{^(Done|Looping over chromosomes\.\.)\. \[$date_re, mem \d+\]$},
      );

    push @skip,
      (qr{^MSG: Transcript contained trans splicing event$}, # XXX: investigate
       qr{^Found a truncated gene$}, # XXX: investigate
#       qr{^(INFO|WARNING):},
#       qr{ complex transfer of non-coding$},
#       qr{ cannot be transferr?ed on },
       qr{^CHANGED OTT...G\d+\.\d+$},
       qr{^-{43}$},
      ) if $opt{'quiet'};

    # More from "--verbose 1"
    push @skip,
      (qr{^Transferring gene OTT\w+ \w+ \d+ \d+$},
       qr{^\tOTT\w+ $biotype_re \d+ \d+ transferred successfully:$},
       qr{^\t{2}transfer to identical DNA$},
       qr{^\t{2}transfer to \d+ cDNA diffs and \d+ protein diffs$},
       qr{^\tSummary for OTT\w+ : \d+ out of \d+ transcripts? transferred$},
       qr{^GENE OTT\w+ \w+ successfully TRANSFERED \(chr[^: ]+:\d+-\d+ => chr[^: ]+:\d+-\d+\)$},
       qr{^========================================$},
      );

    my $skip = join '|', @skip;
    $skip = qr{$skip}o;

    while (<>) {
        if ($_ =~ $skip) {
            print "  |$_" if $opt{'verbose'};
        } else {
# s/\t/\\t/g; s/( +)$/length($1)." spaces"/e;
            print;
        }
    }
}

main();

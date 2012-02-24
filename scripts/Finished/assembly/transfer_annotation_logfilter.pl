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
    my @skip =
      (qr{^\*Changes observed$},
       qr{^\S+ now has changed main key (seq_region_name-seq_region_start-seq_region_end-seq_region_strand-phase-end_phase|start_exon_vega_hashkey-end_exon_vega_hashkey-tl_start-tl_end)$},
       qr{^ (Before| After)-key: \w+[-0-9]+(-\w+[-0-9]+)?$},
       qr{^$},

       # Multiline warning. MSG is not excluded here.
       qr{^(-{20} WARNING -{22}|-{51})$},
       qr{^(FILE|CALLED BY): Bio/EnsEMBL/Transcript\.pm\s+LINE: \d+$},
       qr{^Ensembl API version = \d+$},

       # Boilerplate top&tail
       qr{^Script: \S+:/\S+$},
       qr{^Date: $date_re$},
       qr{^User: \w+$},
       qr{^Parameters:$},
       qr{^ {4}(PARAMETER {12}VALUE *|-{71})$},
       qr{^ {3}[a-z_]{4,14} ?\t{1,2}\S+$},
       qr{^(Done|Looping over chromosomes\.\.)\. \[$date_re, mem \d+\]$},
      );

    push @skip,
      (qr{^MSG: Transcript contained trans splicing event$},
       qr{^CHANGED OTT...G\d+\.\d+$},
       qr{^-{43}$},
      ) if $opt{'quiet'};


    my $skip = join '|', @skip;
    $skip = qr{$skip}o;

    while (<>) {
        if ($_ =~ $skip) {
            print "  |$_" if $opt{'verbose'};
        } else {
            print;
        }
    }
}

main();

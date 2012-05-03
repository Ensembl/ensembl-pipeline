#!/software/bin/perl

use strict;
use warnings;
use Getopt::Long;
use CGI 'escapeHTML';


=head1 NAME

transfer_annotation_logfilter.pl - post-transfer assistance

=head1 SYNOPSIS

 # see serious stuff interactively
 transfer_annotation_logfilter.pl transfer.log | less

 # generate reports
 mkdir html
 transfer_annotation_logfilter.pl --html transfer.log

=head1 DESCRIPTION

This started off as a one-liner to find the serious stuff in the log
output from L<transfer_annotation.pl>, and has grown.

It has ballooned into a mess with complex internal state; but it may
be better to rewrite the generation of its input
(transfer_annotation.pl) than fix this script.

=cut


our $ICONS = '/icons/'; # cheerfully assume
# a) we display under Apache
# b) this directory is visible


sub main {
    my %opt = (verbose => 0, quiet => 0, html => 0, debug => 0);
    die "Syntax: $0 [--quiet] [--verbose] [ <logfile>... or STDIN]\n" unless
      GetOptions(\%opt, 'verbose|v!', 'quiet|q!', 'html|H!', 'debug!');

    my $date_re = qr{[-0-9]{10} [:0-9]{8}};
    my $biotype_re = qr{(?:[35]')?\w+};

    my $gene_div_re = qr{^========================================$}; # --verbose=1

    # These mark start of chromosome, end of header & start of footer, respectively
    my $chr_div_re = qr{^(?: {4}Chromosome (\S+)/(\S+) \.\.|(Looping) over chromosomes\.\.|(Done))\. \[$date_re, mem \d+\]$};

    my @skip =
      (qr{^\*Changes observed$},
       qr{^\S+\ now\ has\ changed\ main\ key\ (
              seq_region_name-seq_region_start-seq_region_end-seq_region_strand-phase-end_phase|
              seq_region_name-seq_region_start-seq_region_end-seq_region_strand-biotype-status-exon_count-description-evidence_count-attrib_string|
              start_exon_vega_hashkey-end_exon_vega_hashkey-tl_start-tl_end)$}x,
       qr{^ (Before| After)-key: \w+[-0-9]+(-\w+)?(-\w+[-0-9]+)?(-\w+=[^=]+)*$},
#      1  Before-key: chr1-03-58468844-58480910--1-nonsense_mediated_decay-UNKNOWN-11--1-cds_end_NF=0-cds_start_NF=0-hidden_remark=frame-shifted-mRNA_end_NF=0-mRNA_start_NF=0-name=AC121091.1-005
       qr{^$},
       qr{^DEBUG:}, # generic but probably not going to be used much in production!

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
       qr{^(Locking \S+ and \S+|Removing \S+ and \S+ Locks)$},
      );

    push @skip,
      (qr{^MSG: Transcript contained trans splicing event$}, # XXX: investigate
       qr{^Found a truncated gene$}, # XXX: investigate
#       qr{^(INFO|WARNING):},
#       qr{ complex transfer of non-coding$},
#       qr{ cannot be transferr?ed on },
       qr{^CHANGED OTT...G\d+\.\d+$}, qr{^-{43}$}, # from Bio::Vega::DBSQL::GeneAdaptor
      ) if $opt{'quiet'};

    # More from "--verbose 1"
    push @skip,
      (qr{^Transferring gene OTT\w+ \w+ \d+ \d+$},
       qr{^\tOTT\w+ $biotype_re \d+ \d+ transferred successfully:$},
       qr{^\t{2}transfer to identical DNA$},
       qr{^\t{2}transfer to 0 cDNA diffs and 0 protein diffs$},
       qr{^\t{2}transfer to 0 cDNA diffs$},
       qr{^\tSummary for OTT\w+ : \d+ out of \d+ transcripts? transferred$},
       qr{^GENE OTT\w+ \w+ successfully TRANSFERED \(chr[^: ]+:\d+-\d+ => chr[^: ]+:\d+-\d+\)$},
      );

    my $skip_re = join '|', @skip;
    $skip_re = qr{$skip_re}o;

    my ($skipped_last, @skipped_text, $showed_any);
    html_for_chr('header') if $opt{'html'};

    my $flush_skipped = sub {
        if ($showed_any) {
            print html_skipped(@skipped_text);
            @skipped_text = ();
        }
    };

    my @footer_info; # collect INFO: for all chromosomes

    while (<>) {
        my @div      =         $_ =~ $gene_div_re;
        my @skipping = @div || $_ =~ $skip_re;
        my @chr      =         $_ =~ $chr_div_re;

        if ($opt{'debug'}) {
            s/\t/\\t/g;
            s/( +)$/"{".length($1)." spaces}"/e;
        }

        if ($opt{'html'}) {
            if (@chr) {
                my ($old, $new, $end_header, $start_footer) = @chr;
                if ($end_header) {
                    # we only need to ensure output is flushed, so it's not on the first chromosome
                } elsif ($start_footer) {
                    $new = 'footer';
                } else {
                    # $new should be set
                    push @footer_info, [ $new ];
                }
                $flush_skipped->();
                html_for_chr($new) if defined $new; # switch output file
            }
            push @{ $footer_info[-1] }, $_ if /^INFO:/; # summary displays twice in total
            if (@skipping) {
                # (no verbose yet, probably too huge to display)
                push @skipped_text, $_;
            }
            $flush_skipped->() if @div;
            print next_item() if @div && $showed_any;
            if (!@skipping) {
                $showed_any = 1;
                $flush_skipped->();
                print html_link_quote($_);
            }

        } else {
            # plaintext
            if (@skipping) {
                print "  |$_" if $opt{'verbose'};
            } else {
                $showed_any = 1;
                print;
            }
        }

        $flush_skipped->() if eof; # not eof()

        # updating state
        if (@div) {
            @skipped_text = ();
            $showed_any = 0;
        }
        $skipped_last = scalar @skipping; # bool
    }

    extra_footer(@footer_info) if $opt{'html'};
}


sub extra_footer {
    my @footer_info = @_;
    foreach my $item (@footer_info) {
        my ($chr, @txt) = @$item;
        print next_item();
        print "Chromosome $chr\n";
        foreach my $ln (@txt) {
            $ln =~ s{^(.*\b)(\d+)/(\d+)\b}{sprintf("%s%d/%d%s", $1,$2,$3, spaced_perc(length($1.$2.$3),$2,$3))}eg;
            chomp $ln;
            my $html = escapeHTML($ln);
            $html = qq{<span class="failinfo">$html</span>} if $ln =~ /miss|split|skip/;
            print "$html\n";
        }
    }
}

sub spaced_perc {
    my ($len_before, $num, $denom) = @_;
    return '' unless $denom;
    my $perc = sprintf('(%.2f%%)', $num/$denom * 100);
    return (' ' x (60 - $len_before - length($perc))).$perc;
}


sub html_head {
    my ($title) = @_;
    my $htitle = escapeHTML($title || 'Extract of transfer log');
    return qq{<html><head>
 <title> $htitle </title>
 <style type="text/css">
  ol.log li { font-family: monospace; white-space: pre-wrap; border-bottom: 1px black dotted }
  .skip     { color: #8a8 }
  li:target { outline: 2px red solid }
  .nav { position: fixed; left: 0; top: 3em; border 1px blue dotted; }
  a.gene { font-weight: bold }
  a.again { text-decoration: line-through }
  .key { float: right; border: 3px solid pink }
  .failinfo { color: #800 }
 </style>
 <script type="text/javascript">
  function at_n() {
    var n = parseInt( location.href.split('#').pop() );
    return (n > 0 ? n : 1);
  }
  function go_n(n) {
    location.replace( location.href.split('#').shift() + '#' + n );
  }
 </script>
</head><body>
<div class="nav">
 <a href="javascript:go_n( at_n() - 1 )"> <img src="$ICONS/up.png"  alt="Previous" /> </a> <br />
 <a href="javascript:go_n( at_n() + 1 )"> <img src="$ICONS/down.png" alt="Next" />     </a>
</div>
<div class="key">
 <h2> Style Key </h2>
 <ol class="log">
  <li>numbered section per gene,
or for other parts of the log </li>
  <li>navigation arrows (left) jump between these
<span class="skip">this text is "not a problem",
it is left in for context</span>
with <a class="gene" href="">link to gene</a> and
the <a class="gene again" href="">same gene again</a></li>
 </ol>
</div>
<h1> $htitle </h1>
};
}
sub html_foot {
    return qq{</li></ol></body></html>\n};
}

{
    my $item_N = 0;
    sub next_item {
        my @out;
        push @out, $item_N ? '</li>' : '<ol class="log">';
        $item_N ++;
        push @out, qq{<li id="$item_N">};
        return @out;
    }

    my $prev_fh;
    sub html_for_chr {
        my ($chr) = @_;
        print html_foot() if $prev_fh;
        my $fn = "html/xfer_filtered.$chr.html";
        open my $fh, '>', $fn or die "Create $fn: $!";
        warn "Created $fn\n";
        select $fh;
        $prev_fh = $fh;
        $item_N = 0;
        print html_head($chr), next_item();
    }
    END {
        print html_foot() if $prev_fh;
    }
}

sub html_skipped {
    my @txt = @_;
    return (q{<span class="skip">},
            (map { html_link_quote($_) } @txt),
            q{</span>});
}

sub html_link_quote {
    my ($txt) = @_;

    my $html = escapeHTML($txt);
    $html =~ s{\b(OTT([A-Z]{3})?G\d+)\b}{annotrack_link($1, $2)}eg;

    return $html;
}

my %species = (qw( MUS mouse SUS pig ));
my %seen_gene; # key = OTTfooGnum
sub annotrack_link {
    my ($ottg, $species_tla) = @_;

    my $species = $species{$species_tla}
      || die "Unknown species TLA '$species_tla' in $ottg";

    my $cls = $seen_gene{$ottg}++ ? "again" : "";

    return qq{<a class="gene $cls" href="http://annotrack.sanger.ac.uk/$species/projects/show/$ottg">$ottg</a>};
}


main();

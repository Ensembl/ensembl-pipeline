#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-PIPELINE/scripts/Finished/pipeline_zap.pl,v $
# $Revision: 1.2 $
use strict;
use warnings;

use Sys::Hostname 'hostname';
use Cwd 'cwd';
use File::Slurp 'write_file';
use YAML 'Dump';

=head1 NAME

pipeline_zap.pl - "Stop!" button for pipeline

=head1 DESCRIPTION

Background process to allow someone else to kill off "my" pipeline
processes, before taking over.

Stop this process with SIGINT.

Stop the processes it sees by deleting the file it writes, to some
(presumably shared) directory.

=head1 CAVEATS

This is a quick hack.

=cut

sub main {
    die "Takes no args yet (QnD)" if @ARGV;
    my $pspat = qr{ruleman|\bperl.*dequeuer};
    my $sleep_for = 180;
    my $putdir = "/nfs/anacode";
# $putdir = "$ENV{HOME}"; $sleep_for=5; # debug

    my $u = (getpwuid($<))[0];
    my $herefile = "$putdir/pipeline_zap.$u.$$.yaml";
    my %info =
      (self => { host => hostname(), prog => $0, pid => $$,
                 cwd => cwd(), herefile => $herefile },
       run_start => scalar localtime($^T),
       sleep_for => $sleep_for,
       pspat => $pspat);

    $SIG{INT} = sub { $info{quit} = 'sigint '.localtime() };
    print "\nHerefile is $herefile\n\n";
    write_file($herefile, 'starting');

    while (!$info{quit}) {
        my @bjob = eval { proclist($pspat, qw( bjobs -w )) };
        @bjob = ("Error: $@", 'wait') unless @bjob;
        my @ps   = eval { proclist($pspat, qw( ps xww )) };
        @ps = ("Error: $@", 'wait') unless @ps;

        $info{next_wake} = localtime(time() + $sleep_for);
        $info{to_kill} = { bjobs => \@bjob, local_ps => \@ps };

        if (! -f $herefile) {
            push @{$info{stop}}, localtime().": $herefile deletion detected ";
            print "Zapping at ".localtime()."\n";
            do_kill(\@ps, qw( kill -INT ));
            do_kill(\@bjob, qw( bkill ));
        }
        if (@bjob < 2 && @ps < 2) {
            # only headers remain
            # (and no error, because that is two elements)
            push @{$info{stop}}, localtime().": Nothing running ";
            $info{quit} = 'procs gone';
        }

        write_file($herefile, Dump(\%info));
        sleep $sleep_for;
    }
    print Dump(\%info);
    write_file($herefile, Dump(\%info));
}

sub proclist {
    my ($pat, @cmd) = @_;
    die "scalar context" unless wantarray;
    open my $fh, '-|', @cmd or die "Pipe from @cmd: $!";
    my @out;
    push @out, scalar <$fh>; # header
    push @out, grep { $_ =~ $pat } <$fh>; # matching processes
    chomp @out;
    close $fh or die "Piped from @cmd: !=$! ?=$? after close";
    return @out;
}

sub do_kill {
    my ($procs, @cmd) = @_;
    my @zap = map { /^\s*(\d+)\s+/ ? ($1) : () } @$procs; # skip header
#unshift @cmd, qw( echo would ); # dry run
    print "do_kill: @cmd @zap\n";
    system(@cmd, @zap) && warn "Zap failed: exit code $?";
}

main();

#!/usr/bin/env perl

### pipe_error_summary.pl

use strict;
use warnings;
use Getopt::Long 'GetOptions';
use DBI;

{
    $| = 1;

    my $dbhost = undef;
    my $dbport = undef;
    my $dbuser = undef;
    my $dbpass = undef;
    my $dbname = undef;
    my $usage = sub{ exec('perldoc', $0) };
    GetOptions(
        'h|help!'   => $usage,
        'dbhost=s'  => \$dbhost,
        'dbport=i'  => \$dbport,
        'dbuser=s'  => \$dbuser,
        'dbpass=s'  => \$dbpass,
        'dbname=s'  => \$dbname,
        ) or $usage->();
    
    my $dbh = DBI->connect(
        "DBI:mysql:database=$dbname;host=$dbhost;port=$dbport",
        $dbuser, $dbpass,
        { AutoCommit => 1, RaiseError => 1 },
        );
    
    my $list_failures = $dbh->prepare(q{
        SELECT a.logic_name
          , js.status
          , j.stderr_file
        FROM analysis a
          , job j
          , job_status js
        WHERE a.analysis_id = j.analysis_id
          AND j.job_id = js.job_id
          AND js.is_current = 'y'
    });
    $list_failures->execute;
    
    my $too_many = 0;
    while (my ($ana_name, $status, $err_file) = $list_failures->fetchrow) {
        my $err_output = slurp_file($err_file);
        if ($err_output =~ /Too many connections/) {
            $too_many++;
            if ($too_many % 100) {
                print ".";
            }
            else {
                print "\nSaw $too_many jobs with 'Too many connections' error\n";
            }
            next;
        }
        else {
            print "\nANALYSIS: $ana_name; STATUS: $status\n$err_output;";
        }
    }
    print "\nSaw $too_many jobs with 'Too many connections' error\n";
}

sub slurp_file {
    my ($file) = @_;
    
    open(my $fh, '<', $file) or die "Can't read '$file'; $!";
    my $contents = join('', <$fh>);
    close $fh or die "Error reading '$file'; $!";
    return $contents;
}


__END__

=head1 NAME - pipe_error_summary.pl

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


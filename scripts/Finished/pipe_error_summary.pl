#!/usr/bin/env perl

### pipe_error_summary.pl

use strict;
use warnings;
use Getopt::Long 'GetOptions';
use DBI;
use Try::Tiny;

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
    

    my @ok_status = qw{
        CREATED           
        READING           
        RUNNING           
        SUBMITTED         
        WAITING           
        WRITING           
    };

    my @error_status = qw{
        AWOL              
        BUS_ERROR         
        FAILED            
        OUT_OF_MEMORY     
        SEGMENTATION_FAULT
    };

    my $error_placeholders = join(',', map { '?' } @error_status);

    my $sql = qq{
        SELECT a.logic_name
          , js.status
          , j.stdout_file
          , j.stderr_file
        FROM analysis a
          , job j
          , job_status js
        WHERE a.analysis_id = j.analysis_id
          AND j.job_id = js.job_id
          AND js.is_current = 'y'
          AND js.time > '2012-01-01'
          AND js.status IN ($error_placeholders)
    };
          # AND a.logic_name = 'Uniprot_raw_test'
          # AND a.logic_name = 'Uniprot_raw'
    warn $sql;

    my $list_failures = $dbh->prepare($sql);
    $list_failures->execute(@error_status);

    my $match_count = 0;
    my $err_count = 0;
    my $sep = 'x' x 40;
    while (my ($ana_name, $status, $out_file, $err_file) = $list_failures->fetchrow) {
        $err_count++;
        my $err_output = slurp_file($err_file);
        # $err_output .= slurp_file($out_file);
        # my $err_output = slurp_file($out_file);
        if (0 and $err_output =~ m{Missing hit_description entry}) {
            $match_count++;
            print STDERR "MATCHED: $ana_name\t$status\n";
        }
        else {
            print "$sep\nANALYSIS: $ana_name; STATUS: $status;\n$err_output\n";
        }
        
        # my $err_output = slurp_file($err_file);
        # if ($err_output =~ /Too many connections/) {
        #     $too_many++;
        #     if ($too_many % 100) {
        #         print ".";
        #     }
        #     else {
        #         print "\nSaw $too_many jobs with 'Too many connections' error\n";
        #     }
        #     next;
        # }
        # else {
        #     print "$sep\nANALYSIS: $ana_name; STATUS: $status\n$err_output;\n";
        # }
    }
    print "\nSaw $match_count matching jobs out of $err_count jobs\n";
}

sub slurp_file {
    my ($file) = @_;

    my $contents;
    try {
        local $/ = undef;
        open(my $fh, '<', $file) or die "Can't read '$file'; $!";
        $contents = <$fh>;
        close $fh or die "Error reading '$file'; $!";        
    }
    catch {
        warn $_;
    };
    return $contents || '';
}


__END__

=head1 NAME - pipe_error_summary.pl

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


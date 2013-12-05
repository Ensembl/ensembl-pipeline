#!/usr/bin/env perl
### fill_hit_desc ###

use warnings ;
use strict;
## CREATED FROM TEMPLATE
use Getopt::Long;
use Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Dfetch;
use Data::Dumper;

{
    my $help   = 0;
    my $dbhost = 'otterpipe1';
    my $dbname = 'pipe_human';
    my $dbuser = 'ottadmin';
    my $dbpass = undef;
    my $dbport = 3322;
    my $table  = 'dna';
    my $query  = undef;
    my $left   = 0;
    my $rev    = 0;
    my $ana_id = 0;
    GetOptions(
           'help|h'   => \$help,
           'dbhost=s' => \$dbhost,
           'dbname=s' => \$dbname,
           'dbuser=s' => \$dbuser,
           'dbpass=s' => \$dbpass,
           'dbport=s' => \$dbport,
           'table=s'  => \$table,
           'ana_id=s' => \$ana_id,
           'reverse'  => \$rev,
           'left'     => \$left,
           ) or useage();
    useage() if $help;

    my $db = Bio::EnsEMBL::Pipeline::DBSQL::Finished::DBAdaptor->new(
        -host   => $dbhost,
        -dbname => $dbname,
        -user   => $dbuser,
        -pass   => $dbpass,
        -port   => $dbport,
    );

    my $ids = [];
    if (! -t STDIN) {
        while (<>) {
            push(@$ids, split);
        }
    }
    elsif (@ARGV) {
        @$ids = @ARGV;
    } else {
        if($left){
            $query = "SELECT DISTINCT(f.hit_name)
                        FROM ${table}_align_feature f LEFT JOIN hit_description h
                          ON f.hit_name = h.hit_name
                       WHERE h.hit_name IS NULL
                       AND f.hit_name NOT LIKE '%\\_%'
                    ORDER BY f.hit_name";
        }
        elsif($ana_id){
            $query = "SELECT DISTINCT(f.hit_name)
                        FROM ${table}_align_feature f LEFT JOIN hit_description h
                          ON f.hit_name = h.hit_name
                       WHERE h.hit_name IS NULL
                       AND f.hit_name NOT LIKE '%\\_%'
                         AND f.analysis_id = '${ana_id}'
                    ORDER BY f.hit_name";
        }else{
            $query = "SELECT DISTINCT(f.hit_name) FROM ${table}_align_feature f ORDER BY f.hit_name";
        }
        $query .= " DESC" if $rev;

        my $sth = $db->prepare($query);
        $sth->execute();
        while (my ($id) = $sth->fetchrow){
            push(@$ids, $id);
        }
    }

    my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Dfetch->new( -type => $table );

    warn scalar(@$ids)." description(s) missing in $dbname\n";

    my $dbh = $db->dbc->db_handle;
    $dbh->begin_work;
    eval {
        my $chunk_size = 10_000;
        my @failed;
        for (my $i = 0; $i < @$ids; $i += $chunk_size) {
            my $j = $i + $chunk_size - 1;
            $j = $#$ids if $j > $#$ids;
            my $chunk = [@$ids[$i..$j]];
            my $fail = $seqfetcher->write_descriptions($db, $chunk, $chunk_size);
            if (@$fail) {
                push(@failed, @$fail)
            }
        }
        if (@failed) {
            warn "Failed to fetch:\n@failed\n";
            printf STDERR "Failed to fill in %d of %d hit descriptions requested\n", scalar(@failed), scalar(@$ids);
        }
        else {
            warn "Filled in hit descriptions for all IDs\n";
        }
        warn "FINISHED OK\n"
    };
    if (my $err = $@) {
        $dbh->rollback;
        throw("No hit_descriptions written: ", $err);
    }
    else {
        $dbh->commit;
    }
}

sub useage{
    exec('perldoc', $0);
    exit;
}
##END TEMPLATE##

=pod

=head1 NAME fill_hit_desc



=head1 AUTHOR



=cut

#!/software/bin/perl

=head1 NAME

fix_overlaps.pl - remove overlapping mappings between two closely related
assemblies

=head1 SYNOPSIS

fix_overlaps.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  --assembly=ASSEMBLY                 assembly version ASSEMBLY
  --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional arguments:

  --chromosomes, --chr=LIST           only process LIST chromosomes
  --altchromosomes, --altchr=LIST     supply alternative chromosome names (the two lists must agree)

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script removes overlapping mappings that were generated by the code in
align_nonident_regions.pl. Mappings are trimmed so that no overlaps are present
in the assembly table, because such overlaps may break the AssemblyMapper when
projecting between the two assemblies.

It also merges adjacent assembly segments which can result from alternating
alignments from clone identity and blastz alignment.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two assemblies
is done by a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

modified by Mustapha Larbaoui <ml6@sanger.ac.uk>

=head1 CONTACT

Please post comments/questions to Anacode
<anacode-people@sanger.ac.uk>

=cut

use strict;
use warnings;
# no warnings 'uninitialized';

use FindBin qw($Bin);
use lib "$Bin";

use Pod::Usage;
use AssemblyMapper::Support;

$| = 1;

my $support = AssemblyMapper::Support->new();

unless ($support->parse_arguments(@_)) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->connect_dbs;
my $dbc = $support->ref_dbc;

my $assembly    = $support->ref_asm;
my $altassembly = $support->alt_asm;

my $select = qq(
  SELECT a.*
  FROM assembly a, seq_region sr1, seq_region sr2,
       coord_system cs1, coord_system cs2
  WHERE a.asm_seq_region_id = sr1.seq_region_id
  AND a.cmp_seq_region_id = sr2.seq_region_id
  AND sr1.coord_system_id = cs1.coord_system_id
  AND sr2.coord_system_id = cs2.coord_system_id
  AND cs1.name = 'chromosome'
  AND cs2.name = 'chromosome'
  AND cs1.version = '$assembly'
  AND cs2.version = '$altassembly'
  AND sr1.name = ?
  AND sr2.name = ?
  ORDER BY a.asm_start
);

my $sth_select = $dbc->prepare($select);

my $sql_insert = qq(INSERT INTO assembly VALUES (?, ?, ?, ?, ?, ?, ?));
my $sth_insert = $dbc->prepare($sql_insert);

my $fmt1 = "%10s %10s %10s %10s %3s\n";

my $ok = $support->iterate_chromosomes(
    prev_stage =>     '30-align_nonident',
    this_stage =>     '40-fix_overlaps',
    worker =>         \&do_fix_overlaps,
    callback_data => {
        dbc        => $dbc,
        select_sth => $sth_select,
        insert_sth => $sth_insert,
    },
    );

$sth_select->finish;
$sth_insert->finish;

$support->finish_log;
exit ($ok ? 0 : 1);

sub do_fix_overlaps {
    my ($asp, $data) = @_;
    my $db_handles = $data;

    my $R_chr = $asp->ref_chr;
    my $A_chr = $asp->alt_chr;

    my @rows  = ();

    $db_handles->{select_sth}->execute( $R_chr, $A_chr );

    # do an initial fetch
    my $last = $db_handles->{select_sth}->fetchrow_hashref;

    # skip seq_regions for which we don't have data
    unless ($last) {
        $support->log( "No $R_chr/$A_chr mappings found. Skipping.\n", 1 );
        return 1;
    }

    push @rows, $last;

    my $i = 0;
    my $j = 0;
    my $k = 0;

    my $r; # declare outside loop or redo won't work as anticipated

  SEGMENT: while ( $last and ( $r = $db_handles->{select_sth}->fetchrow_hashref ) ) {

        # merge adjacent assembly segments (these can result from alternating
        # alignments from clone identity and blastz alignment)
        if (    $last->{'asm_end'} == ( $r->{'asm_start'} - 1 )
            and $last->{'cmp_end'} == ( $r->{'cmp_start'} - 1 )
            and $last->{'ori'}  == $r->{'ori'} )
        {

            $j++;

            # debug warnings
            $support->log_verbose(
                'merging - last:   '
                  . sprintf( $fmt1,
                    map { $last->{$_} }
                      qw(asm_start asm_end cmp_start cmp_end ori) ),
                1
            );
            $support->log_verbose(
                'this: '
                  . sprintf( $fmt1,
                    map { $r->{$_} }
                      qw(asm_start asm_end cmp_start cmp_end ori) ),
                1
            );

            # remove last row
            pop(@rows);

            # merge segments and add new row
            $last->{'asm_end'} = $r->{'asm_end'};
            $last->{'cmp_end'} = $r->{'cmp_end'};
            push @rows, $last;

            next SEGMENT;
        }

      # bridge small gaps (again, these can result from alternating alignments
      # from clone identity and blastz alignment). A maximum gap size of 10bp is
      # allowed
        my $asm_gap = $r->{'asm_start'} - $last->{'asm_end'} - 1;
        my $cmp_gap = $r->{'cmp_start'} - $last->{'cmp_end'} - 1;

        if ( $asm_gap == $cmp_gap and $asm_gap <= 10 and $asm_gap > 0 ) {

            $k++;

            # debug warnings
            $support->log_verbose(
                'bridging - last: '
                  . sprintf( $fmt1,
                    map { $last->{$_} }
                      qw(asm_start asm_end cmp_start cmp_end ori) ),
                1
            );
            $support->log_verbose(
                'this:            '
                  . sprintf( $fmt1,
                    map { $r->{$_} }
                      qw(asm_start asm_end cmp_start cmp_end ori) ),
                1
            );

            # remove last row
            pop(@rows);

            # merge segments and add new row
            $last->{'asm_end'} = $r->{'asm_end'};
            $last->{'cmp_end'} = $r->{'cmp_end'};
            push @rows, $last;

            next SEGMENT;
        }

        # look for overlaps with last segment
        if ( $last->{'asm_end'} >= $r->{'asm_start'} ) {

            $i++;

            # debug warnings
            $support->log_verbose(
                'last:   '
                  . sprintf( $fmt1,
                    map { $last->{$_} }
                      qw(asm_start asm_end cmp_start cmp_end ori) ),
                1
            );
            $support->log_verbose(
                'before: '
                  . sprintf( $fmt1,
                    map { $r->{$_} }
                      qw(asm_start asm_end cmp_start cmp_end ori) ),
                1
            );

            # skip if this segment ends before the last one
            if ( $r->{'asm_end'} <= $last->{'asm_end'} ) {
                $support->log_verbose( "skipped\n\n", 1 );
                next SEGMENT;
            }

            my $overlap = $last->{'asm_end'} - $r->{'asm_start'} + 1;

            $r->{'asm_start'} += $overlap;

            if ( $r->{'ori'} == -1 ) {
                $r->{'cmp_end'} -= $overlap;
            }
            else {
                $r->{'cmp_start'} += $overlap;
            }

            $support->log_verbose(
                'after:  '
                  . sprintf( $fmt1,
                    map { $r->{$_} }
                      qw(asm_start asm_end cmp_start cmp_end ori) )
                  . "\n",
                1
            );

            # Now that the overlap's removed, this segment may be
            # adjacent to, or within bridging distance of,
            # the last segment.
            #
            # Therefore we restart the loop body. Since it doesn't
            # overlap any more, if nothing else, it'll fall through
            # to below anyway...
            #
            redo SEGMENT;

        }

        push @rows, $r;
        $last = $r;
    }

    $support->log( "$R_chr/$A_chr: Fixed $i mappings.\n",  1 );
    $support->log( "$R_chr/$A_chr: Merged $j mappings.\n", 1 );
    $support->log( "$R_chr/$A_chr: Bridged $k gaps.\n",    1 );

    # fix the assembly table
    if ( $support->param('dry_run') ) {
        $support->log("\nNothing else to do for a dry run.\n");
    }
    else {

        # delete all current mappings from the db and insert the corrected ones
        my $c = $db_handles->{dbc}->do(
            qq(
    DELETE a
    FROM assembly a, seq_region sr1, seq_region sr2,
         coord_system cs1, coord_system cs2
    WHERE a.asm_seq_region_id = sr1.seq_region_id
    AND a.cmp_seq_region_id = sr2.seq_region_id
    AND sr1.coord_system_id = cs1.coord_system_id
    AND sr2.coord_system_id = cs2.coord_system_id
    AND cs1.name = 'chromosome'
    AND cs2.name = 'chromosome'
    AND cs1.version = '$assembly'
    AND cs2.version = '$altassembly'
    AND sr1.name = '$R_chr'
    AND sr2.name = '$A_chr' )
        );
        $support->log(
            "\n$R_chr/$A_chr: Deleted $c entries from the assembly table.\n");

        # now insert the fixed entries
        foreach my $r (@rows) {
            $db_handles->{insert_sth}->execute(
                map { $r->{$_} }
                  qw(asm_seq_region_id cmp_seq_region_id
                  asm_start asm_end cmp_start cmp_end ori)
            );
        }

        $support->log( "$R_chr/$A_chr: Added "
              . scalar(@rows)
              . " fixed entries to the assembly table.\n" );
    }

    return 1;
}

# EOF


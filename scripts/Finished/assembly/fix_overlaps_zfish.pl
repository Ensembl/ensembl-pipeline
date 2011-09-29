#!/software/bin/perl

=head1 NAME

fix_overlaps_zfish.pl - remove overlapping mappings between two closely related
set of zfish assemblies, step 3

=head1 SYNOPSIS

fix_overlaps_zfish.pl [arguments]

Required arguments:

    --dbname, db_name=NAME              database name NAME
    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --from_assembly                     old assembly date
    --to_assembly                       new assembly date

Optional arguments:

    --from_cs_version                   coordinate system version, this option will overwrite from_assembly
    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)
    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)

    -v, --verbose=0|1                   verbose logging (default: true)
    -i, --interactive=0|1               run script interactively (default: false)
    -n, --dry=0|1, --dry_run=0|1        don't write results to database (default: false)
    -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script removes overlapping mappings that were generated by the code in
align_nonident_zfish.pl. Mappings are trimmed so that no overlaps are present
in the assembly table, because such overlaps may break the AssemblyMapper when
projecting between the two assemblies.

It also merges adjacent assembly segments which can result from alternating
alignments from clone identity and blastz alignment.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two sets of zfish assemblies
is done by a series of scripts. Please see scripts in

  ensembl-pipeline/scripts/Finished/assembly/

for a high-level description of this process, and POD in the individual scripts
for the details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

modified by Mustapha Larbaoui <ml6@sanger.ac.uk>

=head1 CONTACT

Please post comments/questions to Anacode
<anacode-people@sanger.ac.uk>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);
BEGIN {
    $SERVERROOT = "$Bin/../../../..";
    unshift(@INC, "$Bin");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

$support->param('verbose', 1);
$support->param('interactive', 0);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'from_cs_version=s','from_assembly=s', 'to_assembly=s'
);
$support->allowed_params( $support->get_common_params, 'from_cs_version=s','from_assembly=s', 'to_assembly=s' );

if ( $support->param('help') or $support->error ) {
    warn $support->error if $support->error;
    pod2usage(1);
}

my $from_cs_version = $support->param('from_cs_version');
my $from_assembly = $support->param('from_assembly');
my $to_assembly = $support->param('to_assembly');

throw("Must set from_assembly or from_cs_version and to_assembly parameters!\n") unless(($from_assembly || $from_cs_version) && $to_assembly);

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# database connection
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;
my $sa  = $dba->get_SliceAdaptor;
my $csa  = $dba->get_CoordSystemAdaptor;



my $select = qq(
  SELECT a.*
  FROM assembly a, seq_region sr1, seq_region sr2,
       coord_system cs1, coord_system cs2
  WHERE a.asm_seq_region_id = sr1.seq_region_id
  AND a.cmp_seq_region_id = sr2.seq_region_id
  AND sr1.coord_system_id = cs1.coord_system_id
  AND sr2.coord_system_id = cs2.coord_system_id
  AND cs1.name IN ('chromosome','scaffold')
  AND cs2.name IN ('chromosome','scaffold')
);
$select .= " AND cs1.version = '$from_cs_version' " if $from_cs_version;
$select .= qq(
  AND sr1.name = ?
  AND sr2.name LIKE '%$to_assembly'
  ORDER BY a.asm_start
);

my $sth_select = $dbh->prepare($select);

my $sql_insert = qq(INSERT INTO assembly VALUES (?, ?, ?, ?, ?, ?, ?));
my $sth_insert = $dbh->prepare($sql_insert);

my $fmt1 = "%10s %10s %10s %10s %3s\n";

$support->log_stamped("Looping over chromosomes/scaffold...\n");

my $chr_list = $sa->fetch_all('chromosome');
if($from_cs_version) {
    foreach my $cs (@{$csa->fetch_all}) {
        push @$chr_list, @{$sa->fetch_all($cs->name,$cs->version)} if $cs->version eq $from_cs_version;
    }
}

my $sr_name_to_chr;
foreach my $chr_slice (@$chr_list) {
    my ($chr) = $chr_slice->seq_region_name =~ /chr(.*)_/;
    $chr = $chr_slice->seq_region_name unless $chr;
    $sr_name_to_chr->{$chr_slice->seq_region_name} = $chr;
}

my @from_chrs;
if($from_cs_version) {
    @from_chrs =
    sort { $sr_name_to_chr->{$a->seq_region_name} cmp $sr_name_to_chr->{$b->seq_region_name} }
    grep ( $_->coord_system->version eq $from_cs_version, @$chr_list );
} else {
    @from_chrs =
    sort { $sr_name_to_chr->{$a->seq_region_name} cmp $sr_name_to_chr->{$b->seq_region_name} }
    grep ( $_->seq_region_name =~ /$from_assembly/, @$chr_list );
}

my @to_chrs =
  sort { $sr_name_to_chr->{$a->seq_region_name} cmp $sr_name_to_chr->{$b->seq_region_name} }
  grep ( $_->seq_region_name =~ /$to_assembly/, @$chr_list );

if( !$from_cs_version && scalar(@from_chrs) != scalar(@to_chrs) ) {
    throw("Chromosome lists do not match by length:\n[".
    join(" ",map($_->seq_region_name,@from_chrs))."]\n[".
    join(" ",map($_->seq_region_name,@to_chrs))."]\n");
}

# Check that the chromosome names match
for my $i ( 0 .. scalar(@from_chrs) - 1 ) {
    my $R_sr_name = $from_chrs[$i]->seq_region_name;
    my $A_sr_name = $to_chrs[$i] ? $to_chrs[$i]->seq_region_name : "undef";
    my $R_chr = $sr_name_to_chr->{$R_sr_name};
    my $A_chr = $sr_name_to_chr->{$A_sr_name};

    throw("chromosome names don't match $R_chr != $A_chr\n[".
    join(" , ",map($_->seq_region_name,@from_chrs))."]\n[".
    join(" , ",map($_->seq_region_name,@to_chrs))."]\n") unless
    ( $R_chr eq $A_chr || $from_cs_version);

    $support->log_verbose("$R_sr_name   =>  $A_sr_name\n");
}

my @R_chr_list = map $_->seq_region_name , @from_chrs;
my @A_chr_list = map $_->seq_region_name , @to_chrs;

for my $i ( 0 .. scalar(@R_chr_list) - 1 ) {
    my $R_chr = $R_chr_list[$i];
    my $A_chr = $A_chr_list[$i];
    my @rows  = ();

    $support->log_stamped( "Chromosome/scaffold $R_chr/$to_assembly ...\n", 1 );

    $sth_select->execute( $R_chr );

    # do an initial fetch
    my $last = $sth_select->fetchrow_hashref;

    # skip seq_regions for which we don't have data
    unless ($last) {
        $support->log( "No $R_chr/$to_assembly mappings found. Skipping.\n", 1 );
        next;
    }

    push @rows, $last;

    my $i = 0;
    my $j = 0;
    my $k = 0;

    while ( $last and ( my $r = $sth_select->fetchrow_hashref ) ) {

        # merge adjacent assembly segments (these can result from alternating
        # alignments from clone identity and blastz alignment)
        if (    $last->{'asm_end'} == ( $r->{'asm_start'} - 1 )
            and $last->{'cmp_end'} == ( $r->{'cmp_start'} - 1 )
            and $last->{'ori'}  == $r->{'ori'}
            and $last->{'ori'}  == 1
            and $last->{'cmp_seq_region_id'} == $r->{'cmp_seq_region_id'} )
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

            next;
        }

      # bridge small gaps (again, these can result from alternating alignments
      # from clone identity and blastz alignment). A maximum gap size of 10bp is
      # allowed
        my $asm_gap = $r->{'asm_start'} - $last->{'asm_end'} - 1;
        my $cmp_gap = $r->{'cmp_start'} - $last->{'cmp_end'} - 1;

        if (        $asm_gap == $cmp_gap
                and $asm_gap <= 10
                and $asm_gap > 0
                and $last->{'cmp_seq_region_id'} == $r->{'cmp_seq_region_id'} ) {

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

            next;
        }

        # look for overlaps with last segment
        if ( $last->{'asm_end'} >= $r->{'asm_start'}
             and $last->{'cmp_seq_region_id'} == $r->{'cmp_seq_region_id'} ) {

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
                next;
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
        }

        push @rows, $r;
        $last = $r;
    }

    $support->log( "$R_chr/$to_assembly: Fixed $i mappings.\n",  1 );
    $support->log( "$R_chr/$to_assembly: Merged $j mappings.\n", 1 );
    $support->log( "$R_chr/$to_assembly: Bridged $k gaps.\n",    1 );

    # fix the assembly table
    if ( $support->param('dry_run') ) {
        $support->log("\nNothing else to do for a dry run.\n");
    }
    else {

        # delete all current mappings from the db and insert the corrected ones
        my $delete = qq(
    DELETE a
    FROM assembly a, seq_region sr1, seq_region sr2,
         coord_system cs1, coord_system cs2
    WHERE a.asm_seq_region_id = sr1.seq_region_id
    AND a.cmp_seq_region_id = sr2.seq_region_id
    AND sr1.coord_system_id = cs1.coord_system_id
    AND sr2.coord_system_id = cs2.coord_system_id);
    $delete .= " AND cs1.version = '$from_cs_version' " if $from_cs_version;
    $delete .= qq(
    AND cs1.name IN ('chromosome','scaffold')
    AND cs2.name IN ('chromosome','scaffold')
    AND sr1.name = '$R_chr'
    AND sr2.name LIKE '%$to_assembly' );

        my $c = $dbh->do($delete);

        $support->log(
            "\n$R_chr/$to_assembly: Deleted $c entries from the assembly table.\n");

        # now insert the fixed entries
        foreach my $r (@rows) {
            $sth_insert->execute(
                map { $r->{$_} }
                  qw(asm_seq_region_id cmp_seq_region_id
                  asm_start asm_end cmp_start cmp_end ori)
            );
        }

        $support->log( "$R_chr/$to_assembly: Added "
              . scalar(@rows)
              . " fixed entries to the assembly table.\n" );
    }
}

$sth_select->finish;
$sth_insert->finish;

# finish logfile
$support->finish_log;


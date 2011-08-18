#!/software/bin/perl

=head1 NAME

reverse_existing_mapping.pl - create a new mapping by reversing an existing one.

=head1 SYNOPSIS

reverse_existing_mapping.pl [arguments]

Required arguments:

    --dbname, db_name=NAME              database name NAME
    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --assembly=ASSEMBLY                 assembly version ASSEMBLY
    --altdbname=NAME                    alternative database NAME
    --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY
    --chromosomes, --chr=REF_CHR        chromosome or set name [just one!]
    --altchromosomes, --altchr=ALT_CHR  alternative chromosome or set name [just one!]
    --dstassembly=ASSEMBLY              destination assembly version ASSEMBLY
    --dstchromosome, --dstcr=DST_CHR    destination chromosome or set name [just one!]

Optional arguments:

    --conffile, --conf=FILE             read parameters from FILE

    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)

    -v, --verbose=0|1                   verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script is an optional part of a series of scripts to create a 
mapping between two assemblies. It creates a new mapping by reversing
an existing mapping.

The new mapping will be FROM the altchromosome TO the dstchromosome,
reversing the mapping from the ref chromosome to the alt chromosome.

See "Related files" below for an overview of the whole process.

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

Michael Gray <mg13@sannger/ac/il>

=head1 CONTACT

Please post comments/questions to Anacode
<anacode-people@sanger.ac.uk>

=cut

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin";

use Pod::Usage;

use Bio::EnsEMBL::Utils::Exception qw( throw ); # supply via AssemblyMapper::Support?
use AssemblyMapper::Support;

$| = 1;

my $support = AssemblyMapper::Support->new(
    extra_options     => [
        'dstassembly=s',
        'dstchromosome|dstchr=s',
    ],
    required_params   => [ 'dstassembly', 'dstchromosome', 'chromosomes', 'altchromosomes' ],
    single_chromosome => 1,
    );

unless ($support->parse_arguments(@_)) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->connect_dbs;

my $dest_chr      = $support->s_param('dstchromosome');
my $dest_assembly = $support->s_param('dstassembly');

my $dest_slice = $support->ref_sa->fetch_by_region(
    'chromosome',
    $dest_chr,
    undef,
    undef,
    undef,
    $dest_assembly,
    );
my $dest_seq_region_id = $support->ref_sa->get_seq_region_id($dest_slice);

# Pre-prepare assembly insert
# FIXME: duplication!
#
my $asm_insert_sth = $support->ref_dbc->prepare(qq{
    INSERT IGNORE INTO assembly (asm_seq_region_id, cmp_seq_region_id,
                                 asm_start, asm_end, cmp_start, cmp_end, ori)
    VALUES (?, ?, ?, ?, ?, ?, ? )
});
my $asm_select_sth = $support->ref_dbc->prepare(qq{
    SELECT asm_start, asm_end, cmp_start, cmp_end, ori
      FROM assembly
     WHERE asm_seq_region_id = ?
       AND cmp_seq_region_id = ?
});

# This will only do the single chr, but handles the session stuff for us
#
my $ok = $support->iterate_chromosomes(
    prev_stage    => '40-fix_overlaps',
    this_stage    => '45-reverse_existing_mapping',
    worker        => \&do_reverse_mapping,
    callback_data => {
        dest_chr_tag => "${dest_chr}[${dest_assembly}]",
        dest_sr_id   => $dest_seq_region_id,
        select_sth   => $asm_select_sth,
        insert_sth   => $asm_insert_sth,
        dry_run      => $support->s_param('dry_run'),
    },
    );

$support->log_stamped("-Done.\n\n"); # FIXME: should be in Assembly::Support.

# finish logfile
$support->finish_log;

exit 0;

sub do_reverse_mapping {
    my ($asp, $data) = @_;

    my $ref_sr_id = $asp->ref_seq_region_id;
    my $alt_sr_id = $asp->alt_seq_region_id;
    my $support   = $asp->align_support;

    my $dest_chr_tag   = $data->{dest_chr_tag};
    my $dest_sr_id     = $data->{dest_sr_id};
    my $asm_select_sth = $data->{select_sth};
    my $asm_insert_sth = $data->{insert_sth};
    my $write_db = not($data->{dry_run});

    $asm_select_sth->execute($ref_sr_id, $alt_sr_id);

    while (my $row = $asm_select_sth->fetchrow_hashref) {

        $support->log(sprintf("Adding mapping %s/%d [%9d-%9d] -> %s/%d [%9d-%9d] {%2d}%s\n",
                              $asp->alt_chr, $alt_sr_id,  $row->{cmp_start}, $row->{cmp_end},
                              $dest_chr_tag, $dest_sr_id, $row->{asm_start}, $row->{asm_end},
                              $row->{ori}, $write_db ? "" : " (dry run)",
                      ), 2);
        next unless $write_db;

        $asm_insert_sth->execute(
            $alt_sr_id,
            $dest_sr_id,
            $row->{cmp_start},
            $row->{cmp_end},
            $row->{asm_start},
            $row->{asm_end},
            $row->{ori},
            );
    }

    return 1;
}

# EOF

#!/usr/bin/env perl

=head1 NAME

add_extra_mappings.pl - add extra specified mappings between two assemblies

=head1 SYNOPSIS

add_extra_mappings.pl [arguments]

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
    --mapping                           [list of] mappings of form <ref_start>-<ref_end>:<alt_start>-<alt_end>

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

This script is an optional part of a series of scripts to create a mapping
between two assemblies. It adds extra mappings which will have been generated
by non-standard means.

Multiple mappings can be specified as a comma-separated list to a single
'--mapping' argument, or by specifying multiple '--mapping' arguments.

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
    extra_options     => [ 'mapping=s@' ],
    required_params   => [ 'mapping', 'chromosomes', 'altchromosomes' ],
    single_chromosome => 1,
    );

unless ($support->parse_arguments(@_)) {
    warn $support->error if $support->error;
    pod2usage(1);
}

my @mappings = process_mappings($support);
unless (@mappings) {
    $support->log_warning("No valid mappings so nothing to do.\n");
    exit 1;
}

$support->connect_dbs;

# Pre-prepare assembly insert
    my $sth = $support->ref_dbc->prepare(qq{
    INSERT IGNORE INTO assembly (asm_seq_region_id, cmp_seq_region_id,
                                 asm_start, asm_end, cmp_start, cmp_end, ori)
    VALUES (?, ?, ?, ?, ?, ?, ? )
    });

# This will only do the single chr, but handles the session stuff for us
#
my $ok = $support->iterate_chromosomes(
    prev_stage    => '30-align_nonident',
    this_stage    => '35-add_extra_mappings',
    before_stage  => '40-fix_overlaps',
    worker        => \&do_add_mappings,
    callback_data => {
        mappings => \@mappings,
        sth      => $sth,
        dry_run  => $support->s_param('dry_run'),
    },
    );

$support->log_stamped("-Done.\n\n"); # FIXME: should be in Assembly::Support.

# finish logfile
$support->finish_log;

exit 0;

sub do_add_mappings {
    my ($asp, $data) = @_;

    my $ref_sr_id = $asp->ref_seq_region_id;
    my $alt_sr_id = $asp->alt_seq_region_id;
    my $support   = $asp->align_support;

    my @mappings = @{$data->{mappings}};
    my $sth      = $data->{sth};
    my $write_db = not($data->{dry_run});

    foreach my $m (@mappings) {

        $support->log(sprintf("Adding mapping [%s-%s]->[%s-%s]%s\n",
                              $m->{ref}->{start}, $m->{ref}->{end},
                              $m->{alt}->{start}, $m->{alt}->{end},
                              $write_db ? "" : " (dry run)",
                      ), 2);
        next unless $write_db;

        $sth->execute(
            $ref_sr_id,
            $alt_sr_id,
            $m->{ref}->{start},
            $m->{ref}->{end},
            $m->{alt}->{start},
            $m->{alt}->{end},
            1,
            );
    }

    return 1;
}

sub process_mappings {
    my $support = shift;

    $support->comma_to_list('mapping');
    my @mappings = $support->param('mapping');

    my @results;

  MAPPING: foreach my $m (@mappings) {

      my %mapping;

      my %ref_alt;
      @ref_alt{'ref', 'alt'} = split(':', $m);

      foreach my $ra (keys %ref_alt) {

          my %side;
          @side{'start', 'end'} = split('-', $ref_alt{$ra});

          unless ($side{start} and $side{end}) {
              $support->log_warning("Missing start or end coord for $ra in '$m'\n");
              next MAPPING;
          }

          $side{length} = $side{end} - $side{start} + 1;

          unless ($side{length} > 0) {
              $support->log_warning("Non-positive length for $ra in '$m'\n");
              next MAPPING;
          }

          $mapping{$ra} = \%side;
      }

      unless ($mapping{ref}->{length} == $mapping{alt}->{length}) {
          $support->log_warning("Lengths for ref and alt do not match in '$m'\n");
          next MAPPING;
      }

      push @results, \%mapping;

  } # MAPPING

    return @results;
}

# EOF

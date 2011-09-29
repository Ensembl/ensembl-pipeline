#!/software/bin/perl

=head1 NAME

new_align_session.pl

=head1 SYNOPSIS

new_align_session.pl [arguments]

=head1 DESCRIPTION

Set up align_session entries to allow progress locking of an
alignment session.

Create tables if necessary.

Should be run before align_by_component_identity.

=head1 OPTIONS

[see align_by_component_identity for now]

=head1 CONTACT

Michael Gray B<email> mg13@sanger.ac.uk

=cut

use strict;
use warnings;

use Pod::Usage;

use FindBin qw($Bin);
use lib "$Bin";

use AssemblyMapper::Support;

my $a_support = AssemblyMapper::Support->new(
    extra_options => [
        'author=s',
        'comment=s',
    ],
    required_params => [qw(author)],
    );

unless ($a_support->parse_arguments(@_)) {
    warn $a_support->error if $a_support->error;
    pod2usage(1);
}

$a_support->connect_dbs;

my $ok = check_create_tables($a_support->ref_dbc);
exit(1) unless $ok;

$ok = $a_support->iterate_chromosomes(
    prev_stage =>     undef,
    this_stage =>     '00-new_align_session',
    create_session => 1,
    worker =>         \&check_start_session,
    callback_data =>  'hello',
    );

$a_support->log_warning("\n*** Couldn't create all required sessions ***\n") if not $ok;
$a_support->finish_log;

exit ($ok ? 0 : 2);

sub check_start_session {
    my ($asp, $data) = @_;

    # This is a no-op for new_align_session.
    # All the work is done in AssemblyMapper::Support::iterate_chromosomes().

    return 1;
}

sub check_create_tables {
    my ($dbc) = shift;

    my @tables = $dbc->db_handle->tables(undef, undef, '%', undef);

    unless (grep { /align_session/ } @tables) {
        if ($a_support->param('dry_run')) {
            $a_support->log("[dry_run] Here I would create the 'align_session' table for you\n");
        } else {
            $dbc->do( qq{
                CREATE TABLE IF NOT EXISTS align_session (
                  align_session_id  int(10) unsigned NOT NULL auto_increment,
                  ref_seq_region_id int(10) unsigned NOT NULL,
                  alt_seq_region_id int(10) unsigned NOT NULL,
                  alt_db_name       varchar(50),
                  author            varchar(50),
                  comment           text,

                  PRIMARY KEY (align_session_id)
                ) ENGINE=InnoDB
                }
                );
            $a_support->log_warning("Created align_session table\n");
        }
    }

    unless (grep { /align_stage/ } @tables) {
        if ($a_support->param('dry_run')) {
            $a_support->log("[dry_run] Here I would create the 'align_stage' table for you\n");
        } else {
            $dbc->do( qq{
                CREATE TABLE IF NOT EXISTS align_stage (
                  align_stage_id    int(10) unsigned NOT NULL auto_increment,
                  align_session_id  int(10) unsigned NOT NULL,
                  stage             varchar(50)      NOT NULL,
                  ts                timestamp        DEFAULT current_timestamp,
                  script            varchar(200)     NOT NULL,
                  parameters        text,

                  PRIMARY KEY (align_stage_id),
                  KEY align_session (align_session_id),
                  KEY align_session_stage (align_session_id, stage)
                ) ENGINE=InnoDB;
                }
                );
            $a_support->log_warning("Created align_stage table\n");
        }
    }

    return 1;
}

# EOF

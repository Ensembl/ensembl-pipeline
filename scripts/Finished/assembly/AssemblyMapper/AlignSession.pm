package AssemblyMapper::AlignSession;

use namespace::autoclean;
use Moose;
use MooseX::ClassAttribute;

use AssemblyMapper::AlignStage;

with 'AssemblyMapper::ClassUtils';

has id => (
    is       => 'ro',
    isa      => 'Int',
    writer   => '_set_id',
    );

has ref_seq_region_id => (
    is       => 'ro',
    isa      => 'Int',
    required => 1,
    );

has alt_seq_region_id => (
    is       => 'ro',
    isa      => 'Int',
    required => 1,
    );

has alt_db_name => (
    is       => 'ro',
    isa      => 'Maybe[Str]',
    );

has author => (
    is       => 'ro',
    isa      => 'Maybe[Str]',
    );

has comment => (
    is       => 'ro',
    isa      => 'Maybe[Str]',
    );

class_has Dbc => (
    is       => 'rw',
    isa      => 'Bio::EnsEMBL::DBSQL::DBConnection',
    trigger  => \&_Dbc_set,
    );

# Constructors
#
sub BUILD {
    my $self = shift;

    # Store the object if necessary
    unless ($self->id) {
        # Write it to db
        my $sth = $self->Dbc->prepare(qq(
            INSERT INTO align_session
                (ref_seq_region_id, alt_seq_region_id, alt_db_name, author, comment)
                VALUES (?, ?, ?, ?, ?)
            ));
        $sth->execute($self->ref_seq_region_id,
                      $self->alt_seq_region_id,
                      $self->alt_db_name,
                      $self->author,
                      $self->comment);
        my $id = $sth->{mysql_insertid};
        $self->_set_id($id);
    }

    return $self;
}

sub by_id {
    my ($self, $id) = @_;

    # Recover by id
    my $sth = $self->Dbc->prepare(qq(
        SELECT
          align_session_id AS id,
          ref_seq_region_id,
          alt_seq_region_id,
          alt_db_name,
          author,
          comment
        FROM align_session
        WHERE align_session_id = ?
        ));
    $sth->execute($id);

    my $details = $sth->fetchrow_hashref;
    return unless $details;

    return $self->new($details);
}

sub latest {
    my $self = shift;
    my $params = $self->validate_params(
        \@_,
        ref_seq_region_id => { isa => 'Int' },
        alt_seq_region_id => { isa => 'Int' },
        alt_db_name       => { isa => 'Maybe[Str]', optional => 1 },
        );

    my $alt_db_name_clause = 'AND alt_db_name IS NULL';

    my @bind_params = ($params->{ref_seq_region_id}, $params->{alt_seq_region_id});
    if ($params->{alt_db_name}) {
        $alt_db_name_clause = 'AND alt_db_name = ?';
        push @bind_params, $params->{alt_db_name};
    }

    # Recover latest entry by params
    my $sth = $self->Dbc->prepare(qq(
        SELECT
          align_session_id AS id,
          ref_seq_region_id,
          alt_seq_region_id,
          alt_db_name,
          author,
          comment
        FROM align_session
        WHERE
              ref_seq_region_id = ?
          AND alt_seq_region_id = ?
          ${alt_db_name_clause}
        ORDER BY align_session_id DESC
        ));
    $sth->execute(@bind_params);

    my $details = $sth->fetchrow_hashref;
    return unless $details;

    return $self->new($details);
}

# Methods
#

sub stage_has_run {
    my ($self, $stage) = @_;
    return AssemblyMapper::AlignStage->by_session_stage(
        align_session_id => $self->id,
        stage            => $stage,
        );
}

sub create_stage {
    my $self = shift;
    my $params = $self->validate_params(
        \@_,
        stage      => { isa => 'Str' }, # better val'n, see AlignStage.pm?
        previous   => { isa => 'Maybe[Str]' },
        before     => { isa => 'Maybe[Str]', optional => 1 },
        script     => { isa => 'Str' },
        parameters => { isa => 'Str',  optional => 1 },
        dry_run    => { isa => 'Bool', optional => 1 },
        force_stage=> { isa => 'Bool', optional => 1 },
        );

    my $warn_sub = $params->{force_stage} ? \&Carp::cluck : \&Carp::confess;

    my $stage = $params->{stage};
    &$warn_sub("Already run: '$stage'\n") if $self->stage_has_run($stage);

    my $previous = $params->{previous};
    if ($previous) {
        unless ($self->stage_has_run($previous)) {
            &$warn_sub("Prerequisite '$previous' has not been run for stage '$stage'\n");
        }
    }

    my $before = $params->{before};
    if ($before) {
        if ($self->stage_has_run($before)) {
            &$warn_sub("Following stage '$before' has already been run for stage '$stage'\n");
        }
    }

    return if $params->{dry_run};

    $params->{align_session_id} = $self->id;
    return AssemblyMapper::AlignStage->new($params);
}

sub _Dbc_set {
    my ($self, $dbc, $old_dbc) = @_;
    AssemblyMapper::AlignStage->Dbc($dbc);
}

 __PACKAGE__->meta->make_immutable;

1;

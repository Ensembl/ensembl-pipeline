package AssemblyMapper::AlignStage;

use namespace::autoclean;
use Moose;
use MooseX::ClassAttribute;

with 'AssemblyMapper::ClassUtils';

has id => (
    is       => 'ro',
    isa      => 'Int',
    writer   => '_set_id',
    );

has align_session_id => (
    is       => 'ro',
    isa      => 'Int',
    required => 1,
    );

has stage => (
    is       => 'ro',
    isa      => 'Str',          # FIXME: validate stage syntax
    required => 1,
    );

has ts => (
    is       => 'ro',
    isa      => 'Str',
    );

has script => (
    is       => 'ro',
    isa      => 'Str',
    required => 1,
    );

has parameters => (
    is       => 'ro',
    isa      => 'Str',
    );

class_has Dbc => (
    is       => 'rw',
    isa      => 'Bio::EnsEMBL::DBSQL::DBConnection',
    );

# Constructors
#
sub BUILD {
    my $self = shift;

    # Store the object if necessary
    unless ($self->id) {
        # Write it to db
        my $sth = $self->Dbc->prepare(qq(
            INSERT INTO align_stage
                (align_session_id, stage, script, parameters)
                VALUES (?, ?, ?, ?)
            ));
        $sth->execute($self->align_session_id,
                      $self->stage,
                      $self->script,
                      $self->parameters);
        my $id = $sth->{mysql_insertid};
        $self->_set_id($id);
        # FIXME: what about the timestamp?
    }

    return $self;
}

sub by_id {
    my ($self, $id) = @_;

    # Recover by id
    my $sth = $self->Dbc->prepare(qq(
        SELECT
          align_stage_id AS id,
          align_session_id,
          stage,
          ts,
          script,
          parameters
        FROM align_stage
        WHERE align_stage_id = ?
        ));
    $sth->execute($id);

    my $details = $sth->fetchrow_hashref;
    return unless $details;

    return $self->new($details);
}

sub latest_for_session {
    my ($self, $align_session_id) = @_;

    my $sth = $self->Dbc->prepare(qq(
        SELECT
          align_stage_id AS id,
          align_session_id,
          stage,
          ts,
          script,
          parameters
        FROM align_stage
        WHERE align_session_id = ?
        ORDER BY stage DESC
        ));
    $sth->execute($align_session_id);

    my $details = $sth->fetchrow_hashref;
    return unless $details;

    return $self->new($details);
}

sub by_session_stage {
    my $self = shift;
    my $params = $self->validate_params(
        \@_,
        align_session_id => { isa => 'Int' },
        stage            => { isa => 'Str' }, # better val'n, see above?
        );

    my $sth = $self->Dbc->prepare(qq(
        SELECT
          align_stage_id AS id,
          align_session_id,
          stage,
          ts,
          script,
          parameters
        FROM align_stage
        WHERE
              align_session_id = ?
          AND stage            = ?
        ));
    $sth->execute($params->{align_session_id}, $params->{stage});

    my $details = $sth->fetchrow_hashref;
    return unless $details;

    return $self->new($details);
}

# Methods
#

 __PACKAGE__->meta->make_immutable;

1;

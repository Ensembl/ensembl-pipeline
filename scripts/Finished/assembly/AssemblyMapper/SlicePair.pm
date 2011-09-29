package AssemblyMapper::SlicePair;

use namespace::autoclean;
use Moose;

has align_support => (
    is       => 'ro',
    isa      => 'AssemblyMapper::Support',
    required => 1,
    handles  => [qw(ref_dba ref_sa ref_asm ref_start ref_end alt_dba alt_sa alt_asm alt_start alt_end)],
    );

has iterator_params => (
    is       => 'ro',
    isa      => 'HashRef',
    required => 1,
    );

has session_id => (
    is       => 'rw',
    isa      => 'Int',
    );

has ref_chr => (
    is => 'ro',
    isa => 'Str',
    required => 1,
    );

has ref_slice => (
    is => 'ro',
    isa => 'Bio::EnsEMBL::Slice',
    required => 1,
    );

has ref_seq_region_id => (
    is       => 'ro',
    isa      => 'Int',
    required => 1,
    );

has alt_chr => (
    is => 'ro',
    isa => 'Str',
    required => 1,
    );

has alt_slice => (
    is => 'ro',
    isa => 'Bio::EnsEMBL::Slice',
    required => 1,
    );

has alt_seq_region_id => (
    is       => 'ro',
    isa      => 'Int',
    required => 1,
    );

__PACKAGE__->meta->make_immutable;

1;

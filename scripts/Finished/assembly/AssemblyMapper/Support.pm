package AssemblyMapper::Support;

use namespace::autoclean;
use Moose;

with 'AssemblyMapper::ClassUtils';

use Carp;
use Pod::Usage;
use Readonly;
use Bio::EnsEMBL::Utils::ConversionSupport;
use AssemblyMapper::AlignSession;
use AssemblyMapper::SlicePair;

Readonly my @OTTER_ALIGN_COMMON_OPTIONS => (
    'assembly=s',
    'altassembly=s',
    'altdbname=s',
    'chromosomes|chr=s@',
    'altchromosomes|altchr=s@',
    'force_stage|force-stage!',
    );

Readonly my @OTTER_ALIGN_REQUIRED_PARAMS => (
    'assembly',
    'altassembly',
    );

has support => (
    is       => 'ro',
    isa      => 'Bio::EnsEMBL::Utils::ConversionSupport',
    writer   => '_set_support',
    init_arg => undef,
    handles  => {
        error   => 'error',
        param   => 'param',
        s_param => 's_param',

        ref_asm   => [ s_param => 'assembly'  ],
        ref_start => [ s_param => 'ref_start' ],
        ref_end   => [ s_param => 'ref_end'   ],
        alt_asm   => [ s_param => 'altassembly' ],
        alt_start => [ s_param => 'alt_start' ],
        alt_end   => [ s_param => 'alt_end'   ],

        comma_to_list => 'comma_to_list',
        filehandle    => 'filehandle',
        finish_log    => 'finish_log',
        log           => 'log',
        log_error     => 'log_error',
        log_stamped   => 'log_stamped',
        log_verbose   => 'log_verbose',
        log_warning   => 'log_warning',
    },
    );

has extra_options => (
    is       => 'ro',
    isa      => 'ArrayRef[Str]',
    default  => sub { [] },
    );

has required_params => (
    is       => 'ro',
    isa      => 'ArrayRef[Str]',
    writer   => '_set_required_params',
    );

has single_chromosome => (
    is       => 'ro',
    isa      => 'Bool',
    default  => undef,
    );

has ref_dba => (
    is       => 'ro',
    isa      => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    writer   => '_set_ref_dba',
    init_arg => undef,
    handles  => {
        ref_dbc => 'dbc',
        ref_sa  => 'get_SliceAdaptor',
    },
    );

has alt_dba => (
    is       => 'ro',
    isa      => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    writer   => '_set_alt_dba',
    init_arg => undef,
    handles  => {
        alt_dbc => 'dbc',
        alt_sa  => 'get_SliceAdaptor',
    },
    );

# Constructors
#
sub BUILD {
    my ($self, $args) = @_;

    $self->_set_support(Bio::EnsEMBL::Utils::ConversionSupport->new('/dev/null'));

    # Ensure default required params are present
    #
    my @required_params = @{ $self->required_params || [] };
    if (@required_params) {
        foreach my $default (@OTTER_ALIGN_REQUIRED_PARAMS) {
            push @required_params, $default unless grep { /$default/ } @required_params;
        }
        $self->_set_required_params( [ @required_params ] );
    } else {
        $self->_set_required_params( [ @OTTER_ALIGN_REQUIRED_PARAMS ] );
    }

    return;
}

# Methods
#
sub parse_arguments {
    my ($self, @args) = @_;
    my $support = $self->support;

    $support->parse_common_options(@args);
    $support->parse_extra_options(
        @OTTER_ALIGN_COMMON_OPTIONS,
        @{ $self->extra_options },
        );

    $support->allowed_params(
        $self->get_common_params,
        $self->extra_params,
        );

    return if $support->param('help') or $support->error;

    $support->check_required_params(@{$self->required_params}); # dies if not

    $support->comma_to_list( 'chromosomes', 'altchromosomes' );
    if ($self->single_chromosome) {
        my @ref_chr = $support->param('chromosomes');
        my @alt_chr = $support->param('altchromosomes');
        if ( scalar(@ref_chr) > 1 or scalar(@alt_chr) > 1 ) {
            $support->error('This script does not support chromosome lists');
            return;
        }
    }

    $support->init_log;

    # some scripts need to know whether this was set originally or not
    $support->param(original_altdbname => $support->param('altdbname'));

    # set connection parameters for alternative db.
    # both databases have to be on the same host, so we don't need to configure
    # them separately
    for my $prm (qw(host port user pass dbname)) {
        $support->param( "alt$prm", $support->param($prm) )
            unless ( $support->param("alt$prm") );
    }

    return 1;
}

# List of our common params (inc. ConversionSupport ones)
sub get_common_params {
    my $self = shift;
    return ($self->support->get_common_params,
            $self->opts_to_params(@OTTER_ALIGN_COMMON_OPTIONS),
        );
}

# Turn extra_options into a simple list of params
sub extra_params {
    my $self = shift;
    return $self->opts_to_params(@{$self->extra_options});
}

sub opts_to_params {
    my ($self, @opts) = @_;
    my @params;
    foreach my $opt (@opts) {
        my ($param) = $opt =~ m/^([\w\-]+)[!+=:|]/;
        push @params, $param;
    }
    return @params;
}

sub connect_dbs {
    my $self = shift;
    my $params = $self->validate_params(
        \@_,
        rebless_dba => { isa => 'Str', optional => 1 },
        );
    my $support = $self->support;

    my $ref_dba = $support->get_database('ensembl', '');
    my $alt_dba = $support->get_database('ensembl', 'alt');

    if ($params->{rebless_dba}) {
        bless $ref_dba, $params->{rebless_dba};
        bless $alt_dba, $params->{rebless_dba};
    }

    $self->_set_ref_dba($ref_dba);
    $self->_set_alt_dba($alt_dba);

    AssemblyMapper::AlignSession->Dbc($self->ref_dbc);

    return;
}

sub iterate_chromosomes {
    my $self = shift;
    my $params = $self->validate_params(
        \@_,
        prev_stage     => { isa => 'Maybe[Str]'               }, # better val'n, see AlignStage.pm?
        before_stage   => { isa => 'Str',       optional => 1 }, # better val'n, see AlignStage.pm?
        this_stage     => { isa => 'Str'                      }, # better val'n, see AlignStage.pm?
        create_session => { isa => 'Bool',      optional => 1 },
        worker         => { isa => 'CodeRef'                  },
        callback_data  => { isa => 'Defined',   optional => 1 },
        );
    my $support = $self->support;

    $support->log_stamped("Looping over chromosomes...\n");

    my @ref_chr_list = $support->param('chromosomes');
    if ( !scalar(@ref_chr_list) ) {
        @ref_chr_list = $support->sort_chromosomes;

        if ( scalar( $support->param('altchromosomes') ) ) {
            croak "AltChromosomes list is defined while Chromosomes list is not!";
        }
    }

    my @alt_chr_list = $support->param('altchromosomes');
    if ( !scalar(@alt_chr_list) ) {
        @alt_chr_list = @ref_chr_list;
    }
    elsif ( scalar(@ref_chr_list) != scalar(@alt_chr_list) ) {
        croak "Chromosome lists do not match by length";
    }

    my $ok = 1;

  CHR: for my $i ( 0 .. scalar(@ref_chr_list) - 1 ) {
      my $ref_chr = $ref_chr_list[$i];
      my $alt_chr = $alt_chr_list[$i];

      $support->log_stamped( "Chromosome $ref_chr/$alt_chr ...\n", 1 );

      # fetch chromosome slices
      my $ref_slice =
          $self->ref_sa->fetch_by_region(
              'chromosome',
              $ref_chr,
              $self->ref_start,
              $self->ref_end,
              undef,
              $self->ref_asm,
          );
      my $ref_seq_region_id = $self->ref_sa->get_seq_region_id($ref_slice);

      my $alt_slice =
          $self->alt_sa->fetch_by_region(
              'chromosome',
              $alt_chr,
              $self->alt_start,
              $self->alt_end,
              undef,
              $self->alt_asm,
          );
      my $alt_seq_region_id = $self->alt_sa->get_seq_region_id($alt_slice);

      $support->log("Ref: ".$ref_slice->seq_region_name.", seq_region: ".$ref_seq_region_id."\n", 2);
      $support->log("Alt: ".$alt_slice->seq_region_name.", seq_region: ".$alt_seq_region_id."\n", 2);

      my $align_slice_pair = AssemblyMapper::SlicePair->new(
          align_support      => $self,
          iterator_params    => $params,
          ref_chr            => $ref_chr,
          ref_slice          => $ref_slice,
          ref_seq_region_id  => $ref_seq_region_id,
          alt_chr            => $alt_chr,
          alt_slice          => $alt_slice,
          alt_seq_region_id  => $alt_seq_region_id,
      );

      unless ($self->session_setup($align_slice_pair, $support->param('original_altdbname'))) {
          $ok = 0;
          next CHR;
      }


      $ok &&= &{$params->{worker}}(
          $align_slice_pair,
          $params->{callback_data},
      );

  } # CHR

    return $ok;
}

sub session_setup {
    my ($self, $asp, $alt_db_name) = @_;

    my $params = $asp->iterator_params;

    my $ref_seq_region_id = $asp->ref_seq_region_id;
    my $alt_seq_region_id = $asp->alt_seq_region_id;

    my $session = AssemblyMapper::AlignSession->latest(
        ref_seq_region_id => $ref_seq_region_id,
        alt_seq_region_id => $alt_seq_region_id,
        alt_db_name       => $alt_db_name,
        );
    if ($params->{create_session}) {
        if ($session) {
            $self->support->log_warning(
                "Already a session entry for [${ref_seq_region_id}, ${alt_seq_region_id}]\n", 3);
            return;
        }
        $session = AssemblyMapper::AlignSession->new(
            ref_seq_region_id => $ref_seq_region_id,
            alt_seq_region_id => $alt_seq_region_id,
            alt_db_name       => $alt_db_name,
            author            => $self->support->s_param('author'),
            comment           => $self->support->s_param('comment'),
            );
    }

    $session->create_stage(
        previous   => $params->{prev_stage},
        before     => $params->{before_stage},
        stage      => $params->{this_stage},
        script     => "$0",
        parameters => $self->support->list_all_params,
        dry_run    => $self->support->s_param('dry_run'),
        force_stage=> $self->support->s_param('force_stage'),
        );

    $asp->session_id($session->id);
}

__PACKAGE__->meta->make_immutable;

package Bio::EnsEMBL::Utils::ConversionSupport;

# Ensure parameter is evaluated in scalar context.
# Read-only
#
sub s_param {
    my ($self, $name) = @_;
    return scalar($self->param($name));
}

1;

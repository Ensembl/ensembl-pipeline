package AssemblyMapper::ClassUtils;

use namespace::autoclean;
use Moose::Role;
use MooseX::Params::Validate;

use Carp;
use Scalar::Util;

# This is only necessary because MooseX::Params::Validate only
# introduced hash&hash-ref handling in 0.16, and we have 0.13.
#
sub validate_params {
    my ($class, $arg_ref, %spec) = @_;
    my $arrayref = (   @$arg_ref == 1
                      && ref $arg_ref->[0]
                      && Scalar::Util::reftype( $arg_ref->[0] ) eq 'HASH' )
        ? [ %{$arg_ref->[0]} ]
        : $arg_ref;
    return { validated_hash(
                 $arrayref,
                 %spec,
                 MX_PARAMS_VALIDATE_NO_CACHE => 1
                 ) };
}

1;

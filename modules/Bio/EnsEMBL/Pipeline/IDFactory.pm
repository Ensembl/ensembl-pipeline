# Bio::EnsEMBL::Pipeline::IDFactory
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright (c) EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::IDFactory

=head1 SYNOPSIS

  ...
  my $idfactory = $module->new(-NAME => $name,
                               -CONFIG   => $config);

  while(my $id_set = $idfactory->next(20)) {
    #create some jobs with input ids $id_set
    ...
  }

=head1 DESCRIPTION

This is an abstract base class that provides a common interface over different
methods of retrieving input identifiers.  This class defines a base constructor
and the abstract method next() which can be used to iterate over a collection
of identifiers.  The implementation of how the identifiers are retrieved is
left to the base classes which can sit overtop of databases, file, etc.

Specific information needed by subclasses (e.g. database connection 
information) can be obtained through values defined in the configuration 
object. 

=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods have an underscore prefix.

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::IDFactory;


use Bio::EnsEMBL::Root;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [-CONFIG]: Bio::EnsEMBL::Pipeline::Config 
                 The configuration for the running pipeline
  Arg [-NAME]: string
               The name of this id factory in the configuration
  Example    : $idf = Bio::EnsEMBL::Pipeline::IDFactory->new(
                                           -NAME => $name,
                                           -CONFIG => $conf);
  Description: Instantiates an IDFactory class
  Returntype : Bio::EnsEMBL::Pipeline::IDFactory
  Exceptions : thrown on incorrect arguments
  Caller     : none

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = bless({}, $class);

  my ($config, $name) = $self->_rearrange([qw(CONFIG NAME)], @_);

  $self->throw('config argument is required') 
    if(!ref($config) || !$config->isa('Bio::EnsEMBL::Pipeline::Config'));

  $self->throw('name argument is required') if(!$name);

  $self->{'config'} = $config;
  $self->{'name'} = $name;

  return $self;
}


=head2 name

  Arg [1]    : none 
  Example    : print $idfactory->name();
  Description: Retrieves the name of this IDFactory
  Returntype : string
  Exceptions : none
  Caller     : internal

=cut

sub name {
  my $self = shift;

  return $self->{'name'};
}


=head2 config

  Arg [1]    : none
  Example    : $config = $idfactory->config();
  Description: Retrieves the configuration that this IDFactory is using
  Returntype : Bio::EnsEMBL::Pipeline::Config
  Exceptions : none
  Caller     : internal

=cut

sub config {
  my $self = shift;

  return $self->{'config'};
}


=head2 next

  Arg [1]    : (optional) int $size
               The number of identifiers to retrieve
  Example    : while($ids = $idset->next(20)) { do_something() };
  Description: Retrieves the next bunch of identifiers from this factory
               if $size is not speicifed all of the remaining identifiers are
               retrieved.  If there are no remaining identifiers, undef is 
               returned instead.
  Returntype : Bio::EnsEMBL::Pipeline::IDSet
  Exceptions : none
  Caller     : general

=cut

sub next {
  my $self = shift;

  $self->throw('abstract method next not implmented by subclass');
}





1;

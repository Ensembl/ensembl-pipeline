# Bio::EnsEMBL::Pipeline::IDFactory::Config
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright (c) EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::IDFactory::Config

=head1 SYNOPSIS

  ...
  my $idfactory = Bio::EnsEMBL::Pipeline::IDFactory::Config->new(
                               -NAME   => $name,
                               -CONFIG => $config);

  while(my $id_set = $idfactory->next(20)) {
    #create some jobs with input ids $id_set
    ...
  }

=head1 DESCRIPTION

This is an implementation of the IDFactory interface.  It provides a means to 
retrieve identifiers directly from the configuration file.

=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods have an underscore prefix.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::IDFactory::Config;

use vars qw(@ISA);

use Bio::EnsEMBL::Pipeline::IDFactory;
use Bio::EnsEMBL::Pipeline::IDSet;

@ISA = qw(Bio::EnsEMBL::Pipeline::IDFactory);



=head2 new

  Arg [-CONFIG]: Bio::EnsEMBL::Pipeline::Config
               The configuration for the running pipeline
  Arg [-NAME]: string
               The name of that this ID Factory from the config
  Example    : $idf = Bio::EnsEMBL::Pipeline::IDFactory::Database->new(
                                           -NAME => $name,
                                           -CONFIG => $conf);
  Description: Instantiates a Database IDFactory class.  The configuration 
               object that is passed to this constructor must define an
               id_list parameter that contains a semi-colon delimited 
               list of identifiers

               An example of a configuration file specifying an id factory
               that creates ids from the configuration file follows

               ------
               [SRS_ID_FACTORY]
               module  = Bio::EnsEMBL::Pipeline::IDFactory::Config
               id_list = SPTREMBL,SWISSPROT 

               [SRS_FETCH_TASK]
               id_factory = SRS_ID_FACTORY

               ------

  Returntype : Bio::EnsEMBL::Pipeline::IDFactory::Config
  Exceptions : thrown on incorrect arguments
  Caller     : none

=cut

sub new {
  my $caller = shift;

  #new can be called as either a class or object method
  my $class = ref($caller) || $caller;

  my $self = bless({}, $class);

  $self = $self->SUPER::new(@_);

  my $name = $self->name();
  my $config = $self->config();

  my $id_list = $config->get_parameter($name, 'id_list');

  if(!$id_list) {
    $self->throw("[$name] configuration must define value for the " .
                 "parameter id_list");
  }


  my @ids = split(/;/, $id_list);

  $self->ids(\@ids);

  return $self;
}



sub ids {
  my $self = shift;
  $self->{'idlist'} = shift if(@_);
  return $self->{'idlist'};
}



=head2 next

  Arg [1]    : (optional) int $num_ids
               The number of ids to retrieve
  Example    : while($idset = $id_factory->next(20)) { do something; }; 
  Description: Retrieves an idset containing the next $num_ids from the files
               that this id factory is using. If there are less than 
               $num_ids left in the files, then the remaining ids will
               be returned as part of the idset.  If there are no
               ids left in the files undef will be returned instead.  If 
               the $num_ids argument was not defined, all remaining ids will
               be returned.
  Returntype : Bio::EnsEMBL::Pipeline::IDSet or undef
  Exceptions : none
  Caller     : Tasks

=cut

sub next {
  my $self = shift;
  my $num_ids = shift;

  my $all_ids = $self->ids();
  my @ids;

  while(@$all_ids && scalar(@ids) < $num_ids) {
    push(@ids , shift(@$all_ids));
  }

  #
  # Return undef if ids were requested and none were available
  #
  if((!defined($num_ids) || $num_ids > 0) && !@ids) {
    return undef;
  }

  return Bio::EnsEMBL::Pipeline::IDSet->new(-ID_LIST => \@ids);
}


1;

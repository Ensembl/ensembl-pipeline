# Bio::EnsEMBL::Pipeline::IDFactory::Slice
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright (c) EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::IDFactory::Slice

=head1 SYNOPSIS

  ...
  my $idfactory = Bio::EnsEMBL::Pipeline::IDFactory::Slice->new(
                               -NAME => $name,
                               -CONFIG   => $config);

  while(my $id_set = $idfactory->next(20)) {
    #create some jobs with input ids $id_set
    ...
  }

=head1 DESCRIPTION

This is an implementation of the IDFactory interface.  It provides a means to 
retrieve identifiers which are slice names.  An ensembl database must be
specified in the configuration for this IDFactory so that chromosomal 
information can be obtained. This module is also dependant on the EnsEMBL
perl API, which must be in the PERL5LIB environtment variable for this module
to work.

=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::IDFactory::Slice;

use Bio::EnsEMBL::Pipeline::IDFactory;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::IDSet;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Pipeline::IDFactory);


=head2 new

  Arg [-CONFIG]: Bio::EnsEMBL::Pipeline::Config 
               The configuration for the running pipeline
  Arg [-NAME]: string
               The name of that this ID Factory from the config
  Example    : $idf = Bio::EnsEMBL::Pipeline::IDFactory::Database->new(
                                           -NAME => $name,
                                           -CONFIG => $conf);
  Description: Instantiates a Slice IDFactory class.  The configuration 
               object that is passed to this constructor must define an 
               ensembl database, and the size and overlap of slices created.

               ------
               [ENSEMBL_DATABASE]
               host = ecs1g
               dbname = homo_sapiens_core_15_33
               user = ensro
               pass =

               [SLICE_ID_FACTORY]
               module = Bio::EnsEMBL::Pipeline::IDFactory::Slice
               database = ensembl_database
               overlap = 1000
               size = 1000000

               [GENSCAN_TASK]
               where = LSF:acari
               id_factory = slice_id_factory

               ------

  Returntype : Bio::EnsEMBL::Pipeline::IDFactory::Database
  Exceptions : thrown on incorrect arguments
  Caller     : none

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = bless({}, $class);

  $self = $self->SUPER::new(@_);

  my $name = $self->name();
  my $conf = $self->config();

  my $size    = $conf->get_parameter($name, 'size');
  my $overlap = $conf->get_parameter($name, 'overlap');

  if(!defined($size) || !defined($overlap)) {
    $self->throw("[$name] configuration must define a values for" .
                 " the parameters 'size' and 'overlap'");
  }

  if($overlap >= $size) {
    $self->throw("overlap must be less than size");
  }

  $self->size($size);
  $self->overlap($overlap);

  my $dbheader = $conf->get_parameter($name, 'database');

  if(!$dbheader) {
    $self->throw("[$name] configuration must define a value for" .
                 " the parameter 'database'");
  }

  my $host     = $conf->get_parameter($dbheader, 'host');
  my $user     = $conf->get_parameter($dbheader, 'user');
  my $pass     = $conf->get_parameter($dbheader, 'pass');
  my $dbname   = $conf->get_parameter($dbheader, 'dbname');
  my $port     = $conf->get_parameter($dbheader, 'port');
  my $driver   = $conf->get_parameter($dbheader, 'driver');

  if(!$dbname) {
    $self->throw("Database configuration [$dbheader] must define value for " .
                 "the parameter dbname");
  }

  $self->user($user);
  $self->pass($pass);
  $self->dbname($dbname);
  $self->port($port);
  $self->host($host);
  $self->driver($driver);

  return $self;
}


sub size {
  my $self = shift;
  $self->{'size'} = shift if(@_);
  return $self->{'size'};
}


sub overlap {
  my $self = shift;
  $self->{'overlap'} = shift if(@_);
  return $self->{'overlap'};
}

sub user {
  my $self = shift;
  $self->{'user'} = shift if(@_);
  return $self->{'user'};
}

sub pass {
  my $self = shift;
  $self->{'pass'} = shift if(@_);
  return $self->{'pass'};
}

sub host {
  my $self = shift;
  $self->{'host'} = shift if(@_);
  return $self->{'host'};
}

sub dbname {
  my $self = shift;
  $self->{'dbname'} = shift if(@_);
  return $self->{'dbname'};
}

sub driver {
  my $self = shift;
  $self->{'driver'} = shift if(@_);
  return $self->{'driver'};
}

sub port {
  my $self = shift;
  $self->{'port'} = shift if(@_);
  return $self->{'port'};
}


=head2 next

  Arg [1]    : (optional) int $num
               The number of ids to retrieve
  Example    : while($idset = $id_factory->next(20)) { do something; }; 
  Description: Retrieves an idset containing the next $num of slice 
               identifiers from the ensembl database that this id factory is 
               using. If there are less than $num slices left, then the 
               remaining ids will be returned as part of the idset.  
               If there are no slices left undef will be returned instead.  If 
               the $num argument was not defined, all remaining ids will
               be returned.

               The slice identifiers which are retrieved are formatted as:
               chromosome_name.start-end
               For example "X.1-100000"
  Returntype : Bio::EnsEMBL::Pipeline::IDSet or undef
  Exceptions : none
  Caller     : Tasks

=cut

sub next {
  my $self = shift;
  my $num  = shift;

  my $db = $self->{'db'};

  #
  # Database is loaded on demand the first time next is called
  #
  if(!$db) {
    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $self->dbname,
                                              -user   => $self->user,
                                              -pass   => $self->pass,
                                              -host   => $self->host,
                                              -driver => $self->driver,
                                              -port   => $self->port);
    $self->{'db'} = $db;
    $self->{'chrs'} = $db->get_ChromosomeAdaptor->fetch_all;
    $self->{'cur_start'} = 1;
    $self->{'cur_chr'} = shift(@{$self->{'chrs'}});
  }

  my $cur_chr = $self->{'cur_chr'};
  my $cur_start  = $self->{'cur_start'};
  my $chrs = $self->{'chrs'};
  my $length = $cur_chr->length();
  my $overlap = $self->overlap();
  my $size = $self->size();

  #sanity check for bad data in DB
  $self->throw('Chromosome length is not defined') if(!$length);

  #
  #retrieve as many slices as requested or until no more can be created
  #
  my @ids;
  while((!defined($num) || (scalar(@ids) < $num)) && 
        (scalar(@$chrs) || ($cur_start <= $length))) {

    #
    # move on to next chromosome if this one is done
    #
    if($cur_start > $length) {
      $cur_chr = shift(@$chrs);
      $length = $cur_chr->length();
      $cur_start = 1;
      next;
    }

    my $end = $cur_start + $size - 1;

    if($end > $length) {
      $end = $length;
      #should we try to get a full sized slice if we are at the end?
      #what if it is essential to have a certain overlap?
      #$cur_start = $length - $size + 1;
      #$cur_start = 1 if($cur_start < 1);
    }

    push(@ids, $cur_chr->chr_name . ".$cur_start-$end");

    #
    #if this slice was right to the end than it doesn't make sense to 
    #make another slice with overlap
    #
    if($end < $length) {
      $cur_start += $size - $overlap;
    } else {
      $cur_start = $length +1;
    }
  }

  #save state for next time this method is called
  $self->{'cur_chr'} = $cur_chr;
  $self->{'cur_start'} = $cur_start;

  # return undef if at least one id was requested and there were none
  return undef if((!defined($num) || ($num > 0)) && !@ids);

  return Bio::EnsEMBL::Pipeline::IDSet->new(-id_list => \@ids);
}





1;


# Bio::EnsEMBL::Pipeline::IDFactory::Database
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)
#
# Copyright (c) EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

  Bio::EnsEMBL::Pipeline::IDFactory::Database

=head1 SYNOPSIS

  ...
  my $idfactory = Bio::EnsEMBL::Pipeline::IDFactory::Database->new(
                               -NAME   => $name,
                               -CONFIG => $config);

  while(my $id_set = $idfactory->next(20)) {
    #create some jobs with input ids $id_set
    ...
  }

=head1 DESCRIPTION

This is an implementation of the IDFactory interface.  It provides a means to 
retrieve identifiers stored in a single column of a database table.
It is assumed that the table will have a single unique identifier per row.

=head1 CONTACT

Post general questions to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods have an underscore prefix.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Pipeline::IDFactory::Database;

use vars qw(@ISA);

use DBI;
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
               object that is passed to this constructor must define a
               database, a table, and a column to be used to retrieve
               identifiers.

               Note that a database connection is not actually established
               until the next() method is called for the first time.  The
               database does not need to be present until this time.

               An example of a configuration file specifying an id factory
               that creates ids from the name column of contig table in 
               an ensembl database follows.

               ------
               [ENSEMBL_DATABASE]
               host = ecs1g
               dbname = homo_sapiens_core_15_33
               user = ensro
               pass =

               [CONTIG_ID_FACTORY]
               module = Bio::EnsEMBL::Pipeline::IDFactory::Database
               database = ensembl_database
               table = contig
               column = name

               [CHROMOSOME_ID_FACTORY]
               module = Bio::EnsEMBL::Pipeline::IDFactory::Database
               database = ensembl_database
               table = chromosome
               column = name

               [REPEAT_MASKER_TASK]
               where = LSF:acari
               id_factory = contig_id_factory

               [BLAST_TASK]
               where = LSF:acari
               id_factory = contig_id_factory

               ------

  Returntype : Bio::EnsEMBL::Pipeline::IDFactory::Database
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

  #
  # Read database connection information from the configuration
  #

  my $dbheader = $config->get_parameter($name, 'database');

  if(!$dbheader) {
    $self->throw("[$name] configuration must define a value for " .
                 "the parameter 'database'");
  }

  my $host = $config->get_parameter($dbheader, 'host') || 'localhost';
  my $user = $config->get_parameter($dbheader, 'user') || 'root';
  my $pass = $config->get_parameter($dbheader, 'pass') || '';
  my $dbname = $config->get_parameter($dbheader, 'dbname');
  my $port = $config->get_parameter($dbheader, 'port') || 3306;
  my $driver = $config->get_parameter($dbheader, 'driver') || 'mysql';

  if(!$dbname) {
    $self->throw("Database configuration [$dbheader] must define value for " .
                 "the parameter dbname");
  }

  $self->dsn("DBI:$driver:database=$dbname;host=$host;port=$port");
  $self->user($user);
  $self->pass($pass);


  my $table = $config->get_parameter($name, 'table');
  my $column = $config->get_parameter($name, 'column');

  if(!$table) {
    $self->throw("[$name] configuration must define a value for the parameter"
                 . " table") 
  }
  if(!$column) {
    $self->throw("[$name] configuration must define a vlue for the parameter"
                 . " column");
  }

  $self->column($column);
  $self->table($table);

  return $self;
}

sub dsn {
  my $self = shift;
  $self->{'dsn'} = shift if(@_);
  return $self->{'dsn'};
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


sub table {
  my $self = shift;
  $self->{'table'} = shift if(@_);
  return $self->{'table'};
}

sub column {
  my $self = shift;
  $self->{'column'} = shift if(@_);
  return $self->{'column'};
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

  my $dbh = $self->{'dbh'};
  my $sth = $self->{'sth'};

  #
  # If this is the first time calling this method we need to establish
  # a database connection and execute the query
  #
  if(!$dbh) {
    my $col = $self->column();
    my $tab = $self->table();

    $dbh = DBI->connect($self->dsn(), $self->user(), 
                        $self->pass(), {'RaiseError' => 1});
    $sth = $dbh->prepare("select $col from $tab");
    $sth->execute();

    $self->{'sth'} = $sth;
    $self->{'dbh'} = $dbh;
  }

  #
  # Retrieve ids from the executed statement until there are none left,
  # or we have retrieved the number requested
  #
  my @ids;
  my $id;
  while($sth && (!defined($num_ids) || scalar(@ids) < $num_ids)) {
    ($id) = $sth->fetchrow();
    if(!defined($id)) {
      $sth->finish();
      $self->{'sth'} = undef;
      last;
    }

    push @ids, $id;
  }

  #
  # If > 0 ids were requested and there were none, return undef
  #
  return undef if((!defined($num_ids) || $num_ids > 0) && scalar(@ids) == 0);



  return Bio::EnsEMBL::Pipeline::IDSet->new(-ID_LIST => \@ids);
}


sub DESTROY {
  my $self = shift;

  $self->{'sth'}->finish if($self->{'sth'});
  $self->{'dbh'}->disconnect if($self->{'dbh'});
}


1;

# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=head1 NAME

  Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor - Object representing an instance of an EnsEMBL DB

=head1 SYNOPSIS

  use Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor;

  $db = Bio::EnsEMBL::DBSQL::Protein::DBAdaptor->new ( -host   => 'caldy',
                                                       -user   => 'root',
                                                       -dbname => 'pog',
                                                       -driver => 'mysql',
                                                       -pass   => 'xyz',
                                                     );

=head1 DESCRIPTION

  This object represents a database that is implemented somehow (you shouldn\'t
  care much as long as you can get the object). From the object you can pull
  out other objects by their stable identifier, such Proteins. The Protein gives 
  you a DB::Protein object, from which you can pull out features etc. 

=head1 CONTACT

  Marc Sohrmann (ms2@sanger.ac.uk)

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _.

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::Protein::DBAdaptor;

use vars qw(@ISA);
use strict;
use DBI;

# Object preamble - inherits from Bio::Root::Object

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer;

@ISA = qw(Bio::EnsEMBL::Root);


sub new {
    my($class, @args) = @_;

    my $self = bless {}, $class;

    my ( $dbname,
         $host,
         $driver,
         $user,
         $password,
         $port,
         $sql_big_tables,
       ) = $self->_rearrange([qw(DBNAME
                                 HOST
	                         DRIVER
  	                         USER
	                         PASS
	                         PORT
                                 BIG)],
                              @args);

    $dbname || $self->throw("Database object must have a database name");
    $host || $self->throw("Database object must have a database host");
    $user || $self->throw("Database object must have a user");
  
    if( ! $driver ) {
        $driver = 'mysql';
    }
    if( ! $host ) {
        $host = 'localhost';
    }
    if ( ! $port ) {
        $port = 3306;
    }

    my $dsn = "DBI:$driver:database=$dbname;host=$host;port=$port";
	
    my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 1});
    $dbh || $self->throw("Could not connect to database $dbname as $user using [$dsn] as a locator");

    $self->dbname ($dbname);
    $self->host ($host);
    $self->username ($user);
    $self->password ($password);
    $self->_db_handle ($dbh);
    $self->_analysis_cache ({});

    if ($sql_big_tables) {
        $dbh->do ('SET SQL_BIG_TABLES = 1');
    }

    return $self;
}

# accessor methods

sub host {
    my $self = shift;
    if (@_) {
        $self->{'_host'} = shift;
    }
    return $self->{'_host'};
} 

sub dbname {
    my $self = shift;
    if (@_) {
        $self->{'_dbname'} = shift;
    }
    return $self->{'_dbname'};
} 

sub username {
    my $self = shift;
    if (@_) {
        $self->{'_username'} = shift;
    }
    return $self->{'_username'};
} 

sub password {
    my $self = shift;
    if (@_) {
        $self->{'_password'} = shift;
    }
    return $self->{'_password'};
} 

sub _db_handle{
    my $self = shift;
    if (@_) {
        $self->{'_db_handle'} = shift;
    }
    return $self->{'_db_handle'};
} 


# some methods to get the required Adaptors

=head2 get_JobAdaptor

 Title    : get_JobAdaptor
 Usage    : $db->get_JobAdaptor
 Function : The Adaptor for Job objects in this db
 Example  : 
 Returns  : Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor
 Args     : 
 Throws   :

=cut

sub get_JobAdaptor {
    my ($self) = @_;
    unless ($self->{'_JobAdaptor'}) {
      $self->{'_JobAdaptor'} = Bio::EnsEMBL::Pipeline::DBSQL::JobAdaptor->new ($self);
    }
    return $self->{'_JobAdaptor'};
}

=head2 get_AnalysisAdaptor

 Title    : get_AnalysisAdaptor
 Usage    : $db->get_AnalysisAdaptor
 Function : The Adaptor for Analysis objects in this db
 Example  : 
 Returns  : Bio::EnsEMBL::DBSQL::AnalysisAdaptor
 Args     :
 Throws   : 


=cut

sub get_AnalysisAdaptor {
    my ($self) = @_;
    unless ($self->{'_AnalysisAdaptor'}) {
        $self->{'_AnalysisAdaptor'} = Bio::EnsEMBL::DBSQL::AnalysisAdaptor->new ($self);
    }
    return $self->{'_AnalysisAdaptor'};
}


=head2 get_RuleAdaptor

 Title    : get_RuleAdaptor
 Usage    : $db->get_RuleAdaptor
 Function : The Adaptor for Rule objects in this db
 Example  : 
 Returns  : Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor
 Args     :
 Throws   : 


=cut

sub get_RuleAdaptor {
    my ($self) = @_;
    unless ($self->{'_RuleAdaptor'}) {
        $self->{'_RuleAdaptor'} = Bio::EnsEMBL::Pipeline::DBSQL::RuleAdaptor->new ($self);
    }
    return $self->{'_RuleAdaptor'};
}

=head2 get_StateInfoContainer

 Title    : get_StateInfoContainer
 Usage    : $db->get_StateInfoContainer
 Function : 
 Example  : 
 Returns  : Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer
 Args     :
 Throws   : 


=cut

sub get_StateInfoContainer {
    my ($self) = @_;
    unless ($self->{'_StateInfoContainer'}) {
        $self->{'_StateInfoContainer'} = Bio::EnsEMBL::Pipeline::DBSQL::StateInfoContainer->new ($self);
    }
    return $self->{'_StateInfoContainer'};
}

=head2 get_ProteinAdaptor

 Title    : get_ProteinAdaptor
 Usage    : $db->get_ProteinAdaptor
 Function :
 Example  :
 Returns  : Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinAdaptor 
 Args     :
 Throws   :


=cut

sub get_ProteinAdaptor {
    my ($self) = @_;
    unless ($self->{'_ProteinAdaptor'}) {
        $self->{'_ProteinAdaptor'} = Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinAdaptor->new ($self);
    }
    return $self->{'_ProteinAdaptor'};
}

sub get_Protein_Adaptor {
    my ($self) = @_;
    $self->get_ProteinAdaptor;
}


=head2 get_ProteinFeatureAdaptor

 Title    : get_ProteinFeatureAdaptor
 Usage    : $db->get_ProteinFeatureAdaptor
 Function :
 Example  :
 Returns  : Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinFeatureAdaptor 
 Args     :
 Throws   :


=cut

sub get_ProteinFeatureAdaptor {
    my ($self) = @_;
    unless ($self->{'_ProteinFeatureAdaptor'}) {
        $self->{'_ProteinFeatureAdaptor'} = Bio::EnsEMBL::Pipeline::DBSQL::Protein::ProteinFeatureAdaptor->new ($self);
    }
    return $self->{'_ProteinFeatureAdaptor'};
}

sub get_Protein_Feature_Adaptor {
    my ($self) = @_;
    $self->get_ProteinFeatureAdaptor;
}

sub get_Protfeat_Adaptor {
    my ($self) = @_;
    $self->get_ProteinFeatureAdaptor;
}


=head2 _analysis_cache

 Title   : _analysis_cache
 Usage   : $obj->_analysis_cache()
 Function: 
 Returns : reference to a hash
 Args    : newvalue (optional)


=cut

sub _analysis_cache{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_analysis_cache'} = $value;
    }
    return $obj->{'_analysis_cache'};

}




=head2 prepare

 Title    : prepare
 Usage    : $db->prepare
 Function : prepare sql query
 Example  :
 Returns  : statement handle
 Args     :
 Throws   :


=cut

sub prepare {
    my ($self, $string) = @_;

    unless ($string ) {
        $self->throw ("attempting to prepare an empty SQL query");
    }
    unless (defined $self->_db_handle ) {
       $self->throw("database object has lost its database handle");
    }
    return $self->_db_handle->prepare ($string);
}


=head2 disconnect

 Title    : disconnect
 Usage    : $db->disconnect
 Function : disconnect the database handle
 Example  :
 Returns  : 
 Args     :
 Throws   :


=cut

sub disconnect {
    my ($self) = @_;
    if ($self->_db_handle) {
        $self->_db_handle->disconnect;
    }
}


=head2 DESTROY

 Title    : DESTROY
 Usage    :
 Function : disconnect the database handle
 Example  :
 Returns  : 
 Args     :
 Throws   :


=cut

sub DESTROY {
   my ($self) = @_;
   if ($self->_db_handle) {
       $self->_db_handle->disconnect;
   }
}   

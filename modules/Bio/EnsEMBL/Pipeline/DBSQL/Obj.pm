#
# Bioperl module for TimSQL::Obj
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::TimSQL::Obj - Object representing an instance of an EnsEMBL DB

=head1 SYNOPSIS

    $db = new Bio::TimSQL::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );

    $clone  = $db->get_clone('X45667');

=head1 DESCRIPTION

This object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object). From the object you can pull
out clone objects by their stable identifier, such as Clone (accession number). 
The clone gives you a DB::Clone object, from which you can pull out desired fields. 

=head1 CONTACT

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::TimSQL::Obj;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;

use Bio::TimSQL::Clone;
use DBI;

use Bio::EnsEMBL::DBSQL::DummyStatement;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($db,$host,$driver,$user,$password,$debug) = 
      $self->_rearrange([qw(DBNAME
			    HOST
			    DRIVER
			    USER
			    PASS
			    DEBUG
			    )],@args);

  $db   || $self->throw("Database object must have a database name");
  $user || $self->throw("Database object must have a user");

  $self->{'_lock_table_hash'} = {};

  if( $debug ) {
      $self->_debug($debug);
  } else {
      $self->_debug(0);
  }
  
  if( ! $driver ) {
      $driver = 'mysql';
  }

  if( ! $host ) {
      $host = 'localhost';
  }

  my $dsn = "DBI:$driver:database=$db;host=$host";

  if( $debug && $debug > 10 ) {
      $self->_db_handle("dummy dbh handle in debug mode $debug");
  } else {

      my $dbh = DBI->connect("$dsn","$user",$password);

      $dbh || $self->throw("Could not connect to database $db user $user using [$dsn] as a locator");
      
      if( $self->_debug > 10 ) {
	  $self->warn("Using connection $dbh");
      }
     
      $self->_db_handle($dbh);
  }

  return $make; # success - we hope!

}

=head2 get_Clone

 Title   : get_Clone
 Usage   : $db->get_Clone($disk_id)
 Function: Gets a Clone object with disk_id as its disk_id
 Example : $db->get_Clone($disk_id)
 Returns : Bio::TimSQL::Clone object 
 Args    : disk_id


=cut

sub get_Clone{
   my ($self,$disk_id) = @_;

   $disk_id || $self->throw("Trying to delete a clone without a disk id!");

   my $sth = $self->prepare("select disk_id from clone where disk_id = \"$disk_id\";");
   my $res = $sth ->execute();
   my $rv  = $sth ->rows;

   if( ! $rv ) {
       # make sure we deallocate sth - keeps DBI happy!
       $sth = 0;
       $self->throw("Clone $disk_id does not seem to occur in the database!");
   }

   my $clone = new Bio::TimSQL::Clone( -disk_id    => $disk_id,
					       -dbobj => $self );

   return $clone;
}

=head2 delete_Clone

 Title   : delete_Clone
 Usage   : $db->delete_Clone($disk_id)
 Function: Deletes a clone with a specific disk_id
 Example : $db->delete_Clone($disk_id)
 Returns : nothing
 Args    : disk_id


=cut

sub delete_Clone{
   my ($self,$disk_id) = @_;

   $disk_id || $self->throw("Trying to delete a clone without a disk id!");

   my $sth = $self->prepare("delete from clone where disk_id = \"$disk_id\";");
   my $res = $sth ->execute();
   my $rv  = $sth ->rows;

   if( ! $rv ) {
       # make sure we deallocate sth - keeps DBI happy!
       $sth = 0;
       $self->throw("Clone $disk_id does not seem to occur in the database!");
   }
   return 1;
}

=head2 create_Clone

 Title   : create_Clone
 Usage   : $db->create_Clone($disk_id,$clone_group,$chromosome)
 Function: writes a new clone in the database
 Example : $db->create_Clone('dummy','SU','22')
 Returns : nothing
 Args    : disk_id,clone group,chromosome


=cut

sub create_Clone{
   my ($self,$disk_id,$clone_group,$chromosome) = @_;

   $disk_id || $self->throw("Trying to create a clone without a disk id!\n");
   $clone_group || $self->throw("Trying to create a clone without a clone group!");
   $chromosome || $self->throw("Trying to create a clone without a chromosome id!");
   
   my @sql;

   push(@sql,"lock tables clone write");
   push(@sql,"insert into clone(disk_id,clone_group,chromosome,last_check,created) values('$disk_id','$clone_group','$chromosome',now(),now())");
   push(@sql,"unlock tables");   

   foreach my $sql (@sql) {
     my $sth =  $self->prepare($sql);
     #print STDERR "Executing $sql\n";
     my $rv  =  $sth->execute();
     $self->throw("Failed to insert clone $disk_id") unless $rv;
   }
}



=head2 prepare

 Title   : prepare
 Usage   : $sth = $dbobj->prepare("select seq_start,seq_end from feature where analysis = \" \" ");
 Function: prepares a SQL statement on the DBI handle

           If the debug level is greater than 10, provides information into the
           DummyStatement object
 Example :
 Returns : A DBI statement handle object
 Args    : a SQL string


=cut

sub prepare{
   my ($self,$string) = @_;

   if( ! $string ) {
       $self->throw("Attempting to prepare an empty SQL query!");
   }

   
   if ($self->_diffdump) {
       my $fh=$self->_diff_fh;
       open (FH,"$fh");
       if ($string =~/insert|delete|replace/) {
	   print FH "$string\n";
       }
   }

   if( $self->_debug > 10 ) {
       print STDERR "Prepared statement $string\n";
       my $st = Bio::EnsEMBL::DBSQL::DummyStatement->new();
       $st->_fileh(\*STDERR);
       $st->_statement($string);
       return $st;
   }

   # should we try to verify the string?

   return $self->_db_handle->prepare($string);
}

=head2 diff_fh

 Title   : diff_fh
 Usage   : $obj->diff_fh($newval)
 Function: path and name of the file to use for writing the mysql diff dump
 Example : 
 Returns : value of diff_fh
 Args    : newvalue (optional)


=cut

sub diff_fh{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'diff_fh'} = $value;
    }
    return $self->{'diff_fh'};
    
}


=head2 _diffdump

 Title   : _diffdump
 Usage   : $obj->_diffdump($newval)
 Function: If set to 1 sets $self->prepare to print the diff sql 
           statementents to the filehandle specified by $self->diff_fh
 Example : 
 Returns : value of _diffdump
 Args    : newvalue (optional)


=cut

sub _diffdump{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_diffdump'} = $value;
    }
    return $self->{'_diffdump'};
    
}


=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut

sub _debug{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_debug'} = $value;
    }
    return $self->{'_debug'};
    
}


=head2 _db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function: 
 Example : 
 Returns : value of _db_handle
 Args    : newvalue (optional)


=cut

sub _db_handle{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};

}

=head2 _lock_tables

 Title   : _lock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _lock_tables{
   my ($self,@tables) = @_;
   
   my $state;
   foreach my $table ( @tables ) {
       if( $self->{'_lock_table_hash'}->{$table} == 1 ) {
	   $self->warn("$table already locked. Relock request ignored");
       } else {
	   if( $state ) { $state .= ","; } 
	   $state .= "$table write";
	   $self->{'_lock_table_hash'}->{$table} = 1;
       }
   }

   my $sth = $self->prepare("lock tables $state");
   my $rv = $sth->execute();
   $self->throw("Failed to lock tables $state") unless $rv;

}

=head2 _unlock_tables

 Title   : _unlock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _unlock_tables{
   my ($self,@tables) = @_;

   my $sth = $self->prepare("unlock tables");
   my $rv  = $sth->execute();
   $self->throw("Failed to unlock tables") unless $rv;
   %{$self->{'_lock_table_hash'}} = ();
}


=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub DESTROY{
   my ($obj) = @_;

   $obj->_unlock_tables();

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}

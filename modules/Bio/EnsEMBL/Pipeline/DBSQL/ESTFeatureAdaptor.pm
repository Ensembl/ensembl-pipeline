=head1 NAME

=head1 SYNOPSIS

  
=head1 DESCRIPTION


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::DBSQL::BaseAdaptor
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub get_est_length {
  my ($self,$est_id) = @_;
  my $query = qq(SELECT est_length 
		 FROM est 
		 WHERE est_id = "$est_id"
		);
  my $sth = $self->prepare($query) || $self->throw("can't prepare $query");
  my $res = $sth->execute || $self->throw("can't execute $query");
  
  my @lengths;
  while (my $length = $sth->fetchrow_array) {
    push (@lengths,$length);
  }
  @lengths = sort { $b <=> $a } @lengths;
  my $thislength = $lengths[0];
  return $thislength;
}

sub exists_est {
  my ($self,$est_id) = @_;

  my $query = "select * from est where est_id = '$est_id'";
  
  my $sth = $self->prepare($query);
  my $res = $sth->execute;

  if ($sth->rows > 0) {
    return 1;
  } else {
    return 0;
  }
}

sub create_sql {
  my($self,$dbname) = @_;

  my $create = "create database $dbname;";

  my $est = "create table est ( est_internal_id int(10) unsigned NOT NULL auto_increment,
					  est_id          varchar(40) NOT NULL, 
					  est_length int(10) unsigned NOT NULL,
					  
					  PRIMARY KEY (est_internal_id),
					  UNIQUE  est_id(est_id)
					);";
  
  
  return $create . "\n" . $est . "\n";
}

1;      

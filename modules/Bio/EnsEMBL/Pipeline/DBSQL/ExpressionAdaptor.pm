=head1 NAME

Bio::EnsEMBL::Pipeline::DBSQL::ExpressionAdaptor

=head1 SYNOPSIS

  
=head1 DESCRIPTION

BaseAdaptor for the database from SANBI containing curated expression data of human ESTs

in the main trunk this should inherit from DBConnection!


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Pipeline::DBSQL::ExpressionAdaptor;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


############################################################

sub get_class_ids{
  my ($self, $clonelibrary_id) = @_;
  my $q = qq ( SELECT SequenceClassLink.class_id
	       FROM   SequenceClassLink, Sequence
	       WHERE  SequenceClassLink.sequence_id = Sequence.id
	       AND    Sequence.clonelibrary_id = "$clonelibrary_id"
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  my @class_ids;
  while ( my ($class_id) = $sth->fetchrow_array ){
    push ( @class_ids, $class_id );
  }
  return @class_ids;
}

############################################################

sub get_Clonelib_id_by_name{
  my ($self, $name ) = @_;
  my $q = qq ( SELECT Id 
	       FROM   Clonelib
	       WHERE  Name="$name"
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  my $count = 0;
  my $clone_dbID;
  while ( my ($id) = $sth->fetchrow_array ){
    $clone_dbID = $id;
    $count++;
  }
  if ( $count>1){
    print STDERR "WARNING: 1-to-2 mapping\n";
  }
  return $clone_dbID;
}

############################################################

=head2 get_vocabulary_of_est
  
  Arg        : an EST accession, it will check for version numbers and remove them
  Description: given an EST accession, for each Ontology tree, it finds the full path leading to the 
               node where the EST sits
  Returntype : It returns an array with as many elements as paths the EST have been found in.
               Each element is an arrayref which contains all the terms in the vocabulary,
               starting from the ontology tree name. It does not contain the clone library name, which
               can be retrieved with $self->get_library_by_est($est_id)               
=cut

sub get_vocabulary_of_est{
  my ($self,$est_id) = @_;

  # we get first the clone library for this est:
  my $clone_library = $self->get_library_by_est( $est_id);

  # get the Node ids linking through the clone library id:
  my $q = qq ( SELECT Mapping.NodeId
	       FROM   Mapping, EST
	       WHERE  EST.Accession = "$est_id"
	       AND    EST.ClonelibId = Mapping.ClonelibId
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  
  my @nodes;
  while ( my ($node_id) = $sth->fetchrow_array ){
    push (@nodes, $node_id);
  }

  # now find the whole path to the root of the tree from each of these end nodes
  my @ontologies;
  foreach my $node_id (@nodes){
    
    my @vocabulary;
    my $string  = $self->_find_path($node_id);
    @vocabulary = split '|||',$string;
    push ( @ontologies, \@vocabulary );
  }
  return @ontologies;
}

############################################################

sub _find_path{
  my ($self,$node_id) = @_;
  
  my $parent_node_id = $self->_get_parent_node_id($node_id);
  if ( $parent_node_id == 0 ){
    # we got to the root, get the name of the tree:
    my $ontology = $self->_get_ontology_from_node_id( $node_id);
    return $ontology;
  }
  my $parent_term = $self->find_path($parent_node_id);

  my $term = $self->_get_term_by_node_id($node_id); 
  
  my $string = $parent_term."|||".$term;
  
  return $string;
}

############################################################

sub _get_term_by_node_id{
  my ($self,$node_id) = @_;
  my $q = qq ( SELECT Vocabulary.Term
	         FROM Vocabulary
	        WHERE Vocabulary.NodeId = $node_id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  
  my $term;
  while ( my ($id) = $sth->fetchrow_array ){
    $term = $id;
  }
  return $term;
}

############################################################

sub _get_ontology_from_node_id{
  my ($self,$node_id) = @_;
  my $q = qq ( SELECT Ontology.Name
	         FROM Ontology, Node
	        WHERE Node.Id = $node_id
	          AND Node.OntologyId = Ontology.Id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  
  my $ontology;
  while ( my ($id) = $sth->fetchrow_array ){
    $ontology = $id;
  }
  return $ontology;
}

############################################################

=head2 get_libraryId_by_est
  
  Arg        : an EST Accession ( a string)
  Description: it finds the clone library name in the Clonelib table from which
               this EST was derived. It checks for version numbers and prune them
               as this information is not stored in the expression database
  Returntype : a string hlding the name of the clone library 

=cut

sub get_libraryId_by_est{
  my ($self,$est_id) = @_;
  
  # get rid of version numbers
  if ($est_id =~/(\S+)\.(\d+)/){
    $est_id = $1;
  }

  my $q = qq ( SELECT EST.Accession, EST.ClonelibId
	       FROM   EST
	       WHERE  EST.Accession  = "$est_id"
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  
  my $clone_library;
  while ( my ($id) = $sth->fetchrow_array ){
    $clone_library = $id;
  }
  return $clone_library;
}
############################################################

=head2 get_libraryId_by_estarray
  
  Arg        : an array of EST Accessions
  Description: it finds the clone library id for each est
               this EST was derived. It checks for version numbers and prune them
               as this information is not stored in the expression database
  Returntype : returns an arrray of arrayrefs, each arrayref holds a pair
               with EST.Accession and EST.ClonelibId

=cut

sub get_libraryId_by_estarray{
  my ($self,@est_ids) = @_;
  
  my $est_string = "(";
  foreach my $est_id (@est_ids){
    # get rid of version numbers
    #print STDERR "ExpressionAdaptor: est_id: $est_id\n";
    if ($est_id =~/(\S+)\.(\d+)/){
      $est_id = $1;
    }
    $est_string .= "\"$est_id\"";
    unless ( $est_id eq $est_ids[-1] ){
      $est_string .= ",";
    }
  }
  $est_string .=")";
  
  my $q = qq ( SELECT EST.Accession, EST.ClonelibId
	       FROM   EST
	       WHERE  EST.Accession in $est_string
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  
  my @pairs;
  while ( my @pair = $sth->fetchrow_array ){
    push( @pairs, \@pair );
  }
  return @pairs;
}

############################################################

=head2 _get_parent_node_id
  
 Arg        : node internal id (integer)
 Description: it finds the parent node id for a given node in the tree structure
              of nodes, represented by a table with node and parent node in each row.
              It is used in the recursion to find the complete expression vocabulary for a given
              est.

=cut

sub _get_parent_node_id{
  my ( $self, $node_id ) = @_;
  my $q = qq ( SELECT ParentId
	       FROM   Node
	       WHERE  Id = $node_id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  my $count = 0;
  my $clone_dbID;
  while ( my ($id) = $sth->fetchrow_array ){
    $clone_dbID = $id;
    $count++;
  }
  return $clone_dbID;
}

############################################################

sub get_ests_by_vocabulary{
  my ($self,$term) = @_;
  
}

1;

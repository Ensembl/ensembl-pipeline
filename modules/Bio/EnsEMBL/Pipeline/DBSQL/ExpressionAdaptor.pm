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


#use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnectio);


####################################################################
# 
# Method to get all the ESTs related to an array of Vocabulary Terms
#
#####################################################################

=head2 fetch_ests_by_Term
  
 Arg        : an array of terms in the vocabulary
 Description: it retrieves all the est accessions from the libraries in the ontology tree that are under the Node
              defined by the Term provided in the argument
 Returntype : an array of strings

=cut

sub fetch_ests_by_Term{
  my ($self,@terms) = @_;

  my %lib_arrays;
  foreach my $term (@terms){
    print STDERR "fetching ests for $term\n";
    @{ $lib_arrays{$term} } = $self->get_clonelib_by_Term($term);
  }
  
  # find the intersection
  if ( scalar( @terms ) == 1 ){
    print STDERR "fetching ESTs\n";
    my @est_ids = $self->_get_ests_by_clone_id_array( @{ $lib_arrays{$terms[0]} });
    return @est_ids;
  }
  # this method rely on having no repeated ESTs within each array
  if ( scalar( @terms ) > 1 ){
    print STDERR "finding intersectino of libraries\n";
    my @common_libs;
    my %counter;
    foreach my $term (@terms){
      foreach my $lib ( @{ $lib_arrays{$term} } ){
	$counter{$lib}++;
	if ( $counter{$lib} == scalar(@terms) ){
	  push( @common_libs, $lib );
	}
      }
    }
    unless( @common_libs ){
      print STDERR "No clone libraries with terms: @terms\n";
      exit(0);
    }
    # get all the ests for these clone libraries
    print STDERR "fetching ESTs\n";
    my @est_ids = $self->_get_ests_by_clone_id_array( @common_libs);
    return @est_ids;
  }
}


##############################################################################
# 
# Method to get all the clone library ids related to a given Vocabulary Term
#
##############################################################################

=head2 get_clonelib_by_Term
  
 Arg        : a term in the vocabulary
 Description: it retrieves all the clone library ids in the ontology tree that are under the Node
              defined by the Term provided in the argument
 Returntype : an array of clone library ids (integers)

=cut

sub get_clonelib_by_Term{
  my ($self,$term) = @_;

  #first get the NodeId and the Ontology tree for this term
  my $q = qq ( SELECT Vocabulary.NodeId, Node.OntologyId 
	       FROM   Vocabulary, Node 
	       WHERE  Vocabulary.Term   = "$term" 
	       AND    Vocabulary.NodeId = Node.Id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute      || $self->throw("can't execute: $q");
  my ($node_id, $ontology_id) = $sth->fetchrow_array;
  
  unless ( $ontology_id ){
    print STDERR "term $term is not present in any Tree\n";
  }
  
  # load this ontology tree in memory
  # (this is better than doing recursive SQL queries)
  # but check first whether it has been cached:
  my %tree;
  if ( $self->_tree($ontology_id) ){
    print STDERR "reading a cached tree: $ontology_id\n";
    %tree = %{ $self->_tree($ontology_id) };
  }
  else{
    print STDERR "reading a new tree: $ontology_id\n";
    %tree = $self->get_Tree_by_id($ontology_id);
  }
  
  # get the leaves (clone libraries) of the subtree:
  my @leaves;
  my $ref = \%tree;
  print STDERR "getting libraries with this term\n";
  @leaves = $self->_get_leaves_from_node( $node_id , \@leaves, \%tree);
  
  return @leaves;
}
  

############################################################

=head2 _get_leaves_from_node
  
 Arg        : a node_id (an integer )
 Description: it gets children recursively until it cannot go any further, then it gets all the
              leaves from the last nodes
 Returntype : an array of strings holding the clone library ids

=cut

sub _get_leaves_from_node{
  my ( $self, $node_id, $leaves, $tree ) = @_;
  my %tree = %$tree;
  
  if ( $tree{children}{$node_id}  && @{ $tree{children}{$node_id} } ){
    foreach my $child ( @{ $tree{children}{$node_id} }){
      push( @{ $leaves },  $self->_get_leaves_from_node( $child, $leaves, $tree ) );
    }
  }
  elsif( $tree{libraries}{$node_id} && @{ $tree{libraries}{$node_id} } ){
    push( @{ $leaves }, @{ $tree{libraries}{$node_id} } );
    return @{ $leaves };
  }
}
############################################################

=head2 _get_ests_by_clone_id_array
  
 Arg        : an array of clone library internal ids
 Description: it retrieves all the est accessions from the libraries in this array
 Returntype : an array of strings

=cut

sub _get_ests_by_clone_id_array{
  my ($self, @leaves ) = @_;

  my $string = "(";
  foreach my $leaf (@leaves){
    $string .= "\"$leaf\"";
    unless ( $leaf == $leaves[-1] ){
      $string .= ",";
    }
  }
  $string .=")";
  
  my $q = qq ( SELECT EST.Accession
	       FROM   EST
	       WHERE  EST.ClonelibId in $string
	     ); 
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute      || $self->throw("can't execute: $q");
  my @est_ids;
  while ( my ($est_id) = $sth->fetchrow_array ){
    push( @est_ids, $est_id );
  }
  return @est_ids;
}


########################################################################################
#
# Method to load the entire tree structure for a given ontology given its internal id
#
########################################################################################

=head2 get_Tree_by_id
  
  Arg        : the ontology internal id
 Description : it loads all the information for a tree.
  Returntype : returns a hash which contains all the information necessary to
               reproduce tree
               $tree{term}{$node_id}   is the vocabulary term for that node
               $tree{parent}{$node_id} is the parent node id
               @{ $tree{nodes} }       holds all the node ids in the tree
               @{ $tree{libraries}{$node_id} } holds the leaves ( all the clone libraries hanging from a node )
=cut

sub get_Tree_by_id{
  my ($self,$ontology_id) = @_;

  my %tree;
  
  # get the root of the tree
  my $q = qq ( SELECT Ontology.Name
	       FROM   Ontology
	       WHERE  Ontology.Id = $ontology_id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute      || $self->throw("can't execute: $q");

  $tree{root} = $sth->fetchrow_array;
  
  print STDERR "getting tree structure\n";
  # store terms and parent ids
  my $q = qq ( SELECT Node.Id, Node.ParentId, Vocabulary.Term
	       FROM   Vocabulary, Node
	       WHERE  Vocabulary.NodeId = Node.Id
	       AND    Node.OntologyId = $ontology_id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute      || $self->throw("can't execute: $q");
  while( my ($node_id, $parent_id, $term) = $sth->fetchrow_array){
    $tree{term}{$node_id}   = $term;
    $tree{parent}{$node_id} = $parent_id;
    push( @{ $tree{children}{$parent_id} }, $node_id );  
    push( @{ $tree{nodes} }, $node_id );
  }
  
  print STDERR "getting libraries for this tree\n";
  # get all the leaves in this tree
  my $q = qq ( SELECT Node.Id, Mapping.ClonelibId
	       FROM   Node, Mapping
	       WHERE  Mapping.NodeId  = Node.Id
	       AND    Node.OntologyId = $ontology_id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute      || $self->throw("can't execute: $q");
  my $count = 0;
  while( my ($node_id, $lib_id) = $sth->fetchrow_array){
    $count++;
    push( @{ $tree{libraries}{$node_id} }, $lib_id);
  }
  print STDERR "$count libraries found\n";
  # cache the tree
  $self->_tree($ontology_id,\%tree);
  return %tree;
}

=head2 _tree
  
  Arg        : an ontology_id and optionally the tree ( a hash reference ) created in get_Tree_by_id
 Description : it is a get/set method to cache the tree structure created in get_Tree_by_id
  Returntype : returns a hashref pointing to the tree hash

=cut

sub _tree{
  my ($self,$ontology_id,$tree) = @_;
  unless ( $ontology_id ){
    $self->throw("Cannot retrived a cached tree without an ontology_id");
  }
  if ( $ontology_id && $tree ){
    $self->{_tree}{$ontology_id} = $tree;
  }
  return $self->{_tree}{$ontology_id};
}

########################################################################################
#
# Method to check which ests (from an array) are present in (at least) one of the Vocabulary trees
#
########################################################################################

=head2 get_libraryId_by_estarray
  
  Arg        : an array of EST Accessions
  Description: it finds the clone library id for each est
               this EST was derived. It checks for version numbers and prune them
               as this information is not stored in the expression database
  Returntype : returns an arrray of arrayrefs, each arrayref holds a pair
               with EST.Accession and EST.ClonelibId
  Caller     : It is called in Bio::EnsEMBL::Pipeline::RunnableDB::MapGeneToExpression

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
#
# Method to get all the vocabulary associated to an EST
#
############################################################

=head2 get_Vocabulary_of_est
  
  Arg        : an EST accession (a string), it will check for version numbers and remove them
  Description: given an EST accession, for each Ontology tree where the library can be found, 
               it finds the full path from the root to the node where the EST sits
  Returntype : It returns an array with as many elements as paths the EST have been found in,
               one per ontology tree where it was found in.
               Each element is an arrayref which contains all the terms in the vocabulary,
               starting from the ontology tree name. It does not contain the clone library name, which
               can be retrieved with $self->get_library_by_est($est_id)               
=cut

sub get_Vocabulary_of_est{
  my ($self,$est_id) = @_;

  # remove version numbers
  if ( $est_id =~/(\S+)\.(\d+)/ ){
    $est_id = $1;
  }

  # get the Node ids linking through the clone library id for this est:
  my $q = qq ( SELECT Node.OntologyId, Mapping.NodeId
	       FROM   Mapping, EST, Node
	       WHERE  EST.Accession  = "$est_id"
	       AND    EST.ClonelibId = Mapping.ClonelibId
	       AND    Mapping.NodeId = Node.Id 
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute      || $self->throw("can't execute: $q");
  
  my @nodes;
  my %ontology;
  while ( my ($ontology_id, $node_id) = $sth->fetchrow_array ){
    push (@nodes, $node_id);
    $ontology{$node_id} = $ontology_id;
  }
  unless( @nodes ){
    print STDERR "est $est_id not found in the tree\n";
  }

  # now find the whole path to the root of the tree from each of these end nodes
  my @vocabularies;
  foreach my $node_id (@nodes){
    my @vocabulary;
    my %tree;
    if ( $self->_tree( $ontology{$node_id}) ){
      %tree = %{ $self->_tree( $ontology{$node_id} ) };
    }
    else{
      %tree = $self->get_Tree_by_id($ontology{$node_id}) ;
    }
    $self->_find_path($node_id, \@vocabulary, \%tree);
    push ( @vocabularies, \@vocabulary );
  }
  return @vocabularies;
}

############################################################

=head2 _find_path
  
  Arg        : a node_id (an integer)
  Description: recursive method to find the path from one clonelib (a leave in the tree)
               all the way up to the root of the ontology tree
  Returntype : a string with the vocabulary terms found at a given stage of the recursion. 
               As it finds recursively the terms along the path, it concatenates them into a string.
  Caller     : it is only called from get_vocabulary_of_est

=cut

sub _find_path{
  my ($self,$node_id,$vocabulary, $tree) = @_;

  my %tree = %$tree;

  # recall that
  # $tree{root} is the name of the Onotology tree (the root)
  # $tree{term}{$node_id}   is the vocabulary term for that node
  # $tree{parent}{$node_id} is the parent node id
  # @{ $tree{nodes} }       holds all the node ids in the tree
  # @{ $tree{libraries}{$node_id} } holds the leaves ( all the clone libraries hanging from a node )

  my $parent_node_id = $tree{parent}{$node_id};

  if ( $parent_node_id == 0 || !defined($parent_node_id)){

    # we got to the root, get the name of the tree:
    my $root = $tree{root};
    push ( @{ $vocabulary }, $root );
    
    my $term = $tree{term}{$node_id};
    push ( @{$vocabulary}, $term);
    
    return;
  }
  $self->_find_path($parent_node_id, $vocabulary, $tree);
  
  my $child_term = $tree{term}{$node_id};
  push ( @{$vocabulary}, $child_term );
  
  return;
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

=head2 get_library_Name_by_est
  
  Arg        : an EST Accession ( a string)
  Description: it finds the clone library name in the Clonelib table from which
               this EST was derived. It checks for version numbers and prune them
               as this information is not stored in the expression database
  Returntype : a string holding the name of the clone library 

=cut

sub get_library_Name_by_est{
  my ($self,$est_id) = @_;
  
  # get rid of version numbers
  if ($est_id =~/(\S+)\.(\d+)/){
    $est_id = $1;
  }

  my $q = qq ( SELECT Clonelib.Name
	       FROM   EST, Clonelib
	       WHERE  EST.Accession  = "$est_id"
	         AND  EST.ClonelibId = CLonelib.Id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute      || $self->throw("can't execute: $q");
  
  my $library_name = $sth->fetchrow_array;
  return $library_name;
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
  my $res = $sth->execute      || $self->throw("can't execute: $q");
  my $parent_id;
  while ( my ($id) = $sth->fetchrow_array ){
    $parent_id = $id;
  }
  return $parent_id;
}

############################################################

=head2 get_Ontologies
  
 Arg        : none
 Description: it retrieves all the available ontologies in the SANBI database
 Returntype : an array of strings

=cut

sub get_Ontologies{
  my ($self) = @_;
  my $q = qq ( SELECT Name
	       FROM   Ontology
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  my @ontologies;
  while ( my ($name) = $sth->fetchrow_array ){
    push( @ontologies, $name );
  }
  return @ontologies;
}


############################################################

=head2 get_Vocabulary_by_Ontology
  
 Arg        : a string = one of the available ontologies
              'Anatomical Site', 'Development Stage', 'Pathology', 'Cell Type', 'Preparation'
 Description: it retrieves all the available vocabulary terms within the provided ontology
 Returntype : an array of strings

=cut

sub get_Vocabulary_by_Ontology{
  my ($self,$ontology) = @_;
  # $ontology should be one of these:
  # 'Anatomical Site', 'Development Stage', 'Pathology', 'Cell Type', 'Preparation'

  my $q = qq ( SELECT Vocabulary.Term
	       FROM   Vocabulary, Node, Ontology
	       WHERE  Ontology.Name     = "$ontology"
	         AND  Ontology.Id       = Node.OntologyId
	         AND  Vocabulary.NodeId = Node.Id
	     );
  
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  my $res = $sth->execute || $self->throw("can't execute: $q");
  my @vocabulary;
  while ( my ($term) = $sth->fetchrow_array ){
    push( @vocabulary, $term );
  }
  return @vocabulary;
}

########################################
#
# Methods associated to storage of data
#
########################################

sub store_ensembl_link{
  my ($self, $transcript_id, $est_ids) = @_;
  my @est_ids = @$est_ids;
  my $q = qq(
	     INSERT into EnsemblTranscript ( Id, ESTAccession )
	     VALUES ( ?, ? ) 
	   );
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  foreach my $est_id ( @est_ids ){
    $sth->execute($transcript_id, $est_id) || $self->throw("can't execute: $q");
  }
}

############################################################

sub store_est_transcript_link{
  my ($self, $transcript_id, $est_ids) = @_;
  my @est_ids = @$est_ids;
  my $q = qq(
	     INSERT into EstTranscript ( Id, ESTAccession )
	     VALUES ( ?, ? ) 
	   );
  my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
  foreach my $est_id ( @est_ids ){
    $sth->execute($transcript_id, $est_id) || $self->throw("can't execute: $q");
  }
}
  


1;

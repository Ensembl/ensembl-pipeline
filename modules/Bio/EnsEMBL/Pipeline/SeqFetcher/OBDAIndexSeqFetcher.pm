#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher

=head1 SYNOPSIS

    my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(
										  -db     => $db,
										  -format => $format,
										 );
    my $seq = $seqfetcher->get_Seq_by_acc($acc);

       where $acc is the primary key on which the index has been made (accession or id)

    my $seq = $seqfetcher->get_Seq_by_secondary($name,$acc);

       where $name is the namespace or identifier for the secondary key, and $acc is the accession or id.
       Beware that this method can return multiple sequences as the secondary id is not necessarily unique

=head1 DESCRIPTION

  This is basically a wrapper around the SeqFetcher Bio::DB::Flat::OBDAIndex,
to use it in Ensembl in the same way as other SeqFetcher. 
It reads some configuration info from pipeConf.pl.
 
  Sequences are fetched from a 
database previously formatted with indicate (made by Steve Searle)

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;

use strict;
use Bio::EnsEMBL::Root;
use Bio::DB::RandomAccessI;
use Bio::Seq;
use Bio::DB::Flat::OBDAIndex;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  
  my ($db, $format) = $self->_rearrange(['DB', 'FORMAT'], @args);
  
  # expect an array of dbs
  $self->throw("Expected a reference to an array of db\n") unless ref($db) eq 'ARRAY';

  # $db is a reference to an array of databases:
  if ($db) {
    $self->db($db);
  }
  else{
    $self->throw("Sorry, you must specify a database");
  }
  
  foreach my $database ( $self->db ){
    # we want the form '/path/to/index/dir', so remove the last if there is any '/'
    if ( $database =~/(\S+)\/$/ ){
      $database = $1;
    }

    $self->index_name( $database);

    # get the index name and the index directory out of $database
    my @path = split /\//, $database;

    # the db_name is the last name in the path
    my $db_name = pop( @path );
    if ( $db_name =~/(\S+)\.fa/){
      $db_name = $1;
    }
    if ( !defined( $db_name ) ){
      $self->throw("cannot define db_name");
    }

    # the index_dir is the directory path where the db sits
    my $index_dir = join '/', @path;
    if ( !defined( $index_dir ) ){
      $self->throw("cannot define index_dir");
    }

    # we take as default format = 'FASTA';
    unless ( $format ){
      $format = 'FASTA';
    }


    my $OBDAfetcher = new Bio::DB::Flat::OBDAIndex(-index_dir => $index_dir,
						   -dbname    => $db_name,
						   -format    => $format
						 );

    $self->_seqfetcher($OBDAfetcher);
  }

  return $self;

}

=head2 _seqfetcher

  Function: stores and retrieves an Bio::DB::Flat::OBDAIndex object which is initialised in new()

=cut

sub _seqfetcher{
  my ($self, $fetcher) = @_;
  if ( $fetcher ){
    push( @{ $self->{'_seqfetcher'} }, $fetcher);
  }
  return @{ $self->{'_seqfetcher'} };
}


=head2 db

  Title   : db
  Usage   : $self->db('/data/blastdb/dbname');
  Function: Get/set for dbs to be searched. Checks that the database 
            appropriate files are present, but nothing else.
  Returns : string
  Args    : string

=cut

sub db {

  my ($self, $dbs) = @_;
  if (!defined($self->{'_db'})) {
    $self->{'_db'} = [];
  }
  if (defined $dbs){
    if (ref($dbs) eq 'ARRAY') {
      foreach my $db(@$dbs){
	push (@{$self->{'_db'}},$db);
      }
    }
  }
  return (@{$self->{'_db'}});  
  
}

=head2 get_Seq_by_acc

  Function: Does the sequence retrieval via the OBDAIndex module
  Returns : Bio::Seq

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;

  if (!defined($acc)) {
    $self->throw("No accession input");
  }  
  
  my $seq;
  my @seqfetchers = $self->_seqfetcher;

  foreach my $seqfetcher (@seqfetchers){
    
    eval{
      # note that this only works for OBDAIndex.pm
      $seq = $seqfetcher->get_Seq_by_id($acc);
    };
    if ( $@ ){
      $self->warn("problem fetching sequence for $acc");
    }
    
    if ( $seq ){
      $seq->display_id( $acc );
      $seq->accession_number( $acc );
      $seq->desc("");
      last;
    }
  }

  if(!$seq){
    my ($p, $f, $l) = caller;
    $self->warn("OBDAIndexSeqFetcher: could not find sequence for primary key $acc in index ".$self->index_name." $f:$l\n");
    
  FETCHER:
    foreach my $seqfetcher ( $self->_seqfetcher ){
      
      my @secondary_namespaces = $seqfetcher->secondary_namespaces;
      foreach my $name ( @secondary_namespaces ){
	$self->warn("seqfetcher $seqfetcher is looking in namespace $name for secondary key $acc\n");
	
	my @seqs;
	eval{
	  # this returns potentially an array of Bio::Seq
	  @seqs = $seqfetcher->get_Seq_by_secondary($name,$acc);
	};
	if ( $@ ){
	  $self->warn("problem fetching sequence for secondary key $acc");
	}
	if ( @seqs > 1 ){
	  $self->warn("Multiple sequences (".scalar(@seqs).") for the same secondary accession $acc\n");
	  next;
	}
	
	if ( defined $seqs[0] ){
	  $seqs[0]->display_id( $acc );
	  $seqs[0]->accession_number( $acc );
	  $seqs[0]->desc("");
	  $seq = $seqs[0];
	  last FETCHER;
	}
      }
    }
    unless ($seq){
      $self->warn("could not find sequence for secondary key $acc");
    }
  }

  if ($seq){
    #print STDERR "found sequence!\n";
    #print STDERR "OBDAIndexSeqFetcher: returning sequence:\n";
    #print STDERR "display_id: ".$seq->display_id."\n";
    #print STDERR $seq->seq."\n";
  }
  else{
    print STDERR "sequence not found. Returning undef\n";
  }

  return $seq;
  
}


sub  get_Seq_by_id {
  my $self = @_;
  $self->warn("cannot call  get_Seq_by_id on OBDAIndexSeqFetcher, use get_Seq_by_acc instead");
  return undef;
}

=head2 get_Seq_by_secondary

  Function: Does the sequence retrieval via the OBDAIndex module using the secondary index key
            An index should have been made prior to this on the two keys.
  Returns : Bio::Seq

=cut

sub  get_Seq_by_secondary {
  my ($self, $name, $acc) = @_;

  if (!defined($acc)) {
    $self->throw("No secondary key input");
  }  
  if (!defined($name)){
    $self->throw("No name space for the secondary key");
  }

  my @seqs;
  my @seqfetchers = $self->_seqfetcher;

  foreach my $seqfetcher (@seqfetchers){
    
    eval{
      # this returns potentially an array of Bio::Seq
      @seqs = $seqfetcher->get_Seq_by_secondary($name,$acc);
    };
    if ( $@ ){
      $self->warn("problem fetching sequence for $acc");
    }
    
    if ( @seqs > 1 ){
      $self->warn("Multiple sequences (".scalar(@seqs).") for the same secondary accession $acc\n");
      next;
    }
    if ( defined $seqs[0] ){
      $seqs[0]->display_id( $acc );
      $seqs[0]->accession_number( $acc );
      $seqs[0]->desc("");
      last;
    }
  }
  
  unless (@seqs){
    $self->warn("OBDAIndexSeqFetcher: could not find sequence for $acc");
  }
  
  #print STDERR "OBDAIndexSeqFetcher: returning sequence:\n";
  #print STDERR "display_id: ".$seqs[0]->display_id."\n";
  #print STDERR $seqs[0]->seq."\n";
  return $seqs[0];
}

sub index_name{
  my ($self,$name) = @_;
  if ($name){
    $self->{index_name} = $name;
  }
  return $self->{index_name};
}


1;

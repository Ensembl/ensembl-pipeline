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

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(
									   -db     => $db,
									   -format => $format,
									  );
    my $seq = $obj->get_Seq_by_acc($acc);

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
use Bio::Root::RootI;
use Bio::DB::RandomAccessI;
use Bio::Seq;
use Bio::DB::Flat::OBDAIndex;

use vars qw(@ISA);

@ISA = qw(Bio::Root::RootI Bio::DB::RandomAccessI);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;
  
  my ($db, $format) = $self->_rearrange(['DB', 'FORMAT'], @args);
  
  # expect an array of dbs
  print STDERR "you have passed a ".ref($db)."\n";
$self->throw("Expected a reference to an array of db\n") unless ref($db) eq 'ARRAY';

  # $db is a reference to an array of databases:
  if (defined $db) {
    $self->db($db);
  }
  else{
    $self->throw("Sorry, you must specify a database");
  }
  
  foreach my $database ( $self->db ){

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

    print STDERR "making a OBDAIndex fetcher with db_name = $db_name, index_dir = $index_dir, format = $format\n"; 

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
    
    if ( defined $seq ){
      $seq->display_id( $acc );
      $seq->accession_number( $acc );
      $seq->desc("");
      last;
    }
  }
  if(!defined $seq){
    $self->warn("could not find sequence for $acc");
  }
  
  return $seq;
}


sub  get_Seq_by_id {
  my $self = @_;
  $self->warn("cannot call  get_Seq_by_id on OBDAIndexSeqFetcher, use get_Seq_by_acc instead");
  return undef;
}


1;

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

  Bio::EnsEMBL::Pipeline::SeqFetcher::BioIndex

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::SeqFetcher::BioIndex->new(
								-db    => $db
							     );
    my $seq = $obj->get_Seq_by_acc($acc);

=head1 DESCRIPTION

  Object to retrieve sequences as Bio::Seq, using getseqs by James Cuff. Sequences are fetched from a 
database previously formatted with makeindex

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Pipeline::SeqFetcher::BioIndex;

use strict;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::SeqFetcher::DBIndex;
use Bio::Seq;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($db) = $self->_rearrange(['DB'], @args);

  # expect an array of dbs
  $self->throw("Expected a reference to an array of db [" . ref($db) . "]\n") unless ref($db) eq 'ARRAY';
  if (defined $db) {
    $self->db($db);
  }
  
  return $self; # success - we hope!
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
	$self->throw("are you sure that $db has been formatted with DBIndex? Some files are missing\n") 
	  unless (-d "$db" && -e "$db/config.dat" && -e "$db/fileids.dat");
	
	my $index =  new Bio::EnsEMBL::Pipeline::SeqFetcher::DBIndex(-database => $db);

	$self->index($index);

	push (@{$self->{'_db'}},$db);
      }
    }
  }
  return (@{$self->{'_db'}});  
  
}

=head2 index

 Title   : index
 Usage   : Get/set function for storing the Bio indexes
 Function: Is only called internally to tthe object
 Example : $seqfetcher->index
 Returns : array of Bio::EnsEMBL::Pipeline::SeqFetcher::Index
 Args    : nothing or Bio::EnsEMBL::Pipeline::SeqFetcher::Index


=cut


sub index {
  my ($self,$arg) = @_;


  if (!defined ($self->{_index})) {
    $self->{_index} = [];
  }

  if (defined($arg)) {
    push(@{$self->{_index}},$arg);
  }
  return @{$self->{_index}};
}


=head2 get_Seq_by_acc

  Title   : get_Seq_by_acc
  Usage   : $self->get_Seq_by_acc($accession);
  Function: Does the sequence retrieval via getseqs
  Returns : Bio::Seq
  Args    : 

=cut

sub  get_Seq_by_acc {
  my ($self, $acc) = @_;

  if (!defined($acc)) {
    $self->throw("No accession input");
  }  
  
  my $seq;

  foreach my $index ($self->index) {

    my $entry = $index->get_Seq_by_id($acc);

    if ($entry ne "") {
      print "Entry\n";
      my ($idstr,$seqstr) = split(/\n/,$entry,2);
      my ($id,$desc) = split(/ /,$idstr,2);

      $id =~ s/^>//;
      
      $seqstr =~ s/\n//g;
      
      my $seq =  new Bio::Seq(-id      => $acc,
			      -display_id => $acc,
			      -seq    => $seqstr);
      return $seq;
    }
  }
  return;
}

sub  get_Seq_by_id {
  my $self = @_;
  return undef;
}

1;

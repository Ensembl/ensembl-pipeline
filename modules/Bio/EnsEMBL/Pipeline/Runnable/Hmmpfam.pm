
#
# Ensembl module for Bio::EnsEMBL::Pipeline::Hmmpfam
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Hmmpfam - Runnable for running Hmmpfam

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Runnable::Hmmpfam;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Analysis::Programs qw(hmmpfam); 
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::SeqIO;
use Bio::Tools::HMMER::Results;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($peptide,$database) = 
      $self->_rearrange([qw(PEPTIDE DATABASE)], @args);
  
  if( !defined $peptide || !ref $peptide || !$peptide->isa('Bio::PrimarySeqI') ) {
      $self->throw("No peptide input");
  }
  $self->peptide($peptide);

  if( !defined $database ) {
      $database = 'Pfam';
  }

  $self->database($database);
  $self->{'_sf'} = [];
  return $self;
}


=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run{
   my ($self,@args) = @_;

   my $seqfile = "/tmp/hmmpfam.$$.seq";
   my $seqout = Bio::SeqIO->new( -file => ">$seqfile", -format => 'fasta' );
   $seqout->write_seq($self->peptide);
   $seqout = undef;
   my $db = $self->database;

   open(HMMER,"hmmpfam $db $seqfile |");
   
   my $res = new Bio::Tools::HMMER::Results( -fh => \*HMMER , -type => 'hmmpfam');

   close(HMMER) || $self->throw("Error in running hmmpfam $db $seqfile $!");
   
   foreach my $domain ( $res->each_Domain ) {
       $self->add_SeqFeature($domain);
   }

   unlink($seqfile);

}

=head2 output

 Title   : output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub output{
   my ($self,@args) = @_;
   
   return @{$self->{'_sf'}};
}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_SeqFeature{
   my ($self,$sf) = @_;

   if( !defined $sf ) {
       $self->throw("You have not passed in a sequence feature...");
   }

   push(@{$self->{'_sf'}},$sf);
}


=head2 peptide

 Title   : peptide
 Usage   : $obj->peptide($newval)
 Function: 
 Returns : value of peptide
 Args    : newvalue (optional)


=cut

sub peptide{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'peptide'} = $value;
    }
    return $obj->{'peptide'};

}

=head2 database

 Title   : database
 Usage   : $obj->database($newval)
 Function: 
 Returns : value of database
 Args    : newvalue (optional)


=cut

sub database{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'database'} = $value;
    }
    return $obj->{'database'};

}

1;

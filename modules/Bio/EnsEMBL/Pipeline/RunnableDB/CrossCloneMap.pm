
#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap.pm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap.pm - DESCRIPTION of Object

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


package Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap.pm;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDBI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDBI Bio::Root::RootI);

sub new {
  my($class,@args) = @_;

  my $self = {};

  $self->{'_crossmatch'} = [];
  bless $self,$class;
  
  return $self;
}


=head2 fetch_input

 Title   : fetch_input
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_input{
   my ($self,$input_id) = @_;

   # connects to donor and acceptor database for this clone

   my $new = $self->dbobj->new_dbobj;
   my $old = $self->dbobj->old_dbobj;

   my $newclone = $new->get_Clone($input_id);
   my $oldclone = $old->get_Clone($input_id);

   if( $newclone->version == $oldclone->version ) {
       $self->throw("Input $input_id has identical versioned clones in each database");
   }

   my @newcontigs = $newclone->get_all_Contigs();
   my @oldcontigs = $oldclone->get_all_Contigs();
   
   my @newseq;
   my @oldseq;

   foreach my $contig ( @newcontigs ) {
       my ($acc,$number) = split(/\./,$contig->id);
       my $version = $newclone->version;

       my $seq = Bio::PrimarySeq->new( -display_id => "$acc.$version.$number", -seq => $contig->seq);
       push(@newseq,$seq);
   }

   foreach my $contig ( @oldcontigs ) {
       my ($acc,$number) = split(/\./,$contig->id);
       my $version = $oldclone->version;

       my $seq = Bio::PrimarySeq->new( -display_id => "$acc.$version.$number", -seq => $contig->seq);
       push(@oldseq,$seq);
   }


   foreach my $newseq ( @newseq ) {
       foreach my $oldseq ( @oldseq ) {
	   my $cross = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( 
									  -nocopy => 1,
									  -seq1 => $newseq,
									  -seq2 => $oldseq,
									  );

	   push(@{$self->{'_crossmatch'}},$cross);
       }
   }

   
   return;
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
   my ($self) = @_;

   my @fp;

   
}



=head2 write_output

 Title   : write_output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_output{
   my ($self,@args) = @_;


}


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


package Bio::EnsEMBL::Pipeline::RunnableDB::CrossCloneMap;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  $self->{'_final'} = [];
  $self->{'_crossmatch'} = [];

  my ($cross_db,$score) = $self->_rearrange([qw(CROSSDB SCORE)],@args);

  if ((!$cross_db) || 
      (!$cross_db->isa('Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor'))) {
      $self->throw("You need to provide a CrossMatch database adaptor to run CrossCloneMap!");
  }
  $self->_score($score);

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

   $self->_acc($input_id);
   # connects to donor and acceptor database for this clone
   #This is a crossmatch DB adaptor
   my $new = $self->dbobj->new_dbobj;
   my $old = $self->dbobj->old_dbobj;

   my $newclone = $new->get_Clone($input_id);
   my $oldclone = $old->get_Clone($input_id);

   if( $newclone->embl_version == $oldclone->embl_version ) {
       $self->throw("Input $input_id has identical versioned clones in each database");
   }

   my @newcontigs = $newclone->get_all_Contigs();
   my @oldcontigs = $oldclone->get_all_Contigs();
   
   my @newseq;
   my @oldseq;

   foreach my $contig ( @newcontigs ) {
       $contig->id =~ /\S+\.(\d+)/;
       my $number = $1;
       my $version = $newclone->embl_version;
       my $seq = Bio::PrimarySeq->new( -display_id => "$version\_$number", -seq => $contig->seq);
       push(@newseq,$seq);
   }

   foreach my $contig ( @oldcontigs ) {
       $contig->id =~ /\S+\.(\d+)/;
       my $number = $1;
       my $version = $oldclone->embl_version;
       my $seq = Bio::PrimarySeq->new( -display_id => "$version\_$number", -seq => $contig->seq);
       push(@oldseq,$seq);
   }
   print STDERR "Creating crossmatches for clone ".$self->_acc."\n";
   foreach my $newseq ( @newseq ) {
       foreach my $oldseq ( @oldseq ) {
	   my $cross = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new( 
									  -nocopy => 1,
									  -seq1 => $newseq,
									  -seq2 => $oldseq,
									  -score => $self->_score(),
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
   
   #The feature pair array will represent the matrix of hits between inputs
   my @fp;
   #Run all the crossmatch runnables created in fetch_input
   print STDERR "Running crossmatches for clone ".$self->_acc."\n";
   foreach my $crossmatch (@{$self->{'_crossmatch'}}) {
       #print STDERR "Running crossmatch on clone ".$self->_acc." between new seq id ".$crossmatch->seq1->id." and old seq id ".$crossmatch->seq2->id."\n";
       $crossmatch->run();
       #Push all the feature pairs into the array
       push (@fp,$crossmatch->output);
   }
   
   #Sort the array by score
   my @sorted= sort { $a->score <=> $b->score} @fp;

   my $initial=1;
   my %looks_ok=();

   print STDERR "Analysing output for clone ".$self->_acc."\n";
   #The juicy part: look at the matrix, and sort out the mapping
   FP: foreach my $fp (@sorted) {
       #print STDERR "Got feature pair between contig ".$fp->seqname." (".$fp->start."-".$fp->end.") and contig ".$fp->hseqname." (".$fp->hstart."-".$fp->hend.") with score ".$fp->score."\n";
       #Take the first one as correct, because it hsa the highest score...
       if ($initial) {
	   #print STDERR "Pushing it to final (first fp)\n";
	   push (@{$looks_ok{$fp->seqname}},$fp);
	   $initial=undef;
	   next FP;
       }
       
       #Check if this feature pair is consistent with the rest
       
       #If seqname already in final map, check other matches
       
       if ($looks_ok{$fp->seqname}) {
	   my @ha_match=@{$looks_ok{$fp->seqname}};
	   #print STDERR "Contig already in final map, checking...\n";
	   #sort other matches by start
	   my @s_matches= sort { $a->start <=> $b->start} @ha_match;
	   my $first=$s_matches[0];
	   
	   #Speed up by eliminating the two most obvious cases...
           #If the match is before the first match on this contig, add to final
	   if ($fp->end < $first->start){
	       #print STDERR "Pushing it to final (before first)\n";
	       push (@{$looks_ok{$fp->seqname}},$fp);
	       next FP;
	   }
	   #If the match is after the last match on this contig, add to final 
	   my $size=scalar(@s_matches);
	   if ($fp->start > $s_matches[$size]) {
	       #print STDERR "Pushing it to final (after last)\n";
	       push (@{$looks_ok{$fp->seqname}},$fp);
	       next FP;
	   } 
	   
	   #Loop through all the other matches and check if $fp does not
	   #overlap any of them
	   my $add=1;
	   foreach my $match (@s_matches) {
	       #If fp start or end are contained in any match, overlapping...
	       if ((($fp->start > $match->start) && ($fp->start < $match->end)) || (($fp->end < $match->end) && ($fp->end < $match->start))) {
		   $add=0;
	       }
	   }
	   if ($add) {
	       #print STDERR "Pushing it to final (no overlaps)\n";
	       push (@{$looks_ok{$fp->seqname}},$fp);
	       next FP;
	   } 
	   #In any other case, do not add this feature pair
	   #print STDERR "End of checks for this contig, not added to final\n";
       }
       
       #Else, just add to final map
       else {
	   #print STDERR "This seqname was not found earlier, add to final map\n";
	   push (@{$looks_ok{$fp->seqname}},$fp);
       }
   }
   my @final=();
   foreach my $seqname (keys %looks_ok) {
       foreach my $fp (@{$looks_ok{$seqname}}) {
	   my $sn=$self->_acc.".".$fp->seqname;
	   my $hsn=$self->_acc.".".$fp->hseqname;
	   $sn =~ s/\_/\./g;
	   $hsn =~ s/\_/\./g;
	   $fp->seqname($sn);
	   $fp->hseqname($hsn);
	   push (@final,$fp);
       }
   }
   %looks_ok=();
   
   #foreach my $fp (@final) {
       #print STDERR "In final, got ".$fp->seqname."-".$fp->hseqname."\n";
       
   #}
   push(@{$self->{'_final'}},@final);
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

   print STDERR "Writing output for clone ".$self->_acc."\n";
   my @fp=@{$self->{'_final'}};
   my $fc=$self->dbobj->get_SymmetricContigFeatureContainer;
   $fc->write_FeaturePair_List(@fp);
   return;
}

=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj
 Args    : newvalue (optional)

=head2 score

 Title   : score
 Usage   : $obj->score($newval)
 Function: 
 Example : 
 Returns : value of score
 Args    : newvalue (optional)


=cut

sub _score{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_score'} = $value;
    }
    return $self->{'_score'};

}


=head2 _acc

 Title   : _acc
 Usage   : $obj->_acc($newval)
 Function: 
 Example : 
 Returns : value of _acc
 Args    : newvalue (optional)


=cut

sub _acc{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_acc'} = $value;
    }
    return $self->{'_acc'};

}

1;


# Spangle version of the ensembl Transcript object

# POD documentation - main docs before the code

=head1 NAME

MiniSeq - an artificially constructed cDNA on a genomic sequence

=head1 SYNOPSIS

This module is used when we only want to run an analysis over
a part or multiple parts of a sequence. 

=head1 DESCRIPTION

Contains details of coordinates of all exons that make
up a gene transcript.

Creation:
   
     my $mini = new Bio::EnsEMBL::Pipeline::MiniSeq(-id      => $id,
						    -pairaln => $pairaln);



   $pairaln is a Bio::EnsEMBL::Analysis::PairAlign object containing
   one or more feature pairs that represent the mapping of
   the genomic coords to the cDNA coords

Manipulation:

    my $align   = $mini->get_PairAlign;
    my $cdnaseq = $mini->get_cDNA_sequence;

# We now do some analysis on the cdnaseq that puts sequence
# features on it.  We don't want the coordsin the cDNA frame so
# we now pass them back to the miniseq to convert them into 
# genomic coordinates.

    my  @newfeatures = $mini->convert_FeaturePairs(@featurepairs);


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::MiniSeq;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::PairAlign;
use Bio::PrimarySeq;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($id,$pairalign) = $self->_rearrange([qw(ID 
					      PAIRALIGN)],
					  @args);

  $self->throw("No input id for MiniSeq")        unless defined($id);
  $self->throw("No input pairalign for MiniSeq") unless defined($pairalign);

  $self->throw("[$pairalign] is not a Bio::EnsEMBL::Analysis::PairAlign") 
      unless $pairalign->isa("Bio::EnsEMBL::Analysis::PairAlign");
  
  $self->id($id);
  $self->pairAlign($pairalign);

  return $self;
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {
    my ($self,$arg) = @_;

    if(defined($arg)) {
	$self->{'_id'} = $arg;
    }

    return $self->{'_id'};
}

=head2 pairAlign

 Title   : pairAlign
 Usage   : $self->pairAlign($pair)
 Function: Get/set method for the pairalign object
           that stores the cDNA-genomic exon mapping
 Returns : Bio::EnsEMBL::Analysis::PairAlign
 Args    : 


=cut

sub pairAlign {
   my ($self,$pair) = @_;

   if (defined($pair)) {
       if( ! $pair->isa("Bio::EnsEMBL::Analysis::PairAlign") ) {
	   $self->throw("$pair is not a Bio::EnsEMBL::Analysis::PairAlign!");
       }

       foreach my $p ($pair->eachFeaturePair) {
	   if ($p->strand != 1) {
	       $self->throw("Can't have a PairAlign object where the strand of the first object is reversed");
	   }
       }


       $self->{'_pair'} = $pair;
   }
   
   return $self->{'_pair'};
   
}


=head2 get_cDNA_sequence

 Title   : get_cDNA_sequence
 Usage   : my $seq = $self->get_cDNA_sequence
 Function: Returns the cdna sequence corresponding
           to the cDNA in the pairAlign object
 Example : 
 Returns : Bio::PrimarySeq
 Args    : none


=cut

sub get_cDNA_sequence {
   my ($self) = @_;

   my $seqstr = "";

   my @exons = $self->pairAlign->eachFeaturePair;
   
   return unless (scalar @exons > 0);
   
   ##print STDERR "exons:\n"; 
   foreach my $exon (@exons) {
     ##print STDERR "Have ".$exon."\n";
     #print STDERR $exon->start."-".$exon->end."\n";
     #print STDERR $seqstr."\n";
     $seqstr .= $exon->seq;
   }
   return new Bio::PrimarySeq('-id' => "genomic" ,
                              -seq => $seqstr);
   
}

=head2 convert_FeaturePair

 Title   : convert_FeaturePair
 Usage   : my @newfeatures = $self->convert_FeaturePairs($feature)
 Function: Converts feature pair coordinates on the cDNA sequence
           into an array of feature pairs on the genomic sequence
 Example : 
 Returns : Bio::EnsEMBL::FeaturePair
 Args    : Array of Bio::EnsEMBL::FeaturePair


=cut

sub convert_FeaturePair {
    my ($self,$feature) = @_;

    my @newfeatures;

    my @tmp = $self->pairAlign->convert_FeaturePair($feature);
    push(@newfeatures,@tmp);

    return @newfeatures;
}

=head2 convert_SeqFeature

 Title   : convert_FeaturePair
 Usage   : my @newfeatures = $self->convert_FeaturePairs($feature)
 Function: Converts feature coordinates on the cDNA sequence
           into an array of features on the genomic sequence
 Example : 
 Returns : Bio::EnsEMBL::FeaturePair
 Args    : Array of Bio::EnsEMBL::FeaturePair


=cut

sub convert_SeqFeature {
    my ($self,$feature) = @_;

    my @newfeatures;

    my @tmp = $self->pairAlign->convert_cDNA_feature($feature);
    push(@newfeatures,@tmp);

    return @newfeatures;
}

=head2 convert_PepFeaturePair

 Title   : convert_PepFeaturePair
 Usage   : my @newfeatures = $self->convert_PepFeaturePair($feature)
 Function: Converts feature pair coordinates on the cDNA sequence
           into an array of feature pairs on the genomic sequence
           Peptide coordinates are maintained.
 Example : 
 Returns : Array of Bio::EnsEMBL::FeaturePair
 Args    : Bio::EnsEMBL::FeaturePair


=cut

sub convert_PepFeaturePair {
    my ($self,$feature) = @_;

    my @newfeatures;

    my @tmp = $self->pairAlign->convert_FeaturePair($feature);

    # replace protein coordinates
    $tmp[0]->hstart($feature->hstart);
    $tmp[0]->hend($feature->hend);
    push(@newfeatures,@tmp);

    return @newfeatures;
}

1;

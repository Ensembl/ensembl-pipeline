#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils - 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::DnaPepAlignFeature;

@ISA = qw(Bio::EnsEMBL::Root);

 



###########################################################
#
# METHODS DOING CHECKS
#
###########################################################

sub _check_Transcript{
    my ($self,$transcript) = @_;
    my $slice = $self->query;
    
    my $id = $self->transcript_id( $transcript );
    my $valid = 1;
    
    # check that transcripts are not completely outside the slice
    if ( $transcript->start > $slice->length || $transcript->end < 1 ){
	print STDERR "transcript $id outside the slice\n";
	$valid = 0;
    }
    # allow transcripts that fall partially off the slice only at one end, the 'higher' end of the slice
    elsif ( $transcript->start < 1 && $transcript->end > 1 ){
	print STDERR "transcript $id falls off the slice by its lower end\n";
	$valid = 0;
    }
    
    # sort the exons 
    $transcript->sort;
    my @exons = @{$transcript->get_all_Exons};
    
    if ($#exons > 0) {
	for (my $i = 1; $i <= $#exons; $i++) {
	    
	    # check phase consistency:
	    if ( $exons[$i-1]->end_phase != $exons[$i]->phase  ){
		print STDERR "transcript $id has phase inconsistency\n";
		$valid = 0;
		last;
	    }
	    
	    # check for folded transcripts
	    if ($exons[0]->strand == 1) {
		if ($exons[$i]->start < $exons[$i-1]->end) {
		    print STDERR "transcript $id folds back on itself\n";
		    $valid = 0;
		} 
	    } 
	    elsif ($exons[0]->strand == -1) {
		if ($exons[$i]->end > $exons[$i-1]->start) {
		    print STDERR "transcript $id folds back on itself\n";
		    $valid = 0;
		} 
	    }
	}
    }
    if ($valid == 0 ){
	$self->_print_Transcript($transcript);
    }
    return $valid;
}

############################################################

sub _check_Translation{
  my ($self,$transcript) = @_;
  
  my $id = $self->transcript_id( $transcript );
  
  my $valid = 1;
  
  # check that they have a translation
  my $translation = $transcript->translation;
  my $sequence;
  eval{
    $sequence = $transcript->translate;
  };
  unless ( $sequence ){
    print STDERR "transcript $id has no translation\n";
    return 0;
  }
  if ( $sequence ){
    my $peptide = $sequence->seq;
    if ( $peptide =~ /\*/ ){
      print STDERR "translation of transcript $id has STOP codons\n";
      $valid = 0;
    }
  }
  if ($valid == 0 ){
    $self->_print_Transcript($transcript);
  }
  return $valid;
}

############################################################

sub transcript_id {
  my ( $self, $t ) = @_;
  my $id;
  if ( $t->stable_id ){
    $id = $t->stable_id;
  }
  elsif( $t->dbID ){
    $id = $t->dbID;
  }
  elsif( $t->temporary_id ){
    $id = $t->temporary_id;
  }
  else{
    $id = 'no-id';
  }
  
  if ($t->type){
      $id .= " ".$t->type."\n";
  }
  return $id;
}

############################################################
#
# METHODS DOING THE PRINTING
#
############################################################

sub _print_Transcript{
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my $id;
  if ($transcript->stable_id){
    $id = $transcript->stable_id;
  }
  elsif ( $transcript->dbID ){
    $id = $transcript->dbID;
  }
  else{
    $id = "no id";
  }
  if ( defined( $transcript->type ) ){
    $id .= " ".$transcript->type;
  }
  print STDERR "transcript: ".$id."\n";
  foreach my $exon ( @exons){
    $exon->gffstring;
    #$self->print_Exon($exon);
  }
  if ($transcript->translation){
    $self->_print_Translation($transcript->translation);
  }
}

############################################################

sub _print_Translation{
  my ($self,$translation) = @_;
  
  print STDERR "translation start exon: ".
    $translation->start_Exon->start."-".$translation->start_Exon->end.
      " start: ".$translation->start."\t phase: ".$translation->start_Exon->phase.
	" end_phase: ".$translation->start_Exon->end_phase."\n";
  
  print STDERR "translation end exon: ".
    $translation->end_Exon->start."-".$translation->end_Exon->end.
	" end: ".$translation->end."\t phase: ".$translation->end_Exon->phase.
	  " end_phase: ".$translation->end_Exon->end_phase."\n";
}

############################################################


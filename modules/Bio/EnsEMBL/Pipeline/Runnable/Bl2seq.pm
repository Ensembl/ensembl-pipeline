#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Bl2seq

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new(-seq1 => $seq1,
							    -seq2 => $seq2,
							    -min_eval => $min_eval,
							    -min_score => $min_score);
    or
    
    my $obj = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new()

=head1 DESCRIPTION

This runs bl2seq (executable from the ncbi) and provides feature pair output

-seq1 and -seq2 referred to Bio::PrimarySeq objects

=head2 Methods:

 new,
 seq1,
 seq2,
 min_eval,
 min_score,
 run,
 output.

=head1 CONTACT

ensembl-dev@ebi.ec.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::Runnable::Bl2seq;

use strict;
use vars qw(@ISA);
use Carp;
use Bio::Tools::BPlite;

# Object preamble - inherits from Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;
  
  my ($seq1, $seq2, $min_score, $min_eval) = 
    $self->_rearrange([qw(SEQ1 SEQ2 MIN_SCORE MIN_EVAL)], @args);

  if (! defined $seq1 || ! defined $seq2) {
    $self->throw("Must pass in both seq1 and seq1 args");
  }

  $self->seq1($seq1);
  $self->seq2($seq2);

  if (defined $min_score) {
    $self->min_score($min_score);
  } else {
    $self->min_score(40); 
  }
  if (defined $min_eval) {
    $self->min_eval($min_eval);
  } else {
    $self->min_score(0.01); 
  }
  return $self;
}

=head2 seq1

 Title   : seq1
 Usage   : $obj->seq1($newval)
 Function: 
 Example : 
 Returns : value of seq1
 Args    : newvalue (optional)


=cut

sub seq1 {
  my ($obj,$value) = @_;
  if( defined $value) {
    if (! ref $value || ! $value->isa('Bio::PrimarySeqI')) {
      $obj->throw("$value is not a PrimarySeqI object. Cannot throw");
    }
    my $selfseq = Bio::PrimarySeq->new(-display_id => $value->id,
				       -seq => $value->seq);
    if ($selfseq->length == 0) {
      $obj->throw("attempting to crossmatch seemingly 0 length sequence!");
    }
    $obj->{'seq1'} = $selfseq;
  }
  return $obj->{'seq1'};
}

=head2 seq2

 Title   : seq2
 Usage   : $obj->seq2($newval)
 Function: 
 Example : 
 Returns : value of seq2
 Args    : newvalue (optional)


=cut

sub seq2 {
   my ($obj,$value) = @_;
   if( defined $value) {
      if( !ref $value || !$value->isa('Bio::PrimarySeqI') ) {
	  $obj->throw("$value is not a PrimarySeqI object. Cannot throw");
      }
      my $selfseq = Bio::PrimarySeq->new( -display_id => $value->id , -seq => $value->seq);
      if( $selfseq->length == 0 ) {
	  $obj->throw("attempting to crossmatch seemingly 0 length sequence!");
      }
      $obj->{'seq2'} = $selfseq;

    }
    return $obj->{'seq2'};

}

=head2 min_score

 Title   : min_score
 Usage   : $obj->min_score($newval)
 Function: Get/Set the min_score value
 Example : 
 Returns : value of min_score used to filter bl2seq results
 Args    : newvalue (optional)


=cut

sub min_score {
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'min_score'} = $value;
  }
  return $self->{'min_score'};
}

=head2 min_eval

 Title   : min_eval
 Usage   : $obj->min_eval($newval)
 Function: Get/Set the min_eval value
 Example : 
 Returns : value of -e option used by bl2seq
 Args    : newvalue (optional)


=cut

sub min_eval {
  my ($self,$value) = @_;
  if (defined $value) {
    $self->{'min_eval'} = $value;
  }
  return $self->{'min_eval'};
}

=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run {
  my ($self,@args) = @_;

  # dump sequences to work directory
  my $query = "/tmp/query.".$$;
  my $sbjct = "/tmp/sbjct.".$$;
  
  open(F,">$query") || $self->throw("Cannot make $query $!");
  my $seqout = Bio::SeqIO->new(-fh => \*F,
			       '-format' => 'fasta');
  $seqout->write_seq($self->seq1);
  close(F);
  
  open(F,">$sbjct") || $self->throw("Cannot make $sbjct $!");
  $seqout = Bio::SeqIO->new(-fh => \*F,
			    '-format' => 'fasta');
  $seqout->write_seq($self->seq2);
  close(F);
  
  my $min_eval = $self->min_eval;
  my $min_score = $self->min_score;
  
  my $command = "/usr/local/ensembl/bin/bl2seq -i $query -j $sbjct -p blastn -e $min_eval 2>/dev/null |";
  print STDERR "opening bl2seq pipe\n";
  open (BL2SEQ, $command) || 
    croak "Can't open pipe '$command' : $!";
  
  my $bl2seq_parsing = new Bio::Tools::BPlite (-fh => \*BL2SEQ);

  my @features;

  while (my $sbjct = $bl2seq_parsing->nextSbjct) {
    last unless ($sbjct);
    while (my $hsp = $sbjct->nextHSP) {
      my ($score,$start,$end,$strand,$hstart,$hend,$hstrand) = ($hsp->bits, $hsp->query->start, $hsp->query->end,$hsp->query->strand, $hsp->subject->start,$hsp->subject->end,$hsp->subject->strand);
      
      next if ($score < $min_score);
      
      my $fp = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
      #print STDERR "Processing FP with $start-$end to $hstart-$hend\n";
      
      $fp->start($start);
      $fp->end($end);
      $fp->strand($strand);
      $fp->seqname($self->seq1->id);
      $fp->hstart($hstart);
      $fp->hend($hend);
      $fp->hstrand($hstrand);
      $fp->hseqname($self->seq2->id);
      $fp->score($score);
      
      push @features, $fp;

      $self->_add_fp($fp);
    }
    last; # should be only one subject anyway,
          # and skip the warning message thrown by Bio::Tools::BPlite object
          # Bio::AlignIO seems not work well dealing with bl2seq parsing!!!
  }
  
  unlink($query);
  unlink($sbjct);
}

=head2 _add_fp

 Title   : _add_fp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _add_fp{
   my ($self,@args) = @_;
   push(@{$self->{'_fp_array'}},@args);
}

sub output{
  my ($self,@args) = @_;
  return @{$self->{'_fp_array'}};
}

1;

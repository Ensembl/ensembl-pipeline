#
# Written by Eduardo Eyras
#
# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::AlignWise

=head1 SYNOPSIS

 
=head1 DESCRIPTION


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::AlignWise;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($slice, $seqs, $ids, $options ) = 
	$self->_rearrange([qw(
			      SLICE
			      SEQS
			      IDS
			      OPTIONS
			      )
			   ], @args);
    
    
    $self->{_output} = [];
    # must have a target and a query sequences
    unless( $seqs ){
	$self->throw("AlignWise: need sequences: $seqs");
    }
    $self->seqs(@{$seqs});
    
    # you can pass a sequence object for the target or a database (multiple fasta file);
    if( $ids ){
	$self->ids( @$ids );
    }
    else{
	$self->throw("AlignWise: need ids: $ids");
    }
    if ( $slice ){
	$self->slice( $slice );
    }
    else{
	$self->throw("need a slice");
    }

    my $basic_options = ' -genes ';
    if ( $options ){
	$basic_options .= $options;
    }
    $self->options( $basic_options);
    
    return $self;
}

############################################################
#
# Analysis methods
#
############################################################

=head2 run

Usage   :   $obj->run($workdir, $args)
Function:   Runs exonerate script and puts the results into the file $self->results
            It calls $self->parse_results, and results are stored in $self->output
=cut
  
sub run {
  my ($self) = @_;

  my $verbose = 1;
    
  # set the working directory (usually /tmp)
  $self->workdir('/tmp') unless ($self->workdir());
  my $dir         = $self->workdir();
  
  # results go here:
  $self->results($self->workdir()."/alignwise_results.$$");
  
  my $multi = "$dir/multi_seqs.$$";
  open( MULTI,">$multi") || $self->throw("Could not open $multi $!");
  my @ids  = $self->ids;
  my @seqs = $self->seqs;
  my $focus_id  = shift @ids;
  my $focus_seq = shift @seqs;
  
  print MULTI ">".$focus_id."\n";
  for( my $i=0; $i<length( $focus_seq ); $i+=60 ){
      print MULTI substr( $focus_seq, $i, 60 ),"\n";
  }
  for( my $i = 0; $i<=$#seqs; $i++ ) {
      print MULTI ">".$ids[$i]."\n";
      my $seq = $seqs[$i];
      for( my $j=0; $j<length( $seq ); $j+=60 ){
	  print MULTI substr( $seq, $j, 60 ),"\n";
      }
  }
  close(MULTI);
  
  my $command ="/nfs/acari/birney/bin/alignwise $multi ".$self->options." | ";
  
  print STDERR "running alignwise: $command\n" if $verbose;

  open( ALIGN, $command ) || $self->throw("Error running alignwise $!");
  
  
  ############################################################
  # store each prediction as a transcript with supporting features
  my @transcripts;

  ############################################################
  # parse results - avoid writing to disk the output
  my $in_gene = 0;
  my $transcript = Bio::EnsEMBL::Transcript->new();
  while (<ALIGN>){
    
      #next if $_=~/calculation/;
      #print STDERR $_ if $verbose;

      # gene output is of the form:
      #
      #Gene 1
      #Gene 3502 13989 
      #  Exon 3502 3570 phase -1
      #  Exon 4121 4162 phase -1
      #  Exon 4834 4840 phase -1
      #  Exon 10892 10919 phase -1
      #  Exon 11003 11041 phase -1
      #  Exon 13956 13989 phase -1
      #Gene 2
      #Gene 21246 19480
      #  Exon 21246 21226 phase -1
      #  Exon 19593 19557 phase -1
      #  Exon 19523 19480 phase -1
      #Gene 3
      chomp;
      my @entries = split;
      if ( ( $entries[0] eq 'Gene' && $entries[1] && !$entries[2] ) 
	   ||
	   ( $entries[0] eq '//' )
	   ){
	  if ( $in_gene ){
	      push( @transcripts, $transcript );
	  }
	  $transcript = Bio::EnsEMBL::Transcript->new();
      }
      if ( $entries[0] eq 'Gene' && $entries[2] ){
	  $in_gene = 1;
      }
      if ( $entries[0] eq 'Exon' ){
	  my $exon = Bio::EnsEMBL::Exon->new();
	  my $strand = 1;
	  my ($start, $end) = ( $entries[1], $entries[2] );
	  if ( $start > $end ){
	      $strand = -1;
	      $end    = $entries[1];
	      $start  = $entries[2];
	  }
	  $exon->start( $start );
	  $exon->end( $end );
	  $exon->strand( $strand );
	  $exon->contig( $self->slice );
	  $transcript->add_Exon( $exon );
      }
  
  } # end of while loop
  
  $self->output( @transcripts );
  close(ALIGN);  
    
  ############################################################
  # remove temp file 
  # unlink $multi;

  open( OUT, ">/ecs2/work1/eae/AlignWise/out.gff" ) or die("cannot open file /ecs2/work1/eae/AlignWise/out.gff");

  print STDERR scalar(@transcripts)." transcripts produced:\n";
  my $trans_count = 0;
  my $exon_count = 0;
  foreach my $trans ( @transcripts ){
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript( $trans );
    $trans_count++;
    my $trans_id = $trans_count;
    foreach my $exon ( @{$trans->get_all_Exons} ){
      $exon_count++;
      my $strand = "+";
      if ($exon->strand == -1) {
	$strand = "-";
      }
      my $phase = ".";
      my $score = 100;
      my $g_type = 'alignwise';
      print OUT $exon_count."\t".
	$g_type."\t".
	  "exon\t". 
	    ($exon->start) . "\t" . 
	      ($exon->end) . "\t" .
		$score. "\t" . 
		  $strand . "\t" . 
		    $phase . "\t" . 
		      $trans_id. "\n";
      
    }
  }
  close( OUT );

}

############################################################
#
# get/set methods
#
############################################################

sub seqs {
    my ($self, @seqs) = @_;
    if (@seqs){
	push(@{$self->{_query_seqs}}, @seqs) ;
    }
    return @{$self->{_query_seqs}};
}

############################################################

sub options {
  my ($self, $options) = @_;
  if ($options) {
    $self->{_options} = $options ;
  }
  return $self->{_options};
}

############################################################

sub slice {
  my ($self, $slice) = @_;
  if ($slice) {
    $self->{_slice} = $slice;
  }
  return $self->{_slice};
}

############################################################

sub output {
    my ($self, @output) = @_;
  if (@output) {
      unless( $self->{_output} ){
	  $self->{_output} = [];
      }
      push( @{$self->{_output}}, @output );
  }
  return @{$self->{_output}};
}

############################################################

sub ids {
    my ($self, @ids) = @_;
    if (@ids){
	push(@{$self->{_ids}}, @ids) ;
    }
    return @{$self->{_ids}};
}

############################################################


1;


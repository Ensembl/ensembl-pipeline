#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Tools::TranslationUtils - 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::Tools::TranslationUtils;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;
use Bio::EnsEMBL::Utils::PolyA;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;


@ISA = qw(Bio::EnsEMBL::Root);

 



############################################################

=head2

method to calculate the translation for a given transcript
in slice coordinates. The translation is attached to the transcript
as a translation object.

The method calculates the longest translation starting with Methionine.
If there is none, it will take the longest translation. It is simple
but it seems to give similar results to genomewise (regarding translation length).

=cut

sub compute_translation{
  my ($self,$trans) = @_;
  
  my $verbose = 0;

  my @met_predictions   = $self->run_translate( $trans,1);
  my @nomet_predictions = $self->run_translate( $trans );
  
  my $count = 0;
  while ( $count < 2 && $met_predictions[$count] ){
    my @entry = @{$met_predictions[$count]};
    #print STDERR "MET length:$entry[0] start:$entry[1] end:$entry[2]\n";
    $count++;
  }
  $count = 0;
  while ( $count < 2 && $nomet_predictions[$count] ){
    my @entry = @{$nomet_predictions[$count]};
    #print STDERR "NO_MET length:$entry[0] start:$entry[1] end:$entry[2]\n";
	$count++;
    }
    my $translation = $trans->translation;
    my $length = $trans->seq->length;
    my $best;
    if ( @met_predictions && @nomet_predictions ){
	my $met_best   = $met_predictions[0];
	my $nomet_best = $nomet_predictions[0];
	if ( $nomet_best->[0] > 2*$met_best->[0] ){
	    $best = $nomet_best;
	}
	else{
	    $best = $met_best;
	}
    }
    elsif( @met_predictions ){
	$best = $met_predictions[0];
    }
    elsif( @nomet_predictions ){
	$best = $nomet_predictions[0];
    }
    my @entry = @{$best};
    my $orf_start = $entry[1];
  my $orf_end   = $entry[2];
  print STDERR "BEST length:$entry[0] start:$entry[1] end:$entry[2]\n" if $verbose;
    my @exons;
    my $strand = $trans->start_Exon->strand;
    if ( $strand == 1 ){
	@exons = sort { $a->start <=> $b->start } @{$trans->get_all_Exons};
    }
    else{
	@exons = sort { $b->start <=> $a->start } @{$trans->get_all_Exons};
    }
    my $transl_start;
    my $transl_end;
    my $transl_start_Exon;
    my $transl_end_Exon;
    my $exon_count = 0;
    print STDERR "transcript length: $length\n" if $verbose;
    my $pos = 1;
    foreach my $exon ( @exons ){
	$exon_count++;
	print STDERR "exon:$exon_count exon_length:".$exon->length." pos:$pos orf_start:$orf_start orf_end:$orf_end pos+:".($pos + $exon->length - 1)."\n" if $verbose;
	if ( $orf_start >= $pos && $orf_start <= $pos + $exon->length - 1 ){
	    $transl_start_Exon = $exon;
	    $transl_start      = $orf_start - $pos + 1;
	    print STDERR "start found\n" if $verbose;
	    $exon->phase(0);
	}
	if ( $orf_end >= $pos && $orf_end <= $pos + $exon->length - 1 ){
	    $transl_end_Exon   = $exon;
	    $transl_end        = $orf_end - $pos + 1;
	    print STDERR "end found\n" if $verbose;
	}
	$pos += $exon->length;
	#last if ( $pos > $length );
    }
  #print STDERR "original translation:\n";
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($trans);

    my $newtranslation;
    if ( $transl_start && $transl_end &&  $transl_start_Exon && $transl_end_Exon ){
	$newtranslation = Bio::EnsEMBL::Translation->new();
	$newtranslation->start( $transl_start );
	$newtranslation->end( $transl_end );
	$newtranslation->start_Exon( $transl_start_Exon );
	$newtranslation->end_Exon( $transl_end_Exon );
	$trans->translation($newtranslation);
	#print STDERR "modified translation:\n";
	#Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($trans);
      }
  else{
      print STDERR "problem making the translation\n";
    }
  $trans->flush_Exons;
  foreach my $exon ( @exons ){
    $trans->add_Exon($exon);
  }
  return $trans;
}

############################################################

sub run_translate{
    my ($self,$trans,$met) = @_;
    
    my $verbose = 0;

    my $trans_id = $trans->stable_id || $trans->dbID;
    unless ( $trans_id ){
      if ( $trans->type ){
	$trans_id = $trans->type;
      }
      else{
	$trans_id = "transcript_".$$;
      }
    }
    my $seq = $trans->seq;
    unless ( $seq->display_id ){
	$seq->display_id( $trans_id );
    }
    my $length = $seq->length;
    
    print STDERR "display_id = ".$seq->display_id."\n";
    ############################################################
    # create file
    my $file = "/tmp/"."cdna_".$$.".fa";
    open ( SEQ, ">$file" ) || die("could not open file $!");
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
				 '-fh'     => \*SEQ);
    $seqout->write_seq($seq);
    close(SEQ);
    
    my $command;
    if ( $met){
	$command = "/usr/local/ensembl/bin/translate -m $file |";
    }
    else{
	$command = "/usr/local/ensembl/bin/translate $file |";
    } 
    open ( ORF, $command ) || die( "Error running translate" );
    ############################################################
    # out put is of the form:
    #> gi|20070124|ref|NM_000918.2|.44    length 62, nt 2236..2051
    #AHDRRRSPGLREGEGPGLCRAPGLAATSSSSRHGGHPDRIRKSPFTQKCKSHDQSWRHCRRY
    #> gi|20070124|ref|NM_000918.2|.45    length 34, nt 2047..1946
    #VTMSSPAPSLPHGGQASPRRPGQGGTNTLMSKNV
    #> gi|20070124|ref|NM_000918.2|.46    length 34, nt 1942..1841
    #KSHRRNFQKEEKPPAGGRQRDSEHGSKHSGQTHV
    
    my @orf_predictions;
  ORF:
    while ( <ORF> ){
	chomp;
	next ORF unless /\>/;
	print STDERR "$_\n" if $verbose;
	my @entries = split;
	next ORF unless ( $entries[3] && $entries[5] );
	my $id = $entries[1];
	my $orf_length = $entries[3];
	$orf_length =~s/\,//;
	$entries[5] =~/(\d+)\.\.(\d+)/;
	my $orf_start = $1;
	my $orf_end   = $2;
	next ORF if $orf_start>=$orf_end;
	print STDERR "id:$id\torf_length:$orf_length\tstart:$orf_start\tend:$orf_end\n" if $verbose;
	my @prediction = ($orf_length,$orf_start,$orf_end);
	push( @orf_predictions, \@prediction );
    }
    my @sorted_predictions = 
	map { $_->[1] } sort { $b->[0] <=> $a->[0] } map { [$_->[0], $_] } @orf_predictions;
    return @sorted_predictions;
}


############################################################

 1;

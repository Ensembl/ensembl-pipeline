#
# Ensembl module for Bio::EnsEMBL::Pipeline::Runnable::Genomewise.pm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Genomewise.pm - Runs Genomewise and makes a set of Transcripts

=head1 SYNOPSIS

my $genomewise = Bio::EnsEMBL::Pipeline::Runnable::Genomewise->new( -seq              => $primary_seq,
                                                                    -switch           => $switch,
                                                                    -smell            => $smell,
                                                                    -cds              => $cds_transcript
								    -skip_small_exons => $min_exon_size,
                                                                  );

where 
-seq    is the sequence where we are going to run genomewise on
-switch is the cost to switch evidence. It is irrelevant if one only passes 1 transcript to Genomewise.
        With more than one transcript, it is the cost of jumping from one to the other to find the CDS.
-smell  is the space allowed to move off the splice site of the evidence. The minimum is 0.
-cds    if this option is included, it will feed a cds evidence to genomewise. This is passed in
        in the form of a transcript with exons and each exon must have phases. This option
        has not been tested in genomewise itself yet.
-skip_small_exons   
        minimum exon size which is accepted to be run in genomewise, i.e. smaller exons are not
        passed as evidence for genomewise. This could be useful (maybe with cross-species alignments)
        as genomewise sometimes shifts around very small exons.

add evidence in the form of a transcript:
   $genomewise->add_Transcript($transcript);

   $genomewise->run;
   my @transcripts = $genomewise->output;


=head1 DESCRIPTION

This class is a wrapper to around the genomewise binary.
You pass in a piece of sequence and evidence in the form of exons in a transcript. It is also possible to pass in
cds and framshifts (see genomewise -help for details). The evidence must be in the coordinate system
of the PrimarySeq passed in. 

Genomewise then will try to find the most likely ORF across the sequence and it will produce
a transcript structure with CDS and UTRs. It can even accept multiple transcripts. It will then choose
the best across the exons. It aims for longest ORF and it tests for start/stop codons and
splice site consensus sequences.

At the moment you can pass in the options -switch and -smell.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=cut

# Let the code begin...


package Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    my( $seq, $switch, $smell, $cds, $skip_small_exons ) = $self->_rearrange([qw(SEQ
										 SWITCH
										 SMELL
										 CDS
										 SKIP_SMALL_EXONS
										 )],
									     @args);
    
    unless ( $seq ){
	$self->throw->("must provide a PrimarySeq with option -seq");
    }
    $self->seq($seq);
    
    # read 'skip' argument
    if ( $skip_small_exons ){
	$self->_skip_small_exons( $skip_small_exons );
    }
    
    # read smell argument (can be 0 ) or default it to 8
    if (defined($smell) ){
	$self->smell($smell);
    }
    else{
	$self->smell(8);
    }
    # read the switch argument or default it to 10000
    if (defined($switch)){
      $self->switch($switch);
    }
    else{
      $self->switch(10000);
    }
    # any cds?
    if ( $cds ){
      $self->cds($cds);
    }


    $self->{'_transcript_array'} = [];
    $self->{'_output_array'}     = [];
    return $self;
}

=head2 run
 
      Arg: no argument
Function : runs genomewise using as evidence the transcript stored in $self->each_Transcript and as sequence
  the one stored in $self->seq. These are populated when the Runnable is created (see for instance EST_GeneBuilder
  where MiniGenomewise runnables are created, or in Minigenomewise where Genomewise runnables are created in a 
  miniseq.)

Returntype: it does not return anything, but instead it populates the get/set container $self->output with
  a transcript with translation and UTRs (and exons, of course).
  Exception : throws if in the reverse strand, genomewise cannot handle this. Also throws if the output
  of does not make sense or if it calls an empty container.

=cut

sub run{
  my ($self) = @_;

#  $self->analysis($self->_make_Analysis) unless $self->analysis;

  if( !defined $self->workdir ) {
    $self->workdir('/tmp');
  }

  my $dir = $self->workdir();
  
  # file with the genomic sequence
  my $genome_file = "$dir/gowise_seq.$$";
  
  # file with the evidence
  my $evi_file    = "$dir/gowise_evi.$$";
  
  open(F,">$genome_file") || $self->throw("Could not open $genome_file $!");
  open(E,">$evi_file")    || $self->throw("Could not open $evi_file $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'fasta','-fh' => \*F);
  $seqout->write_seq($self->seq);
  
  $seqout = undef;
  close(F);
  
  #### collect info from transcripts and put it in evidence file
  my %supp_evidence;
  
  # can pass multiple transcripts to genomewise
  foreach my $transcript ( @{$self->get_all_Transcripts}) {
      
      #print STDERR "In Genomewise, before going through genomewise\n";
      my $excount = 0;
      
    EXON:
      foreach my $exon ( @{$transcript->get_all_Exons} ) {
	  $excount++;
	  
	  # store supporting evidence
	  my @evi = @{$exon->get_all_supporting_features};
	  $supp_evidence{ $exon } = \@evi;
	  
	  if( $exon->strand == -1 ) {
	      $self->throw("Genomewise cannot handle reverse strand exons. Must flip outside");
	  }
	  # skip small exons if we have chosen to do so
	  if ( $self->_skip_small_exons ){
	      if ( ($exon->end - $exon->start + 1 ) < $self->_skip_small_exons ){
		  next EXON;
	      }
	  }
	  print E "exon ",$exon->start," ",$exon->end,"\n";
      }
      print E "//\n";
  }
  
  if ( $self->cds ){
      foreach my $exon ( @{$self->cds->get_all_Exons} ){
	  print E $exon->start." ".$exon->end." ".$exon->phase."\n";
      }
      print E "//\n";
  }
  close(E);
  
  my $switch = $self->switch;
  my $smell  = $self->smell;
  print STDERR "running genomewise with smell: $smell and switch: $switch\n";
  if ( $self->_skip_small_exons ){
      print STDERR "skipping exons smaller than ".$self->_skip_small_exons."\n";
  }

  #### Steve's version (fixed)
  #  open(GW,"/nfs/acari/searle/progs/ensembl-trunk/wise2/src/models/genomewise -switch 10000 -silent -nogff -smell 8 -notrans -nogenes -geneutr $genome_file $evi_file |");
  
  ### 16th August 2002, the version in the blades is the latest, the version in the slates is old.
  open(GW,"/usr/local/ensembl/bin/genomewise -switch $switch -silent -nogff -smell $smell -notrans -nogenes -geneutr $genome_file $evi_file |");
  
  # parse gff output for start, stop, strand, phase
  my $genename = '';
  
 GENE:
  while( <GW> ) {
    /\/\// && last;
    if( /Gene/ ) {
      
      my $t = Bio::EnsEMBL::Transcript->new();
      my $trans = Bio::EnsEMBL::Translation->new();
      $t->translation($trans);
      $self->output($t);
      
      my $seen_cds = 0;
      my $seen_utr3 = 0;
      my $prev = undef;
      while( <GW> ) {
	
	## this will print out the structure foudn by genomewise
	#print STDERR "$_";
	#
	chomp;
	#print "Seen $_\n";
	  if( /End/ ) {
	    if( $seen_utr3 == 0 ) {
	      #  print STDERR "Seen utr - setting end to prev\n";
	      #  ended on a cds 
	      $trans->end_Exon($prev);
	      $trans->end($prev->length);
	    }
	    next GENE;
	  }
	
	if( /utr5\s+(\d+)\s+(\d+)/) {
	  my $start = $1;
	  my $end   = $2;
	  my $strand = 1;
	  
	  my $exon = Bio::EnsEMBL::Exon->new();
	  $exon->start($start);
	  $exon->end  ($end);
	  $exon->strand($strand);
	  # there isn't really a phase in the UTR exon, but we set it to -1 otherwise later stages fail
	  $exon->phase(-1);
	  $exon->end_phase(-1);
	  $t->add_Exon($exon);
	  $prev = $exon;
	  next;
	}
	if( /cds\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/ ) {
	  my $start = $1;
	  my $end   = $2;
	  my $strand = 1;
	  my $phase = $3;
	  
	  my $this_start  = $1;
	  my $this_end    = $2;
	  my $this_length = $this_end - $this_start + 1;
	  
	  my $cds_start = $start;
	  my $exon;
	  
	  if( $seen_cds == 0 && defined $prev && $prev->end+1 == $start ) {
	    # we have an abutting utr5
	    $start = $prev->start();
	    $exon  = $prev;
	    # end_phase is number of bases at the end of the exon which do not fall in a codon
	    my $end_phase = $this_length%3;
	    $exon->start($start);
	    $exon->end  ($end);
	    $exon->strand($strand);
	    $exon->phase($prev->phase);
	    #print STDERR "Genomewise: setting end_phase to $end_phase\n";
	    $exon->end_phase($end_phase);
	  } 
	  else {
	    # make a new Exon
	    $exon  = Bio::EnsEMBL::Exon->new();
	    $t->add_Exon($exon);
	    $exon->phase($phase);
	    $exon->start($start);
	    $exon->end  ($end);
	    $exon->end_phase( ( $end - $start + 1 + $phase ) %3 );
	    $exon->strand($strand);
	  }
	  if( $seen_cds == 0 ) {
	    $trans->start_Exon($exon);
	    $trans->start($cds_start-$start+1);
	  }
	  $seen_cds = 1;
	  $prev = $exon;
	  next;
	}
	
	if( /utr3\s+(\d+)\s+(\d+)/ ) {
	  #	print "Found utr3\n";
	  my $start = $1;
	  my $end   = $2;
	  my $strand = 1;
	  
	  my $orig_end;
	  my $exon;
	  if( $seen_utr3 == 0 && defined $prev && $prev->end+1 == $start ) {
	    
	    #   abutting 3utr
	    #      _____________
	    #...--|______|______|--...
	    #   cds($prev) utr3
	    
	    $orig_end = $prev->end;
	    $exon     = $prev;
	    $start    = $prev->start;
	    $exon->start($start);
	    $exon->end  ($end);
	    $exon->strand($strand);
	    $exon->phase($prev->phase);
	    $exon->end_phase(-1);
	    
	  } 
	  else {
	    
	    #   not abutting; should be fine.
	    $exon =  Bio::EnsEMBL::Exon->new();
	    $t->add_Exon($exon);
	    $exon->start($start);
	    $exon->end  ($end);
	    $exon->strand($strand);	    
	    $exon->phase(-1);
	    $exon->end_phase(-1);
	  }
	    
	  # there isn't really a phase in the UTR exon, but we set it to -1 otherwise later stages fail
	  if( $seen_utr3 == 0 ) {
	    if( !defined $orig_end ) {
	      $self->throw("This should not be possible. Got a switch from cds to utr3 in a non-abutting exon");
	    } 
	    else {
	      $trans->end_Exon($exon);
	      
	      # position of the end of translation in exon-local coordinates
	      $trans->end($orig_end - $start + 1);
	    }
	  }
	  $seen_utr3 = 1;
	  next;
	}
	
	# else - worrying 
	
	chomp;
	$self->throw("Should not able to happen - unparsable line in geneutr $_");
      }
    }
    chomp;
    #print STDERR "genomic file: $genome_file, evidence file: $evi_file\n";
    $self->throw("Should not able to happen - unparsable in between gene line $_");
  }
  #print STDERR "genomic file: $genome_file, evidence file: $evi_file\n";
  

  ############################################################
  ##### need to put back the supporting evidence
  ############################################################

  # since exons have been created anew, need to check overlaps
  my @trans_out = $self->output;
  my @trans_in  = @{$self->get_all_Transcripts};
  
  # let's try to make it fast:
  if ( scalar( @trans_out) == 1 && scalar( @trans_in ) == 1 ){
      my ( $consecutive_overlap, $mismatches ) = 
	  $self->_check_consecutive_overlap( $trans_in[0], $trans_out[0] );
      
      if ( $consecutive_overlap == 1 ){
	  my @exons_in  = sort { $a->start <=> $b->start } @{$trans_in[0] ->get_all_Exons};
	  my @exons_out = sort { $a->start <=> $b->start } @{$trans_out[0]->get_all_Exons};
	  
	EXON_2_EXON:
	  while ( scalar ( @exons_in ) > 0 ){
	      my $exon_in  = shift( @exons_in  );
	      my $exon_out = shift( @exons_out );  
	      
	      # check just in case
	      if ( $exon_in->overlaps( $exon_out ) ){
		  foreach my $feature ( @{ $supp_evidence{ $exon_in } } ){
		      $exon_out->add_supporting_features( $feature );
		  }
	      }
	  } # end of EXON_2_EXON
      }
      else{
	  
	  my @exons_in  = sort { $a->start <=> $b->start } @{$trans_in[0] ->get_all_Exons};
	  my @exons_out = sort { $a->start <=> $b->start } @{$trans_out[0]->get_all_Exons};
	  
	  print STDERR "passing evi info between 2 transcripts with non-consecutive overlaping exons\n";
	  foreach my $exon_in ( @exons_in ){
	      foreach my $exon_out( @exons_out ){
		  if ( $exon_out->overlaps($exon_in) ){
		      foreach my $feature ( @{ $supp_evidence{ $exon_in } } ){
			  $exon_out->add_supporting_features( $feature );
		      }
		  }
	      }
	  }
	  # test:
	  print STDERR "before genomewise\n:";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence( $trans_in[0] );
	  print STDERR "after genomewise\n:";
	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence( $trans_out[0] );
      }

  }
  else{
      # if we have more than one transcript at one or both sides, we also have to check them all
      print STDERR "passing evi info between more than 2 transcripts\n";
      foreach my $tran_in ( @trans_in){
	  my @exons_in  = sort { $a->start <=> $b->start } @{$trans_in[0] ->get_all_Exons};
	  
	  foreach my $tran_out ( @trans_out ){
	      my @exons_out =  sort { $a->start <=> $b->start } @{$trans_out[0]->get_all_Exons};
	      
	      foreach my $exon_in ( @exons_in ){
		  
		  foreach my $exon_out( @exons_out ){
		      if ( $exon_out->overlaps($exon_in) ){
			  foreach my $feature ( @{ $supp_evidence{ $exon_in } } ){
			      $exon_out->add_supporting_features( $feature );
			  }
		      }
		  }
	      }
	  }
      }   
  }
  


  unlink $genome_file;
  unlink $evi_file;
}

############################################################

sub output{
  my ($self,@args) = @_;
  unless ( $self->{'_output_array'} ){
    $self->{'_output_array'} = [];
  }
  if ( @args ){
    push ( @{$self->{'_output_array'}}, @args );
  }
  return @{$self->{'_output_array'}};
}




=head2 get_all_Transcripts

 Title   : get_all_Transcripts
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Transcripts {
   my ($self) = @_;

   return $self->{'_transcript_array'};
}


=head2 add_Transcript

 Title   : add_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Transcript{
   my ($self,@args) = @_;

   foreach my $val ( @args ) {
     $self->throw("[$val] is not a Transcript...") unless $val->isa("Bio::EnsEMBL::Transcript");
       push(@{$self->{'_transcript_array'}},$val);
   }
   return $self->{_transcript_array};
}

############################################################

sub seq{
  my ($self,$seq) = @_;
  if( $seq ) {
    if( !ref $seq || !$seq->isa("Bio::PrimarySeqI") ) {
      $self->throw("No primary sequence provided...");
    }
    
    $self->{_seq} = $seq;
  }
  return $self->{_seq};

}

############################################################

sub switch{
  my ($self,$switch) = @_;
  if( defined($switch) ) {
    $self->{_switch} = $switch;
  }
  return $self->{_switch};
}

############################################################

sub cds{
  my ($self,$cds) = @_;
  if ( $cds ){
    $self->{_cds} = $cds;
  }
  return $self->{_cds};
}

############################################################

sub smell{
  my ($self,$smell) = @_;
  if( defined($smell) ) {
    $self->{_smell} = $smell;
  }
  return $self->{_smell};
}

############################################################

sub _skip_small_exons{
    my ($self,$skip_small_exons) =@_;
    if ( $skip_small_exons ){
	$self->{_skip_small_exons} = $skip_small_exons;
    }
    return $self->{_skip_small_exons};
}

############################################################

sub _transfer_transcript_supporting_evidence{
    my ($self,$trans_source,$trans_target) = @_;
    
    my @exons_source  = @{$trans_source->get_all_Exons};
    my @exons_target  = @{$trans_target->get_all_Exons};
    
    # most of the time the numbers of exons doesn't vary
    if ( scalar( @exons_source ) == scalar ( @exons_target ) ){
	#print STDERR "passing evi info between 2 transcripts with same number of exons\n";
	while ( scalar ( @exons_source ) > 0 ){
	    my $exon_in  = shift( @exons_source  );
	    my $exon_out = shift( @exons_target );  
	    
	    # check just in case
	    if ( $exon_in->overlaps( $exon_out ) ){
	      Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence( $exon_in,$exon_out);
	    }
	    else{
		$self->warn("Trying to pass evidence between exons that do not overlap, this won't work!");
	    }
	}
    }
    else{
	# if not the same number of exons, we cannot know how the split happened
	print STDERR "passing evi info between 2 transcripts with different number of exons\n";
	foreach my $exon_in ( @exons_source ){
	    foreach my $exon_out( @exons_target ){
		if ( $exon_out->overlaps($exon_in) ){
		  Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence( $exon_in,$exon_out);
		}
	    }
	}
    }
}

############################################################

sub _check_consecutive_overlap{
    my ( $self, $t1, $t2 ) = @_;
    my @e1 = sort { $a->start <=> $b->start } @{$t1->get_all_Exons};
    my @e2 = sort { $a->start <=> $b->start } @{$t2->get_all_Exons};
    my $consecutive_overlap = 1;
    my $mismatches = 0;
    
  EXON_2_EXON:
    while ( scalar ( @e1 ) > 0 ){
	my $exon_in  = shift( @e1 );
	my $exon_out = shift( @e2 );  
       	unless ( $exon_in->overlaps( $exon_out ) ){
	    $consecutive_overlap = 0;
	    $mismatches++;
	}
    }
    return ( $consecutive_overlap, $mismatches );
}

############################################################



1;

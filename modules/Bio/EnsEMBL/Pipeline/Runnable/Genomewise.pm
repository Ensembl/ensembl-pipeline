
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

   $run = Bio::EnsEMBL::Pipeline::Runnable::Genomewise->new();

   $trans = Bio::EnsEMBL::Transcript->new();
   $exon  = Bio::EnsEMBL::Exon->new();
   $exon->start(123);
   $exon->end(345);
 
   $trans->add_Exon($exon);
   # add more Exons, and more Transcripts as desired

   $run->add_Transcript($trans);
   $run->seq($genomic_seq);

   $run->run;


=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Transcript;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);

sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
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
  foreach my $transcript ( @{$self->get_all_Transcripts}) {
    
    #print STDERR "In Genomewise, before going through genomewise\n";
    my $excount = 0;
    foreach my $exon ( @{$transcript->get_all_Exons} ) {
      $excount++;

      # store supporting evidence
      my @evi = @{$exon->get_all_supporting_features};
      print STDERR "Exon $excount: ".scalar(@evi)." features\n";

      $supp_evidence{ $exon } = \@evi;

      ## test
      if (@evi){
	foreach my $evi ( @evi ){
	  print STDERR $evi." - ".$evi->gffstring."\n";
	}
      }
      else{
	;#print STDERR "No supporting evidence\n";
      }
      
      if( $exon->strand == -1 ) {
	$self->throw("Genomewise cannot handle reverse strand exons. Must flip outside");
      }
      print E "exon ",$exon->start," ",$exon->end,"\n";
    }
    print E "//\n";
  }
  close(E);
  

#
#   open(GW,"genomewise $genome_file $evi_file |");
#  system "/nfs/acari/birney/prog/wise2/src/models/genomewise -silent -nogff -notrans -nogenes -geneutr $genome_file $evi_file > /tmp/test.out";
#  open(GW,"</tmp/test.out");
  

# in acari, genomewise is in '/nfs/acari/birney/prog/wise2/src/bin/'
# or in '/usr/local/ensembl/bin/genomewise'

#### Ewan's version
#
#  open(GW,"/nfs/acari/birney/prog/wise2/src/bin/genomewise -silent -nogff -smell 4 -notrans -nogenes -geneutr $genome_file $evi_file |");

#### Steve's version (fixed)

  open(GW,"/nfs/acari/searle/progs/ensembl-trunk/wise2/src/models/genomewise -switch 10000 -silent -nogff -smell 8 -notrans -nogenes -geneutr $genome_file $evi_file |");

  
  # parse gff output for start, stop, strand, phase
  my $genename = '';

 GENE:
  while( <GW> ) {
    #print STDERR $_;

      /\/\// && last;
      if( /Gene/ ) {
	
	my $t = Bio::EnsEMBL::Transcript->new();
	my $trans = Bio::EnsEMBL::Translation->new();
	$t->translation($trans);
	$self->output($t);
	#push(@{$self->{'_output_array'}},$t);
	
	my $seen_cds = 0;
	my $seen_utr3 = 0;
	my $prev = undef;
	while( <GW> ) {
	  print STDERR "$_";
	  chomp;
	  #print "Seen $_\n";
	  if( /End/ ) {
	    if( $seen_utr3 == 0 ) {
	      #  print STDERR "Seen utr - setting end to prev\n";
	      # ended on a cds 
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
	    # there isn't really a phase in teh UTR exon, but we set it to -1 otherwise later stages fail
	    $exon->phase(-1); 
	    $t->add_Exon($exon);
	    $prev = $exon;
	    next;
	  }
	  if( /cds\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/ ) {
	    my $start = $1;
	    my $end   = $2;
	    my $strand = 1;
	    my $phase = $3;
	    
	    my $cds_start = $start;
	    my $exon;
	    if( $seen_cds == 0 && defined $prev && $prev->end+1 == $start ) {
	      # we have an abutting utr5
	      $start = $prev->start();
	      $exon  = $prev;
	    } else {
	      # make a new Exon
	      $exon  = Bio::EnsEMBL::Exon->new();
	      $t->add_Exon($exon);
	    }
	    $exon->start($start);
	    $exon->end  ($end);
	    $exon->strand($strand);
	    $exon->phase($phase);
	    
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
	      
	    } else {
	      
	      #   not abutting; should be fine.
	      $exon =  Bio::EnsEMBL::Exon->new();
	      $t->add_Exon($exon);
	      $exon->phase(-1); 
	    }
	    
	    # there isn't really a phase in the UTR exon, but we set it to -1 otherwise later stages fail
	    
	    $exon->start($start);
	    $exon->end  ($end);
	    $exon->strand($strand);
	    
	    if( $seen_utr3 == 0 ) {
	      if( !defined $orig_end ) {
		$self->throw("This should not be possible. Got a switch from cds to utr3 in a non-abutting exon");
	      } else {
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
  

  ##### need to put back the supporting evidence
  # since exons have been created anew, need to check overlaps
  my @trans_out = $self->output;
  my @trans_in  = @{$self->get_all_Transcripts};
  
  # let's try to make it fast:
  if ( scalar( @trans_out) == 1 && scalar( @trans_in ) == 1 ){
    my @exons_in  = @{$trans_in[0] ->get_all_Exons};
    my @exons_out = @{$trans_out[0]->get_all_Exons};
    
    # most of the time the numbers of exons doesn't vary
    if ( scalar( @exons_in ) == scalar ( @exons_out ) ){
      #print STDERR "passing evi info between 2 transcripts with same number of exons\n";
      while ( scalar ( @exons_in ) > 0 ){
	my $exon_in  = shift( @exons_in  );
	my $exon_out = shift( @exons_out );  
	
	# check just in case
	if ( $exon_in->overlaps( $exon_out ) ){
	  foreach my $feature ( @{ $supp_evidence{ $exon_in } } ){
	    $exon_out->add_Supporting_Feature( $feature );
	  }
	}
	else{
	  $self->throw("Try to pass evidence between exons that do not overlap, this is not good!");
	}
      }
    }
    else{
      # if not the same number of exons, we cannot know how the split happened
      print STDERR "passing evi info between 2 transcripts with different number of exons\n";
      foreach my $exon_in ( @exons_in ){
	foreach my $exon_out( @exons_out ){
	  if ( $exon_out->overlaps($exon_in) ){
	    foreach my $feature ( @{ $supp_evidence{ $exon_in } } ){
	      $exon_out->add_Supporting_Feature( $feature );
	    }
	  }
	}
      }
    }

  }
  else{
    # if we have more than one transcript at one or both sides, we also have to check them all
    print STDERR "passing evi info between more than 2 transcripts\n";
    foreach my $tran_in ( @trans_in){
      my @exons_in  = @{$trans_in[0] ->get_all_Exons};
      
      foreach my $tran_out ( @trans_out ){
	my @exons_out = @{$trans_out[0]->get_all_Exons};
	
	foreach my $exon_in ( @exons_in ){
	  
	  foreach my $exon_out( @exons_out ){
	    if ( $exon_out->overlaps($exon_in) ){
	      foreach my $feature ( @{ $supp_evidence{ $exon_in } } ){
		$exon_out->add_Supporting_Feature( $feature );
	      }
	    }
	  }
	}
      }
    }   
  }
    
  ## now we check that all the evi info has been correctly passed
  foreach my $tran_out ( @trans_out ){
    #print STDERR "In Genomewise, AFTER going through genomewise\n";
    my $count = 0;
    foreach my $exon_out ( @{$tran_out->get_all_Exons} ) {
      $count++;
      my @evi = @{$exon_out->get_all_supporting_features};
      #print "Exon $count: ".scalar(@evi)." features\n";
      if (@evi){
  	#foreach my $evi ( @evi ){
        #print STDERR $evi." - ".$evi->gffstring."\n";
  	#}
      }
      else{
  	;#print STDERR "No supporting evidence\n";
      }
    }
  }

#    print STDERR " Transcript  : ".$t."\n";
#    print STDERR " Translation : ".$t->translation."\n";
#    print STDERR " translation starts: ".$t->translation->start."\n";
#    print STDERR " translation end   : ".$t->translation->end."\n";
#    print STDERR " start exon        : ".$t->translation->start_Exon."\n";
#    print STDERR " end exon          : ".$t->translation->end_Exon."\n";
#    foreach my $exon (@{$t->get_all_Exons}){
#      print STDERR "     Exon          : " . $exon . " ".$exon->phase  . " " . $exon->end_phase ."\tstarts: ".$exon->start."\tends: ".$exon->end."\n";
#    }
#  }
  # tidy up output files.
   #unlink $genome_file;
   #unlink $evi_file;

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



=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: 
 Returns : value of seq
 Args    : newvalue (optional)


=cut

sub seq{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      if( !ref $value || !$value->isa("Bio::PrimarySeqI") ) {
	  $obj->throw("No primary sequence provided...");
      }

      $obj->{'_seq'} = $value;
    }
    return $obj->{'_seq'};

}

1;

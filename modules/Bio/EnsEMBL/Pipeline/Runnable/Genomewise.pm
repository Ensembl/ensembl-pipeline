
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

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;
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

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run{
  my ($self) = @_;
  
  if( !defined $self->workdir ) {
    $self->workdir('/tmp');
  }
  
  my $dir = $self->workdir();
  my $genome_file = "$dir/gowise_seq.$$";
  my $evi_file    = "$dir/gowise_evi.$$";
  
  open(F,">$genome_file") || $self->throw("Could not open $genome_file $!");
  open(E,">$evi_file") || $self->throw("Could not open $evi_file $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'fasta','-fh' => \*F);
  $seqout->write_seq($self->seq);
  
  $seqout = undef;
  close(F);
  
  foreach my $t ( $self->each_Transcript ) {
    foreach my $e ( $t->get_all_Exons ) {

      if( $e->strand == -1 ) {
	$self->throw("Genomewise cannot handle reverse strand exons. Must flip outside");
      }
      print E "exon ",$e->start," ",$e->end,"\n";
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
  
  # we try now without -smell 4

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
	  push(@{$self->{'_output_array'}},$t);

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
		      $trans->end_exon($prev);
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
		    $trans->start_exon($exon);
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
			  $trans->end_exon($exon);
			  
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
  print STDERR "genomic file: $genome_file, evidence file: $evi_file\n";
# foreach my $t ( @{ $self->{'_output_array'} } ){
#    print STDERR "\nIn Genomewise.run\n";
#    print STDERR " Transcript  : ".$t."\n";
#    print STDERR " Translation : ".$t->translation."\n";
#    print STDERR " translation starts: ".$t->translation->start."\n";
#    print STDERR " translation end   : ".$t->translation->end."\n";
#    print STDERR " start exon        : ".$t->translation->start_exon."\n";
#    print STDERR " end exon          : ".$t->translation->end_exon."\n";
#    foreach my $exon ($t->get_all_Exons){
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

   return @{$self->{'_output_array'}};
}




=head2 each_Transcript

 Title   : each_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_Transcript{
   my ($self) = @_;

   @{$self->{'_transcript_array'}};
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
#     if (!defined $val || !$val->isa('Bio::EnsEMBL::Transcript') ) {
#       if( !ref $val || !$val->isa('Bio::EnsEMBL::Transcript') ) {

 #      }
       push(@{$self->{'_transcript_array'}},$val);
   }
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


#
# Ensembl module for Bio::EnsEMBL::Pipeline::Runnable::Genomewise.pm
#
# Cared for by EWan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Genomewise.pm - DESCRIPTION of Object

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


package Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Transcript;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI);

sub new {
    my($class,@args) = @_;
    
    my $self = {};
    
    bless $self,$class;
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
       foreach my $e ( $t->each_Exon ) {
	   if( $e->strand == -1 ) {
	       $self->throw("Genomewise cannot handle reverse strand exons. Must flip outside");
	   }
	   print E "exon ",$e->start," ",$e->end,"\n";
       }
       print E "//\n";
   }
   close(E);
       
   open(GW,"genomewise $genome_file $evi_file |");
   
   while( <GW> ) {
       /Gene\s+\d+\s+\d+/ && do {
	   my $t = Bio::EnsEMBL::Transcript->new();
	   push(@{$self->{'_output_array'}},$t);
	   while( <GW> ) {
	       /Gene/ && last;
	       /Exon\s+(\d+)\s+(\d+)/ || next;
	       my $start = $1;
	       my $end   = $2;
	       my $strand = 1;

	       my $exon = Bio::EnsEMBL::Exon->new();
	       $exon->start($start);
	       $exon->end  ($end);
	       $exon->strand($strand);
	       $t->add_Exon($exon);
	   }
       };
   }

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
       if( !ref $val || !$val->isa('Bio::EnsEMBL::Transcript') ) {
	   $self->throw("[$val] is not a Transcript...");
       }
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

      $obj->{'seq'} = $value;
    }
    return $obj->{'seq'};

}

#
# Ensembl module for Bio::EnsEMBL::Pipeline::Runnable::ORF.pm
#
# Cared for by EWan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::ORF.pm - DESCRIPTION of Object

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


package Bio::EnsEMBL::Pipeline::Runnable::ORF;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::SeqFeature;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI Bio::Root::RootI);

sub new {
    my($class,@args) = @_;
    
    my $self = {};
    
    bless $self,$class;

    my ($seq,$length) = $self->_rearrange([qw(SEQ LENGTH)],@args);

    if( !defined $seq || !ref $seq || !$seq->isa('Bio::PrimarySeqI') ) {
	$self->throw("no sequence passed into runnable");
    }
    
    if( !defined $length ) {
	$length = 20;
    }

    $self->seq($seq);
    $self->length($length);

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

   my $seq = $self->seq;

   if( !defined $seq ) {
       $self->throw("No sequence in run routine");
   }

   # generate 3 frame translations in the forward frame

   my $seq0 = $seq->translate(undef,undef,0)->seq;
   my $seq1 = $seq->translate(undef,undef,1)->seq;
   my $seq2 = $seq->translate(undef,undef,2)->seq;
   my $olen = length($seq0);

   my $start0 = 0;
   my $start1 = 0;
   my $start2 = 0;

   while ($seq0 =~ /([^\*]{1,})[\*\w{$olen}]/g) {
       
       my $match = $1;
       
       my $prestring = $&;

       $prestring =~ s/\*//;
       
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($start0*3+1);
       $sf->end((length($prestring)+$start0)*3+1);
       $sf->strand(1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       $start0 = $start0 + length($prestring);

       #print STDERR $sf->start."\n";
       #print STDERR $sf->end."\n";
       $seq0 =~ s/$prestring\*//g;
       
   }
	  
   while ($seq1 =~ /([^\*]{1,})[\*\w{$olen}]/g) {
       
       my $match = $1;
       
       my $prestring = $&;

       $prestring =~ s/\*//;
       
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($start1*3+1);
       $sf->end((length($prestring)+$start1)*3+1);
       $sf->strand(1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       $start1 = $start1 + length($prestring);

       #print STDERR $sf->start."\n";
       #print STDERR $sf->end."\n";
       $seq1 =~ s/$prestring\*//g;
       
   }
	  
   while ($seq2 =~ /([^\*]{1,})[\*\w{$olen}]/g) {
       
       my $match = $1;
       
       my $prestring = $&;

       $prestring =~ s/\*//;
       
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($start2*3+1);
       $sf->end((length($prestring)+$start2)*3+1);
       $sf->strand(1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       $start2 = $start2 + length($prestring);

       #print STDERR $sf->start."\n";
       #print STDERR $sf->end."\n";
       $seq2 =~ s/$prestring\*//g;
       
   } 

   my $rev = $seq->revcom();

   $seq0 = $rev->translate(undef,undef,0)->seq;
   $seq1 = $rev->translate(undef,undef,1)->seq;
   $seq2 = $rev->translate(undef,undef,2)->seq;

   $start0 = 0;
   $start1 = 0;
   $start2 = 0;
   
   #print STDERR "SEQS: $seq0\t$seq1\t$seq2\n"; 

   my $len = $seq->length();
 while ($seq0 =~ /([^\*]{1,})[\*\w{$olen}]/g) {
       
       my $match = $1;
       
       my $prestring = $&;

       $prestring =~ s/\*//;
       
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($len - $start0*3+1);
       $sf->end($len - (length($prestring)+$start0)*3+1);
       $sf->strand(-1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       $start0 = $start0 + length($prestring);

       #print STDERR $sf->start."\n";
       #print STDERR $sf->end."\n";
       $seq0 =~ s/$prestring\*//g;
       
   }
	  
   while ($seq1 =~ /([^\*]{1,})[\*\w{$olen}]/g) {
       
       my $match = $1;
       
       my $prestring = $&;

       $prestring =~ s/\*//;
       
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($len - $start1*3+1);
       $sf->end($len - (length($prestring)+$start1)*3+1);
       $sf->strand(-1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       $start1 = $start1 + length($prestring);

       #print STDERR $sf->start."\n";
       #print STDERR $sf->end."\n";
       $seq1 =~ s/$prestring\*//g;
       
   }
	  
   while ($seq2 =~ /([^\*]{1,})[\*\w{$olen}]/g) {
       
       my $match = $1;
       
       my $prestring = $&;

       $prestring =~ s/\*//;
       
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($len - $start2*3+1);
       $sf->end(($len - length($prestring)+$start2)*3+1);
       $sf->strand(-1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       $start2 = $start2 + length($prestring);

       #print STDERR $sf->start."\n";
       #print STDERR $sf->end."\n";
       $seq2 =~ s/$prestring\*//g;
       
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

=head2 length

 Title   : length
 Usage   : $obj->length($newval)
 Function: 
 Example : 
 Returns : value of length
 Args    : newvalue (optional)


=cut

sub length{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'length'} = $value;
    }
    return $obj->{'length'};

}


=head2 add_output

 Title   : add_output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_output{
   my ($self,$val) = @_;

   push(@{$self->{'_output_array'}},$val);
}



=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: 
 Example : 
 Returns : value of seq
 Args    : newvalue (optional)


=cut

sub seq{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'seq'} = $value;
    }
    return $obj->{'seq'};

}










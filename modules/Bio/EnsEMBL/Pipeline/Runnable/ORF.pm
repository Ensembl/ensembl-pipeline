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
   my $olen = $self->length;


  
   while( $seq0 =~ /([^\*]{$olen,})[\*\s]/g ) {
       my $match = $1;
       # print STDERR "Got a $match ",length($match)," $olen \n";
       # deal with coordinate to the left - $`
       my $prestring = $`;
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start(length($prestring)*3+1);
       $sf->end((length($prestring)+length($match))*3+1);
       $sf->strand(1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       
   }
   
   while( $seq1 =~ /([^\*]{$olen,})[\*\s]/g ) {
       my $match = $1;
       # deal with coordinate to the left - $`
       my $prestring = $`;
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start(length($prestring)*3+2);
       $sf->end((length($prestring)+length($match))*3+2);
       $sf->strand(1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
   }
	  
   while( $seq2 =~ /([^\*]{$olen,})[\*\s]/g ) {
       my $match = $1;
       # deal with coordinate to the left - $`
       my $prestring = $`;
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start(length($prestring)*3+3);
       $sf->end((length($prestring)+length($match))*3+3);
       $sf->strand(1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
       
   }

   my $rev = $seq->revcom();

   $seq0 = $rev->translate(undef,undef,0)->seq;
   $seq1 = $rev->translate(undef,undef,1)->seq;
   $seq2 = $rev->translate(undef,undef,2)->seq;
   
   my $len = $seq->length();

   while( $seq0 =~ /([^\*]{$olen,})[\*s\s]/g ) {
       my $match = $1;
       # deal with coordinate to the left - $`
       my $prestring = $`;
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->end($len - length($prestring)*3);
       $sf->start($len - (length($prestring)+length($match))*3);
       $sf->strand(-1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
   }

   while( $seq1 =~ /([^\*]{$olen,})[\*\s]/g ) {
       my $match = $1;
       # deal with coordinate to the left - $`
       my $prestring = $`;
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->end($len - length($prestring)*3);
       $sf->start($len - (length($prestring)+length($match))*3);
       $sf->strand(-1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
   }
  
   while( $seq2 =~ /([^\*]{$olen,})[\*\s]/g ) {
       my $match = $1;
       # deal with coordinate to the left - $`
       my $prestring = $`;
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->end($len - length($prestring)*3);
       $sf->start($len - (length($prestring)+length($match))*3);
       $sf->strand(-1);
       $sf->{'_peptide'} = $match;
       $self->add_output($sf);
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












#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::DnaHmm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::DnaHmm - DnaHmm

=head1 SYNOPSIS

  # standard runnabledb useage

=head1 DESCRIPTION

This runnabledb is really quite involved. It takes a RawContig,
predicts all ORFs and then removes those which overlap with BLAST hits
> 150bits and all repeats. The results peptides are then screened
against seg. The final datase is then searched against a library of
HMMs (this is the first expensive operation), using hmmpfam. The HMMs
which get greater than 25 bits in that round are then passed into
genewise against the original DNA sequence. The final gene fragment
predictions are written out into the gene table.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::DnaHmm;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::ORF;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg;
use Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam;
use Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

# new() is written here 

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::EPCR object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::EPCR object
    Args    :   -dbobj:     A Bio::EnsEMBL::DB::Obj, 
                input_id:   Contig input id , 
                -analysis:  A Bio::EnsEMBL::Analysis 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    
    $self->{'_blast'}   = [];
    $self->{'_repeat'}  = [];
        
    my ( $dbobj, $input_id, $analysis, $threshold) = 
            $self->_rearrange (['DBOBJ', 'INPUT_ID', 'ANALYSIS'], @args);
    
    $self->throw('Need database handle') unless ($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DB::ObjI")  
                unless ($dbobj->isa ('Bio::EnsEMBL::DB::ObjI'));
    $self->dbobj($dbobj);
    
    $self->throw("No input id provided") unless ($input_id);
    $self->input_id($input_id);
    
    #$self->throw("Analysis object required") unless ($analysis);
    #$self->throw("Analysis object is not Bio::EnsEMBL::Analysis")
     #           unless ($analysis->isa("Bio::EnsEMBL::Analysis"));
    #$self->analysis($analysis);
    

    return $self;
}


=head2 fetch_input

 Title   : fetch_input
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_input{
   my ($self) = @_;

   my $rc = $self->dbobj->get_Contig($self->input_id);

   my $seq = Bio::PrimarySeq->new( -id => $self->input_id, -seq => $rc->seq);
   $self->seq($seq);

   push(@{$self->{'_blast'}},$rc->get_all_SimilarityFeatures);
   push(@{$self->{'_repeat'}},$rc->get_all_RepeatFeatures);



   return;

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

   my $orf = Bio::EnsEMBL::Pipeline::Runnable::ORF->new( -seq => $self->seq);
   $orf->run;
   my @orf = $orf->output;
   
   $orf = undef;

   # build up a string of seq length, setting 1 when there is 
   # repeat or Blast score > 150 bits

   my $str = '0' x $self->seq->length;

 

   foreach my $bl ( @{$self->{'_blast'}} ) {
       substr($str,$bl->start,$bl->length) = '1' x $bl->length;
   }

   foreach my $rep ( @{$self->{'_repeat'}} ) {
       substr($str,$rep->start,$rep->length) = '1' x $rep->length;
   }

   # now loop over ORF hits - if any '1s' in the string, bug out

   my @finalorf;

   foreach my $orf ( @orf ) {
       if( substr($str,$orf->start,$orf->length) =~ /1/ ) {
	   next;
       } else {
	   push(@finalorf,$orf);
       }
   }

   # loop over finalorf set, write out the files to a tempfile, seg it,
   # run hmmpfam, store the HMMs that hit in a hash

   my $count;
   my %hmmhash;

   foreach my $orf ( @finalorf ) {
       $count++;
       my $tempseq = Bio::PrimarySeq->new( -id => $self->input_id."_orf_".$count , -seq => $orf->{'_peptide'} );

       print STDERR $orf->{'_peptide'}."\n";
       
       my $seganalysis = $self->dbobj->get_AnalysisAdaptor->fetch_by_newest_logic_name('Seg');

       print STDERR "TEMPSEQ: $seganalysis\n";

       
       my $seg = Bio::EnsEMBL::Pipeline::Runnable::Protein::Seg->new(-clone => $tempseq, -analysis => $seganalysis);
       
       print STDERR "Runing Seg\n";
       $seg->run;
    
       
       my @out = $seg->output;
       
       my $count;

       foreach my $o(@out) {
	   
	   my $low_length = $o->end - $o->start + 1;
	   
	   substr($orf->{'_peptide'},$o->start,$low_length) = 'X' x $low_length;
	   
	   print STDERR $orf->{'_peptide'}."\n";
       }
   
       my $tempseq2 = Bio::PrimarySeq->new( -id => $self->input_id."seg".$count , -seq => $orf->{'_peptide'});
       
       
       my $hmmanalysis = $self->dbobj->get_AnalysisAdaptor->fetch_by_newest_logic_name("Pfam");
       
       my $pfam = Bio::EnsEMBL::Pipeline::Runnable::Protein::Hmmpfam->new(-clone => $tempseq2, -analysis=> $hmmanalysis);
       
       $pfam->run;
       
       foreach my $domain ( $pfam->output ) {
	   
	   if($domain->feature1->score > 25) {
	       
	       #Get the the genomic sequence corresponding to the orf + dowstream and upstream region (length to be defined)
	       #my $subseq_start = $orf->start - 100;
	       #if ($subseq_start < 1) {
	#	   $subseq_start = 1;
	       #}
	       
	       #my $subseq_end  = $orf->end + 100;
	       #if ($subseq_end > $self->seq->length) {
		#   $subseq_end = $self->seq->length;
	       #}
	       
	       #my $subseq = $self->seq->subseq($subseq_start,$subseq_end);
	       #$hmmhash{$domain->hmmname} = 1;
	       
	       #Push each genomic sequence corresponding to a given Hmm
	       #push (@{$hmmhash{$domain->feature2->seqname}},$subseq);
	       
	       my $hmmnames = $domain->feature2->seqname;

	       my $hmmtemp = "/tmp/hmmtemp".$$;
	       my $hmmfile = $hmmanalysis->db_file();
	       system("hmmfetch $hmmfile $hmmnames > $hmmtemp ");

	       my $gw = Bio::EnsEMBL::Pipeline::Runnable::GenewiseHmm->new( -hmmfile => $hmmtemp, -genomic => $self->seq);
	       $gw->run;

	       #my @f = $gw->output;

	       #foreach my $f (@f) {
		#   print STDERR ("Aligned output is \t" . $f->start . "\t" . $f->end . "\t" . $f->score . "\n");
		 #  print $f;
	       #}
	      
	   }
       }
   }
   return;
}

=head2 write_output

 Title   : write_output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_output{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function yet");

}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}


sub analysis {
    my ($self, $analysis) = @_;
    
    if ($analysis)
    {
        $self->throw("Not a Bio::EnsEMBL::Analysis object")
            unless ($analysis->isa("Bio::EnsEMBL::Analysis"));
        $self->{'_analysis'} = $analysis;
        $self->parameters($analysis->parameters);
    }
    return $self->{'_analysis'}
}

=head2 input_id

 Title   : input_id
 Usage   : $obj->input_id($newval)
 Function: 
 Example : 
 Returns : value of input_id
 Args    : newvalue (optional)


=cut

sub input_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'input_id'} = $value;
    }
    return $obj->{'input_id'};

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

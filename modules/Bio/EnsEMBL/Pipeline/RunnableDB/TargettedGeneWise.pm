
#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise.pm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise.pm - DESCRIPTION of Object

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


package Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );



# Object preamble - inheriets from Bio::Root::RootI


use Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;

BEGIN { print STDERR "\n\n***I'm here!***\n"; };

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::Pipeline::RunnableDBI );

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@_);    
           
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $dbobj, $input_id ) = $self->_rearrange(['DBOBJ',
						 'INPUT_ID'], @args);
       
    $self->throw("No database handle input")           unless defined($dbobj);
    $self->throw("[$dbobj] is not a Bio::EnsEMBL::Pipeline::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::Pipeline::DB::ObjI");
    $self->dbobj($dbobj);
    $dbobj->static_golden_path_type('UCSC');
    $self->throw("No input id input") unless defined($input_id);
    $self->input_id($input_id);
    
    return $self; # success - we hope!
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
   my ($self,@args) = @_;


   my $input = $self->input_id;

   if( !($input =~ /(\S+):(\d+),(\d+);(\S+)/) ) {
       $self->throw("Not a valid input id... $input");
   }

   my $fpc    = $1;
   my $start  = $2;
   my $end    = $3;
   my $pid    = $4;

   my @fp;

   push(@fp,$pid);

   my $vc = $self->dbobj->get_StaticGoldenPathAdaptor->fetch_VirtualContig_by_fpc_name($fpc);
   $self->vc($vc);

   my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( -genomic => $vc->primary_seq,-ids => \@fp);

   $self->runnable($r);

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
   my ($self,@args) = @_;

   $self->runnable->run();
   $self->convert_output

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

   return @{$self->{'_output'}};
}



=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut


sub write_output {
    my($self,@genes) = @_;

    my $db = $self->dbobj();
    my $gene_obj = $db->gene_Obj;

    my $vc = $self->vc;
    my @newgenes;
    my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
    foreach my $gene (@genes) { 
	my $newgene = $vc->convert_Gene_to_raw_contig($gene);
	push(@newgenes,$newgene);
	$newgene->_dump(\*STDERR);
	
	foreach my $trans ( $newgene->each_Transcript ) {
	    $seqio->write_seq($trans->translate);
	}
    }

    

    $self->throw("exiting before real write");

    eval {
	my $gene_obj = $self->dbobj->gene_Obj();

        foreach my $gene (@newgenes) {	    
	    $gene->type('pruned');

	    my ($geneid) = $gene_obj->get_New_external_id('gene',$GENE_ID_SUBSCRIPT,1);

	    $gene->id($geneid);
	    print (STDERR "Writing gene " . $gene->id . "\n");

            # Convert all exon ids and save in a hash
            my %namehash;
	    my @exons = $gene->each_unique_Exon();
	    my @exonids = $gene_obj->get_New_external_id('exon',$EXON_ID_SUBSCRIPT,scalar(@exons));
	    my $count = 0;
            foreach my $ex (@exons) {
		$namehash{$ex->id} = $exonids[$count];
		$ex->id($exonids[$count]);
		print STDERR "Exon id is ".$ex->id."\n";
		$count++;
            }

	    my @transcripts = $gene->each_Transcript;
	    my @transcript_ids = $gene_obj->get_New_external_id('transcript',$TRANSCRIPT_ID_SUBSCRIPT,scalar(@transcripts));
	    my @translation_ids = $gene_obj->get_New_external_id('translation',$PROTEIN_ID_SUBSCRIPT,scalar(@transcripts));
	    $count = 0;
	    foreach my $tran (@transcripts) {
		$tran->id             ($transcript_ids[$count]);
		$tran->translation->id($translation_ids[$count]);
		$count++;

		my $translation = $tran->translation;

		print (STDERR "Transcript  " . $tran->id . "\n");
		print (STDERR "Translation " . $tran->translation->id . "\n");

		foreach my $ex ($tran->each_Exon) {
                    my @sf = $ex->each_Supporting_Feature;
                    print STDERR "Supporting features are " . scalar(@sf) . "\n";

                    if ($namehash{$translation->start_exon_id} ne "") {
		      $translation->start_exon_id($namehash{$translation->start_exon_id});
                    }
                    if ($namehash{$translation->end_exon_id} ne "") {
		      $translation->end_exon_id  ($namehash{$translation->end_exon_id});
                    }
		    print(STDERR "Exon         " . $ex->id . "\n");
		}
		
	    }

	$gene_obj->write($gene);
       }
    };
    if ($@) {
	
	$self->throw("Error writing gene for " . $self->input_id . " [$@]\n");
    } else {
	# nothing
    }



}


sub convert_output {
    my ($self) =@_;
    my @tmpf = $self->runnable->output;
    my @genes;
    my $count = 1;
    my $time  = time; chomp($time);
    my $vc = $self->vc();

    foreach my $tmpf (@tmpf) {
	my $gene   = new Bio::EnsEMBL::Gene;
	my $tran   = new Bio::EnsEMBL::Transcript;
	my $transl = new Bio::EnsEMBL::Translation;

	$gene->id($self->input_id . ".genewise.$count");
	$gene->created($time);
	$gene->modified($time);
	$gene->version(1);

	$tran->id($self->input_id . ".genewise.$count");
	$tran->created($time);
	$tran->modified($time);
	$tran->version(1);

	$transl->id($self->input_id . ".genewise.$count");
	$transl->version(1);



	$gene->add_Transcript($tran);
	$tran->translation($transl);

	push(@genes,$gene);

	my $excount = 1;
	my @exons;

	foreach my $subf ($tmpf->sub_SeqFeature) {
	    my $exon = new Bio::EnsEMBL::Exon;
	    $exon->id($self->input_id . ".genewise.$count.$excount");
	    $exon->created($time);
	    $exon->modified($time);
	    $exon->version(1);
	    $exon->contig_id($vc->id);
	    $exon->seqname($vc->id);
	    $exon->start($subf->start);
	    $exon->end  ($subf->end);
	    $exon->strand($subf->strand);
	    
	    # This is dummy phase
	    $exon->phase(0);
	    push(@exons,$exon);

	    $excount++;
	}
	$count++;
	if ($exons[0]->strand == -1) {
	    @exons = sort {$b->start <=> $a->start} @exons;
	} else {
	    @exons = sort {$a->start <=> $b->start} @exons;
	}

	foreach my $exon (@exons) {
	    $tran->add_Exon($exon);
	}

	$transl->start_exon_id($exons[0]->id);
	if( $exons[0]->strand == 1 ) {
	    $transl->start($exons[0]->start);
	} else {
	    $transl->start($exons[0]->end);
	} 
	

	$transl->end_exon_id  ($exons[$#exons]->id);

	if( $exons[$#exons]->strand == 1 ) {
	    $transl->end($exons[$#exons]->end);
	} else {
	    $transl->end($exons[$#exons]->start);
	} 

    }

    $self->{'_output'} = \@genes;

}

=head2 input_id

 Title   : input_id
 Usage   : $obj->input_id($newval)
 Function: 
 Returns : value of input_id
 Args    : newvalue (optional)


=cut

sub input_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'input_id'} = $value;
    }
    return $obj->{'input_id'};

}

=head2 runnable

 Title   : runnable
 Usage   : $obj->runnable($newval)
 Function: 
 Returns : value of runnable
 Args    : newvalue (optional)


=cut

sub runnable{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'runnable'} = $value;
    }
    return $obj->{'runnable'};

}

=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)


=cut

sub vc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'vc'} = $value;
    }
    return $obj->{'vc'};

}

1;

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

Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneDNA

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


package Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneDNA;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );

# Object preamble - inheriets from Bio::Root::RootI


use Bio::Root::Object;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

           
  my( $dbobj, $input_id ) = $self->_rearrange(['DBOBJ',
					       'INPUT_ID'], @args);
  
  $self->throw("No database handle input")           unless defined($dbobj);
#    $self->throw("[$dbobj] is not a Bio::EnsEMBL::Pipeline::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::Pipeline::DB::ObjI");
  # this not right designwise
  
  $self->throw("[$dbobj] is not a Bio::EnsEMBL::DBSQL::Obj") unless $dbobj->isa("Bio::EnsEMBL::DBSQL::Obj");
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

  my $entry = $self->input_id;
  my @fps;
  my $fpc;
  my $pid; 
  
  my $start;
  my $end;
  my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;


  
  # input format: chr22:10602496,10603128;Q9UGV6:AC00012
  # changed to use : as field separator - for bsub
  if( !($entry =~ /(\S+):(\d+),(\d+):(\S+):(\S+)/) ) {
      $self->throw("Not a valid input id... $entry");
  }
  
  print STDERR "input: ".$entry . "\n";
  
  if ($fpc) { $self->throw("mixed fpc contigs") unless $fpc = $1; }
  if ($pid) { $self->throw("mixed protein hits") unless $pid = $4; }
 
  $fpc = $1;
  $pid = $4;
  my $dnaid = $5;

  my $fpcstart   = $2;
  my $fpcend     = $3;
  my $fpcstrand  = 1;
  
  if ($2 > $3) { # let blast sort it out
      $fpcstart  = $3;
      $fpcend    = $2;
      $fpcstrand = -1;
  }
  
  $start = $fpcstart unless (defined $start && $start < $fpcstart);
  $end = $fpcend unless (defined $end && $end > $fpcend);

    
  
  my $sgpa = $self->dbobj->get_StaticGoldenPathAdaptor();

  print STDERR "$fpc $start $end\n";

  my ($chrname,$chrstart,$chrend) = $sgpa->convert_fpc_to_chromosome($fpc,$start-500,$end+500);
  # I *think* the chr start and end will be handled by sgpa->fetch_rawcontigs_by_chr_start_end
  # can make sure $start > 0, but what to compare $end with?
  #  my ($chrname,$chrstart,$chrend) = $sgpa->convert_fpc_to_chromosome($fpc,$start-250000,$end+250000);
  print STDERR "$chrname $chrstart $chrend\n";
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$chrstart,$chrend);
  
  $self->vc($vc);
  
  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( -genomic => $vc->primary_seq,
								    -ids => [ $pid ] );
  
  my $cdna = $seqfetcher->run_efetch($dnaid);

  my $ex = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new( -genomic => $vc->primary_seq,
							     -est => [$cdna]
							     );



    
  $self->runnable($r);
  $self->exonerate_runnable($ex);

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

   $self->exonerate_runnable->run();
   $self->runnable->run();

   # temporary solution before michele gets 
   # into sorting out this mess...
   push(@{$self->{'_output'}},$self->runnable->output);
   push(@{$self->{'_output'}},$self->exonerate_runnable->output);


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



=head2 exonerate_runnable

 Title   : exonerate_runnable
 Usage   : $obj->exonerate_runnable($newval)
 Function: 
 Returns : value of exonerate_runnable
 Args    : newvalue (optional)


=cut

sub exonerate_runnable{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'exonerate_runnable'} = $value;
    }
    return $obj->{'exonerate_runnable'};

}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut


sub write_output {
  my ($self) = @_;

  $self->throw("Not implemented yet!");
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
	    

	    $exon->phase($subf->feature1->{_phase});
	    $exon->attach_seq($self->vc->primary_seq);
	    
	    # fix source tag and primary tag for $subf - this isn;t the right place to do this.
	    $subf->source_tag('TGW');
	    $subf->hsource_tag('TGW');
	    $subf->primary_tag('TGW');
	    $subf->hprimary_tag('TGW');
	    $subf->score(100);
	    $subf->hscore(100);

	    $exon->add_Supporting_Feature($subf);

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
	if( $exons[0]->phase == 0 ) {
	    $transl->start(1);
	} elsif( $exons[0]->phase == 1 ) {
	    $transl->start(3);
	} else {
	    $transl->start(2);
	}
	
	
	$transl->end_exon_id  ($exons[$#exons]->id);
	my $endexon = $exons[$#exons];

	if( $endexon->end_phase == 1 ) {
	    $transl->end($endexon->length -1 );
	} elsif ( $endexon->end_phase == 2 ) {
	    $transl->end($endexon->length -2 );
	} else {
	    $transl->end($endexon->length);
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

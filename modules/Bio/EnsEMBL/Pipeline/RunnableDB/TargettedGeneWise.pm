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
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Gene;

BEGIN { print STDERR "\n\n***I'm here!***\n"; };

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::Pipeline::RunnableDB );

sub new {
  my ($class,@args) = @_;
  my $self = bless {}, $class;

    $self->{'_fplist'} = []; #create key to an array of feature pairs
    
    my( $dbobj, $input_id ) = $self->_rearrange(['DBOBJ',
						 'INPUT_ID'], @args);
       
    $self->throw("No database handle input")           unless defined($dbobj);
#    $self->throw("[$dbobj] is not a Bio::EnsEMBL::Pipeline::DB::ObjI") unless $dbobj->isa("Bio::EnsEMBL::Pipeline::DB::ObjI");
    # this not right designwise

    $self->throw("[$dbobj] is not a Bio::EnsEMBL::DBSQL::Obj") unless $dbobj->isa("Bio::EnsEMBL::DBSQL::Obj");
    $self->dbobj($dbobj);
    $dbobj->static_golden_path_type('UCSC');
  
    $self->throw("No input id input") unless defined($input_id);
#    $self->throw("$input_id is not an array ref") unless ref($input_id) eq "ARRAY";
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

  my $input_id = $self->input_id;
  my $input;

  # is it an array ref or just a string?
  if (ref($input_id) eq "ARRAY") {
    $input = $input_id;
  }
  else {
    my @tmp;
    push (@tmp,$input_id);
    $input = \@tmp;
  }

  my @fps;
  my $fpc;
  my $pid; 
  
  my $start;
  my $end;
  
  foreach my $entry(@$input) {
    # input format: chr22:10602496,10603128;Q9UGV6:1,105
    # changed to use : as field separator - for bsub
#    if( !($entry =~ /(\S+):(\d+),(\d+);(\S+):(\d+),(\d+)/) ) {
    if( !($entry =~ /(\S+):(\d+),(\d+):(\S+):(\d+),(\d+)/) ) {
      $self->throw("Not a valid input id... $entry");
    }
    
    print STDERR "input: ".$entry . "\n";

    if ($fpc) { $self->throw("mixed fpc contigs") unless $fpc = $1; }
    if ($pid) { $self->throw("mixed protein hits") unless $pid = $4; }
    $fpc = $1;
    $pid = $4;
    my $fpcstart   = $2;
    my $fpcend     = $3;
    my $fpcstrand  = 1;
    
    if ($2 > $3) { # let blast sort it out
      $fpcstart  = $3;
      $fpcend    = $2;
      $fpcstrand = -1;
    }
    
    if ($5 > $6) { 
      # there's something seriously amiss
      $self->throw("Proteins only have one strand!!!\n");
    }
    
    $start = $fpcstart unless (defined $start && $start < $fpcstart);
    $end = $fpcend unless (defined $end && $end > $fpcend);

  }
  
  push (@fps,$pid);
  
  my $sgpa = $self->dbobj->get_StaticGoldenPathAdaptor();

  print STDERR "$fpc $start $end\n";


# 500bp flanking seems awfully short. Trouble is, upping it to, say, 250000 creates problems
# with "extra" blast hits. These should be got by similarity genewises later on?
#  my ($chrname,$chrstart,$chrend) = $sgpa->convert_fpc_to_chromosome($fpc,$start-500,$end+500);
  my ($chrname,$chrstart,$chrend) = $sgpa->convert_fpc_to_chromosome($fpc,$start-10000,$end+10000);
  print STDERR "$chrname $chrstart $chrend\n";
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$chrstart,$chrend);
  
  $self->vc($vc);
  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( -genomic => $vc->primary_seq,
								    -ids => \@fps);
  
  
  
  $self->runnable($r);
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
   $self->convert_output();

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
  my ($self) = @_;
  my @genes = $self->output();
  my $db = $self->dbobj();
  my $gene_obj = $db->gene_Obj;
  
  my $vc = $self->vc;
  my @newgenes;
  my $seqio = Bio::SeqIO->new(-fh => \*STDERR);
  foreach my $gene (@genes) {
    print STDERR "**** remapping\n";
    # bug in here ?
    my $newgene = $vc->convert_Gene_to_raw_contig($gene);
    push(@newgenes,$newgene);

    # VC phase not reset when remapped to raw contig???
    #    print STDERR "checking remapped translation\n";
    #    foreach my $trans ( $newgene->each_Transcript ) {
    #      $seqio->write_seq($trans->translate);
    #    }
      
    print STDERR "checking non-remapped translation\n";
    foreach my $trans ( $gene->each_Transcript ) {
      $seqio->write_seq($trans->translate);
    }

    foreach my $t ( $newgene->each_Transcript ) {
      print STDERR "  Trans ", $t->id(), " :\n";
      foreach my $e ( $t->each_Exon ) {
	print STDERR " ",$e->id(),"," , $e->start(), "-", $e->end(), "\n";
      }
      print STDERR "\n";
    }
  }
  
#  $self->throw("exiting before real write");
  
  eval {
    foreach my $gene (@newgenes) {	    
      $gene->type('pruned_TGW');
      
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

#
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer->new(-dbobj => $db
								     -input_id => $id,
								     -alnprog => 'crossmatch',
								     -min_score => $score);
    $obj->fetch_input();
    $obj->run();

    my @newfeatures = $obj->output();
    
    $obj->write_output();

=head1 DESCRIPTION

    This runnabledb creates the crossmatch runnable needed to 
    compare two contigs from databases from different organisms

    This runnabledb is more complex than usual in that it needs
    to be connected to several dbs:

    alnprog: string defining the alginment runnable used. Can be 'crossmatch', 
             'bl2seq' or 'exonarate' (the latter not implemented at the moment).
             'crossmath' is the default.

    dbobj: this is where the output is finally written, 
    storing the hits between two Ensembl databases. 
    It is actually a Bio::EnsEMBL::Comparer::DBSQL::DBAdaptor 
    database, not a normal pipeline EnsEMBL database

    query_db and target_db: these are the two databases from which the 
    input dna is fetched

    The input id has the following format: 
    db_name1:contig_id1::db_name2:contig_id2

    contig_id1 being considered as the reference contig to define the 
    reference AlignBlockSet

    min_score is an argument for the 'crossmatch' runnable to set the 
    minimum score used in the crossmatch run. In the case of 'bl2seq' runnable 
    it sets the minimum Eval used in bl2seq run.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::CrossComparer;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::CrossMatch;
use Bio::EnsEMBL::Pipeline::Runnable::Exonerate;
use Bio::EnsEMBL::Pipeline::Runnable::Bl2seq;
use Bio::PrimarySeq;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::FeaturePair;
use Bio::SeqIO;
use Bio::Root::RootI;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($score, $alnprog) = $self->_rearrange([qw(MIN_SCORE ALNPROG)],@args);
  bless $self, $class;
  $score ||= 100;
  $self->min_score($score);
  $alnprog ||= 'crossmatch';
  $self->alnprog($alnprog);
  return $self; 
}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input dna for crossmatch from the databases
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my($self) = @_; 

    $self->throw("No input id") unless defined($self->input_id);
    $self->throw("alnprog should be 'crossmatch' or 'bl2seq'or 'exonerate'")
      unless (defined $self->alnprog &&
	      ($self->alnprog eq 'crossmatch' ||
	       $self->alnprog eq 'bl2seq' ||
	       $self->alnprog eq 'exonerate'));

    my $input_id  = $self->input_id;

    my ($db1,$db2,$c1,$c2);
    if ($input_id =~ /(\S+):(\S+)::(\S+):(\S+)/ ) {
	$db1 = $1;
	$c1 = $2;
	$db2 = $3;
	$c2 = $4;
    }
    else {
	$self->throw("Input id not in correct format: got $input_id, should be parsable by (\w+)\:(\S+)\/(\w+)\:(\S+)");
    }

    $self->_c1_id($c1);
    $self->_c2_id($c2);

    my $gadp = $self->dbobj->get_GenomeDBAdaptor();

    $db1 = $gadp->fetch_by_species_tag($db1)->ensembl_db();
    $db2 = $gadp->fetch_by_species_tag($db2)->ensembl_db();

    my $contig1 = $db1->get_Contig($c1);
    my $contig2 = $db2->get_Contig($c2);

    my $seq1 = Bio::PrimarySeq->new( -display_id => 'seq1', -seq => $contig1->seq);
    my $seq2 = Bio::PrimarySeq->new( -display_id => 'seq2', -seq => $contig2->seq);
    
    my $alnrunnable;

    if ($self->alnprog eq 'crossmatch') {
      $alnrunnable = Bio::EnsEMBL::Pipeline::Runnable::CrossMatch->new(-nocopy => 1,
								       -seq1 => $seq1,
								       -seq2 => $seq2,
								       -score => $self->min_score,
								       -minmatch => 14,
								       -masklevel => 80);
    } elsif ($self->alnprog eq 'bl2seq') {
      $alnrunnable = Bio::EnsEMBL::Pipeline::Runnable::Bl2seq->new(-seq1 => $seq1,
								   -seq2 => $seq2,
								   -min_score => 40,
								   -min_eval => $self->min_score);
    } elsif ($self->alnprog eq 'exonerate') {
      $self->throw("exonarate aligment runnable not implemented yet");
#      $alnrunnable = Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(-genomic => $seq2,
#								      -est => [$seq1]);
    } 
    
    $self->runnable($alnrunnable);
}

=head2 run

    Title   :   run
    Usage   :   $self->run()
    Function:   Runs the crossmatch analysis, producing feature pairs
    Returns :   Nothing, but $self->output is filled 
                contains the crossmatch feature pairs.
    Args    :   None

=cut

sub run {
    my ($self) = @_;
    $self->throw("Can't run - no runnable objects") unless $self->runnable;
    $self->runnable->run();
    $self->output($self->runnable->output);
}

=head2 output

    Title   :   output
    Usage   :   $self->output() or $self->output(@)
    Function:   Push an array of Bio::EnsEMBL::FeaturePairs if Args is given and
                Returns all output feature pairs
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   Array of Bio::EnsEMBL::FeaturePairs (optional)

=cut

sub output {
   my ($self,@args) = @_;
   if(@args) {
     push @{$self->{'output'}}, @args;
   }
   return @{$self->{'output'}};
}

=head2 _greedy_filter

    Title   :   _greedy_filter
    Usage   :   _greedy_filter(@)
    Function:   Clean the Array of Bio::EnsEMBL::FeaturePairs in two step, 
                First, determines the highest scored hit, and fix the expected strand hit
                Second, hits on expected strand are kept if they do not overlap, 
                either on query or subject sequence, previous strored, higher scored hits.
    Returns :   Array of Bio::EnsEMBL::FeaturePairs
    Args    :   Array of Bio::EnsEMBL::FeaturePairs

=cut

sub _greedy_filter {
  my (@features) = @_;

  @features = sort {$b->score <=> $a->score} @features;

  my @features_filtered;
  my $ref_strand;
  foreach my $fp (@features) {
    if (! scalar @features_filtered) {
        push @features_filtered, $fp;
	$ref_strand = $fp->hstrand;
        next;
    }
    next if ($fp->hstrand != $ref_strand);
    my $add_fp = 1;
    foreach my $feature_filtered (@features_filtered) {
      my ($start,$end,$hstart,$hend) = ($feature_filtered->start,$feature_filtered->end,$feature_filtered->hstart,$feature_filtered->hend);
      if (($fp->start >= $start && $fp->start <= $end) ||
	  ($fp->end >= $start && $fp->end <= $end) ||
	  ($fp->hstart >= $hstart && $fp->hstart <= $hend) ||
	  ($fp->hend >= $hstart && $fp->hend <= $hend)) {
	$add_fp = 0;
	last;
      }
    }
    push @features_filtered, $fp if ($add_fp);
  }

  return @features_filtered;
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output()
    Function:   Writes contents of $self->output into $self->dbobj
    Returns :   1
    Args    :   None

=cut

sub write_output {
  my ($self) = @_;
  
  my @features = _greedy_filter($self->output);

  my $db = $self->dbobj();
  my $gadb = $db->get_GenomeDBAdaptor();
  $db->get_DnaFragAdaptor();

  my $input_id = $self->input_id();
  my ($species_tag1,$contig_id1,$species_tag2,$contig_id2);

  if ($input_id =~ /^(\S+):(\S+)::(\S+):(\S+)$/) {
    ($species_tag1,$contig_id1,$species_tag2,$contig_id2) = ($1,$2,$3,$4);
  } else {
    die "\$input_id should be dbname1:contig_id1::dbname2:contig_id2 in CrossComparer.pm\n";
  }
  
  # Using $contig_id1 as reference sequence for the reference AlignBlockSet

  my $gdb = $gadb->fetch_by_species_tag($species_tag1);
  my $dnafrag = Bio::EnsEMBL::Compara::DnaFrag->new();
  $dnafrag->name($contig_id1);
  $dnafrag->genomedb($gdb);

  # Sorting @feature from maxend to minend of reference sequence
  # and defining the offset_max for the reference AlignBlockSet
  
  @features = sort {$b->end <=> $a->end} @features;
  my $offset_max = $features[0]->end;

  # sorting @feature from minstart to maxstart of reference sequence
  # and defining the offset_min for the reference AlignBlockSet
  
  @features = sort {$a->start <=> $b->start} @features;
  my $offset_min = $features[0]->start;

  # Defining the align_row_id for the reference sequence;

  my $current_align_row_id = 1;

  # We know that reference AlignBlockSet only has one AlignBlock
  
  my $abs = Bio::EnsEMBL::Compara::AlignBlockSet->new();
  my $ab = Bio::EnsEMBL::Compara::AlignBlock->new();
  my $align_start = 1;
  my $align_end = $offset_max - $offset_min + 1;
  $ab->align_start($align_start);
  $ab->align_end($align_end);
  
  $ab->start($align_start);
  $ab->end($align_end);
  $ab->strand(1);
  $ab->dnafrag($dnafrag);
    
  $abs->add_AlignBlock($ab);
  
  # Defining an alignement and adding the reference AlignBlockSet

  my $aln = Bio::EnsEMBL::Compara::GenomicAlign->new();
  $aln->add_AlignBlockSet($current_align_row_id,$abs);
  
  # Defining the align_row_id for the query sequence;

  $current_align_row_id++;
  
  # Using $contig_id2 as query sequence for the query AlignBlockSet

  $gdb = $gadb->fetch_by_species_tag($species_tag2);
  $dnafrag = Bio::EnsEMBL::Compara::DnaFrag->new();
  $dnafrag->name($contig_id2);
  $dnafrag->genomedb($gdb);
  
  $abs = Bio::EnsEMBL::Compara::AlignBlockSet->new();
  
  foreach my $f (@features) {
    my $ab = Bio::EnsEMBL::Compara::AlignBlock->new();
    my $align_start = $f->start - $offset_min + 1;
    my $align_end = $f->end - $offset_min + 1;
    $ab->align_start($align_start);
    $ab->align_end($align_end);
    $ab->start($f->hstart);
    $ab->end($f->hend);
    if ($f->strand == 1) {
      $ab->strand($f->hstrand);
    } elsif ($f->strand == -1) {
      $ab->strand(- $f->hstrand);
    }
    $ab->dnafrag($dnafrag);
    
    $abs->add_AlignBlock($ab);
  }

  # Adding the reference AlignBlockSet to the alignment data

  $aln->add_AlignBlockSet($current_align_row_id,$abs);

  # Storing alignment in the corresponding database

  my $galnad = $db->get_GenomicAlignAdaptor();
  $galnad->store($aln);

  return 1;
}

=head2 runnable

 Title   : runnable
 Usage   : $obj->runnable($runnable)
 Function: get/set for runnable
 Returns : path to runnable
 Args    : runnable (optional)


=cut

sub runnable {
   my ($self, $runnable) = @_;

   if( defined $runnable ) {
       $self->{'_runnable'} = $runnable;
   }
   return $self->{'_runnable'};
}


=head2 min_score

 Title   : min_score
 Usage   : $obj->min_score($newval)
 Function: Getset for min_score value
 Returns : value of min_score
 Args    : newvalue (optional)


=cut

sub min_score{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'min_score'} = $value;
    }
    return $obj->{'min_score'};

}

=head2 alnprog

 Title   : alnprog
 Usage   : $obj->alnprog($string)
 Function: Get/set the alnprog value
 Returns : value of alnprog
 Args    : string, could be 'crossmatch', 'exonerate' or 'bl2seq' (optional)


=cut

sub alnprog {
   my ($self,@args) = @_;
   if (@args) {
      my ($value) = @args;
      $self->{'alnprog'} = $value;
    }
    return $self->{'alnprog'};
}

=head2 _c1_id

 Title   : _c1_id
 Usage   : $obj->_c1_id($newval)
 Function: Getset for _c1_id value
 Returns : value of _c1_id
 Args    : newvalue (optional)


=cut

sub _c1_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_c1_id'} = $value;
    }
    return $obj->{'_c1_id'};

}

=head2 _c2_id

 Title   : _c2_id
 Usage   : $obj->_c2_id($newval)
 Function: Getset for _c2_id value
 Returns : value of _c2_id
 Args    : newvalue (optional)


=cut

sub _c2_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_c2_id'} = $value;
    }
    return $obj->{'_c2_id'};

}

1;

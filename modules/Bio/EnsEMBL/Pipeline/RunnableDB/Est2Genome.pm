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

Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome->new(
					     -dbobj     => $db,
					     -input_id  => $id
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome;
use Bio::EnsEMBL::Pipeline::GeneConf qw (
					 GB_EST_DATABASES
					 GB_EST_GENETYPE
					);


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
			   -SEQFETCHER  => $sf);
                           
    Function:   creates a Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::Est2Genome object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Bio::DB::RandomAccessI Object (required)
=cut


sub new {
    my ($new,@args) = @_;
    my $self = $new->SUPER::new(@args);    
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)
    my ($type, $threshold) = $self->_rearrange([qw(TYPE THRESHOLD)], @args);
    $self->{'_fplist'} = []; #create key to an array of feature pairs
    $self->throw("no protein source databases defined in GeneConf::GB_EST_DATABASES\n") 
    unless scalar(@{$GB_EST_DATABASES});
	 
    foreach my $db(@{$GB_EST_DATABASES}){
      my $seqfetcher = $self->make_seqfetcher($db->{'index'});
      $self->add_seqfetcher_by_type($db->{'type'}, $seqfetcher);
    }
    
   return $self;
}


sub type {
  my ($self,$type) = @_;

  if (defined($type)) {
    $self->{_type} = $type;
  }
  return $self->{_type};
}

sub threshold {
  my ($self,$threshold) = @_;

  if (defined($threshold)) {
    $self->{_threshold} = $threshold;
  }
  return $self->{_threshold};
}

=head2 dbobj

    Title   :   dbobj
    Usage   :   $self->dbobj($obj);
    Function:   Gets or sets the value of dbobj
    Returns :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object
                (which extends Bio::EnsEMBL::DB::ObjI)
    Args    :   A Bio::EnsEMBL::Pipeline::DB::ObjI compliant object

=head2 input_id

    Title   :   input_id
    Usage   :   $self->input_id($input_id);
    Function:   Gets or sets the value of input_id
    Returns :   valid input id for this analysis (if set) 
    Args    :   input id for this analysis 

=head2 seqfetcher

    Title   :   seqfetcher
    Usage   :   $self->seqfetcher($seqfetcher)
    Function:   Get/set method for SeqFetcher
    Returns :   Bio::DB::RandomAccessI object
    Args    :   Bio::DB::RandomAccessI object

=head2 fetch_output

    Title   :   fetch_output
    Usage   :   $self->fetch_output
    Function:   Fetchs output data from a frozen perl object
    Returns :   array of exons (with start and end)
    Args    :   none

=cut




=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetchs input data for est2genome from the database
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my( $self) = @_;
  #print "running fetch input\n";  
  $self->throw("No input id") unless defined($self->input_id);
  
  
  my $contig    = $self->dbobj->get_RawContigAdaptor->fetch_by_name($self->input_id);
  my $genseq   = $contig->get_repeatmasked_seq;

  
 
  my $alignadaptor = $self->dbobj->get_DnaAlignFeatureAdaptor();
  
  foreach my $database(@{$GB_EST_DATABASES}){
    my %ests;
    my @estseqs;
    print STDERR "Fetching features for " . $database->{'type'} . 
      " with score above " . $database->{'threshold'}. "\n";
    
    my @features = 
      $alignadaptor->fetch_by_Contig_and_score($contig, 
					       $database->{'threshold'}, 
					       $database->{'type'});
    print STDERR "have ".scalar(@features)." to run against est2genome\n\n";
    foreach my $f (@features) {
      my $hid = $f->hseqname;
      if(!defined $ests{$hid}){
	$ests{$hid}= 1;
      }
      
    }
    my $seqfetcher = $self->get_seqfetcher_by_type($database->{'type'});
    foreach my $id(keys %ests){
      my $est = $seqfetcher->get_Seq_by_acc($id);
      
      push(@estseqs, $est);
      
    }
    my $runnable  = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniEst2Genome->new('-genomic'     => $genseq, 
									       '-queryseq' => \@estseqs,
									       '-seqfetcher'  => $seqfetcher);
    
    $self->runnable($runnable);
  
  }    

}
      
  
    
    
    
    



sub runnable {
    my ($self,$arg) = @_;
 
    if(!defined($self->{'_seqfetchers'})) {
      $self->{'_seqfetchers'} = [];
    }
    
    if (defined($arg)) {
      $self->throw("[$arg] is not a Bio::EnsEMBL::Pipeline::RunnableI") unless $arg->isa("Bio::EnsEMBL::Pipeline::RunnableI");
	
      push(@{$self->{_runnable}}, $arg);
    }

    return @{$self->{_runnable}};
}

sub run {
    my ($self) = @_;
 
    foreach my $runnable ($self->runnable) {
      $runnable || $self->throw("Can't run - no runnable object");
      #print STDERR "using ".$runnable."\n";
      #my $estseqs =  $runnable->queryseq;
      #print STDERR "there are ".scalar(@$estseqs)." to run against\n";
      #foreach my $est(@$estseqs){
#	print STDERR $est->accession_number()."\n";
 #     }
      $runnable->run;
  #    print STDERR "have finshed running \n\n";
    }
   
   
    $self->_convert_output();
    #print "have run est2genome\n";
}

sub output {
    my ($self) = @_;
    return @{$self->{_output}};
}

sub write_output {

  my($self) = @_;
  
  my @genes    = $self->output();
  
  my $db       = $self->dbobj();


  my $gene_adaptor= $self->dbobj->get_GeneAdaptor;

 
    
    
  #print "writting output there should be ".scalar(@genes)." written to the database\n"; 
 GENE: foreach my $gene (@genes) {	
    # do a per gene eval...
    eval {
      #print "gene = ".$gene."\n";
      my $dbid = $gene_adaptor->store($gene);
      #print "gene has been store as ".$dbid."\n";
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
    
  }
  
  return 1;
}


sub _convert_output {
  my ($self) = @_;
  my $count  = 1;
  my $time   = time; 
  chomp($time);
  my @genes;
  my $genetype = $self->analysis->logic_name;
  my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;
  my $analysis = $anaAdaptor->fetch_by_logic_name($genetype);
  if(!$analysis){
    # make a new analysis object
    $analysis = new Bio::EnsEMBL::Analysis
      (-db              => 'dbEST',
       -db_version      => '1',
       -program         => 'Est2Genome',
       -program_version => 1,
       -gff_source      => 'Est2Genome',
       -gff_feature     => 'gene',
       -logic_name      => 'Est2Genome',
       -module          => 'Est2Genome',
      );
  } 
  # make an array of genes for each runnable
  #print $analysis."\n";
  foreach my $runnable ($self->runnable) {
    my @g = $self->_make_genes($count, $time, $runnable, $self->analysis);
    push(@genes, @g);
  }
  #print STDOUT "genes = @genes\n";

    
  if (!defined($self->{_output})) {
    $self->{_output} = [];
  }
  
  push(@{$self->{_output}},@genes);
}


sub _make_genes {
  my ($self, $count, $time, $runnable, $analysis) = @_;
  my $contig = $self->dbobj->get_RawContigAdaptor->fetch_by_name($self->input_id);
  #print STDERR "contig = ".$contig." dbId = ".$contig->dbID."\n";
  my $genetype = $GB_EST_GENETYPE;
  my @tmpf = $runnable->output; # an array of SeqFeaturesm one per gene prediction, with subseqfeatures
  #print STDERR "we'll have " . scalar(@tmpf) . " genes\n";
  my @genes;
  #print "making genes\n";
  foreach my $tmpf(@tmpf) {
    my $gene   = Bio::EnsEMBL::Gene->new();
    my $tran   = Bio::EnsEMBL::Transcript->new();
    my $transl = Bio::EnsEMBL::Translation->new();
    #print "gene analysis = ".$analysis."\n";
    $gene->type($genetype);
    $gene->analysis($analysis);
    
  
    
    $count++;
    
    $gene->add_Transcript($tran);
    $tran->translation($transl);
    
    my $excount = 1;
    my @exons;
    
    foreach my $exon_pred ($tmpf->sub_SeqFeature) {
      # make an exon
      my $exon = new Bio::EnsEMBL::Exon;
      
     
      $exon->contig_id($contig->dbID);
    
      $exon->seqname($self->input_id);
      $exon->start($exon_pred->start);
      $exon->end  ($exon_pred->end);
      $exon->strand($exon_pred->strand);
      
      #	$exon->phase($subf->feature1->{_phase});
      
      $exon->phase($exon_pred->phase);
      $exon->attach_seq($contig->primary_seq);
      # fix source tag and primary tag for $exon_pred - this isn;t the right place to do this.
      

      $exon->score(100); # ooooooohhhhhh
      $exon->analysis($gene->analysis);
      #print "num subfeatures: " . scalar($exon_pred->sub_SeqFeature) . "\n";

      # sort out supporting evidence for this exon prediction
      foreach my $subf($exon_pred->sub_SeqFeature){
	$subf->seqname($contig->dbID);
	$subf->feature1->seqname($contig->dbID);	
	$subf->feature1->analysis($exon->analysis);
	#print " supporting feature analysis ".$exon->analysis."\n";
	$subf->feature2->analysis($exon->analysis);
	$subf->analysis($exon->analysis);
	#print $subf->analysis."\n";
	$exon->add_Supporting_Feature($subf);
      }
      
      push(@exons,$exon);
      
      $excount++;
    }
    
    if ($#exons < 0) {
      print STDERR "Odd.  No exons found\n";
    } else {
      #print "gene = ".$gene."\n";
      push(@genes,$gene);
      
      if ($exons[0]->strand == -1) {
	@exons = sort {$b->start <=> $a->start} @exons;
      } else {
	@exons = sort {$a->start <=> $b->start} @exons;
      }
      
      foreach my $exon (@exons) {
	$tran->add_Exon($exon);
      }
      
      $transl->start_exon($exons[0]);
      $transl->end_exon($exons[$#exons]);
      
      if ($exons[0]->phase == 0) {
	$transl->start(1);
      } elsif ($exons[0]->phase == 1) {
	$transl->start(3);
      } elsif ($exons[0]->phase == 2) {
	$transl->start(2);
      }
      
      my $endexon = $exons[$#exons];
      
      if( $endexon->end_phase == 1 ) {
	$transl->end($endexon->length -1 );
      } elsif ( $endexon->end_phase == 2 ) {
	$transl->end($endexon->length -2 );
      } else {
	$transl->end($endexon->length);
      }
      

    }
  }
  #print "_make genes made genes @genes\n";
  return @genes


}

sub make_seqfetcher {
  my ( $self, $index ) = @_;
 
  my $seqfetcher;
  #print STDERR "making seqfetch with index ".$index."\n";
  if($index){
    my @db = ( $index );
    $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::GetSeqs->new('-db' => \@db,);
  }
  else{
    # default to Pfetch
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  }
  
  return $seqfetcher;

}

sub add_seqfetcher_by_type{
  my ($self, $type, $seqfetcher) = @_;
  $self->throw("no type specified\n") unless defined ($type); 
  $self->throw("no suitable seqfetcher specified: [$seqfetcher]\n") 
    unless defined ($seqfetcher) && $seqfetcher->isa("Bio::DB::RandomAccessI"); 
  $self->{'_seqfetchers'}{$type} = $seqfetcher;
}


sub get_seqfetcher_by_type{
  my ($self, $type) = @_;
  my %seqfetchers = $self->each_seqfetcher_by_type;
  foreach my $dbtype(keys %seqfetchers){
    if ($dbtype eq $type){
      return $seqfetchers{$dbtype};
    }
  }
}

sub each_seqfetcher_by_type {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB hash of seqfetchers
   return %{$self->{'_seqfetchers'}};
}

sub each_seqfetcher {
  my ($self) = @_;
  
  if(!defined($self->{'_seqfetchers'})) {
     $self->{'_seqfetchers'} = {};
   }
    
  # NB array of seqfetchers
   return values(%{$self->{'_seqfetchers'}});
}

1;

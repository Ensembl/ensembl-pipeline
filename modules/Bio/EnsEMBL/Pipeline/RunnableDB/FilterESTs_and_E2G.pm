#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G->new(
                                                      -db          => $db,
                                                      -input_id    => $id,
                                                      -seq_index   => $index
                                                      );
    $obj->fetch_input
    $obj->run

    my @genes = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

ensembl-dev@ebi.ac.uk
Dan Andrews <dta@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G;

use vars qw(@ISA);
use strict;
use POSIX;  #(used for ceil in the latter part of fetch_input)

# Object preamble
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;
use Bio::EnsEMBL::Pipeline::Runnable::ESTFeatureFilter;
use Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Pipeline::Tools::BPlite;

use Bio::EnsEMBL::Pipeline::ESTConf qw (
                                        EST_REFDBHOST
                                        EST_REFDBNAME
                                        EST_REFDBUSER
                                        EST_DBNAME
                                        EST_DBHOST
                                        EST_DBUSER 
                                        EST_DBPASS
                                        EST_SOURCE
                                        EST_INDEX
                                        EST_MIN_PERCENT_ID
                                        EST_MIN_COVERAGE
                                        EST_INPUTID_REGEX
                                        EST_GENETYPE
                                        EST_FEATFILT_COVERAGE
                                        EST_FEATFILT_MINSCORE
					EST_REPEAT_MASKING
                                       );


@ISA = qw( Bio::EnsEMBL::Pipeline::RunnableDB );


=head2 new

  Args [1]   : -db - A Bio::EnsEMBL::DBSQL::DBAdaptor (required),
                     The database which output is written to 
  Args [2]   : -input_id - Contig input id (required), 
  Args [3]   : -seqfetcher - A Sequence Fetcher Object (required),
  Args [4]   : -analysis - A Bio::EnsEMBL::Analysis (optional) ;
  Example    : $self->new(-DB          => $db
                          -INPUT_ID    => $id
                          -ANALYSIS    => $analysis);
  Description: creates a Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTsobject
  Returntype : A Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTsobject
  Exceptions : None
  Caller     : General

=cut


sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
           
  # we force it to use BioIndex SeqFetcher
  my $seqfetcher = $self->make_seqfetcher();
  $self->seqfetcher($seqfetcher);
    
  # Pull config info from ESTConf.pl      
  my $refdbname = $EST_REFDBNAME;
  my $refdbuser = $EST_REFDBUSER;
  my $refdbhost = $EST_REFDBHOST;

  #print STDERR "refdb: $refdbname $refdbhost $refdbuser\n";
  my $estdbname = $EST_DBNAME;
  my $estdbuser = $EST_DBUSER;
  my $estdbhost = $EST_DBHOST;
  my $estpass   = $EST_DBPASS;

  #print STDERR "estdb: $estdbname $estdbhost $estdbuser $estpass\n";
         
  # database with the dna:
  my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $refdbhost,         
                                                 -user   => $refdbuser,
                                                 -dbname => $refdbname);
         
  # database where the exonerate est/cdna features are:
  my $estdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $estdbhost,         
                                                 -user   => $estdbuser,
                                                 -dbname => $estdbname,
                                                 -pass   => $estpass);
         
  $self->estdb($estdb);
  $self->estdb->dnadb($refdb);
         
  # need to have an ordinary adaptor to the est database for gene writes
  $self->db->dnadb($refdb);
         
  # get the appropriate analysis from the AnalysisAdaptor
  unless (defined $self->analysis){

    my $genetype = $EST_GENETYPE;

    my $anaAdaptor = $self->db->get_AnalysisAdaptor;
    my @analyses = $anaAdaptor->fetch_by_logic_name($genetype);
    if ((scalar @analyses) > 1){
      $self->throw("More than one analysis found.");
    }
    
    $self->analysis($analyses[0]);    
  }

         
  return $self;
}



=head2 estdb

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor $value
  Example    : $self->estdb($obj);
  Description: Gets or sets the value of estdb. ESTs are read from this db.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : thrown if $value arg is not a DBAdaptor
  Caller     : general

=cut

sub estdb {
    my( $self, $value ) = @_;
    
    if ($value) 
    {
        $value->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
            || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
        $self->{'_est_db'} = $value;
    }
    return $self->{'_est_db'};
}


=head2 write_output

  Arg [1]    : none 
  Example    : $self->write_output
  Description: Writes genes to db, and also writes out exons as features with an 
               appropriate analysis type
  Returntype : none
  Exceptions : thrown if the db is not available
  Caller     : run_RunnableDB

=cut

sub write_output {
    my($self) = @_;
    
    my $estdb = $self->db;

    if( !defined $estdb ) {
      $self->throw("unable to make write db");
    }
    
    $self->write_genes();
}


=head2 write_genes

  Arg [1]    : none
  Example    : $self->write_genes; 
  Description: Writes genes to db
  Returntype : none
  Exceptions : none
  Caller     : write_output

=cut

sub write_genes {
  my ($self) = @_;
  my $gene_adaptor = $self->db->get_GeneAdaptor;

 GENE: foreach my $gene ($self->output) {       
    eval {
      $gene_adaptor->store($gene);
      print STDERR "Wrote gene " . $gene->dbID . "\n";
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n" .
        $@ . "Skipping this gene\n";
    }
    
  }
}


=head2 fetch_input

  Arg [1]    : none
  Example    : $runnable->fetch_input
  Description: Fetches input databa for ExonerateESTs and makes runnable
  Returntype : none
  Exceptions : thrown if $self->input_id is not defined
  Caller     : run_RunnableDB

=cut

sub fetch_input {
  my ($self) = @_;
  
  $self->throw("No input id") unless defined($self->input_id);

  # get Slice of input region
  $self->input_id  =~ /$EST_INPUTID_REGEX/;

  my $chrid     = $1;
  my $chrstart  = $2;
  my $chrend    = $3;

  my $slice_adaptor = $self->estdb->get_SliceAdaptor();
  my $slice         = $slice_adaptor->fetch_by_chr_start_end($chrid,$chrstart,$chrend);

  $self->query($slice);

  # find exonerate features amongst all the other features  
  my $allfeatures = $self->estdb->get_DnaAlignFeatureAdaptor->fetch_all_by_Slice($slice);

  my @exonerate_features;
  my %exonerate_ests;
  my $est_source = $EST_SOURCE;

  foreach my $feat(@$allfeatures){
    unless(defined($feat->analysis) && 
           defined($feat->score) && 
           defined($feat->analysis->db) && 
           $feat->analysis->db eq $est_source) {
      $self->warn( "FilterESTs_and_E2G: something went wrong:\n" .
                   "analysis: ".$feat->analysis." analysis_db: " .
                   $feat->analysis->db." =? est_source: ".$est_source."\n");
      next;
    }      
    
    # only take high scoring ests
    if($feat->percent_id >= $EST_MIN_PERCENT_ID){
      if(!defined $exonerate_ests{$feat->hseqname}){
        push (@{$exonerate_ests{$feat->hseqname}}, $feat);
      }
      push (@exonerate_features, $feat);
    }
  }

  # empty out massive arrays
  $allfeatures = undef;

  #print STDERR "exonerate features left with percent_id >= $EST_MIN_PERCENT_ID : " 
  # . scalar(@exonerate_features) . "\n";
  #print STDERR "num ests " . scalar(keys %exonerate_ests) . "\n\n";
  
  unless( @exonerate_features ){
    print STDERR "No exonerate features left, exiting...\n";
    exit(0);
  }
  
  # filter features, current depth of coverage 10, and group successful ones by est id
  my %filtered_ests;
  

  # Apply a feature filter
  my $coverage  = $EST_FEATFILT_COVERAGE;
  my $min_score = $EST_FEATFILT_MINSCORE;

  my $filter = Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter->new('-coverage' => $coverage,
                                                                    '-minscore' => $min_score,
                                                                    '-prune'    => 1,
                                                                   );
  my @filteredfeats = $filter->run(@exonerate_features);
  
  # empty out massive arrays
  @exonerate_features = ();

  foreach my $f(@filteredfeats){
    push(@{$filtered_ests{$f->hseqname}}, $f);
  }
  #print STDERR "num filtered features ". scalar( @filteredfeats) . "\n";  

  # empty out massive arrays
  @filteredfeats = ();

  #print STDERR "num filtered ests " . scalar(keys %filtered_ests) . "\n";

  # reinstate blast
  my @blast_features = $self->blast(keys %filtered_ests);
  #print STDERR "back from blast with " . scalar(@blast_features) . " features\n";
  
  unless (@blast_features) {
    $self->warn("Odd - no exonerate features, cannot make runnables\n");
    return;
  }

  my %final_ests;

  foreach my $feat(@blast_features) {
    my $id = $feat->hseqname;

    # very annoying white space nonsense
    $id =~ s/\s//;
    $feat->hseqname($id);
    push(@{$final_ests{$id}}, $feat);
  }

  # make one runnable per EST set
  my $rcount = 0;
  my $single = 0;
  my $multi  = 0;
  
  my $efa = new Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor($self->db);
  
  # only fetch this once for the whole set or it's SLOW!
  my $genomic  = $self->query->get_repeatmasked_seq($EST_REPEAT_MASKING);
  
 ID:    
  foreach my $id(keys %final_ests) {
    # length coverage check for every EST
    
    my $hitlength;
    my $hitstart;
    my $hitend;

    foreach my $f (@{$final_ests{$id}}){
      if(!defined $hitstart || (defined $hitstart && $f->hstart < $hitstart)){
        $hitstart = $f->hstart;
      }

      if(!defined $hitend || (defined $hitend && $f->hend > $hitend)){
        $hitend = $f->hend;
      }
    }
    
    $hitlength = $hitend - $hitstart + 1;
    my $estlength = $efa->get_est_length($id);
    if(!defined $estlength || $estlength < 1){
      print STDERR "problem getting length for [$id]\n";
      next ID;
    }
    
    my $coverage = ceil(100 * ($hitlength/($estlength)));
    if($coverage < $EST_MIN_COVERAGE){
      print STDERR "rejecting $id for insufficient coverage ( < $EST_MIN_COVERAGE ): " 
        ."$coverage %\n";
      if(scalar(@{$final_ests{$id}}) == 1){
        $single++;
      }
      else{
        $multi++;
      }
      next ID;
    }
  
    # make MiniEst2Genome runnables
    # to repmask or not to repmask?    
    my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome(
                                               '-genomic'    => $genomic,
                                               '-features'   => \@{$final_ests{$id}},
                                               '-seqfetcher' => $self->seqfetcher,
                                               '-analysis'   => $self->analysis
                                                                  );
    $self->runnable($e2g);
    $rcount++;
  }

  print STDERR "number of e2gs: $rcount\n";  
  print STDERR "rejected $single single feature ests\n";
  print STDERR "rejected $multi multi feature ests\n";
}


=head2 run

  Arg [1]    : none 
  Example    : $runnable->run
  Description: runs the list of est2genome runnables generated in fetch_input and
               the converts output to remapped genes.
  Returntype : none
  Exceptions : Thrown if there are no runnables to run.
  Caller     : run_RunnableDB

=cut

sub run {
  my ($self) = @_;

  $self->throw("Can't run - no runnable objects") unless defined($self->runnable);
  
  foreach my $runnable($self->runnable) {
    $runnable->run;
  }

  $self->convert_output;

}




=head2 convert_output

  Arg [1]    : none
  Example    : $self->convert_output()
  Description: Converts est2genome output into an array of genes remapped
               into genomic coordinates
  Returntype : Nothing, but $self->{_output} contains remapped genes
  Exceptions : none
  Caller     : run

=cut

sub convert_output {
  my ($self) = @_;
  my $count  = 1;
  my @genes;

  # make an array of genes for each runnable
  foreach my $runnable ($self->runnable) {
    my @results = $runnable->output;
    foreach my $result (@results){
print STDERR "RESULT is a " . $result . "\n";
    }
    #print STDERR "runnable produced ".@results." results\n";
    my @g = $self->make_genes($count, \@results);
    #print STDERR "have made ".@g." genes\n";
    $count++;
    push(@genes, @g);
  }
###
#foreach my $gene (@genes){
#  my $transcripts = $gene->get_all_Transcripts;
#  foreach my $transcript (@$transcripts){
#    my $exons = $transcript->get_all_Exons;
#    foreach my $exon (@$exons){
#      my $supporting_evidence = $exon->get_all_supporting_features;
#      foreach my $supp_feat (@$supporting_evidence){
#       print STDERR $supp_feat . "\n";
#      }
#    }
#  }
#}
###
  my @remapped = $self->remap_genes(@genes);    
  $self->output(@remapped);
}

=head2 make_genes

  Arg [1]    : int $count: integer, 
  Arg [2]    : listref of Bio::EnsEMBL::SeqFeatures with exon sub_SeqFeatures
  Example    : $self->make_genes($count, $genetype, \@results) 
  Description: converts the output from MiniEst2Genome into Bio::EnsEMBL::Genes in
               Slice coordinates. The genes have type exonerate_e2g, 
               and have $analysis_obj attached. Each Gene has a single Transcript, 
               which in turn has Exons(with supporting features) and a Translation
  Returntype : array of Bio::EnsEMBL::Gene
  Exceptions : None
  Caller     : Internal

=cut

sub make_genes {
  my ($self, $count, $results) = @_;
  my $slice = $self->query;
  my $genetype = 'exonerate_e2g';
  my @genes;
  
  foreach my $tmpf(@$results) {
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->temporary_id($self->input_id . ".$genetype.$count");

    my $transcript = $self->make_transcript($tmpf, $self->query, $genetype, $count);

    $gene->analysis($self->analysis);
    $gene->add_Transcript($transcript);
    $count++;

    # and store it
    push(@genes,$gene);
  }
  return @genes;

}


=head2 make_transcript

  Arg [1]    : Bio::EnsEMBL::SeqFeature - a gene
  Arg [2]    : Bio::EnsEMBL::Slice - the slice in question
  Arg [3]    : String - logic name/genetype
  Arg [4]    : Int  - transcript count (optional)
  Example    : $self->make_transcript($tmpf, $self->query, $genetype, $count);
  Description: Generates transcripts from genes returned from est2genome.
  Returntype : A single Bio::EnsEMBL::Transcript
  Exceptions : Thrown if gene is not a SeqFeatureI
  Caller     : make_genes

=cut



sub make_transcript{
  my ($self, $gene, $slice, $genetype, $count) = @_;
  $genetype = 'unspecified' unless defined ($genetype);
  $count = 1 unless defined ($count);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {$self->throw("$gene must be Bio::EnsEMBL::SeqFeatureI\n");}
  

  my $transcript   = new Bio::EnsEMBL::Transcript;
  $transcript->temporary_id($slice->id . ".$genetype.$count");

  my $translation  = new Bio::EnsEMBL::Translation;    
  $translation->temporary_id($slice->id . ".$genetype.$count");

  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
     
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->temporary_id($slice->id . ".$genetype.$count.$excount");
    $exon->contig_id   ($slice->id);
    $exon->start       ($exon_pred->start);
    $exon->end         ($exon_pred->end);
    $exon->strand      ($exon_pred->strand);
    
    $exon->phase      ($exon_pred->phase);
    $exon->end_phase  ($exon_pred->end_phase );
    $exon->contig     ($slice);
    $exon->score      ($exon_pred->score);
    $exon->adaptor    ($self->estdb->get_ExonAdaptor);



    my @sfs = $exon_pred->sub_SeqFeature;
    # sort out supporting evidence for this exon prediction
    if(@sfs){
      my $align = new Bio::EnsEMBL::DnaDnaAlignFeature(-features => \@sfs);
      $align->seqname($self->input_id);
      $align->contig($slice);
      $align->score(100);
      $align->analysis($self->analysis);
      $exon->add_supporting_features($align);
    } 
    
    push(@exons,$exon);
    
    $excount++;
  }
  
  unless (@exons) {
    $self->warn("Odd.  No exons found.  Creating transcript with no exons.");
    return $transcript;
  } 

  if ($exons[0]->strand == -1) {
    @exons = sort {$b->start <=> $a->start} @exons;
  } else {
    @exons = sort {$a->start <=> $b->start} @exons;
  }
  
  foreach my $exon (@exons) {
    $transcript->add_Exon($exon);
  }
  
  $translation->start_Exon($exons[0]);
  $translation->end_Exon  ($exons[$#exons]);
  
  if ($exons[0]->phase == 0) {
    $translation->start(1);
  } elsif ($exons[0]->phase == 1) {
    $translation->start(3);
  } elsif ($exons[0]->phase == 2) {
    $translation->start(2);
  }
  
  $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
  
  return $transcript;
}


=head2 remap_genes

  Arg [1]    : A list of Bio::EnsEMBL::Gene objects
  Example    : $self->remap_genes(@genes);
  Description: Re-maps predicted genes to genomic coordinates.
  Returntype : array of Bio::EnsEMBL::Gene
  Exceptions : Warns if gene cannot be re-mapped
  Caller     : run

=cut


sub remap_genes {
  my ($self, @genes) = @_;
  my $slice = $self->query;
  my @remapped;
  
 GENEMAP:
  foreach my $gene(@genes) {
    #     print STDERR "about to remap " . $gene->temporary_id . "\n";
    my @t = @{$gene->get_all_Transcripts};
    my $tran = $t[0];
    eval {
      $gene->transform;
      # need to explicitly add back genetype and analysis.
      $gene->type($gene->type);
      $gene->analysis($gene->analysis);
      
      # temporary transfer of exon scores. Cannot deal with stickies so don't try
      
      my @oldtrans  = @{$gene->get_all_Transcripts};
      my @oldexons  = @{$oldtrans[0]->get_all_Exons};
      
      my @newtrans  = @{$gene->get_all_Transcripts};
      my @newexons  = @{$newtrans[0]->get_all_Exons};
      
      if($#oldexons == $#newexons){
        # 1:1 mapping; get_all_Exons gives ordered array of exons
        foreach( my $i = 0; $i <= $#oldexons; $i++){
          $newexons[$i]->score($oldexons[$i]->score);
        }
      }
      
      else{
        $self->warn("cannot transfer exon scores for " . $gene->id . "\n");
      }
      
      push(@remapped,$gene);
      
    };
    if ($@) {
      $self->warn("Couldn't reverse map gene " . $gene->temporary_id . " [$@]\n");
    }
   }

  return @remapped;
}

=head2 output

  Arg [1]    : A list of remapped Bio::EnsEMBL::Gene 
  Example    : $self->output(@remapped);
  Description: Get/set output list of re-mapped genes.
  Returntype : array of Bio::EnsEMBL::Gene
  Exceptions : None
  Caller     : convert_output

=cut

sub output {
   my ($self,@feat) = @_;

   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(@feat){
     push(@{$self->{'_output'}},@feat);
   }

   return @{$self->{'_output'}};
}

=head2 blast

  Arg [1]    : A list of EST id strings. 
  Example    : $self->blast(@ids)
  Description: Creates a db of ESTs and calls $self->run_blast
  Returntype : array of blast features
  Exceptions : Warns is no EST ids are passed in.
  Caller     : fetch_input

=cut

sub blast{
   my ($self, @allids) = @_;

   #print STDERR "retrieving ".scalar(@allids)." EST sequences\n";
      
   my @estseq = $self->get_Sequences(\@allids);

   if ( !scalar(@estseq) ){
     $self->warn("Odd - no ESTs retrieved\n");
     return ();
   }

   #print STDERR scalar(@estseq) . " ests retrieved\n";

   my $numests = scalar(@estseq);

   my $blastdb = $self->make_blast_db(@estseq);

   my @features = $self->run_blast($blastdb, $numests);

   unlink $blastdb;
   unlink $blastdb.".csq";
   unlink $blastdb.".nhd";
   unlink $blastdb.".ntb";
 
   return @features;
 }

=head2 get_Sequences

  Arg [1]    : A list of id strings. 
  Example    : $self->get_Sequences(@ids)
  Description: Creates a db of ESTs and calls $self->run_blast
  Returntype : array of est sequences
  Exceptions : Warns if sequences cant be fetched.
  Caller     : blast

=cut

sub get_Sequences {
  my ($self, $allids) = @_;
  my @estseq;

 ACC:
  foreach my $acc(@$allids) {
    my $seq;

    #print STDERR "getting sequence for $acc\n";
    eval{
      $seq = $self->seqfetcher->get_Seq_by_acc($acc);
    };

    if(!defined $seq){
      my $msg = "Problem fetching sequence for $acc\n";
      if(defined $@){ $msg .= "$@\n"; }
      $self->warn($msg);
    }
    else {
      push(@estseq, $seq);
    }

  }

  return (@estseq);

}

=head2 make_blast_db

  Arg [1]    : List of sequences - String list.
  Example    : $self->make_blast_db(@seqs)
  Description: Creates and formats a BLAST db 
  Returntype : a filename
  Exceptions : None
  Caller     : blast

=cut

sub make_blast_db {
    my ($self, @seq) = @_;

    my $blastfile = '/tmp/FEE_blast.' . $$ . '.fa';
    my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
                                '-file'   => ">$blastfile");

    foreach my $seq (@seq) {

      $seqio->write_seq($seq);
    }
    
    close($seqio->_filehandle);
    
    my $status = system("pressdb $blastfile");
    
    return $blastfile;
  }

=head2 run_blast

  Arg [1]    : blast database filename - String 
  Arg [2]    : number of sequences - Int
  Example    : $self->run_blast($db, $numests)
  Description: runs blast between $self->query and $db, allowing a max of 
               $numests alignments. parses output.
  Returntype : List of Bio::EnsEMBL::FeaturePair each representing a BLAST hit
  Exceptions : None
  Caller     : blast

=cut

sub run_blast {
  my ($self, $estdb, $numests) = @_;
  my @results;
  
  # prepare genomic seq
  my $seqfile  = "/tmp/FEE_genseq." . $$ . ".fa";
  my $blastout = "/tmp/FEE_blastout." . $$ . ".fa";;
  my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
                              -file   => ">$seqfile");
  $seqio->write_seq($self->query);
  close($seqio->_filehandle);

  # set B here to make sure we can show an alignment for every EST
  my $command   = "wublastn $estdb $seqfile B=" . $numests . 
    " -hspmax 1000  2> /dev/null >  $blastout";
  #print STDERR "Running BLAST:\n";
  print STDERR "$command\n";
  my $status = system( $command );
  
  my $blast_report = new Bio::EnsEMBL::Pipeline::Tools::BPlite(-file=>$blastout);

 HIT:
  while(my $hit = $blast_report->nextSbjct) {
    my $estname;

    while (my $hsp = $hit->nextHSP) {
      if(defined $estname && $estname ne $hsp->subject->seqname){
        $self->warn( "trying to switch querynames halfway through a blast hit" .
                     " for $estname - big problem!\n");
        next HIT;
      }
      else{
        $estname = $hsp->subject->seqname;
      }

      # if both genomic and est strands are the same, convention is to set both to be 1
      # if they differ, convention is to set genomic strand to -1, est strand to 1
      my $strand;
      my $hstrand = 1;      
      if($hsp->query->strand == $hsp->subject->strand) {
        $strand = 1;

      } else {
        $strand = -1;
      }
 
      my $fp = new Bio::EnsEMBL::FeaturePair();
      $fp->start  ($hsp->query->start  );
      $fp->end    ($hsp->query->end    );
      $fp->seqname($hsp->query->seqname);
      $fp->score  ($hsp->query->score  );
      $fp->strand ($strand             );

      $fp->hstart ($hsp->subject->start);
      $fp->hend   ($hsp->subject->end  );
      $fp->hstrand($hstrand            );
      $fp->hseqname($hsp->subject->seqname);

      
      #print STDERR $fp->gffstring."\n";
      if ($fp) {
        push (@results, $fp);
      }
    }
  }
  
  unlink $blastout;
  unlink $seqfile;
  
  return @results; 
    
}

=head2 make_seqfetcher

  Arg [1]    : None
  Example    : $self->make_seqfetcher
  Description: makes a Bio::EnsEMBL::SeqFetcher to be used for fetching EST sequences. If 
               $est_genome_conf{'est_index'} is specified in EST_conf.pl, then a Getseqs 
               fetcher is made, otherwise it will be Pfetch. NB for analysing large numbers 
               of ESTs eg all human ESTs, pfetch is far too slow ...
  Returntype : A seqfetcher object of some kind.
  Exceptions : Warns if index is not defined in the EST_conf
  Caller     : new

=cut


sub make_seqfetcher {
  my ( $self ) = @_;
  my $index   = $EST_INDEX;

  my $seqfetcher;
  if(defined $index && $index ne ''){
    my @db = ( $index );
  
    ## SeqFetcher to be used with 'indicate' indexing:
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher('-db' => \@db, );
    
  }
  else{
    $self->throw( "cannot create a seqfetcher from $index");
  }

  return $seqfetcher;

}

1;

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

Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G->new(
									  -dbobj     => $db,
									  -input_id  => $id,
									 );
    $obj->fetch_input
    $obj->run

    my @genes = $obj->output;


=head1 DESCRIPTION
Reads in all files with given subscript from given directory, parses out features, refilters and
passes them to MiniEst2Genome for alignment.
NB This relies on you having previously run ExonerateESTs OVER THE SAME GENOMIC CHUNK.

All a bit of a kludge at the moment :-(

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::FilterESTs_and_E2G;

use vars qw(@ISA);
use strict;
use POSIX;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome;
use Bio::EnsEMBL::Pipeline::SeqFetcher::getseqs;
use Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::FeatureAdaptor;
use Bio::Tools::BPlite;
use Bio::EnsEMBL::Pipeline::GeneConf qw (EXON_ID_SUBSCRIPT
					 TRANSCRIPT_ID_SUBSCRIPT
					 GENE_ID_SUBSCRIPT
					 PROTEIN_ID_SUBSCRIPT
					 );

@ISA = qw( Bio::EnsEMBL::Pipeline::RunnableDB );

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS    => $analysis);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object (required),
                -analysis:   A Bio::EnsEMBL::Pipeline::Analysis (optional) 
=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)

    my( $estdbname, $estdbhost, $estdbuser) = $self->_rearrange([qw(ESTDBNAME
								    ESTDBHOST
   							            ESTDBUSER)],
								 @args);

    if(!defined $self->seqfetcher) {
      # hard code to humanest blast2 databases for the moment - v naughty.

     my @dbs = qw ( /data/blastdb/humanest );
     my $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::getseqs(
								'-db'    => \@dbs,
							       );
      $self->seqfetcher($seqfetcher);
    }

    # if we have all the parameters for a estdb, make one
    # otherwise, assume the estdb must be the same as the dbobj
    if(defined $estdbname & defined $estdbhost && defined $estdbuser){
      my $estdb = new Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor(-host   => $estdbhost,		
								    -user   => $estdbuser,
								    -dbname => $estdbname,
								   );
      my $est_ext_feature_factory = $estdb->get_EstAdaptor();
      $self->dbobj->add_ExternalFeatureFactory($est_ext_feature_factory);

      # need to have an ordinary adaptor to the est database for gene writes
      my $edba = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $estdbhost,		
						    -user   => $estdbuser,
						    -dbname => $estdbname,
						   );
      $self->estdb($edba);

    }
    else { $self->throw("expecting exonerate data in an external feature factory\n"); };

    return $self;
}

=head2 estdb

    Title   :   estdb
    Usage   :   $self->estdb($obj);
    Function:   Gets or sets the value of estdb
    Returns :   A Bio::EnsEMBL::DBSQL::DBAdaptor compliant object
    Args    :   A Bio::EnsEMBL::DBSQL::DBAdaptor compliant object

=cut

sub estdb {
    my( $self, $value ) = @_;
    
    if ($value) 
    {
        $value->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
            || $self->throw("Input [$value] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
        $self->{'_estdb'} = $value;
    }
    return $self->{'_estdb'};
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes genes to db, and also writes out exons as features with an appropriate analysis type
    Returns :   
    Args    :   none

=cut

sub write_output {

    my($self) = @_;
    
    #    $self->throw("exiting before write");
    
    my $estdb = $self->estdb;
    my $refdb = $self->dbobj;

    if( !defined $estdb ) {
      $self->throw("unable to make write db");
    }
    
    if( !defined $refdb ) {
      $self->throw("unable to make ref db");
    }

    $self->write_genes();
    $self->write_exons_as_features();
}

=head2 write_genes

    Title   :   write_genes
    Usage   :   $self->write_genes
    Function:   Writes genes to db
    Returns :   nothing
    Args    :   none

=cut

sub write_genes {
  my ($self) = @_;
  my $gene_obj = $self->estdb->gene_Obj;

  my @newgenes = $self->output;
  print STDERR "genes: " . scalar(@newgenes) . "\n";
  return unless ($#newgenes >= 0);
  
  # get new ids
  eval {
    
    my $genecount  = 0;
    my $transcount = 0;
    my $translcount = 0;
    my $exoncount  = 0;
    
    # get counts of each type of ID we need.
    
    foreach my $gene ( @newgenes ) {
      $genecount++;
      
      foreach my $trans ( $gene->each_Transcript ) {
	$transcount++;
	$translcount++;
      }
      
      foreach my $exon ( $gene->each_unique_Exon() ) {
	$exoncount++;
      }
    }
    
    #	$self->throw("exiting bfore write");
    
    # get that number of ids. This locks the database
    
    my @geneids  =  $gene_obj->get_New_external_id('gene',$GENE_ID_SUBSCRIPT,$genecount);
    my @transids =  $gene_obj->get_New_external_id('transcript',$TRANSCRIPT_ID_SUBSCRIPT,$transcount);
    my @translids=  $gene_obj->get_New_external_id('translation',$PROTEIN_ID_SUBSCRIPT,$translcount);
    my @exonsid  =  $gene_obj->get_New_external_id('exon',$EXON_ID_SUBSCRIPT,$exoncount);
    
    # database locks are over.
    
    # now assign ids. gene and transcripts are easy. Exons are harder.
    # the code currently assummes that there is one Exon object per unique
    # exon id. This might not always be the case.
    
    foreach my $gene ( @newgenes ) {
      $gene->id(shift(@geneids));
      my %exonhash;
      foreach my $exon ( $gene->each_unique_Exon() ) {
	my $tempid = $exon->id;
	$exon->id(shift(@exonsid));
	$exonhash{$tempid} = $exon->id;
      }
      foreach my $trans ( $gene->each_Transcript ) {
	$trans->id(shift(@transids));
	$trans->translation->id(shift(@translids));
	$trans->translation->start_exon_id($exonhash{$trans->translation->start_exon_id});
	$trans->translation->end_exon_id($exonhash{$trans->translation->end_exon_id});
      }
      
    }
    
    # paranoia!
    if( scalar(@geneids)  != 0 || scalar(@exonsid)   != 0 || 
	scalar(@transids) != 0 || scalar(@translids) != 0 ) {
      $self->throw("In id assignment, left with unassigned ids ".
		   scalar(@geneids)  . " " .
		   scalar(@transids) . " " .
		   scalar(@translids)." " .
		   scalar(@exonsid));
    }
    
  };
  if( $@ ) {
    $self->throw("Exception in getting new ids. Exiting befor write\n\n$@" );
  }
  
  
  # this now assummes that we are building on a single VC.
  
  #    $self->throw("Bailing before real write\n");
    
 GENE: foreach my $gene (@newgenes) {	
    # do a per gene eval...
    eval {
      print STDERR $gene->id . "\n";
      # need to do EVIL things to gene_obj to get this to work	  
      $gene_obj->write($gene);
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
    
  }
}

=head2 write_exons_as_features

    Title   :   write_exons_as_features
    Usage   :   $self->write_exons_as_features
    Function:   Converts the exons into features and writes to the feature table
    Returns :   nothing
    Args    :   none

=cut


sub write_exons_as_features {
  my ($self) = @_;
  
  # for writing features
  my $feat_adaptor = $self->estdb->get_FeatureAdaptor;
  my %contig_features;
  my $source_tag   = 'est';
  my $primary_tag  = 'est';
  my $analysis     = $self->get_exon_analysis;

  $self->throw("no analysis\n") unless defined $analysis;

  # process genes
  my @genes = $self->output;
  return unless ($#genes >= 0);

  # convert exons to features
 GENE:
  foreach my $gene(@genes){
    foreach my $transcript($gene->each_Transcript){
    EXON:
      foreach my $exon($transcript->each_Exon){
	my $hstart;
	my $hend;
	my $hid;

	foreach my $sf($exon->each_Supporting_Feature){
	  if(defined $hid){
	    if ($hid ne $sf->hseqname){
	      $self->warn("trying to change hid between supporting features for same exon: " . $exon->id . "\n");
	      next EXON;
	    }
	  }
	  else{
	    $hid    = $sf->hseqname;
	  }

	  if(!defined $hstart || (defined $hstart && $hstart > $sf->hstart)){
	    $hstart = $sf->hstart;
	  }

	  if(!defined $hend   || (defined $hend   && $hend   < $sf->hend)){
	    $hend   = $sf->hend;
	  }

	}

	# score and percent_id are effectively the same for est_genome
	my $genomic = new Bio::EnsEMBL::SeqFeature  (-start       =>   $exon->start,
						     -end         =>   $exon->end,
						     -seqname     =>   $exon->contig_id,
						     -strand      =>   $exon->strand,
						     -score       =>   $exon->score,
						     -percent_id  =>   $exon->score, 
						     -phase       =>   $exon->phase,
						     -end_phase   =>   $exon->end_phase,
						     -source_tag  =>   $source_tag,
						     -primary_tag =>   $primary_tag,
						     -analysis    =>   $analysis );
	
	my $est     = new Bio::EnsEMBL::SeqFeature  (-start       =>   $hstart,
						     -end         =>   $hend,
						     -seqname     =>   $hid,
						     -strand      =>   '1',
						     -score       =>   $exon->score,
						     -percent_id  =>   $exon->score, 
						     -source_tag  =>   $source_tag,
						     -primary_tag =>   $primary_tag,
						     -analysis    =>   $analysis );
	
	my $fp      = new Bio::EnsEMBL::FeaturePair (-feature1 => $genomic,
						     -feature2 => $est) ;
	
	push(@{$contig_features{$exon->contig_id}}, $fp);
      }
    }
  }
  
  # write the features
  foreach my $contig_id( keys %contig_features){
    my $contig;
    eval{
      $contig =   $self->dbobj->get_Contig($contig_id);
    };
    if($@){
      print STDERR "No contig for $contig_id part 1\n$@\n";
    }

    # lock db only once per contig - may still be too slow.
    my @features = @{$contig_features{$contig_id}};
    print STDERR "writing exon features for $contig_id\n";
    $feat_adaptor->store($contig, @features);
  }

}

=head2 get_exon_analysis

  Title   : get_exon_analysis
  Usage   : get_exon_analysis
  Function: checks estdb for a pre-existing analysis to attach to exon features, and 
            makes a new one if necessary
  Returns : Bio::EnsEMBL::Analysis
  Args    : none

=cut


sub get_exon_analysis{

  my ($self) = @_;

  my $logicname  = 'ex_e2g_feat';
  my $anaAdaptor = $self->estdb->get_AnalysisAdaptor;
  my @analyses   = $anaAdaptor->fetch_by_logic_name($logicname);
  my $analysis;
  
  if(scalar(@analyses) > 1){
    $self->throw("panic! > 1 analysis for $logicname\n");
  }
  elsif(scalar(@analyses) == 1){
    $analysis = $analyses[0];
  }
  else{
    # only need to insert ONCE.
    $analysis = new Bio::EnsEMBL::Analysis(
					   -db              => 'dbEST',
					   -db_version      => 1,
					   -program         => 'exonerate_e2g',
					   -program_version => 3,
					   -gff_source      => 'exonerate_e2g',
					   -gff_feature     => 'similarity',
					   -logic_name      => $logicname,
					   -module          => 'Filter_ESTs_and_E2G',
					  );
  }

  return $analysis;

}

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for ExonerateESTs and makes runnable
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my ($self) = @_;
  
  $self->throw("No input id") unless defined($self->input_id);

  # get virtual contig of input region
  my $chrid     = $self->input_id;
     $chrid     =~ s/\.(.*)-(.*)//;
  my $chrstart  = $1;
  my $chrend    = $2;
  $self->dbobj->static_golden_path_type('UCSC');
  my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
  my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);
  $contig->_chr_name($chrid);
  $self->vc($contig);
  print STDERR "contig: " . $contig . "\n";

  # find exonerate features amongst all the other features  
  my @allfeatures = $contig->get_all_ExternalFeatures();
#  my @allfeatures = $self->cheat_with_feature_load();

  print STDERR "got " . scalar(@allfeatures) . " external features\n";

  my @exonerate_features;
  my %exonerate_ests;

  foreach my $feat(@allfeatures){
    if (defined($feat->analysis)      && defined($feat->score) && 
	defined($feat->analysis->db)  && $feat->analysis->db eq "dbEST") {
      # percent_id cutoff 90% for exonerate ... take all features for a sequence as long as 
      # one of them gets over this threshold
      if(($feat->percent_id > 89 || defined $exonerate_ests{$feat->hseqname})){ 
	  push (@{$exonerate_ests{$feat->hseqname}}, $feat);
	  push (@exonerate_features, $feat);
      }
    }
  }

# temporary
#  foreach (my $i = 0; $i < 40; $i++){
#    my $feat = $allfeatures[$i];
#    if(($feat->percent_id > 89 || defined $exonerate_ests{$feat->hseqname})){ 
#      push (@{$exonerate_ests{$feat->hseqname}}, $feat);
#      push (@exonerate_features, $feat);
#    }
#    else {print STDERR "rejecting " . $feat->hseqname . " percid: " . $feat->percent_id . "\n"; }
#  }
  
  print STDERR "exonerate features: " . scalar(@exonerate_features) . "\n";
  print STDERR "num ests " . scalar(keys %exonerate_ests) . "\n";
  
  # filter features, current depth of coverage 10, and group successful ones by est id
  my %filtered_ests;
  my $filter = new Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter( '-coverage' => 10,
								    '-minscore' => 500);
  my @filteredfeats = $filter->run(@exonerate_features);
  
  foreach my $f(@filteredfeats){
    push(@{$filtered_ests{$f->hseqname}}, $f);
  }
  
  print STDERR "num filtered ests " . scalar(keys %filtered_ests) . "\n";
  
  # run blast for each of the ests retained so far - previously ran exonerate on raw 
  # contigs, now need to run on virtual contig ready for making miniseq; need to be 
  # sure we have even the very small features
  # group new features by est
  my @ids = keys %filtered_ests;

  my @blast_features = $self->blast(@ids);

  print STDERR "back from blast with " . scalar(@blast_features) . " features\n";

  # make sure we can go on before we try to dosomething stupid
  if(!defined @blast_features) {
    $self->warn("Odd - no exonerate features, cannot make runnables\n");
    return;
  }

  my %final_ests;
  foreach my $feat(@blast_features) {
    push(@{$final_ests{$feat->hseqname}}, $feat);
  }

  # make one runnable per EST set
 ID:    

  foreach my $id(keys %final_ests) {
    # I have thought about doing a strand split here - sending plus and minus features separately - but
    # currently have decided to just use al the features to make a MiniSeq and let Est2Genome find the 
    # best alignment over the whole region. Otherwise we will end up with squillions of extra genes.

    my @features = @{$final_ests{$id}};

    # reject ESTs with only 1 blast hit unless that hit covers ?90% of the length of the EST
    if (scalar(@features) == 1 ){
      $id =~ s/\s//;
      my $est = $self->{'_seq_cache'}{$id};
      my $hitlength = $features[0]->hend - $features[0]->hstart + 1;
      my $coverage = ceil(100 * ($hitlength/($est->length)));
      if($coverage < 90){
	$self->warn("rejecting $id for insufficient coverage: $coverage %\n");
	next ID;
      }
    }
    
    # make MiniEst2Genome runnables
# to repmask or not to repmask?    
    my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome('-genomic'  => $self->vc->get_repeatmasked_seq,
								   '-features' => \@features,
								   '-seqfetcher' => $self->seqfetcher);
    $self->runnable($e2g);

  }
  
}

=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   Calls run method of each runnable, & converts output into remapped genes
    Returns :   Nothing
    Args    :   None

=cut

sub run {
  my ($self) = @_;

  $self->throw("Can't run - no runnable objects") unless defined($self->{_runnables});
  
  foreach my $runnable($self->runnable) {
    $runnable->run;

  }

  $self->convert_output;

}

=head2 convert_output

    Title   :   convert_output
    Usage   :   $self->convert_output()
    Function:   Converts est2genome output into an array of genes remapped into genomic coordinates
    Returns :   Nothing, but $self->{_output} contains remapped genes
    Args    :   None
=cut

# get merged features into a form where they can be stored in the database.
sub convert_output {
  my ($self) = @_;
  my $count  = 1;
  my $time   = time; chomp($time);
  my $genetype = 'exonerate_e2g';
  my @genes;

# get the appropriate analysis from the AnalysisAdaptor
  my $anaAdaptor = $self->estdb->get_AnalysisAdaptor;
  my @analyses = $anaAdaptor->fetch_by_logic_name($genetype);

  my $analysis_obj;
  if(scalar(@analyses) > 1){
    $self->throw("panic! > 1 analysis for $genetype\n");
  }
  elsif(scalar(@analyses) == 1){
    $analysis_obj = $analyses[0];
  }
  else{
    # make a new analysis object
    $analysis_obj = new Bio::EnsEMBL::Analysis
      (-db              => 'dbEST',
       -db_version      => 1,
       -program         => $genetype,
       -program_version => 1,
       -gff_source      => $genetype,
       -gff_feature     => 'gene',
       -logic_name      => $genetype,
       -module          => 'FilterESTs_and_E2G',
      );
  }

  # make an array of genes for each runnable
  foreach my $runnable ($self->runnable) {
    my @results = $runnable->output;
    my @g = $self->make_genes($count, $genetype, $analysis_obj, \@results);
    $count++;
    push(@genes, @g);
  }

  my @remapped = $self->remap_genes(@genes);	
  $self->output(@remapped);
}

=head2 make_genes

    Title   :   make_genes
    Usage   :   $self->make_genes($count, $genetype, \@results)
    Function:   converts the output from $runnable into Bio::EnsEMBL::Genes in
           $contig(VirtualContig) coordinates. The genes have type $genetype, 
           and have $analysis_obj attached. Each Gene has a single Transcript, 
           which in turn has Exons(with supporting features) and a Translation
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   $count: integer, $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis, 
           $runnable: Bio::EnsEMBL::Pipeline::RunnableI

=cut

sub make_genes {
  my ($self, $count, $genetype, $analysis_obj, $results) = @_;
  my $contig = $self->vc;
  my @genes;
  
  foreach my $tmpf(@$results) {

    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->id($self->input_id . ".$genetype.$count");
    $gene->version(1);

    my $transcript = $self->make_transcript($tmpf,$self->vc,$genetype,$count);
    $gene->analysis($analysis_obj);
    $gene->add_Transcript($transcript);
    $count++;

    # and store it
    push(@genes,$gene);
  }
  return @genes;

}

=head2 make_transcript

 Title   : make_transcript
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub make_transcript{
  my ($self, $gene, $contig, $genetype, $count) = @_;
  $genetype = 'unspecified' unless defined ($genetype);
  $count = 1 unless defined ($count);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}
  unless ($contig->isa ("Bio::EnsEMBL::DB::ContigI"))
    {print "$contig must be Bio::EnsEMBL::DB::ContigI\n";}

  my $time  = time; 
  chomp($time);

  my $transcript   = new Bio::EnsEMBL::Transcript;
  $transcript->id($contig->id . ".$genetype.$count");
  $transcript->version(1);

  my $translation  = new Bio::EnsEMBL::Translation;    
  $translation->id($contig->id . ".$genetype.$count");
  $translation->version(1);

  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->id($contig->id . ".$genetype.$count.$excount");
    $exon->contig_id($contig->id);
    $exon->created($time);
    $exon->modified($time);
    $exon->version(1);
      
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase($exon_pred->phase);
    $exon->attach_seq($contig);
    $exon->score($exon_pred->score);

    # sort out supporting evidence for this exon prediction
    foreach my $subf($exon_pred->sub_SeqFeature){
      $subf->feature1->source_tag($genetype);
      $subf->feature1->primary_tag('similarity');
      $subf->feature1->analysis($exon_pred->analysis);
	
      $subf->feature2->source_tag($genetype);
      $subf->feature2->primary_tag('similarity');
      $subf->feature2->analysis($exon_pred->analysis);
      
      $exon->add_Supporting_Feature($subf);
    }
    
    push(@exons,$exon);
    
    $excount++;
  }
  
  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
  } 
  else {
    
    print STDERR "num exons: " . scalar(@exons) . "\n";

    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }
    
    foreach my $exon (@exons) {
      $transcript->add_Exon($exon);
    }
    
    $translation->start_exon_id($exons[0]->id);
    $translation->end_exon_id  ($exons[$#exons]->id);
    
    if ($exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($exons[0]->phase == 2) {
      $translation->start(2);
    }
    
    $translation->end  ($exons[$#exons]->end - $exons[$#exons]->start + 1);
  }
  
  return $transcript;
}


=head2 remap_genes

    Title   :   remap_genes
    Usage   :   $self->remap_genes(@genes)
    Function:   Remaps predicted genes into genomic coordinates
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, array of Bio::EnsEMBL::Gene

=cut

sub remap_genes {
  my ($self, @genes) = @_;
  my $contig = $self->vc;
  my @remapped;

 GENEMAP:
  foreach my $gene(@genes) {
     print STDERR "about to remap " . $gene->id . "\n";
    my @t = $gene->each_Transcript;
    my $tran = $t[0];
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);

      # temporary transfer of exon scores. Cannot deal with stickies so don't try

      my @oldtrans = $gene->each_Transcript;
      my @oldexons  = $oldtrans[0]->each_Exon;

      my @newtrans = $newgene->each_Transcript;
      my @newexons  = $newtrans[0]->each_Exon;

      if($#oldexons == $#newexons){
	# 1:1 mapping; each_Exon gives ordered array of exons
	foreach( my $i = 0; $i <= $#oldexons; $i++){
	  $newexons[$i]->score($oldexons[$i]->score);
	}
      }

      else{
	$self->warn("cannot transfer exon scores for " . $newgene->id . "\n");
      }

      push(@remapped,$newgene);
      
    };
     if ($@) {
       print STDERR "Couldn't reverse map gene " . $gene->id . " [$@]\n";
     }
   }

  return @remapped;
}


=head2 _print_FeaturePair

    Title   :   print_FeaturePair
    Usage   :   $self->_print_FeaturePair($pair)
    Function:   Prints attributes of a Bio::EnsEMBL::FeaturePair
    Returns :   Nothing
    Args    :   A Bio::EnsEMBL::FeaturePair

=cut

sub _print_FeaturePair {
  my ($self,$pair) = @_;
  
  print $pair->seqname . "\t" . $pair->start . "\t" . $pair->end . "\t" . 
    $pair->score . "\t" . $pair->strand . "\t" . $pair->hseqname . "\t" . 
      $pair->hstart . "\t" . $pair->hend . "\t" . $pair->hstrand . "\n";
}

=head2 output

    Title   :   output
    Usage   :   $self->output
    Function:   Returns output from this RunnableDB
    Returns :   Array of Bio::EnsEMBL::Gene
    Args    :   None

=cut

sub output {
   my ($self,@feat) = @_;

   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(defined @feat){
     push(@{$self->{'_output'}},@feat);
   }

   return @{$self->{'_output'}};
}

=head2 vc

 Title   : vc
 Usage   : $obj->vc($newval)
 Function: 
 Returns : value of vc
 Args    : newvalue (optional)

=head2 estfile

 Title   : estfile
 Usage   : $obj->estfile($newval)
 Function: 
 Returns : value of estfile
 Args    : newvalue (optional)


=cut

sub estfile {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_estfile'} = $value;
    }
    return $obj->{'_estfile'};

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
   my ($self,@genes) = @_;

   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(defined @genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
}

=head2 blast

 Title   : blast
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub blast{
   my ($self, @allids) = @_;

   my @estseq = $self->get_Sequences(\@allids);
   if ( !scalar(@estseq) ){
     $self->warn("Odd - no ESTs retrieved\n");
     return ();
   }

   print STDERR scalar(@estseq) . " ests retrieved\n";

   my $numests = scalar(@estseq);

   my $blastdb = $self->make_blast_db(@estseq);

   my @features = $self->run_blast($blastdb, $numests);


   unlink $blastdb;
   unlink $blastdb.".csq";
   unlink $blastdb.".nhd";
   unlink $blastdb.".ntb";
   
   return @features;
 }

=head2 

 Title   : get_Sequences
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Sequences {
  my ($self, $allids) = @_;
  my @estseq;

 ACC:
  foreach my $acc(@$allids) {
    my $seq;

    if (defined($self->{'_seq_cache'}{$acc})){
      push (@estseq, $seq);
      next ACC;
    }

    eval{
      $seq = $self->seqfetcher->get_Seq_by_acc($acc);
    };
    if(!defined $seq){
      my $msg = "Problem fetching sequence for $acc\n";
      if(defined $@){ $msg .= "$@\n"; }
      $self->warn($msg);
    }
    else {
      $self->{'_seq_cache'}{$acc} = $seq;
      push(@estseq, $seq);
    }
  }

  return (@estseq);

}

=head2 

 Title   : make_blast_db
 Usage   : $self->make_blast_db(@seq)
 Function: creates a wublastn formatted database from @seq
 Example :
 Returns : name of blast dbfile
 Args    : @seq: Array of Bio::Seq


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


=head2 

 Title   : run_blast
 Usage   : $self->run_blast($db, $numests)
 Function: runs blast between $self->vc and $db, allowing a max of $numests alignments. parses output
 Example :
 Returns : array of Bio:EnsEMBL::FeaturePair representing blast hits
 Args    : $estdb: name of wublast formatted database; $numests: number of ests in the database


=cut

sub run_blast {
  my ($self, $estdb, $numests) = @_;
  my @results;
  
  # prepare genomic seq
  my $seqfile  = "/tmp/FEE_genseq." . $$ . ".fa";
  my $blastout = "/tmp/FEE_blastout." . $$ . ".fa";;
  my $seqio = Bio::SeqIO->new('-format' => 'Fasta',
			      -file   => ">$seqfile");
  $seqio->write_seq($self->vc);
  close($seqio->_filehandle);

  # set B here to make sure we can show an alignment for every EST
  my $command   = "wublastn $estdb $seqfile B=" . $numests . " -hspmax 1000  2> /dev/null >  $blastout";
  my $status = system( $command );
  
  my $blast_report = new Bio::Tools::BPlite(-file=>$blastout);

 HIT:
  while(my $hit = $blast_report->nextSbjct) {
    my $estname;

    while (my $hsp = $hit->nextHSP) {
      if(defined $estname && $estname ne $hsp->subject->seqname){
	$self->warn( "trying to switch querynames halfway through a blast hit for $estname - big problem!\n");
	next HIT;
      }
      else{
	$estname = $hsp->subject->seqname;
      }

      my $genomic = new Bio::EnsEMBL::SeqFeature (
						 -start       => $hsp->query->start,
						 -end         => $hsp->query->end,
						 -seqname     => $hsp->query->seqname,
						 -strand      => $hsp->query->strand,
						 -score       => $hsp->query->score,
						 -source_tag  => 'blast',
						 -primary_tag => 'blast',
						);
      
      my $est = new Bio::EnsEMBL::SeqFeature  ( -start       => $hsp->subject->start,
						-end         => $hsp->subject->end,
						-seqname     => $hsp->subject->seqname,
						-strand      => $hsp->subject->strand,
						-score       => $hsp->subject->score,
						-source_tag  => 'blast',
						-primary_tag => 'blast',
					      );

      # if both genomic and est strands are the same, convention is to set both to be 1
      # if they differ, convention is to set genomic strand to -1, est strand to 1
      if($genomic->strand == $est->strand){
	$genomic->strand(1);
	$est->strand(1);
      }
      else{
	$genomic->strand(-1);
	$est->strand(1);
      }
      #create featurepair
      my $fp = new Bio::EnsEMBL::FeaturePair  (-feature1 => $genomic,
					       -feature2 => $est) ;
      if ($fp) {
	push (@results, $fp);
      }
    }
  }
  
  unlink $blastout;
  unlink $seqfile;
  
  return @results; 
    
}

# temporary
sub cheat_with_feature_load{
  my ($self) = @_;
  my @features;

  my $featfile = "469909_est_feat";
#  my $featfile = "all_est_feat";
  open(FEAT, "<$featfile");

  while(<FEAT>){
    my ($id, $contig, $seq_start, $seq_end, $score, $strand, $analysis, 
	$name, $hstart, $hend, $hid, $evalue, $perc_id, $phase, $end_phase) = split;

    my $f1 = new Bio::EnsEMBL::SeqFeature(
					  -start       => $seq_start,
					  -end         => $seq_end,
					  -seqname     => $contig,
					  -strand      => $strand,
					  -score       => $score,
					  -source_tag  => 'est',
					  -primary_tag => 'est',
					 );
    my $f2 = new Bio::EnsEMBL::SeqFeature(
					  -start       => $hstart,
					  -end         => $hend,
					  -seqname     => $hid,
					  -strand      => '1',
					  -score       => $score,
					  -source_tag  => 'est',
					  -primary_tag => 'est',);
    my $fp = new Bio::EnsEMBL::FeaturePair(
					   -feature1 => $f1,
					   -feature2 => $f2,
					  );

    $fp->p_value($evalue);
    $fp->phase($phase);
    $fp->end_phase($end_phase);
    $fp->percent_id($perc_id);
    push (@features, $fp);
  }

  close FEAT;

  return @features;
}

1;



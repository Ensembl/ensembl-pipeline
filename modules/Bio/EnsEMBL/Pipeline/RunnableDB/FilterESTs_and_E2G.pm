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
									  -dbobj       => $db,
									  -input_id    => $id,
									  -seq_index   => $index,
									  -golden_path => $gp,
									 );
    $obj->fetch_input
    $obj->run

    my @genes = $obj->output;


=head1 DESCRIPTION

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
use Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::FeatureAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::BioIndex;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::Tools::BPlite;
use FileHandle;
#use diagnostics;

use Bio::EnsEMBL::Pipeline::ESTConf qw (
					EST_GOLDEN_PATH
					EST_REFDBHOST
					EST_REFDBNAME
					EST_REFDBUSER
					EST_DBNAME
					EST_DBHOST
					EST_DBUSER 
					EST_DBPASS
					EST_SOURCE
					EST_INDEX
				       );

					#EST_REFDBPASS #not needed, it is a reference db anyway

@ISA = qw( Bio::EnsEMBL::Pipeline::RunnableDB );

=head2 new

    Title   :   new
    Usage   :   $self->new(-DBOBJ       => $db
                           -INPUT_ID    => $id
                           -ANALYSIS      => $analysis
			   -REFDBNAME     => $refdbname
			   -REFDBHOST     => $refdbhost
			   -REFDBUSER     => $refdbuser
			   -SEQ_INDEX     => $seq_index
);
                           
    Function:   creates a 
                Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Returns :   A Bio::EnsEMBL::Pipeline::RunnableDB::ExonerateESTs
                object
    Args    :   -dbobj:      A Bio::EnsEMBL::DB::Obj (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object (required),
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
           
    # dbobj, input_id, seqfetcher, and analysis objects are all set in
    # in superclass constructor (RunnableDB.pm)

     my( $refdbname, $refdbhost, $refdbuser, $refpass, $path ) = $self->_rearrange([qw(REFDBNAME
										       REFDBHOST
										       REFDBUSER
										       REFPASS
										       GOLDEN_PATH)],
										   @args);

	 # we force it to use BioIndex SeqFetcher
	 my $seqfetcher = $self->make_seqfetcher();
	 $self->seqfetcher($seqfetcher);

    # check options in EST_conf.pl 
    #if(!defined $self->seqfetcher) {
    #  my $seqfetcher = $self->make_seqfetcher();
    #  $self->seqfetcher($seqfetcher);
    #  
    #}

    $path = $EST_GOLDEN_PATH;
    $path = 'UCSC' unless (defined $path && $path ne '');
    #print STDERR "path: $path\n";
    $self->dbobj->static_golden_path_type($path);

#print STDERR "refdb: $refdbname $refdbhost $refdbuser (refdb pass no needed really)\n";

    $refdbname = $EST_REFDBNAME unless (defined $refdbname && $refdbname ne '');
    $refdbuser = $EST_REFDBUSER unless (defined $refdbuser && $refdbuser ne '');
    $refdbhost = $EST_REFDBHOST unless (defined $refdbhost && $refdbhost ne '');
    #$refpass   = $EST_REFDBPASS unless (defined $refpass   && $refpass   ne '');

#print STDERR "refdb: $refdbname $refdbhost $refdbuser $refpass\n";

    my $estdbname = $EST_DBNAME;
    my $estdbuser = $EST_DBUSER;
    my $estdbhost = $EST_DBHOST;
    my $estpass   = $EST_DBPASS;

#print STDERR "estdb: $estdbname $estdbhost $estdbuser $estpass\n";
    # if we have all the parameters for a refdb, make one

    # otherwise, assume the refdb must be the same as the dbobj
    if(defined $refdbname & defined $refdbhost && defined $refdbuser){
      my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $refdbhost,		
						     -user   => $refdbuser,
						     -dbname => $refdbname,
						     -pass   => $refpass,
						    );

   $refdb->static_golden_path_type($path);
      my $estdb = new Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor(-host   => $estdbhost,		
								    -user   => $estdbuser,
								    -dbname => $estdbname,
								    -pass   => $estpass,
								   );
     
my $est_ext_feature_factory = $estdb->get_EstAdaptor();

      print "exff: $est_ext_feature_factory\n";
      
      $refdb->add_ExternalFeatureFactory($est_ext_feature_factory);
      $self->estdb($refdb);
      $self->estdb->static_golden_path_type($path);

      # need to have an ordinary adaptor to the est database for gene writes
      $self->dbobj->dnadb($refdb);

    }
    else { $self->throw("expecting exonerate data in an external feature factory\n"); };

    if(!defined $self->analysis){ $self->make_analysis; }

    return $self;
}

=head2 estdb

    Title   :   estdb
    Usage   :   $self->estdb($obj);
    Function:   Gets or sets the value of estdb. This is a handle to a database 
                containing dna (contig, sequence) information with the database containing 
                exonerate features as an ExternalFeatureFactory.
    Returns :   A Bio::EnsEMBL::DBSQL::DBAdaptor compliant object
    Args    :   A Bio::EnsEMBL::DBSQL::DBAdaptor compliant object

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

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes genes to db, and also writes out exons as features with an appropriate analysis type
    Returns :   
    Args    :   none

=cut

sub write_output {

    my($self) = @_;
    
    #    $self->throw("exiting before write");
    
    my $estdb = $self->dbobj;

    if( !defined $estdb ) {
      $self->throw("unable to make write db");
    }
    
    $self->write_genes();
#    $self->write_exons_as_features();
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
  my $gene_adaptor = $self->dbobj->get_GeneAdaptor;

 GENE: foreach my $gene ($self->output) {	
    eval {
      $gene_adaptor->store($gene);
      print STDERR "wrote gene " . $gene->dbID . "\n";
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
    
  }
}

# not yet ported
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
  my $feat_adaptor = $self->dbobj->get_FeatureAdaptor;
  my $contig_adaptor = $self->dbobj->get_RawContigAdaptor;
  my %contig_cache; # keep track of which contig internal_ids we have looked up so far
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
      foreach my $exon($transcript->get_all_Exons){
	my $hstart;
	my $hend;
	my $hid;

	foreach my $sf($exon->each_Supporting_Feature){
	  if(defined $hid){
	    if ($hid ne $sf->hseqname){
	      $self->warn("trying to change hid between supporting features for same exon: " . $exon->temporary_id . "\n");
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
	
	my $cid = $contig_cache{$exon->contig_id};
	if(!defined $cid){
	  $cid = $contig_adaptor->get_id_by_internal_id($exon->contig_id);
	  if (defined $cid){ $contig_cache{$exon->contig_id} = $cid; }
	}
	if(defined $cid){ push(@{$contig_features{$cid}}, $fp); }
      }
    }
  }
  
  # write the features
CONTIG:  foreach my $contig_id( keys %contig_features){
    my $contig;
    eval{
      $contig =   $self->dbobj->get_Contig($contig_id);
    };
    if($@){
      print STDERR "No contig for $contig_id part 1\n$@\n";
      next CONTIG;
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

  my $logicname  = 'est';
  my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;
  my @analyses   = $anaAdaptor->fetch_by_logic_name($logicname);
  my $analysis;
  my $est_source = $EST_SOURCE;

  if(scalar(@analyses) > 1){
    $self->throw("panic! > 1 analysis for $logicname\n");
  }
  elsif(scalar(@analyses) == 1){
    $analysis = $analyses[0];
  }
  else{
    # only need to insert ONCE.
    $analysis = new Bio::EnsEMBL::Analysis(
					   -db              => $est_source,
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
  
  print STDERR "In fetch_input\n";
  $self->throw("No input id") unless defined($self->input_id);

  # get virtual contig of input region
  my $chrid     = $self->input_id;
     $chrid     =~ s/\.(.*)-(.*)//;
  my $chrstart  = $1;
  my $chrend    = $2;

  my $stadaptor = $self->estdb->get_StaticGoldenPathAdaptor();
  my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);
  $contig->_chr_name($chrid);
  $self->vc($contig);

  # find exonerate features amongst all the other features  
  my @allfeatures = $contig->get_all_ExternalFeatures();

  print STDERR "got " . scalar(@allfeatures) . " external features\n";
  my @exonerate_features;
  my %exonerate_ests;
  my $est_source = $EST_SOURCE;

  foreach my $feat(@allfeatures){
    if (defined($feat->analysis)      && defined($feat->score) && 
	defined($feat->analysis->db)  && $feat->analysis->db eq $est_source) {
      # only take high scoring ests
      if($feat->percent_id >= 97){
	if(!defined $exonerate_ests{$feat->hseqname}){
	  push (@{$exonerate_ests{$feat->hseqname}}, $feat);
	}
	push (@exonerate_features, $feat);
      }
    }
  }

  # empty out massive arrays
  @allfeatures = ();

  print STDERR "exonerate features left with percent_id > 97 : " . scalar(@exonerate_features) . "\n";
  print STDERR "num ests " . scalar(keys %exonerate_ests) . "\n\n";
  
  # filter features, current depth of coverage 10, and group successful ones by est id
  my %filtered_ests;

  #my @time1 = times();
  # use coverage 5 for now.
  my $filter = new Bio::EnsEMBL::Pipeline::Runnable::FeatureFilter( '-coverage' => 10,
								    '-minscore' => 500,
								    '-prune'    => 1,
								  );
  my @filteredfeats = $filter->run(@exonerate_features);
  #my @time2 = times();
  #print STDERR "Filter time: user = ".($time2[0] - $time1[0])."\tsystem = ".($time2[1] - $time1[1])."\n";
  
  # empty out massive arrays
  @exonerate_features = ();

  foreach my $f(@filteredfeats){
    push(@{$filtered_ests{$f->hseqname}}, $f);
  }
  print STDERR "num filtered features ". scalar( @filteredfeats) . "\n";  

  # empty out massive arrays
  @filteredfeats = ();

  print STDERR "num filtered ests " . scalar(keys %filtered_ests) . "\n";

# reinstate blast
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
    my $id = $feat->hseqname;
    # print STDERR "id-$id-\n";
    # very annoying white space nonsense
    $id =~ s/\s//;
    $feat->hseqname($id);
    push(@{$final_ests{$id}}, $feat);
    my @fe = @{$final_ests{$id}};
  }

  # make one runnable per EST set
  my $rcount = 0;
  my $single = 0;
  my $multi  = 0;
  
  my $efa = new Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor($self->dbobj);

  # only fetch this once for the whole set or it's SLOW!
  my $genomic  = $self->vc->get_repeatmasked_seq;
  
  # keep track of those ESTs who make it into a MiniEst2genome
  my %accepted_ests;
  
 ID:    
  foreach my $id(keys %final_ests) {
    # length coverage check for every EST
    
    my $hitlength;
    my $hitstart;
    my $hitend;
    foreach my $f(@{$final_ests{$id}}){
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
      print STDERR "problem getting length for $id\n";
      next ID;
    }


    my $coverage = ceil(100 * ($hitlength/($estlength)));
    if($coverage < 90){
      print STDERR "rejecting $id for insufficient coverage: $coverage %\n";
      if(scalar(@{$final_ests{$id}}) == 1){
	$single++;
      }
      else{
	$multi++;
      }
      next ID;
    }
  
    # before making a MiniEst2Genome, check that the one we're about to create
    # is not redundant with any one we have created before
    my $do_comparison_stuff = 0;
    if ( $do_comparison_stuff == 1 ){
    
      foreach my $id2 ( keys( %accepted_ests ) ){
	
	# compare $id with each $id2
	# if $id is redundant, skip it
	my @feat1 = sort{ $a->start <=> $b->start } @{$final_ests{$id}};
	my @feat2 = sort{ $a->start <=> $b->start } @{$accepted_ests{$id2} };
	#print STDERR "comparing ".$id."(".scalar(@feat1).") with ".$id2." (".scalar(@feat2).")\n";    
	
	if ( scalar( @feat1 ) == scalar( @feat2 ) ){
	  print STDERR "$id and $id2 have the same number of features\n";
	  
	  # first, let's make a straightforward check for exac matches:
	  my $label = 0;
	  while ( $label < scalar( @feat1 )                      &&
		  $feat1[$label]->start == $feat2[$label]->start &&
		  $feat1[$label]->end   == $feat2[$label]->end   ){	        
	    print STDERR ($label+1)." == ".($label+1)."\n";
	    $label++;
	  }
	  if ( $label == scalar( @feat1 ) ){
	    print STDERR "EXACT MATCH between $id and $id2 features, skipping $id\n";
	  }
	  # make also a test for overlaps
	  $label = 0;
	  while ( $label < scalar( @feat1 ) && $feat1[$label]->overlaps( $feat2[$label] )  ){	        
	    print STDERR ($label+1)." overlaps ".($label+1)."\t";
	    print STDERR $feat1[$label]->start.":".$feat1[$label]->end."   ".
	      $feat2[$label]->start.":".$feat2[$label]->end."\n";
	    $label++;
	  }
	  if ( $label == scalar( @feat1 ) ){
	    print STDERR "approximate MATCH between $id and $id2 features, skipping $id\n";
	  }		
	}
      }
      
    }

    # make MiniEst2Genome runnables
    # to repmask or not to repmask?    
    my $e2g = new Bio::EnsEMBL::Pipeline::Runnable::MiniEst2Genome(
								   '-genomic'  => $genomic,
								   '-features' => \@{$final_ests{$id}},
								   '-seqfetcher' => $self->seqfetcher,
								   '-analysis' => $self->analysis
								  );
    $self->runnable($e2g);
    $rcount++;
  
    # store in a hash of arrays the features put in a MiniEst2Genome
    $accepted_ests{$id} = $final_ests{$id};
  }

  print STDERR "number of e2gs: $rcount\n";  
  print STDERR "rejected $single single feature ests\n";
  print STDERR "rejected $multi multi feature ests\n";
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

  $self->throw("Can't run - no runnable objects") unless defined($self->runnable);
  
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
  my @genes;

  # make an array of genes for each runnable
  foreach my $runnable ($self->runnable) {
    my @results = $runnable->output;
    my @g = $self->make_genes($count, \@results);
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
    Args    :   $count: integer, $runnable: Bio::EnsEMBL::Pipeline::RunnableI

=cut

sub make_genes {
  my ($self, $count, $results) = @_;
  my $contig = $self->vc;
  my $genetype = 'exonerate_e2g';
  my @genes;
  
  foreach my $tmpf(@$results) {

    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->temporary_id($self->input_id . ".$genetype.$count");

    my $transcript = $self->make_transcript($tmpf, $self->vc, $genetype, $count);
    $gene->analysis($self->analysis);
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

  my $transcript   = new Bio::EnsEMBL::Transcript;
  $transcript->temporary_id($contig->id . ".$genetype.$count");

  my $translation  = new Bio::EnsEMBL::Translation;    
  $translation->temporary_id($contig->id . ".$genetype.$count");

  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
     
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->temporary_id($contig->id . ".$genetype.$count.$excount");
    $exon->contig_id($contig->id);
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
      $subf->feature1->analysis($self->analysis);
	
      $subf->feature2->source_tag($genetype);
      $subf->feature2->primary_tag('similarity');
      $subf->feature2->analysis($self->analysis);
      
      $exon->add_Supporting_Feature($subf);
    }
    
    push(@exons,$exon);
    
    $excount++;
  }
  
  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
  } 
  else {
    
#    print STDERR "num exons: " . scalar(@exons) . "\n";

    if ($exons[0]->strand == -1) {
      @exons = sort {$b->start <=> $a->start} @exons;
    } else {
      @exons = sort {$a->start <=> $b->start} @exons;
    }
    
    foreach my $exon (@exons) {
      $transcript->add_Exon($exon);
    }
    
    $translation->start_exon($exons[0]);
    $translation->end_exon  ($exons[$#exons]);
    
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
    #     print STDERR "about to remap " . $gene->temporary_id . "\n";
    my @t = $gene->each_Transcript;
    my $tran = $t[0];
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);
      
      # temporary transfer of exon scores. Cannot deal with stickies so don't try
      
      my @oldtrans = $gene->each_Transcript;
      my @oldexons  = $oldtrans[0]->get_all_Exons;
      
      my @newtrans = $newgene->each_Transcript;
      my @newexons  = $newtrans[0]->get_all_Exons;
      
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
      print STDERR "Couldn't reverse map gene " . $gene->temporary_id . " [$@]\n";
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

   print STDERR "retrieving ".scalar(@allids)." EST sequences\n";
   
   my $time1 = time();
   my @estseq = $self->get_Sequences(\@allids);
   my $time2 = time();
   print STDERR "SeqFetcher time: user = ".($time2 - $time1)."\n";
   #print STDERR "SeqFetcher time: user = ".($time2[0] - $time1[0])."\tsystem = ".($time2[1] - $time1[1])."\n";

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
   # empty seq array
   @estseq = ();

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

#    if (defined($self->{'_seq_cache'}{$acc})){
#      push (@estseq, $seq);
#      next ACC;
#    }

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
#      $self->{'_seq_cache'}{$acc} = $seq;
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
  #print STDERR "Running BLAST:\n";
  #print STDERR "$command\n";
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

 Title   : make_seqfetcher
 Usage   :
 Function: makes a Bio::EnsEMBL::SeqFetcher to be used for fetching EST sequences. If 
           $est_genome_conf{'est_index'} is specified in EST_conf.pl, then a Getseqs 
           fetcher is made, otherwise it will be Pfetch. NB for analysing large numbers 
           of ESTs eg all human ESTs, pfetch is far too slow ...
 Example :
 Returns : Bio::EnsEMBL::SeqFetcher
 Args    :


=cut

sub make_seqfetcher {
  print STDERR "making a seqfetcher\n";
  my ( $self ) = @_;
  my $index   = $EST_INDEX;

  my $seqfetcher;
  if(defined $index && $index ne ''){
    #my $index = '/data/blastdb/Ensembl/human.ests';
    my @db = ( $index );
    
    #$seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs('-db' => \@db,);
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::BioIndex('-db' => \@db,);
  }
  #if(defined $index && $index ne ''){
  #  my @db = ( $index );
  #  $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs('-db' => \@db,);
  #}
  #else{
  #  # default to Pfetch
  #  $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  #}
  else{
    $self->throw( "cannot create a seqfetcher BioIndex from $index");
  }

  return $seqfetcher;

}

sub make_analysis {
  my ($self) = @_;
  
  # get the appropriate analysis from the AnalysisAdaptor
  my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;
  my @analyses = $anaAdaptor->fetch_by_logic_name($self->genetype);
  
  my $analysis_obj;
  if(scalar(@analyses) > 1){
    $self->throw("panic! > 1 analysis for " . $self->genetype . "\n");
  }
  elsif(scalar(@analyses) == 1){
    $analysis_obj = $analyses[0];
  }
  else{
    # make a new analysis object
    $analysis_obj = new Bio::EnsEMBL::Analysis
      (-db              => 'dbEST',
       -db_version      => 1,
       -program         => $self->genetype,
       -program_version => 1,
       -gff_source      => $self->genetype,
       -gff_feature     => 'gene',
       -logic_name      => $self->genetype,
       -module          => 'FilterESTs_and_E2G',
      );
  }

  $self->analysis($analysis_obj);

}

sub genetype {
  my ($self) = @_;
  return 'exonerate_e2g';
}

1;

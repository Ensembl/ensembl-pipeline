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

Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder->new(
								       -dbobj     => $db,
								       -input_id  => $id
								      );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

New version of EST_GeneBuilder to process est2genome gene predictions and feed them
to genomewise to create transcripts with translations and UTRs.
The algorithm is different from the original EST_GeneBuilder. This one is more straightforward.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder;

use diagnostics;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise;
use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::TranscriptCluster;

# config file; parameters searched for here if not passed in as @args
# do FILE == the file FILE gets executed as script
do "Bio/EnsEMBL/Pipeline/EST_conf.pl";

# use new Adaptor to get some extra info from the ESTs
use Bio::EnsEMBL::Pipeline::DBSQL::ESTFeatureAdaptor;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    # dbobj input_id mandatory and read in by BlastableDB
    if (!defined $self->seqfetcher) {
      my $seqfetcher = $self->make_seqfetcher();
      $self->seqfetcher($seqfetcher);
    }
    # dbobj needs a reference dna database
    my $refdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
						  -host             => $::db_conf{'refdbhost'},
						  -user             => $::db_conf{'refdbuser'},
						  -dbname           => $::db_conf{'refdbname'},
						  -pass             => $::db_conf{'refdbpass'},
						  -perlonlyfeatures => 0,
						 );

    $self->dbobj->dnadb($refdb);

    my $path = $::db_conf{'golden_path'};
    $path    = 'UCSC' unless (defined $path && $path ne '');
    $self->dbobj->static_golden_path_type($path);

    return $self; 
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
  my ( $self ) = @_;
  my $index   = $::est_genome_conf{'est_index'};

  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs('-db' => \@db,);
  }
  else{
    # default to Pfetch
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
  }

  return $seqfetcher;

}


=head2 RunnableDB methods

=head2 analysis

    Title   :   analysis
    Usage   :   $self->analysis($analysis);
    Function:   Gets or sets the stored Analusis object
    Returns :   Bio::EnsEMBLAnalysis object
    Args    :   Bio::EnsEMBL::Analysis object

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

############################################################

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

sub write_output {
  my ($self) = @_;
  
  my $gene_adaptor = $self->dbobj->get_GeneAdaptor;
  
 GENE: 
  foreach my $gene ($self->output) {	
    # test
    #my @transcripts = $gene->each_Transcript;
    #my $tran = $transcripts[0];
    #my $strand = $tran->start_exon->strand;
    #print STDERR "EST_GeneBuilder. you would be writting a gene on strand $strand\n";
      
    eval {
     $gene_adaptor->store($gene);
      print STDERR "wrote gene " . $gene->dbID . "\n";
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
  }
}

############################################################

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :    string: chr1.1-10000

=cut

sub fetch_input {
    my( $self) = @_;
    my $strand;

    print STDERR "IN THE 121 BRANCH\n";

    # the type of the genes being read is specified in Bio/EnsEMBL/Pipeline/EST_conf.pl 
    my $genetype =  $::est_genome_conf{'genetype'};

    #print STDERR "Fetching input: " . $self->input_id. " \n";
    $self->throw("No input id") unless defined($self->input_id);

    # get genomic region 
    my $chrid    = $self->input_id;
    $chrid       =~ s/\.(.*)-(.*)//;
    my $chrstart = $1;
    my $chrend   = $2;

    print STDERR "Chromosome id = $chrid , range $chrstart $chrend\n";

    my $stadaptor = $self->dbobj->get_StaticGoldenPathAdaptor();
    my $contig    = $stadaptor->fetch_VirtualContig_by_chr_start_end($chrid,$chrstart,$chrend);
    $contig->_chr_name($chrid);
    $self->vc($contig);
    print STDERR "got vc\n";
    print STDERR "length ".$contig->length."\n";
    
    # forward strand
    $strand = 1;
    print STDERR "\n****** forward strand ******\n\n";

    # get genes
    my @genes  = $contig->get_Genes_by_Type($genetype);
    
    print STDERR "Number of genes = " . scalar(@genes) . "\n";
    if(!scalar(@genes)){
      $self->warn("No forward strand genes found");
    }
    my @plus_transcripts;
    my $single = 0;
    # split by strand
GENE:    
    foreach my $gene (@genes) {
      my @transcripts = $gene->each_Transcript;
      
      # skip genes with more than one transcript
      if( scalar(@transcripts) > 1 ){
	$self->warn($gene->temporary_id . " has more than one transcript - skipping it\n");
	next GENE;
      }

      # skip genes with one exon
      my @exons = $transcripts[0]->get_all_Exons;
      if(scalar(@exons) == 1){
	$single++;
	next GENE;
      }

      # keep only genes in the forward strand
      if($exons[0]->strand == 1){
	push (@plus_transcripts, $transcripts[0]);
	next GENE;
      }
    }  # end of GENE

    print STDERR "In EST_GeneBuilder.fetch_input(): ".scalar(@plus_transcripts) . " forward strand genes\n";
    print STDERR "($single single exon genes thrown away)\n";

    # process transcripts in the foread strand
    if( scalar(@plus_transcripts) ){
      
      my @transcripts  = $self->_process_Transcripts(\@plus_transcripts,$strand);
      
      # make a genomewise runnable for each cluster of transcripts
      foreach my $tran (@transcripts){
	
	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
									    -genomic => $contig->primary_seq,
      								   );
#	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise();

	$self->add_runnable($runnable,$strand);
#	$runnable->seq($contig->primary_seq);
	$runnable->add_Transcript($tran);
      }
    }
    
    
    # minus strand - flip the vc and hope it copes ...
    # but SLOW - get the same genes twice ...
    
    print STDERR "\n****** reverse strand ******\n\n";

    $strand = -1;
    my $revcontig = $contig->invert;
    my @revgenes  = $revcontig->get_Genes_by_Type($genetype);
    my @minus_transcripts;
    
    print STDERR "Number of genes = " . scalar(@revgenes) . "\n";
    if(!scalar(@genes)){
      $self->warn("No reverse strand genes found");
    }

    $single=0;
  REVGENE:    
    foreach my $gene (@revgenes) {
      my @transcripts = $gene->each_Transcript;
      
      # throw away genes with more than one transcript
      if(scalar(@transcripts) > 1 ){
	$self->warn($gene->temporary_id . " has more than one transcript - skipping it\n");
	next REVGENE;
      }
      
      my @exons = $transcripts[0]->get_all_Exons;
      
      # throw away single-exon genes
      if(scalar(@exons) == 1){
	$single++;
	next REVGENE;
      }
      
      # these are really - strand, but the VC is reversed, so they are realtively + strand
      elsif($exons[0]->strand == 1){
	push (@minus_transcripts, $transcripts[0]);
	next REVGENE;
      }
    }
    print STDERR "In EST_GeneBuilfer.fetch_input(): ".scalar(@minus_transcripts) . " reverse strand genes\n";
    print STDERR "($single single exon genes thrown away)\n";
    
    if(scalar(@minus_transcripts)){
      
      my @transcripts = $self->_process_Transcripts(\@minus_transcripts,$strand);  
      
      foreach my $tran (@transcripts) {
	
	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
									    -genomic => $revcontig->primary_seq,
      								   );
#	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise();

	$self->add_runnable($runnable, $strand);
#	$runnable->seq($contig->primary_seq);
	$runnable->add_Transcript($tran);
      }
    }
}

############################################################

=head2 _flush_Transcripts

    Title   :   _flush_Transcripts
    Usage   :   $self->_flush_Transcripts
    Function:   it empties out the array of transcripts $self->{'_transcripts'}
    
=cut  

sub _flush_Transcripts {
  my ($self) = @_;
  $self->{'_transcripts'} = [];
  return;
}

############################################################

=head2 _process_Transcripts

    Title   :   _process_Transcripts
    Usage   :   @new_transcripts= $self->_cluster_Transcripts(@read_transcripts)
    Function:   main magic and witchcraft on the transcripts. 
                It checks, clusters and  merges an input array of transcripts
    Returns :   @Bio::EnsEMBL::Transcript
    Args    :   @Bio::EnsEMBL::Transcript

=cut

sub _process_Transcripts {
  my ($self, $alltranscripts, $strand) = @_;

  print STDERR "EST_GeneBuilder: processing input transcripts...\n";

  # first check transcripts and hold info about est_evidence, etc...
  my @transcripts = $self->_check_Transcripts($alltranscripts,$strand);
  print STDERR scalar(@transcripts)." transcripts returned from _check_Transcripts\n";

  if ( scalar(@transcripts) == 0 ){
    print STDERR "No transcripts created, process stopped\n";
    return;
  }
 
  # cluster the transcripts according to exon-overlap
  my @transcript_clusters = $self->_cluster_Transcripts(\@transcripts);
  print STDERR scalar(@transcript_clusters)." clusters returned from _cluster_Transcripts\n";

  # merge the transcripts in each cluster according to consecutive exon overlap
  my @merged_transcripts  = $self->_merge_Transcripts(\@transcript_clusters,$strand);
  print STDERR scalar(@merged_transcripts)." transcripts returned from _merge_Transcripts\n";

  # check the splice_sites
  my @resulting_transcripts = $self->_check_splice_Sites(\@merged_transcripts,$strand);
  print STDERR scalar(@resulting_transcripts)." transcripts returned from _check_splice_Sites\n";

  return @resulting_transcripts;
}

############################################################

=head2 _check_Transcripts

    Title   :   _check_Transcripts
    Usage   :   @transcripts = $self->_check_Transcripts(@transcripts)
    Function:   checks transcripts obtained from EST2Genome for consistency among exons
                in strand, hit_name (hid), exon content,
                and also checks that the hits associated to consecutive exons do not have a
                discontinuity (in hit-coordinates) larger than than a certain limit
    Returns :   @Bio::EnsEMBL::Transcript
    Args    :   @Bio::EnsEMBL::Transcript, ref to hash for linking hid to transcript, ref to hash for linking 
                exon to hid, 

=cut

sub _check_Transcripts {
  my ($self, $ref_transcripts, $strand) = @_;

  # the source_tag of the supporting evidence is specified in Bio/EnsEMBL/Pipeline/EST_conf.pl
  my $evidence_tag = $::evidence_conf{'evidence_tag'};
  
  # the minimum allowed perc. identity of the evidence per exon 
  my $min_similarity  = $::evidence_conf{'min_evidence_similarity'};

  # the maximum allowed discontinuity in EST hits
  my $max_est_gap  = $::evidence_conf{'max_evidence_discontinuity'};

  print STDERR "EST_GeneBuilder: checking consistency of transcripts...\n";

  my @allexons;       # here we'll put all exons that pass the check
  my @alltranscripts; # here we'll put all the transcripts that pass the check
  my %hid_trans;
  my $exon_adaptor    = $self->dbobj->get_ExonAdaptor;
  my $feature_adaptor = $self->dbobj->get_FeatureAdaptor;
  my $total_rejected        = 0;
  
 TRANSCRIPT: 
  foreach my $transcript (@$ref_transcripts){
    my @exons = $transcript->get_all_Exons;
    #print STDERR "Transcript with ".scalar(@exons)." exons\n";
    my $hid;
    my $this_strand;
    my @accepted_exons; # here we hold the good exons in this transcript
    my $rejected = 0;
    my $exon_count = 0;
    my $seqname;
    
  EXON:
    foreach my $exon (@exons){
      my $hstart;
      my $hend;
      #print STDERR " --- Exon $exon_count ---\n";

      ## check strand consistency
      if(!defined $this_strand) { 
	$this_strand = $exon->strand; 
      }
      if($this_strand ne $exon->strand){
	$self->warn("strand not consistent among exons for " . $transcript->temporary_id . " - skipping it\n");
	next TRANSCRIPT;
      }
      
      # check contig consistency
      unless ($seqname){
	$seqname = $exon->seqname;
      }
      if ( !( $exon->seqname eq $seqname ) ){
	print STDERR "transcript ".$transcript->stable_id." (".$transcript->dbID.") is partly".
		    " outside the contig, skipping it...\n";
	next TRANSCRIPT;
      }

      # get the supporting_evidence for each exon
      $exon_adaptor->fetch_evidence_by_Exon($exon);
      my @nonsorted_sf = $exon->each_Supporting_Feature;      
      my @sf = sort { $a->hstart <=> $b->hstart } @nonsorted_sf;
            
      # check that you get suporting_evidence at all
      if ( scalar( @sf ) == 0 ){
	$self->warn("exon $exon with no supporting evidence, possible sticky exon, ".
		    "exon_id=".$exon->dbID."\n");
	next EXON;
      }

#      ## print out all the info ##
#      print STDERR "exon_id: ".$exon->dbID."\n";
#      print STDERR "contig: ".$exon->contig->internal_id." start: ".$exon->start." end: ".$exon->end.
#                   " exon length: ".$exon->length."\n";
#      print STDERR "supporting evidence:\n";
#      foreach my $f (@sf){
#	print STDERR $f->hseqname." score: ".$f->score." hstart: ".$f->hstart." hend: ".$f->hend.
#	  " seq_start: ".$f->start." seq_end: ".$f->end."\n";
#      }
            
####### check consitency of supporting evidence
      my $total_length;
      my $numerator;
      my $total_evidence_length;
      my $gap  = 0;
      my $total_hgap   = 0;
      my $exon_length  = $exon->length;
      my $feature_count = 0;
    
    EVIDENCE:
      foreach my $feature ( @sf ) {
	$feature_count++;

	# because we get all the supporting features indiscriminately...
	next EVIDENCE unless $feature->source_tag eq $evidence_tag;
	
	# this might sound weird, but sometimes features are repeated, don't know yet why
	if ( $feature_count != $#sf+1 ){
	  if ( $sf[$feature_count]->hstart == $feature->hstart &&
	       $sf[$feature_count]->hend   == $feature->hend ){
	    #print STDERR "feature is repeated!, skipping it\n";
	    next EVIDENCE;
	  }
	}
	# check that all exons have the same feature
	if(!defined $hid) { 
	  $hid = $feature->hseqname; 
	}
	if($hid ne $feature->hseqname){
	  $self->warn("hid not consistent among exons for " . $transcript->temporary_id . " - skipping it\n");
	  next TRANSCRIPT;
	}
	
######### calculate percentage identity with respect to the alignment length
	#
	#  ESTs -->   ACGTACGT-ACGT --> 3 EST pieces
	#  exon -->   ACG-ACGTGACGT
	#             ^^^ ^^^^ ^^^^--> evidence_length
	#  alignment length = evidence_length + gaps
	#
	
	my $length              = $feature->end - $feature->start + 1;
	my $score               = $feature->score;
	my $hgap                = 0;
	$numerator             += $length * $score;
	$total_evidence_length += $length;

	if ( $feature_count != $#sf+1 ){
	  if ( $sf[$feature_count]->hstart > $feature->hend ){
	    $hgap = $sf[$feature_count]->hstart - $feature->hend - 1;
	  }
	  elsif 
	    ( $sf[$feature_count]->hstart <= $feature->hend && !($sf[$feature_count]->overlaps($feature))){
	    $hgap = $feature->hstart - $sf[$feature_count]->hend-1;
	  }
	  else{
	    $hgap = 0;
	  }
	  if ( $sf[$feature_count]->overlaps( $feature )){
	    print STDERR "features are overlapping!\n";
	  }
	  $total_hgap   += $hgap;
	}
	      
      } # end of EVIDENCE 

###### calculate alignment length and the similarity score (exon_score)
	
      $gap           = $exon_length - $total_evidence_length;
      $total_length  = $total_evidence_length + $gap + $total_hgap;
      my $exon_score = $numerator / $total_length;
      
      # round it up
      my $number   = 100 * $exon_score;
      my $decimals = $number % 100;
      my $addition = int( $decimals / 50 );  # gives 1 ( or 0 ) if $decimals > ( or < ) 50
      $exon_score  = int( $exon_score ) + $addition;
      
      # reject EXONS with perc_identity below $min_similarity = $::evidence_conf{'min_evidence_similarity'}
      if ($min_similarity > $exon_score ){
	print STDERR "rejecting Exon due to low percentage identity = $exon_score\n";
	$rejected++;
	$total_rejected++;

	# if we have rejected all exons in this transcript, forget about the transcript
	if ( $rejected == scalar(@exons) ){
	  print STDERR "all exons rejected in this transcript, rejecting transcript\n";
	  next TRANSCRIPT;
	}
	
	next EXON;
      }
      else{
	push (@accepted_exons, $exon);
      }
 
####### alternatively we can get the percentage identity from the feature table

#      # get the percentage identity for the est_hit on this exon from the feature table
#      my $exon_contig = $exon->contig;
      
#      # get the features by hid
#      my @hid_features = $feature_adaptor->fetch_by_hid($hid); 

#      # keep only the features corresponding to the current exon_contig
#      my @selected_features;
#      foreach my $f ( @hid_features ){
	
#	# feature->start/end is returned in vc coordinates
#	  print STDERR "from feature percentage identity: ".$f->percent_id.
#	    " seq_start: ".$f->start." seq_end: ".$f->end."\n";
#	if ($f->start == $exon->start && $f->end == $exon->end ){
#	#  print STDERR "from feature percentage identity: ".$f->percent_id.
#	#    " seq_start: ".$f->start." seq_end: ".$f->end."\n";
#	  print STDERR "chosen!\n";
#	}
#      }

      
####### check the gap with the evidence of the next exon
      # if the ESTs are of good quality, this should not reject any
      if ( $exon_count != $#exons ){
	my $est_gap = 0;
	my $exon2 = $exons[$exon_count+1];
	
	$exon_adaptor->fetch_evidence_by_Exon( $exon2 );
	my @nonsorted_sf2 = $exon2->each_Supporting_Feature; 
	my @sf2 = sort { $a->hstart <=> $b->hend } @nonsorted_sf2;

	if ( scalar( @sf2 ) != 0 ){

	  # if the hstart increases per exon, the EST runs in the same direction of the gene 
	  if ( $sf[0]->hstart < $sf2[0]->hstart ){
	    $est_gap = abs( $sf2[0]->hstart  - $sf[$#sf]->hend );
	  }
	  # if hstart decreases that means that the EST runs in the opposite direction
	  elsif (  $sf[0]->hstart > $sf2[0]->hstart ){
	    $est_gap = abs( $sf[0]->hstart - $sf2[$#sf2]->hend);
	  }
	  # else, same EST piece is hitting two exons, not good!
	  else{
	    $self->warn( "same bit of evidence is hitting two exons!\n");
	  }
	  
	  # skip EXONS that have too large gaps in the EST hit
	  if ( $est_gap > $max_est_gap ){
	    print STDERR "EST with too large gaps skipping it\n";
	    next TRANSCRIPT;
	  }
	}
      }
      
    $exon_count++;
    
    }  # end of EXON

    # put the accepted exons into this transcript
    $transcript->flush_Exon;
    foreach my $exon (@accepted_exons){
      $transcript->add_Exon($exon);
    }

    push(@allexons, @accepted_exons);
    
    # finally check the one-to-one correspondance between ESTs and input_transcripts
    if(defined $hid && defined $hid_trans{$hid}) { 
      $self->warn("$hid is being used by more than one transcript!\n"); 
    }
    # we hold which transcript is associated with each hit
    $hid_trans{$hid} = $transcript;

    # if the transcript made it to this point, keep it
    push (@alltranscripts, $transcript);
  }
  print STDERR $total_rejected." exons rejected due to low similarity score\n";

  return @alltranscripts;
}

############################################################

=head2 _cluster_Transcripts

    Title   :   _cluster_Transcripts
    Usage   :   @clusters = $self->_cluster_Transcripts(\@transcripts)
    Function:   it clusters transcripts, it is very fast but it may miss some transcripts because
                it only compares with the current transcript. It could be improved by comparing also
                with one or two previous clusters
    
=cut 
  
sub _cluster_Transcripts{
  my ($self,$ref_transcripts) = @_;
  my @transcripts = @{ $ref_transcripts };
  my @clusters;
  print STDERR "EST_GeneBuilder: clustering transcripts...\n";
			 
  # first sort the transcripts by their start position coordinate
  my %start_table;
  my $i=0;
  foreach my $transcript (@transcripts){
    $start_table{$i} = $transcript->start_exon->start;
    $i++;
  }
  my @sorted_transcripts=();
  foreach my $pos ( sort { $start_table{$a} <=> $start_table{$b} } keys %start_table ){
    push (@sorted_transcripts, $transcripts[$pos]);
  }
  @transcripts = @sorted_transcripts;
  foreach my $tran (@transcripts){
    print STDERR "$tran, start: ".$tran->start_exon->start."\n";
  }
  # create a new cluster 
  my $cluster = Bio::EnsEMBL::Utils::TranscriptCluster->new();
  my $cluster_count = 1;

  # put the first transcript into these cluster
  $cluster->put_Transcripts( $sorted_transcripts[0] );
  push( @clusters, $cluster );
    
  # loop over the rest of the genes
 LOOP1:
  for (my $c=1; $c<=$#sorted_transcripts; $c++){
    my $found=0;

    # compare with the transcripts in this cluster
  LOOP2:
    foreach my $t_in_cluster ( $cluster->get_Transcripts ){       
      if ( $self->_compare_Transcripts( $sorted_transcripts[$c], $t_in_cluster ) ){	
	$cluster->put_Transcripts( $sorted_transcripts[$c] );                       
	$found=1;
	next LOOP1;
      }
    }
    # if not in this cluster compare to the previous clusters:

    # to restrict this to the ($limit) previous clusters
    # set my $limit = 6; (for example) and include in the while the following condition
    # while ( !(...)  && !($lookup > $limit) )

    if ( $found == 0 && $cluster_count > 1 ) {
      my $lookup = 1;
      while ( !($cluster_count <= $lookup ) ){ 
	print STDERR "cluster_count: $cluster_count, looking at ".($cluster_count - $lookup)."\n";
	my $previous_cluster = $clusters[ $cluster_count - 1 - $lookup ];
	foreach my $t_in_cluster ( $previous_cluster->get_Transcripts ){
	  if ( $self->_compare_Transcripts( $sorted_transcripts[$c], $t_in_cluster ) ){	
	    $previous_cluster->put_Transcripts( $sorted_transcripts[$c] );                       
	    $found=1;
	    next LOOP1;
	  }
	}
	$lookup++;
      }
    }
    # if not-clustered create a new TranscriptCluster
    if ( $found == 0 ){  
      $cluster = new Bio::EnsEMBL::Utils::TranscriptCluster; 
      $cluster->put_Transcripts( $sorted_transcripts[$c] );
      push( @clusters, $cluster );
      $cluster_count++;
    }
  }
  
  # print out the clusters
  my $number  = 1;
  foreach my $cluster (@clusters){
    my $count = 1;
    print STDERR "cluster $number :\n";
    foreach my $tran ($cluster->get_Transcripts){
      print STDERR "$count:\n";
      foreach my $exon ( $tran->get_all_Exons ){
	print STDERR $exon->start.":".$exon->end." ";
      }
      print STDERR "\n";
      $count++;
    }
    $number++;
  }		
  
  return @clusters;
}

############################################################

=head2 _compare_Transcripts()

 Title: _compare_Transcripts()
 Usage: compares the exons of two transcripts according to overlap and returns the number of overlaps

=cut

sub _compare_Transcripts {        
  my ($self,$transcript1,$transcript2) = @_;
  my @exons1   = $transcript1->get_all_Exons;
  my @exons2   = $transcript2->get_all_Exons;
  my $overlaps = 0;
  
  foreach my $exon1 (@exons1){
    foreach my $exon2 (@exons2){
      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	return 1;
      }
    }
  }
  return 0;
}    

#########################################################################

=head2

 Title   : _merge_Transcripts
 Function: reads all the est2genome transcripts that have been clustered and merges those that
           are entirely embedded in each other, producing brand new transcripts from that,
           see below the description of the algorithm
=cut

sub _merge_Transcripts{
  my ($self,$ref_transcript_clusters,$strand) = @_;

  print STDERR "EST_GeneBuilder: merging transcripts...\n";

  my @total_merged_transcripts;
  
  # look in each cluster
  my $count =0;
 CLUSTER:
  foreach my $cluster ( @{ $ref_transcript_clusters } ){
    
    $count++;
    my @merged_transcripts = ();

    # get the transcripts in this cluster
    my @transcripts = $cluster->get_Transcripts;
    
    # sort the transcripts by the number of exons
    @transcripts = sort { scalar( $b->get_all_Exons ) <=> scalar( $a->get_all_Exons ) } @transcripts;
    
    # test
    #print STDERR "\nNew Cluster:\n";
    #foreach my $tran (@transcripts){
    #  print STDERR "transcript: $tran\n";
    #  foreach my $exon ( $tran->get_all_Exons ){
    #	print STDERR $exon->start.":".$exon->end."  ";
    #  }
    #  print STDERR "\n";
    #}
    
    # now we loop over the transcripts, 
    my %is_merged;
    
    # the algorithm goes roughly as follows
    #
    # FOREACH transcript
    #   loop over the rest
    #     if it merges, add it to the current 'cluster'
    #     loop over the rest and only add to the cluster if: at least merges to one in the cluster and 
    #         with those that doesn't merge it should not overlap at all (to avoid merging things that shouldn't)
    # 
    #     allow for transcripts to be used twice
    #     once we've gone through the whole lot, merge the chosen transcripts and put the resulting transcript
    #     in an array.
    #
    # proceed as above with the next transcript in the list, but once the merged transcript is produced,
    # before accepting, check that it does not merge with transcripts previously produced, if it does, reject 
    # (since by construction it is necessarily embedded in that one)
    #
    # iterate until you've gone through all transcripts
    
    # for each transcript
  TRAN1:
    for (my $u=0; $u<scalar(@transcripts); $u++){
      my @current_list = ();
      push (@current_list, $transcripts[$u]);

      # go over the rest
    TRAN2:
      for (my $v=$u+1; $v<scalar(@transcripts); $v++){
	
	my $overlap_ifnot_merged   = 0;
	my $merge_to_current       = 0;
	
	for (my $w=0; $w<scalar(@current_list); $w++){
	  
	  # in order to merge a transcript...
	  #print STDERR "comparing $current_list[$w] ($w) and $transcripts[$v] ($v)\n";
	  my ($merge,$overlaps) = $self->_test_for_Merge( $current_list[$w], $transcripts[$v] );
	  #print STDERR "merge = $merge, overlaps = $overlaps\n";
	  
	  # ...there must be at least one in @current_list to which $transcripts[$v] merges
	  if ( 1 == $merge ){
	    $merge_to_current = 1;
	  }

	  # ...and it should not overlap with those to which it does not merge
	  unless (1 == $merge){
	    $overlap_ifnot_merged += $overlaps;
	  }
	}
	
	# ... and then it can merge to the list in @current_list
	if ( 1 == $merge_to_current && 0 == $overlap_ifnot_merged ){
	  push (@current_list, $transcripts[$v]);
	}
	
      } # end of TRAN2
      
      # create a new transcript with those merged above
      my $new_transcript = $self->_produce_Transcript( \@current_list, $strand );
      
      my $found  = 0;

      # for subsequent transcritps iterate through the merged_transcripts to see if it merges to any previous one
      if ( scalar(@merged_transcripts) >= 1){
      MERGED:
	foreach my $merged_tran (@merged_transcripts){
	  my ($merge,$overlaps) = $self->_test_for_Merge( $new_transcript, $merged_tran );
	  if ( 1 == $merge ){
	    $found = 1;
	    last MERGED;
	  }
	}
      }
      # if it doesn't merge, then add to @merged_transcripts
      # this will add as well the first new_transcript created
      if ($found == 0 ){
	push( @merged_transcripts, $new_transcript); 
      }
      else{
	# we don't add it
      }
    }   # end of TRAN1
    
    # test
    print STDERR "Resulting merged transcripts:".scalar(@merged_transcripts)."\n";
    foreach my $tran (@merged_transcripts){
      print STDERR "transcript: $tran\n";
      foreach my $exon ($tran->get_all_Exons){
    	print STDERR $exon->start.":".$exon->end."  ";
      }
      print STDERR "\n";
    }
    
  push (@total_merged_transcripts, @merged_transcripts);
  }     # end of CLUSTER		       
 
return @total_merged_transcripts;
		      
}

#########################################################################

=head2

 Title   : _test_for_Merge
 Function: this function is called from _merge_Transcripts and actually checks whether two transcripts
           inputs merge.
 Returns : It returns two numbers ($merge,$overlaps), where
           $merge = 1 (0) when they do (do not) merge,
           and $overlaps is the number of exon-overlaps.

=cut

sub _test_for_Merge{
  my ($self,$tran1,$tran2) = @_;
  my @exons1 = $tran1->get_all_Exons;
  my @exons2 = $tran2->get_all_Exons;	
 
  my $foundlink = 0; # flag that gets set when starting to link exons
  my $start     = 0; # start looking at the first one
  my $overlaps  = 0; # independently if they merge or not, we compute the number of exon overlaps
  my $merge     = 0; # =1 if they merge

EXON1:
  for (my $j=0; $j<=$#exons1; $j++) {
  
  EXON2:
    for (my $k=$start; $k<=$#exons2; $k++){
    #print STDERR "comparing ".($j+1)." and ".($k+1)."\n";
	    
      # if exon 1 is not the first, check first whether it matches the previous exon2 as well, i.e.
      #                        ____     ____        
      #              exons1 --|____|---|____|------ etc... $j
      #                        ____________  
      #              exons2 --|____________|------ etc...  $k
      #
      if ($foundlink == 1 && $j != 0){
	if ( $k!= 0 && $exons1[$j]->overlaps($exons2[$k-1]) ){
	  #print STDERR ($j+1)." <--> ".($k)."\n";
	  $overlaps++;
          next EXON1;
	}
      }
      
      # if texons1[$j] and exons2[$k] overlap go to the next exon1 and  next $exon2
      if ( $exons1[$j]->overlaps($exons2[$k]) ){
	#print STDERR ($j+1)." <--> ".($k+1)."\n";
        $overlaps++;
	$foundlink = 1;
      }          
      else {  
	# look at the next exon if you haven't found an overlap yet
	if ( $foundlink == 0 ){
	  next EXON2;
	}
	# leave if we stop finding links between exons before the end of transcripts
	if ( $foundlink == 1 ){
	  $merge = 0;
	  last EXON1;
	}
      }
      
      # if foundlink = 1 and we get to the end of either transcript, we merge them!
      if ( $foundlink == 1 && ( $j == $#exons1 || $k == $#exons2 ) ){
	
	### some prints for test ###
	#print STDERR "\nmerging transcripts $tran1 and $tran2:\n";
	#foreach my $e1 ( @exons1 ){
	#  print STDERR $e1->start.":".$e1->end."\t";
	#}
	#print STDERR "\n";
	#foreach my $e2 ( @exons2 ){
	#  print STDERR $e2->start.":".$e2->end."\t";
	#}
	#print STDERR "\n\n";
	
	# and we can leave
        $merge = 1;
	last EXON2;
      }
      # if foundlink = 1 but we're not yet at the end, go to the next exon 
      if ( $foundlink == 1 ){
	
	# but first check whether in exons2 there are further exons overlapping exon1, i.e.
        #                       ____________        
	#             exons1 --|____________|------ etc...
	#                       ____     ___  
	#             exons2 --|____|---|___|------ etc...
	# 
	my $addition = 0;
	while ( $k+1+$addition < scalar(@exons2) && $exons1[$j]->overlaps($exons2[$k+1+$addition]) ){
	  #print STDERR ($j+1)." <--> ".($k+2+$addition)."\n";
	  $overlaps++;
          $addition++;
	}      
	$start = $k+1+$addition;
	next EXON1;
      }    
    } # end of EXON2 
  }   # end of EXON1      

  # if we haven't returned at this point, they don't merge, thus
  return ($merge,$overlaps);
}
  


############################################################

=head2

 Title   : _produce_Transcript
 Function: reads all the est2genome transcript that can be merged and make a single transcript
           out of them
=cut

sub _produce_Transcript{
  my ($self,$merged,$strand) = @_;

  my @allexons;
  my %exon2transcript;			
  my %is_first;
  my %is_last;			
	       
  # collect all exons
  foreach my $tran (@{ $merged }){
    my @exons = $tran->get_all_Exons;
    @exons    = sort { $a <=> $b } @exons;
    
    push ( @allexons, @exons );
    
    # keep track of whether the exons is first or last and the transcript it belongs to
    for (my $i = 0; $i< scalar( @exons ); $i++){
      if ( 0 == $i ){
	$is_first{$exons[$i]} = 1;
      }
      else{
	$is_first{$exons[$i]} = 0;
      }
      if ( $#exons == $i ){
	$is_last{$exons[$i]} = 1;
      }
      else{
	$is_last{$exons[$i]} = 0;
      }
      $exon2transcript{$exons[$i]} = $tran;
    }
  }

  # cluster them
  my $cluster_list = $self->_cluster_Exons( @allexons );

  # set start and end in the exons
  $cluster_list = $self->_set_splice_Ends($cluster_list,\%exon2transcript,\%is_first,\%is_last,$strand);

  # create a new transcript with these exons
  my $transcript    = Bio::EnsEMBL::Transcript->new();
  my @exon_clusters = $cluster_list->sub_SeqFeature;
  
  foreach my $exon_cluster (@exon_clusters){
    my $new_exon = Bio::EnsEMBL::Exon->new();
    $new_exon->start ($exon_cluster->start);
    $new_exon->end   ($exon_cluster->end);
    
    ###  dont't set strand yet, genomewise cannot handle that ###
    
    foreach my $exon ( $exon_cluster->sub_SeqFeature ){
      foreach my $evidence ( $exon->each_Supporting_Feature ){
	$new_exon->add_Supporting_Feature($evidence);
      }
    }
    $transcript->add_Exon($new_exon);
  }
  
  return $transcript;
}

############################################################

=head2

 Title   : _cluster_Exons
 Function: it cluster exons according to exon overlap

=cut

sub _cluster_Exons{
  my ($self, @exons) = @_;

  #print STDERR "EST_GeneBuilder: clustering exons...\n";
  
  # no point if there are no exons!
  return unless ( scalar( @exons) > 0 );   

  # keep track about in which cluster is each exon
  my %exon2cluster;
  
  # main cluster feature - holds all clusters
  my $cluster_list = new Bio::EnsEMBL::SeqFeature; 
  
  # sort exons by start coordinate
  @exons = sort { $a->start <=> $b->start } @exons;

  # Create the first exon_cluster
  my $exon_cluster = new Bio::EnsEMBL::SeqFeature;
  
  # Start off the cluster with the first exon
  $exon_cluster->add_sub_SeqFeature($exons[0],'EXPAND');
  $exon_cluster->strand($exons[0]->strand);    
  $cluster_list->add_sub_SeqFeature($exon_cluster,'EXPAND');
  
  # Loop over the rest of the exons
  my $count = 0;

 EXON:
  foreach my $exon (@exons) {
    if ($count > 0) {
      my @overlap = $self->match($exon, $exon_cluster);    
      
      # Add to cluster if overlap AND if strand matches
      if ( $overlap[0] && ( $exon->strand == $exon_cluster->strand) ) { 
	$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
      }  
      else {
	# Start a new cluster
	$exon_cluster = new Bio::EnsEMBL::SeqFeature;
	$exon_cluster->add_sub_SeqFeature($exon,'EXPAND');
	$exon_cluster->strand($exon->strand);
		
	# and add it to the main_cluster feature
	$cluster_list->add_sub_SeqFeature($exon_cluster,'EXPAND');	
      }
    }
    $count++;
  }
  return $cluster_list;
}

############################################################


=head2 _set_splice_Ends

    Title   :   _set_splice_Ends
    Usage   :   $cluster_list = $self->_set_splice_Ends($cluster_list)
    Function:   Resets the ends of the clusters to the most frequent coordinate
   
=cut

sub _set_splice_Ends {
  my ($self, $cluster_list, $ref_exon2transcript_hash, $ref_is_first, $ref_is_last, $strand) = @_;

  # hash having exons as keys and mother-transcript as value
  my %exon2transcript = %$ref_exon2transcript_hash;

  # keep track of whether the exon is first or last in the transcript
  my %is_first = %$ref_is_first;
  my %is_last  = %$ref_is_last;

  #print STDERR "EST_GeneBuilder: setting common ends...\n";

  # get the exon clusters
  my @exon_clusters = $cluster_list->sub_SeqFeature;

  # sort clusters according to their start coord.
  @exon_clusters = sort { $a->start <=> $b->start  } @exon_clusters;
  my $count =  0;
  
  #### set first the start-coordinates ####
  my $position = 0;
 CLUSTER:		     
  foreach my $cluster (@exon_clusters) {
    $position++;
    my %start;
    my $new_start;
    my $max_start = 0;
    
  EXON:
    foreach my $exon ($cluster->sub_SeqFeature){

      # for a start-coord in the middle don't use the exon if it is the first of a transcript
      if ( $position > 1 && $position <= scalar( @exon_clusters ) ){
	if ( $is_first{ $exon } == 1 ){
	  next EXON;
	}
      }
      $start{$exon->start}++;
    }
    # take the most common start (note that we do not resolve ambiguities here)
    while ( my ($key,$value) = each %start ){
      if ($value > $max_start){
	$new_start = $key;
	$max_start = $value;
      }
    }
    # if we have too little exons to obtain the start, take the original value
    if ( $max_start == 0 ){
      $new_start = $cluster->start;
    }

    # the first cluster is a special case - potential UTRs, take the longest one.
    if( $position == 1) {
      $new_start = $cluster->start;  
    }   
    # reset the starts of all the exons in this cluster
    foreach my $exon ($cluster->sub_SeqFeature) {
      $exon->start($new_start);
    } 
  }

  #### now set the end coordinate ####
  $position = 0;
 CLUSTER:		     
  foreach my $cluster (@exon_clusters) {
    $position++;
    my %end;
    my $new_end;
    my $max_end = 0;
    
  EXON:
    foreach my $exon ($cluster->sub_SeqFeature){

      # for an end-coord in the middle don't use the exon if it is the last of a transcript
      if ( $position >= 1 && $position < scalar( @exon_clusters ) ){
	if ( $is_last{ $exon } == 1 ){
	  next EXON;
	}
      }
      $end{$exon->end}++;
    }
    # take the most common end (note that we do not resolve ambiguities here)
    while ( my ($key,$value) = each %end ){
      if ($value > $max_end){
	$new_end = $key;
	$max_end = $value;
      }
    }
    # if we have too little exons to obtain the end, take the original value
    if ( $max_start == 0 ){
      $new_end = $cluster->end;
    }

    # the last cluster is a special case - potential UTRs, take the longest one.
    if( $position == 1) {
      $new_end = $cluster->end;  
    }
    # reset the ends of all the exons in this cluster
    foreach my $exon ($cluster->sub_SeqFeature) {
      $exon->end($new_end);
    } 
  }

  return $cluster_list;
}

############################################################

 
sub _check_splice_Sites{
  my ($self, $ref_transcripts, $strand) = @_;

  print STDERR "EST_GeneBuilder: checking splice sites in strand $strand...\n";

  # get the contig being analysed
  my $contig = $self->vc;
  
  # for reverse strand,  invert the contig, since the exons were retrieved in revcontig
  my $revcontig = $contig->invert;   
  
  my ($upstream_correct, $downstream_correct) = (0,0);
  my $site_count = 0;
  
 TRANSCRIPT:
  foreach my $transcript (@$ref_transcripts){
    
    # count the number of exons
    my $count = 0; 
    my @exons = $transcript->get_all_Exons;
    
  EXON:
    foreach my $exon ( @exons ){

      # take the 2 bases right before the exon and after the exon
      my ($upstream, $downstream);
        
      # forward strand
      if ($strand == 1){
	my $seq = $contig->primary_seq; # a Bio::PrimarySeq object
	
	# catch possible exceptions in gettting the sequence (it might not be there!)
	eval{
	  $upstream   = $seq->subseq( ($exon->start)-2 , ($exon->start)-1 ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->start)-2).",".(($exon->start)-1).")\n";
	  print STDERR $@;
	  $upstream = 'NN';
	}
	eval{
	  $downstream = $seq->subseq( ($exon->end)+1   , ($exon->end)+2   ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->end)+1).",".(($exon->end)+2).")\n";
	  print STDERR $@;
	  $downstream = 'NN';
	}
	# print-outs to test it
	#      if ( $count ==0 ){
	#	print STDERR "FIRST EXON-->".$downstream;
	#      }
	#      if ( $count != 0 && $count != $#clusters ){
	#	print STDERR $upstream."-->EXON-->".$downstream;
	#      }
	#      if ( $count == $#clusters ){
	#	print STDERR $upstream."-->LAST EXON";
	#      }
	#      print "\n";
	
	# the first and last exon are not checked - potential UTR's
	if ( $count != 0 && $upstream eq 'AG') {       
	  $upstream_correct++;
	}
	if ( $count !=$#exons && $downstream eq 'GT') { 
	  $downstream_correct++;
	}
	$count++;
      } # end of forward strand
      
      # reverse strand
      if ($strand == -1 ){
	my $seq = $revcontig->primary_seq; # a Bio::PrimarySeq object
	
	# catch possible exceptions in gettting the sequence (it might not be there!)
	eval{
	  $upstream = $seq->subseq( ($exon->start)-2 , ($exon->start)-1 ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->start)-2).",".(($exon->start)-1).")\n";
	  print STDERR $@;
	  $upstream ='NN';
	}
	eval{
	  $downstream   = $seq->subseq( ($exon->end)+1   , ($exon->end)+2   ); 
	};
	if ($@){
	  print STDERR "Unable to get subsequence (".(($exon->end)+1).",".(($exon->end)+2).")\n";
	  print STDERR $@;
	  $downstream = 'NN';
	}
	#  in the reverse strand we're looking at coordinates in the reversed-complement contig:
	#
	#        $contig : --------------------TG----GA---------->  forward strand
	#                                      AC    CT             reverse strand 
	#                           downstream   EXON   upstream
	#
	#
	#                   upstream   EXON   downstream              
	#     $revcontig : ----------AG----GT-------------------->    forward strand
	#                                                             reverse strand
	#
	# and take the reverse complement
	( $downstream = reverse( $downstream) ) =~ tr/ACGTacgt/TGCAtgca/;
	( $upstream   = reverse( $upstream )  ) =~ tr/ACGTacgt/TGCAtgca/;
	
	# so the conserved sequence should be AC<-EXON<-CT printed as in the reverse strand
	
	#      if ( $count ==0 ){
	#	print STDERR "LAST EXON<--".$upstream;
	#      }
	#      if ( $count != 0 && $count != $#clusters ){
	#	print STDERR $downstream."<--EXON<--".$upstream;
	#      }
	#      if ( $count == $#clusters ){
	#	print STDERR $downstream."<--FIRST EXON";
	#      }
	#      print "\n";
	
	# the first and last exon are not checked - potential UTR's
	if ( $count != $#exons && $downstream eq 'AC') {       
	  $upstream_correct++;
	}
	if ( $count != 0 && $upstream eq 'CT') { 
	  $downstream_correct++;
	}
	$count++;
      } # end of reverse strand
            
    }   # end of EXON
    
    $site_count += $#exons; # count the number of splice sites
    
  }   # end of TRANSCRIPT
    
  print STDERR "upstream splice-sites correct: ".$upstream_correct. 
    " out of ".($site_count)." splice-sites\n";
  print STDERR "downstream splice-sites correct: ".$downstream_correct. 
    " out of ".($site_count)." splice-sites\n";
  
  # eventually, we may modify start/end of exons according to these checks
  return @$ref_transcripts;    
}

############################################################

=head2 _put_Transcript

 Title   : _put_Transcript
 Usage   : $self->add_Transcript
 Function: method to add transcripts into the array $self->{'_transcripts'} 
 Returns : nothing
 Args    : Bio::EnsEMBL::Transcript

=cut

sub _put_Transcript {
  my ($self,$transcript) = @_;
  $self->throw("No transcript input") unless defined($transcript);
  $self->throw("Input must be Bio::EnsEMBL::Transcript") unless $transcript->isa("Bio::EnsEMBL::Transcript");
  if ( !defined( $self->{'_transcripts'} ) ){
    @{ $self->{'_transcripts'} } = ();
  }
  push( @{ $self->{'_transcripts'} }, $transcript );
}

############################################################

=head2 _get_all_Transcripts

 Title   : _get_all_Transcripts
 Usage   : my @transcripts = $self->_get_all_Transcripts;
 Function: method to get all the transcripts stored in @{ $self->{'_transcripts'} } 
 Example : 
 Returns : Bio::EnsEMBL::Transcript
 Args    : nothing

=cut

sub _get_all_Transcripts {
  my ($self) = @_;
  if ( !defined( $self->{'_transcripts'} ) ) {
    @{ $self->{'_transcripts'} } = ();
    print STDERR "The transcript array you're trying to get is empty\n";
  }
  my @trans = @{ $self->{'_transcripts'} };
  return @trans;
}

############################################################

sub add_runnable{
  my ($self, $value, $strand) = @_;

  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 
  if (defined($value)) {
    
    if ($value->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
      
      if( $strand == -1 ){
	push(@{$self->{'_reverse_runnables'}},$value);
      }
      elsif( $strand == 1){
	push(@{$self->{'_forward_runnables'}},$value);
      }
      else{
	$self->throw( "Cannot add a runnable with strand = $strand" );
      }
    
    } 
    else {
      $self->throw("[$value] is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
  }
}

sub each_runnable{
  my ($self, $strand) = @_;
  
  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 

  if( $strand == -1 ){
    return @{$self->{'_reverse_runnables'}};
  }
  elsif ($strand == 1){
    return @{$self->{'_forward_runnables'}};
  }
  else{
    $self->throw( "there are no runnables with strand = $strand" );
  }
  
}

sub run {
  my ($self) = @_;
  my $strand;

  my $cc=1;
  # run genomewise 

  # genes get written in the database with the type specified in Bio/EnsEMBL/Pipeline/EST_conf.pl
  my $genetype = $::genomewise_conf{'genetype'};
  
  # sort out analysis here or we will get into trouble with duplicate analyses
  my $anaAdaptor = $self->dbobj->get_AnalysisAdaptor;
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
	(-db              => 'NULL',
	 -db_version      => 1,
	 -program         => $genetype,
	 -program_version => 1,
	 -gff_source      => $genetype,
	 -gff_feature     => 'gene',
	 -logic_name      => $genetype,
	 -module          => 'EST_GeneBuilder',
      );
  }

  $self->genetype($genetype);
  $self->analysis($analysis_obj);

  # plus strand
  $strand = 1;
  my $tcount=0;
  foreach my $gw_runnable( $self->each_runnable($strand) ){
    $tcount++;
    $gw_runnable->run;
    
    # convert_output
    $self->convert_output($gw_runnable, $strand);
  }
  print STDERR $tcount." transcripts run in genomewise in the forward strand\n";

  # minus strand
  $strand = -1;
  my $tcount2=0;
  foreach my $gw_runnable( $self->each_runnable($strand)) {

#    print STDERR "In EST_GeneBuilder...reverse_strand...about to run genomewise";
#    my $ttcount2=0;
#    foreach my $t ($gw_runnable->get_all_Transcripts){
#      $ttcount2++;
#      $tcount2++;
#    }
#    print STDERR " with $ttcount2 transcript(s)\n";

    $tcount2++;
    $gw_runnable->run;
    
    # convert_output
    $self->convert_output($gw_runnable, $strand);
  }
  print STDERR $tcount2." transcripts run in genomewise in the reverse strand\n";
}

sub genetype {
  my ($self, $genetype) = @_;

  if(defined $genetype){
    $self->{'_genetype'} = $genetype;
  }

  return $self->{'_genetype'};
}

# override method from RunnableDB.pm
sub analysis {
  my ($self, $analysis) = @_;

  if(defined $analysis){
    $self->throw("$analysis is not a Bio::EnsEMBL::Analysis") unless $analysis->isa("Bio::EnsEMBL::Analysis");
    $self->{'_analysis'} = $analysis;
  }

  return $self->{'_analysis'};
}

# convert genomewise output into genes

sub convert_output {
  my ($self, $gwr, $strand) = @_;

  my @genes = $self->make_genes($gwr, $strand);

  my @remapped = $self->remap_genes(\@genes, $strand);

  # store genes
  $self->output(@remapped);
 
}

############################################################

sub make_genes {

  my ($self, $runnable, $strand) = @_;
  my $genetype = $self->genetype;
  my $analysis_obj = $self->analysis;
  my $count = 0;
  my @genes;

  my $time  = time; chomp($time);
  my $contig = $self->vc;

  # are we working on the reverse strand?
  if( $strand == -1 ){
    $contig = $contig->invert;
  }

  # transcripts come with a translation, that's sorted out in MiniGenomwise and Genomewise
  my @trans = $runnable->output;
  
  foreach my $transcript (@trans) {
  
   ##test
  #  print STDERR "\nIn EST_GeneBuilder.make_genes\n";
  #  print STDERR " Transcript        : ".$transcript."\n";
  #  my $ecount=1;
  #  foreach my $exon ( $transcript->get_all_Exons ){
  #    print STDERR "Exon $ecount: ".$exon->start." ".$exon->end."\n";
  #    $ecount++;
  #  }

#    print STDERR " Translation       : ".$transcript->translation."\n";
#    print STDERR " translation starts: ".$transcript->translation->start."\n";
#    print STDERR " translation ends  : ".$transcript->translation->end."\n";
#    print STDERR " start exon: ".$transcript->translation->start_exon
#      ." starts: ".$transcript->translation->start_exon->start
#      ." ends: ".$transcript->translation->start_exon->end."\n";
#    print STDERR " end  exon : ".$transcript->translation->end_exon
#      ." starts: ".$transcript->translation->end_exon->start
#      ." ends: ".$transcript->translation->end_exon->end."\n";
    $count++;
    my $gene   = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->temporary_id($contig->id . ".$genetype.$count");
    $gene->analysis($analysis_obj);

    # add transcript to gene
    $transcript->temporary_id($contig->id . ".$genetype.$count");
    $gene->add_Transcript($transcript);

    # and store it
    push(@genes,$gene);

    # sort the exons 
    $transcript->sort;
    my $excount = 1;
    my @exons = $transcript->get_all_Exons;

    foreach my $exon(@exons){

      #print STDERR "  Exon ".$excount." : ".$exon."  start: ".$exon->start." end: ".$exon->end."\n";
      $exon->temporary_id($contig->id . ".$genetype.$count.$excount");
      $exon->contig_id($contig->id);
      $exon->attach_seq($contig->primary_seq);
      
      # if strand = -1 we have inverted the contig, thus
      $exon->strand(1);
      $excount++;
      # when the gene gets stored, the strand is flipped automatically
    }
    #put temporary_id (this hasn't been put in neither Genomewise nor MiniGenomewise
    my $translation = $transcript->translation;
    $translation->temporary_id($contig->id . ".$genetype.$count");
  }
  return @genes;
}

=head2 remap_genes

    Title   :   remap_genes
    Usage   :   $self->remap_genes(@genes)
    Function:   Remaps predicted genes into genomic coordinates
    Returns :   array of Bio::EnsEMBL::Gene
    Args    :   Bio::EnsEMBL::Virtual::Contig, array of Bio::EnsEMBL::Gene

=cut

sub remap_genes {
  my ($self, $genes, $strand) = @_;
  my $contig = $self->vc;
  if ( $strand == -1 ){
    $contig = $contig->invert;
  }

  #print STDERR "genes before remap: " . scalar(@$genes) . "\n";
  
  my @newf;
  my $trancount=1;
  foreach my $gene (@$genes) {
    #test2
    my @trans = $gene->each_Transcript;

# test
#    
#foreach my $tran (@trans){
#      print STDERR "In EST_GeneBuilder.remap_genes: \n";
#      print STDERR "  transcript        : ".$tran."\n";
#      print STDERR "  strand            : ".$tran->start_exon->strand."\n";
#      print STDERR "  translation       : ".$tran->translation."\n";
#      print STDERR "  translation starts: ".$tran->translation->start."\n";
#      print STDERR "  translation ends  : ".$tran->translation->end."\n";
#      print STDERR "  start exon        : ".$tran->translation->start_exon
#                  ."\tstart: ".$tran->translation->start_exon->start."\tends: "
#		  .$tran->translation->start_exon->end."\n";
#      print STDERR "  end  exon         : ".$tran->translation->end_exon
#	          ."\tstart: ".$tran->translation->end_exon->start."\tends: "
#		  .$tran->translation->end_exon->end."\n";
#      print STDERR "  temporary_id      : ".$tran->translation->temporary_id."\n";
#   }


    
    eval {
      my $newgene = $contig->convert_Gene_to_raw_contig($gene);
      # need to explicitly add back genetype and analysis.
      $newgene->type($gene->type);
      $newgene->analysis($gene->analysis);

      push(@newf,$newgene);
      
    };
    if ($@) {
      print STDERR "Couldn't reverse map gene " . $gene->temporary_id . " [$@]\n";
    }
    
  }
  
  return @newf
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

# Adapted from Spangle

=head2 match

 Title   : match
 Usage   : @overlap = $feat->match($f)
 Function: Returns an array of 3 numbers detailing
           the overlap between the two features.
           $overlap[0]  = 1 if any overlap, 0 otherwise
           The remaining elements of the array are only set if $overlap[0] = 1
           $overlap[1]  = left hand overlap (-ve if $f starts within $self, +ve if outside)
           $overlap[2]  = right hand overlap (-ve if $f ends within $self, $+ve if outside)
 Returns : @int
 Args    : Bio::SeqFeature::Generic

=cut

sub match {
  my ($self, $f1,$f2) = @_;
  
  my ($start1,
      $start2,
      $end1,
      $end2,
      $rev1,
      $rev2,
     );

  # Swap the coords round if necessary
  if ($f1->start > $f1->end) {
    $start1 = $f1->end;
    $end1   = $f1->start;
    $rev1   = 1;
  } else {
    $start1 = $f1->start;
    $end1   = $f1->end;
  }

  if ($f2->start > $f2->end) {
    $start2 = $f2->end;
    $end2   = $f2->start;
    $rev2   = 1;
  } else {
    $start2 = $f2->start;
    $end2   = $f2->end;
  }

  # Now check for an overlap
  if (($end2 > $start1 && $start2 < $end1) ) {
	
	#  we have an overlap so we now need to return 
	#  two numbers reflecting how accurate the span 
	#  is. 
	#  0,0 means an exact match with the exon
	# a positive number means an over match to the exon
	# a negative number means not all the exon bases were matched

	my $left  = ($start2 - $start1);
	my $right = ($end1 - $end2);
	
	if ($rev1) {
	    my $tmp = $left;
	    $left = $right;
	    $right = $tmp;
	}
	
	my @overlap;

	push (@overlap,1);

	push (@overlap,$left);
	push (@overlap,$right);

	return @overlap;
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
   my ($self,@genes) = @_;

   if (!defined($self->{'_output'})) {
     $self->{'_output'} = [];
   }
    
   if(defined @genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
}

############################################################




1;



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

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Pipeline::RunnableDB::EST_GeneBuilder;

#use diagnostics;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise;
use Bio::EnsEMBL::Pipeline::Runnable::Genomewise;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


# config file; parameters searched for here if not passed in as @args
require "Bio/EnsEMBL/Pipeline/EST_conf.pl";
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

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   nothing
    Args    :    string: chr1.1-10000

=cut

sub fetch_input {
    my( $self) = @_;

    print STDERR "IN THE 120 BRANCH\n";

    # the type of the genes being read is specified in Bio/EnsEMBL/Pipeline/EST_conf.pl 
    my $genetype =  $::est_genome_conf{'genetype'};

    print STDERR "Fetching input: " . $self->input_id. " \n";
    $self->throw("No input id") unless defined($self->input_id);

    # get genomic 
    my $chrid  = $self->input_id;
       $chrid =~ s/\.(.*)-(.*)//;

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

      my @exons = $transcripts[0]->get_all_Exons;

      # skip genes with one exon
      if(scalar(@exons) == 1){
	$single++;
	next GENE;
      }

      # keep only genes in the forward strand
      if($exons[0]->strand == 1){
	push (@plus_transcripts, $transcripts[0]);
	next GENE;
      }
    }
    print STDERR "In EST_GeneBuilder.fetch_input(): ".scalar(@plus_transcripts) . " forward strand genes\n";
    print STDERR "($single single exon genes thrown away)\n";

    # cluster each strand
    if(scalar(@plus_transcripts)){
      
      $self->_flush_Transcripts;  # this empties out the array $self->{'_transcripts'}
      
      my @transcripts  = $self->cluster_transcripts(\@plus_transcripts);
      
      #my $num1=1;
#      foreach my $tran (@transcripts){
#	my @exons = $tran->get_all_Exons;

#	print STDERR "In EST_GeneBuilder.fetch_input()\n";
#	print STDERR "transcript ".$num1." with ".scalar(@exons)." exons\n";
#	my $ccount=1;
#	foreach my $exon (@exons){
#	  print STDERR "Exon $ccount : ".$exon->start." ".$exon->end."\n";
#	  $ccount++;
#	}
		
#	$num1++;
#      }
      
      # make a genomewise runnable for each cluster of transcripts
      foreach my $tran (@transcripts){
	#print STDERR "new genomewise with one ".ref($tran). "\n";

	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
									    -genomic => $contig->primary_seq,
      								   );
#	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise();

	$self->add_runnable($runnable);
#	$runnable->seq($contig->primary_seq);
	$runnable->add_Transcript($tran);
      }
    }
    
    #my @trans = $self->_get_all_Transcripts;
#    print STDERR "In EST_GeneBuilder.fetch_input() after setting all runnables, forward_strand!!\n";
#    foreach my $t (@trans){
#      my $ccount=1;
#      foreach my $exon ($t->get_all_Exons){
#	print STDERR "Exon $ccount : ".$exon->start." ".$exon->end."\n";
#	  $ccount++;
#	}
#    }
	


    
    # minus strand - flip the vc and hope it copes ...
    # but SLOW - get the same genes twice ...
    
    print STDERR "\n****** reverse strand ******\n\n";

    my $revcontig = $contig->invert;
    my @revgenes  = $revcontig->get_Genes_by_Type($genetype, 'evidence');
    my @minus_transcripts;
    
    print STDERR "Number of genes = " . scalar(@revgenes) . "\n";
    if(!scalar(@genes)){
      $self->warn("No reverse strand genes found");
    }

    $single=0;
    REVGENE:    
    foreach my $gene (@revgenes) {
      my @transcripts = $gene->each_Transcript;
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
      # owwwww my brain hurts
      elsif($exons[0]->strand == 1){
	push (@minus_transcripts, $transcripts[0]);
	next REVGENE;
      }
    }
    print STDERR "In EST_GeneBuilfer.fetch_input(): ".scalar(@minus_transcripts) . " reverse strand genes\n";
    print STDERR "($single single exon genes thrown away)\n";
    
    if(scalar(@minus_transcripts)){
      
      $self->_flush_Transcripts;  # this empties out the array $self->{'_transcripts'}
      
#      # test
#      my $num3=1;
#      foreach my $tran (@minus_transcripts){
#	my @exons = $tran->get_all_Exons;
	
#	print STDERR "In EST_GeneBuilder.fetch_input() before cluster transcripts\n";
#	print STDERR "transcript ".$num3." : ".$tran."\n";
#	my $ccount=1;
#	foreach my $exon (@exons){
#	  print STDERR "Exon $ccount : ".$exon->start." ".$exon->end."\n";
#	  $ccount++;
#	}

#	$num3++;
#      } 

      my $reverse = 1;
      my @transcripts = $self->cluster_transcripts(\@minus_transcripts,$reverse);  
      
      # minus transcripts need some fancy schmancy feature and contig inversion. 
      
#      my $num2=1;
#      foreach my $tran (@transcripts){
#	my @exons = $tran->get_all_Exons;
	
#	print STDERR "In EST_GeneBuilder.fetch_input() after cluster transcripts\n";
#	print STDERR "transcript ".$num2." : ".$tran." exons\n";
#	my $ccount=1;
#	foreach my $exon (@exons){
#	  print STDERR "Exon $ccount : ".$exon->start." ".$exon->end."\n";
#	  $ccount++;
#	}

#	$num2++;
#      } 

      foreach my $tran (@transcripts) {
	#print STDERR "new genomewise with ".$tran."\n";
	
	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise(
									    -genomic => $revcontig->primary_seq,
      								   );
#	my $runnable = new Bio::EnsEMBL::Pipeline::Runnable::Genomewise();

	$self->add_runnable($runnable, $reverse);
#	$runnable->seq($contig->primary_seq);
	$runnable->add_Transcript($tran);
      }
    }
}

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

=head2 cluster_transcripts

    Title   :   cluster_transcripts
    Usage   :   @new_transcripts= $self->cluster_transcripts(@read_transcripts)
    Function:   clusters input array of transcripts
    Returns :   @Bio::EnsEMBL::Transcript
    Args    :   @Bio::EnsEMBL::Transcript

=cut

sub cluster_transcripts {
  my ($self, $alltranscripts, $reverse) = @_;

  # need all the exons - we're going to cluster them; while we're at it, check consistency of hid & strand
  my %est2transcript;  # relate transcript to est id
  my %exon2est;        # relate exon to est id
  
  my ($ref_transcripts,$ref_exons) = $self->check_transcripts(\%est2transcript, \%exon2est, @$alltranscripts);
  
  my @exons = @{ $ref_exons };
  $self->_filtered_Transcripts( @{ $ref_transcripts } );

  # store references to the exon and transcript hashes
  $self->{'_est2transcript'} = \%est2transcript;
  $self->{'_exon2est'}       = \%exon2est;
 
  #print STDERR scalar(@exons) . " exons found\n";
  
  # cluster the exons
  my $main_cluster    = $self->make_clusters(@exons);
 
  # set start and end in the exons
  $main_cluster       = $self->find_common_ends($main_cluster,$reverse);

  # form cluster-pairs
  my @linked_clusters = $self->link_clusters($main_cluster);
  #print STDERR "linked_clusters are ".ref( $linked_clusters[0] )."\n";

  # link recursively the cluster-pairs to build transcripts
  my @transcripts     = $self->process_clusters(\@linked_clusters);

  return @transcripts;

}

=head2
  
   Title   : _filtered_Transcripts
   Function: this method gets/sets an array of Bio::EnsEMBL::Transcript objects which have been filtered
             previously in the method 'check_Transcripts'

=cut

sub _filtered_Transcripts {
   my ($self, @filtered_transcripts) = @_;
   unless ( $self->{'_filtered_transcripts'} ){
     $self->{'_filtered_transcripts'} = [];
   }
   if ( @filtered_transcripts ){
     unless ( $filtered_transcripts[0]->isa('Bio::EnsEMBL::Transcript') ){
       $self->throw('Must pass a Bio::EnsEMBL::Transcript object');
     }
     push( @{ $self->{'_filtered_transcripts'} }, @filtered_transcripts );
   }
   return @{ $self->{'_filtered_transcripts'} }
}

=head2 check_transcripts

    Title   :   check_transcripts
    Usage   :   $self->check_transcripts
    Function:   checks transcripts obtained from EST2Genome for consistency among exons
                in strand, hid, and also checks for exon content
    Returns :   @Bio::EnsEMBL::Exon
    Args    :   @Bio::EnsEMBL::Transcript, ref to hash for linking hid to transcript, ref to hash for linking exon to hid, 

=cut

sub check_transcripts {
  my ($self, $hid_trans, $exon_hid, @transcripts) = @_;
  my @allexons;       # here we'll put all exons that pass the check
  my @alltranscripts; # here we'll put all the transcripts that pass the check
  my $exon_adaptor = $self->dbobj->get_ExonAdaptor;

  TRANS: 
  foreach my $transcript (@transcripts){
    my @exons = $transcript->get_all_Exons;
    my $hid;
    my $strand;
    my $num4=0;
    
    foreach my $exon(@exons){
      my $hstart;
      my $hend;
      $num4++;
    
      # check strand consistency
      if(!defined $strand) { 
	$strand = $exon->strand; 
      }
      if($strand ne $exon->strand){
	$self->warn("strand not consistent among exons for " . $transcript->temporary_id . " - skipping it\n");
	next TRANS;
      }
      
      # check supporting_feature consistency
      # print STDERR "In EST_GeneBuilder.check_transcripts()...";
      $exon_adaptor->fetch_evidence_by_Exon($exon);
      my @sf = $exon->each_Supporting_Feature;      

      # the source_tag of the supporting evidence is specified in Bio/EnsEMBL/Pipeline/EST_conf.pl
      my $evidence_tag = $::evidence_conf{'evidence_tag'};
    SF:
      foreach my $feature ( @sf ) {
	# because we get all the supporting features indiscriminately
	next SF unless $feature->source_tag eq $evidence_tag;
	
	if(!defined $hid) { 
	  $hid = $feature->hseqname; 
	}
	
	if($hid ne $feature->hseqname){
	  $self->warn("hid not consistent among exons for " . $transcript->temporary_id . " - skipping it\n");
	  next TRANS;
	}
      }

      if ( defined( @sf ) && scalar( @sf ) != 0 ) {
	$hstart = $sf[0]->hstart;
	$hend   = $sf[$#sf]->hend;
      }

      $$exon_hid{$exon}{"hid"}    = $hid;
      $$exon_hid{$exon}{"hstart"} = $hstart;
      $$exon_hid{$exon}{"hend"}   = $hend;
    }
    push(@allexons, @exons);
    
    if(defined $hid && defined $$hid_trans{$hid}) { 
      $self->warn("$hid is being used by more than one transcript!\n"); 
    }
#    print STDERR "hid: ".$hid."\ttranscript: ".$transcript."\n";
    $$hid_trans{$hid} = $transcript;
    push (@alltranscripts, $transcript);
  }

  return (\@alltranscripts,\@allexons);
}



sub make_clusters {
  my ($self, @exons) = @_;

  my $main_cluster = new Bio::EnsEMBL::SeqFeature; # main cluster feature - holds all subclusters
  @exons = sort { $a->start <=> $b->start } @exons;

  return unless ( scalar( @exons) > 0 ); # no point if there are no exons!
  
  # Create the first cluster object
  my $subcluster = new Bio::EnsEMBL::SeqFeature;
  
  # Start off the cluster with the first exon
  $subcluster->add_sub_SeqFeature($exons[0],'EXPAND');
  $subcluster->strand($exons[0]->strand);  

  $main_cluster->add_sub_SeqFeature($subcluster,'EXPAND');
  
  # Loop over the rest of the exons
  my $count = 0;

 EXON:
  foreach my $e (@exons) {
    
    if ($count > 0) {
      my @overlap = $self->match($e, $subcluster);    
      
      # Add to cluster if overlap AND if strand matches
      if ( $overlap[0] && ( $e->strand == $subcluster->strand) ) { # strand is not checked in 'match'
	$subcluster->add_sub_SeqFeature($e,'EXPAND');
      }  
      else {
	# Start a new cluster
	$subcluster = new Bio::EnsEMBL::SeqFeature;
	$subcluster->add_sub_SeqFeature($e,'EXPAND');
	$subcluster->strand($e->strand);
	
	# and add it to the main_cluster feature
	$main_cluster->add_sub_SeqFeature($subcluster,'EXPAND');	
      }
    }
    $count++;
  }
    
  # If we're building genes from ESTs, we remove any clusters that have only 1 exon 
  # On the other hand, this may cause fragmentation problems ... if only 1 EST spans a gap 
  # between 2 clusters that really are related
  
  # If we're building genes from cDNAs, we need to take clusters with 1 exon as well, since
  # we have less hits, but these are of better quality than the ESTs

  my @exon_clusters = $main_cluster->sub_SeqFeature;
  print STDERR "In EST_GeneBuilder.make_clusters(), ".scalar( @exon_clusters )." exon clusters found\n";
  
## we allow clusters with cDNA (or mRNA) evidence if they have 
## one or more exons, hence in this case 'strict_lower_bound' = 0, but we don't allow 
## single-exon clusters with only ESTs as supporting evidence, hence in this case 'strict_lower_bound' = 1

  my $strict_lower_bound = $::evidence_conf{'strict_lower_bound'};
  
  # empty out the list of clusters
  $main_cluster->flush_sub_SeqFeature;
  print STDERR "keeping clusters with ".($strict_lower_bound+1)." or more exons...\n";

  # store cluster and transcript in which each exon lives:
  my %exon_in_cluster;
  my %exon_in_transcript;
  
  foreach my $exon_cluster ( @exon_clusters ){
    foreach my $exon ( $exon_cluster->sub_SeqFeature ){
      $exon_in_cluster{ $exon } = $exon_cluster;
    }
  }
  my @transcripts   = $self->_filtered_Transcripts; # transcripts filtered in 'check_Transcripts'
  foreach my $tran ( @transcripts ){
    foreach my $exon ($tran->get_all_Exons ){
      $exon_in_transcript{$exon} = $tran;
    }
  }

  # prune the exon clusters
 EXON_CLUSTER:
  foreach my $exon_cluster (@exon_clusters) {
    my @exons = $exon_cluster->sub_SeqFeature;
    
    # throw away empty clusters
    if ( scalar(@exons) == 0 ){
      next EXON_CLUSTER;
    }

    # in drosophila ESTs and cDNAs are together, so we need to differentiate
   # my $found_est  = 0;
#    my $found_cdna = 0;
#    foreach my $exon ( @exons ){
#      foreach my $feature ( $exon->each_Supporting_Feature ){
#	my $hid = $feature->hseqname;
#	if ( $hid =~ /\./g ){  # ESTs have labels like AT19106.5prime or like LD23056.contig.5_3prime
#	  $found_est = 1;
#	}
#	else{                  # whereas cDNAs have labels like AI944764 (without dots)
#	  $found_cdna = 1;
#	}
#      }
#    }
#    if ( $found_est == 1 && $found_cdna != 1){   # keep clusters with one exon only if it has a cDNA hit
#       $strict_lower_bound = 1;
#    }
#    else{
#       $strict_lower_bound = 0;
#     }

    # throw away clusters with too little exons (this should be set in Bio::EnsEMBL::Pipeline::EST_conf.pl)
    unless ( scalar(@exons) > $strict_lower_bound ){

      # but allow for sister clusters, i.e. if the exon is part of a transcript...
      my $transcript = $exon_in_transcript{ $exons[0] };

      # ...with exons which are in clusters with more than one exon
      my @clusters;
      foreach my $exon ( $transcript->get_all_Exons ){
	push ( @clusters, $exon_in_cluster{$exon} );
      }
      # my threshold
      my $howmany=1;

      # count how many clusters there are with more than $howmany exon
      my $count = 0;
      foreach my $cluster ( @clusters ){
	my @exons = $cluster->sub_SeqFeature;
	if ( scalar(@exons) > $howmany ){
	  $count++;
	}
      }
      
      # throw away this cluster unless there is at least $somany sister-cluster with more than one exon
      my $somany=1;
      unless ( $count >= $somany ){
	next EXON_CLUSTER; 
      }
    }
    
    # keep the rest
    $main_cluster->add_sub_SeqFeature($exon_cluster,'EXPAND');
  }
  
  @exon_clusters = $main_cluster->sub_SeqFeature;
  print STDERR "In EST_GeneBuilder.make_clusters(), ".scalar( @exon_clusters )." exon clusters kept\n";
      
  return $main_cluster;
}


 

=head2 find_common_ends

    Title   :   find_common_ends
    Usage   :   $cluster = $self->find_common_ends($precluster)
    Function:   Resets the ends of the clusters to the most frequent coordinate
    Returns :   
    Args    :   

=cut

# should we introduce some sort of score weighting?

sub find_common_ends {
  my ($self, $main_cluster,$reverse) = @_;
  my $exon_hid = $self->{'exon2est'};
  my @clusters = $main_cluster->sub_SeqFeature;
  @clusters    = sort { $a->start <=> $b->start  } @clusters;
  my $count    =  0 ;
  my ($upstream_correct, $downstream_correct) = (0,0);
  
  print STDERR "finding common ends and checking for splice sites...\n";
  foreach my $cluster(@clusters) {
    my %start;
    my %end;
    my $newstart;
    my $newend;
    my $maxend = 0;
    my $maxstart = 0;
    
    # count frequency of starts & ends
    foreach my $exon ($cluster->sub_SeqFeature) {
      my $est = $$exon_hid{$exon}{'hid'};
      $start{$exon->start}++;
      $end{$exon->end}++;
      
    }

    #print STDERR "cluster $count\nstarts: ";
    while ( my ($key, $value) = each %start) {
      if ($value > $maxstart){
	$newstart = $key;
	$maxstart = $value;
      } elsif ($value == $maxstart){
	#print STDERR "$key is just as common as $newstart - ignoring it\n";
      }
      #print STDERR "$key => $value\n";
    }
    #print STDERR "ends: ";

    while ( my ($key, $value) = each %end) {
      if ($value > $maxend){
	$newend = $key;
	$maxend = $value;
      } elsif ($value == $maxend){
	#print STDERR "$key is just as common as $newend - ignoring it\n";
      }
      #print STDERR "$key => $value\n";
    }

    # if we haven't got a clear winner, we might as well stick with what we had
    if( $maxstart <= 1 ) {
      $newstart = $cluster->start;
    }
    if( $maxend <= 1 ) {
      $newend = $cluster->end;
    }

    # first and last clusters are special cases - potential UTRs, take the longest one.
    if( $count == 0) {
      #print STDERR "first cluster\n";
      $newstart = $cluster->start;
    } elsif ( $count == $#clusters ) {
      #print STDERR "last cluster\n";
      $newend = $cluster->end;
    }
    #print STDERR "we used to have: " . $cluster->start . " - " . $cluster->end . "\n";
    #print STDERR "and the winners are ... $newstart - $newend\n";


    # reset the ends of all the exons in this cluster to $newstart - $newend
    foreach my $exon ($cluster->sub_SeqFeature) {
      $exon->start($newstart);
      $exon->end($newend);
    }
    
    # check the splice sites for the whole cluster
    my $contig = $self->vc;
    my ($upstream, $downstream);
        
    # take the 2 bases right before the exon and after the exon
    
    # forward strand
    if ($cluster->strand == 1 && !defined( $reverse) ){
      my $seq = $contig->primary_seq; # a Bio::PrimarySeq object
      
      # catch possible exceptions in gettting the sequence (it might not be there!)
      eval{
	$upstream   = $seq->subseq( ($cluster->start)-2 , ($cluster->start)-1 ); 
      };
      if ($@){
	print STDERR "Unable to get subsequence (".(($cluster->start)-2).",".(($cluster->start)-1).")\n";
	print STDERR $@;
	$upstream = 'NN';
      }
      eval{
	$downstream = $seq->subseq( ($cluster->end)+1   , ($cluster->end)+2   ); 
      };
      if ($@){
	print STDERR "Unable to get subsequence (".(($cluster->end)+1).",".(($cluster->end)+2).")\n";
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
      if ( $count !=$#clusters && $downstream eq 'GT') { 
	$downstream_correct++;
      }
      $count++;
    }
    
    # reverse strand
    if ($reverse){
      my $revcontig = $contig->invert;   # invert the contig, since the exons were retrieved in revcontig
      my $seq = $revcontig->primary_seq; # a Bio::PrimarySeq object
      
      # catch possible exceptions in gettting the sequence (it might not be there!)
      eval{
	$upstream = $seq->subseq( ($cluster->start)-2 , ($cluster->start)-1 ); 
      };
      if ($@){
	print STDERR "Unable to get subsequence (".(($cluster->start)-2).",".(($cluster->start)-1).")\n";
	print STDERR $@;
	$upstream ='NN';
      }
      eval{
	$downstream   = $seq->subseq( ($cluster->end)+1   , ($cluster->end)+2   ); 
      };
      if ($@){
	print STDERR "Unable to get subsequence (".(($cluster->end)+1).",".(($cluster->end)+2).")\n";
	print STDERR $@;
	$downstream = 'NN';
      }
      #  the reverse strand we're looking at coordinates in the reversed-complement contig:
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
      if ( $count != $#clusters && $downstream eq 'AC') {       
	$upstream_correct++;
      }
      if ( $count != 0 && $upstream eq 'CT') { 
	$downstream_correct++;
      }
      $count++;
    }
  }
  print STDERR "upstream splice-sites correct: ".$upstream_correct. 
               " out of ".($#clusters)." potential splice-sites\n";
  print STDERR "downstream splice-sites correct: ".$downstream_correct. 
               " out of ".($#clusters)." potential splice-sites\n";

  return $main_cluster;
}

=head2 is_in_cluster

  Title   : is_in_cluster
  Usage   : my $in = $self->is_in_cluster($id, $cluster)
  Function: Checks whether a feature with id $id is in a cluster
  Returns : 0,1
  Args    : String,Bio::SeqFeature::Generic

=cut

sub is_in_cluster {
  my ($self, $id, $cluster) = @_;

  my $exon2est = $self->{'_exon2est'};

  foreach my $exon ( $cluster->sub_SeqFeature ) {
    return (1) if $$exon2est{$exon}{'hid'} eq $id;
  }
  return 0;
}

=head2 link_clusters

    Title   :   link_clusters
    Usage   :   $self->link_clusters
    Function:   Links clusters into pairs if any of their exons share the same EST as supporting feature
    Returns :   an array of exon-clusters (Bio::EnsEMBL::SeqFeature) with the links set
    Args    :   an array of exon-clusters (Bio::EnsEMBL::SeqFeature) 

=cut

sub link_clusters{
  my ($self, $main_cluster, $tol) = @_;
  my $exon_hid = $self->{'_exon2est'};
  $tol = 10 unless $tol; # max allowed discontinuity in the EST sequence
  
  my @clusters = $main_cluster->sub_SeqFeature; # exon clusters
  print STDERR "Created " . scalar (@clusters) . " clusters\n";

  # Sort the clusters by start position
  @clusters = sort { $a->start <=> $b->start } @clusters;

  # Loop over all clusters
 CLUSTER1:  
  for (my $c1 = 0; $c1 < $#clusters; $c1++) {
    #print "\nFinding links in cluster number $c1... \n";
    my $cluster1 = $clusters[$c1];             
    my @exons1   = $cluster1->sub_SeqFeature; 
    my $c2       = $c1 + 1;
    my $limit    = 2;  # how many exons away we look to find a link, 2 is the limit used by sub make_ExonPairs 
                       # in GeneBuilder.pm
    
    # NOTE!!!!
    # Maybe put this definition inside CLUSTER2 loop ?????
    my %used_evidence; # we keep track of the evidence that have been used to make links
    
    # Now look (in each exon ) in the next $limit clusters to find a link
  CLUSTER2:    
    while ( $c2 <= ( $c1 + $limit ) && $c2 <= $#clusters ) {
      # print STDERR "...trying to link clusters $c1 and $c2...\n";
      my $cluster2 = $clusters[$c2];      
      my @exons2   = $cluster2->sub_SeqFeature;
      
      # Only look at a cluster if it starts after the end of the current cluster 
      if ($cluster2->start > $cluster1->end) {
	my $foundlink = 0;

	# Now loop over the supporting features in the first linking cluster
      EXON1:
	foreach my $exon1 (@exons1) { 
	  # Does the second cluster contain this database id?
	  my $hid1 = $$exon_hid{$exon1}{'hid'}; 
	  my $exon_link;
	  
	EXON2:
	  foreach my $exon2 (@exons2){ 

	    # if we have the same id in two clusters then...
	    if ($$exon_hid{$exon2}{'hid'} eq $hid1){
	      $exon_link = $exon2;
	      if ( $used_evidence{$hid1} && $used_evidence{$hid1} == 1 ){
		# print STDERR $hid1." used already, searching further...\n";
		next EXON2;
	      }
	      
	      # Check the strand is the same and only allow a $tol discontinuity in the h-sequence
	      my $tol1 = abs($$exon_hid{$exon1}{'hend'}   - $$exon_hid{$exon_link}{'hstart'});
	      my $tol2 = abs($$exon_hid{$exon1}{'hstart'} - $$exon_hid{$exon_link}{'hend'});

	      if ( $exon1->strand == $exon_link->strand && ( $tol1 < $tol || $tol2 < $tol ) ) {
		push(@{ $cluster1->{'_forward'} } ,$cluster2);
		print STDERR "*** linking cluster $c1 to cluster $c2 ***";
		# print STDERR "hid1 : ".$hid1." hid2 : ".$$exon_hid{$exon_link}{'hid'}."\n";
		
		# just announcing alternative splicing
		if ( $c2 >= $c1+$limit && scalar( @{ $cluster1->{'_forward'} } ) ==2 ){
		  print STDERR " -------------- Alternative Splicing ---------------\n";
		}

		# we've found a link, we can stop searching in this cluster
		$used_evidence{$hid1} = 1;
		$c2++;
		next CLUSTER2;		
	      }
	      # else we search further in the cluster
	      else{
		print STDERR "cannot make link: tol1 = $tol1, tol2 = $tol2\n";
		next EXON2;
	      }
	    }

	  } # end of EXON2

	}   # end of EXON1

      }     # end of if ($cluster2->start > $cluster1->end)

      # if we are here that means that no link was found with cluster2
      # print STDERR "couldn't find link between cluster $c1 and cluster $c2\n";
      $c2++;
    }       # end of CLUSTER2

  }         # end of CLUSTER1

  return @clusters;

}

=head2 process_clusters

 Title   : process_clusters
 Usage   : my @transcripts = $self->process_clusters(@clusters);
 Function: Links Cluster Pairs into Transcripts
 Example : 
 Returns : Array of Bio::EnsEMBL::Transcript
 Args    : Array of exon-clusters, which here are Bio::EnsEMBL::SeqFeature

=cut

sub process_clusters{
  my ($self, $clusters) = @_;
  my @clusters = @{ $clusters };

  # comment: this link all the exon clusters (pairs) into a whole transcript. It starts from the first exon 
  # (cluster). Is there a possibility to start in the second, maybe third exon? this would provide further 
  # alternative transcripts. One reason for doing this is that the second (or third) exon pair may have been 
  # paired due to a different est. Could we push this to say that any exon where an EST is first used could be a 
  # potential start of an alternative transcript? In any case, this would seem more logical with
  # cDNA's but a bit risky with EST's since the latter are just one end of the RNA, so different EST's could in 
  # fact be different pieces of the same RNA. Does this make any sense?

  foreach my $cluster (@clusters){
    
    if ( $self->_isHead($cluster,\@clusters) == 1){
      # only start a transcript on this cluster if it is not linked to another cluster from the left 
      
      # get exon from cluster (this is a fresh new object, see _get_RepresentativeExon method)
      my $exon = $self->_get_RepresentativeExon($cluster);
      
      # create a new transcript
      my $transcript = new Bio::EnsEMBL::Transcript;
      
      # with the new representative exon in it
      $transcript->add_Exon($exon);
      $self->_put_Transcript($transcript);
      
      # cluster the exons into transcripts according to the cluster-pairs made before
      $self->_recurseTranscript($cluster,$transcript);
    }
  }
  # at this point we could check/set transcript ids, order of exons, etc...

  print STDERR "In EST_GeneBuilder.process_clusters: "
    .scalar( $self->_get_all_Transcripts )." transcripts have been created\n";
  return $self->_get_all_Transcripts;
}

=head2 _recurseTranscript

 Title   : _recurseTranscript
 Usage   : $self->_recurseTranscript($cluster,\@clusters);
 Function: borrows from _recurseTranscript in GeneBuilder.pm, it links all cluster pairs
           recursively to form transcripts, allowing alternative transcripts, and it puts them into
           $self->{'_transcripts'}
 Returns : nothing
 Args    : it reads an exon-cluster and the transcript it sits in

=cut

sub _recurseTranscript { 
  my ($self,$cluster,$current_transcript) = @_; 
  
  # copy first the information about which clusters the current transcript holds into a temporary transcript
  my $temp_transcript = new Bio::EnsEMBL::Transcript;
  foreach my $ex ( $current_transcript->get_all_Exons ){
    $temp_transcript->add_Exon($ex);  
  }

  if ( defined( $cluster->{'_forward'} ) ){
    my @cluster_links = @{ $cluster->{'_forward'} }; # this gives the clusters to which $cluster links
    my $count = 0;
  LINK:
    foreach my $cluster_link ( @cluster_links ){
      my $exon_link = $self->_get_RepresentativeExon($cluster_link);
      
      if ($count>0){ # this occurs for a second link at a given cluster
	
	#print STDERR "\n STARTING AN ALTERNATIVE TRANSCRIPT \n\n";

	# Start a new transcript
	my $new_transcript = new Bio::EnsEMBL::Transcript;
	
	# put in it everything we have at this level
	foreach my $temp_exon ( $temp_transcript->get_all_Exons ){
	  
	  # we put a new one to make sure that exons are not shared between transcripts
	  my $exon_copy = new Bio::EnsEMBL::Exon;
	  $exon_copy->start($temp_exon->start);
	  $exon_copy->end($temp_exon->end);
	  $exon_copy->strand($temp_exon->strand);
	  foreach my $sf ( $temp_exon->each_Supporting_Feature ){
	    $exon_copy->add_Supporting_Feature($sf);
	  }
	  $new_transcript->add_Exon($exon_copy);
	}
	
	# add the next exon in the link
	$new_transcript->add_Exon($exon_link);
	
	# and add the new transcript to the array of transcripts $self->{'_transcripts'}
	$self->_put_Transcript($new_transcript);
	
	# continue linking from this one
	$self->_recurseTranscript($cluster_link, $new_transcript);
      }
      else{ # this occurs for the first time in (the first link) thus fills up the current transcript
	
	# add the next exon in the link to the current transcript
	$current_transcript->add_Exon($exon_link);
	
	# continue linking from this one
	$self->_recurseTranscript($cluster_link, $current_transcript);
      }
      $count++; 
    }         # end of LINK:
  }           # end of if-scope
  return;     # eventually, @cluster_links will be empty, and we return
}
 
=head2 _get_RepresentativeExon

 Title   : _get_RepresentativeExon
 Usage   : my $exon = $self->_get_RepresentativeExon($cluster);
 Function: it gives you one exon from the exon clusters, since find_common_ends should have occurred
           before hand (check taht if you have modified anything), we should get somethings close to
           a proper exon with the start and end of the $cluster. Note that a new exon object is created,
           since alternative transcript may share exons, adn we do not want them to share them as objects as well,
           
 Example : 
 Returns : an Bio::EnsEMBL::Exon object
 Args    : one exon-cluster, which in this context is a Bio::EnsEMBL::SeqFeature

=cut

sub _get_RepresentativeExon {
  my ($self,$cluster) = @_;
  
  # get the exons in the cluster
  my @exons = $cluster->sub_SeqFeature; 

  # simply take the first one (any one will do, I guess)
  my $exon = $exons[0];

  if ( defined($exon) ){
    unless ( $exon->isa("Bio::EnsEMBL::Exon") ){
      $self->throw("[$exon] is a ".ref($exon)." but should be a Bio::EnsEMBL::Exon");
    }
  }
  else{
    $self->throw("exon-cluster [$cluster] does not contain any exon");
  }

  # make a new exon out of this one, otherwise alternative transcripts will share exons objects
  # and they will get modified twice (in two different transcripts) by MiniGenomewise
  my $exon_copy = new Bio::EnsEMBL::Exon;
  $exon_copy->start($exon->start);
  $exon_copy->end($exon->end);
  $exon_copy->strand($exon->strand);
  foreach my $sf ( $exon->each_Supporting_Feature ){
    $exon_copy->add_Supporting_Feature($sf);
  }

  # at this point we could add all the EST-supporting-features of the rest of the exons in the cluster
  # to the supporting features of this exon, something like:
  shift @exons;
  # my $source_tag = $::evidence_conf{'evidence_tag'};
  foreach my $ex ( @exons ){
    foreach my $sf ($ex->each_Supporting_Feature){
      # add only those that relate to EST's?
      #next unless ( $sf->source_tag eq $source_tag;
      $exon_copy->add_Supporting_Feature($sf);
    }
  }

  return $exon_copy;
}
  
=head2 _isHead

 Title   : _isHead
 Usage   : my $is_first_cluster = $self->_isHead($cluster,\@clusters); 
 Function: checks through all clusters in @clusters to see whether
           $cluster is linked to a preceeding cluster.
 Example : 
 Returns : 0 (if the cluster is not the first), 1 (if the cluster is the first, i.e. the head)
 Args    : one exon-cluster and one Array-ref to exon-clusters, these are Bio::EnsEMBL::SeqFeature

=cut

sub _isHead {
  my ($self,$cluster,$clusters) = @_;
  foreach my $cluster2 (@{ $clusters }){
    foreach my $forward_link ( @{ $cluster2->{'_forward'} } ){
     
      my @overlap = $self->match($cluster,$forward_link); 
      if ( $cluster->start == $forward_link->start ){
	return 0;
      }
     ## this thing below is just an alternative way to see whether two clusters are the same

#      if ( $overlap[0] == 1 && $cluster->strand == $forward_link->strand ){
#	# if they overlap and share strand, we can infer they are the same
#	# this method should therefore only be used after make_clusters and find_common_ends have occurred
#	print STDERR "returning from _isHead\n";
#	return 0; 
#      }     
####################################
      
    }
  }
  return 1; 
}

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



sub add_runnable{
  my ($self, $value, $reverse) = @_;

  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 
  if (defined($value)) {
    if ($value->isa("Bio::EnsEMBL::Pipeline::RunnableI")) {
      if(defined $reverse){
	push(@{$self->{'_reverse_runnables'}},$value);
      }
      else {
	push(@{$self->{'_forward_runnables'}},$value);
      }
    } else {
      $self->throw("[$value] is not a Bio::EnsEMBL::Pipeline::RunnableI");
    }
  }
  
} 

sub each_runnable{
  my ($self, $reverse) = @_;
  
  if (!defined($self->{'_forward_runnables'})) {
    $self->{'_forward_runnables'} = [];
  }
  
  if (!defined($self->{'_reverse_runnables'})) {
    $self->{'_reverse_runnables'} = [];
  } 

  if(defined $reverse){
    return @{$self->{'_reverse_runnables'}};
  }
  else {
    return @{$self->{'_forward_runnables'}};
  }
  
}

sub run {
  my ($self) = @_;

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
  
#  ### test
#  my $tcount=0;
#  foreach my $gw_runnable( $self->each_runnable) {
    
#    print STDERR "In EST_GeneBuilder...forward_strand...the third, uf!\n";
#    print STDERR "reading transcripts in ".$gw_runnable."\n";
#    foreach my $t ($gw_runnable->get_all_Transcripts){
#      print STDERR "Transcript ".($tcount+1).": ".$t."\n";
#      my $count=1;
#      foreach my $exon ($t->get_all_Exons){
#	print STDERR "Exon $count : ".$exon->start." ".$exon->end."\n";
#      }
#      $tcount++;
#    }
#  }
#  print STDERR "In EST_GeneBuilder...from the runnables we have: ".$tcount." transcripts\n";
#
#  ### end of test

  my $tcount=0;
  foreach my $gw_runnable( $self->each_runnable) {

#    print STDERR "In EST_GeneBuilder...forward_strand...about to run genomewise";
#    my $ttcount=0;
#    foreach my $t ($gw_runnable->get_all_Transcripts){
#      $ttcount++;
#      $tcount++;
#    }
#    print STDERR " with $ttcount transcript(s)\n";
   
    $tcount++;
    $gw_runnable->run;

    # convert_output
    $self->convert_output($gw_runnable);
  }
  print STDERR $tcount." transcripts run in genomewise in the forward strand\n";

  # minus strand
  my $reverse = 1;
  my $tcount2=0;
  foreach my $gw_runnable( $self->each_runnable($reverse)) {

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
    $self->convert_output($gw_runnable, $reverse);
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
  my ($self, $gwr, $reverse) = @_;

  my @genes = $self->make_genes($gwr, $reverse);

  my @remapped = $self->remap_genes(\@genes, $reverse);

  # store genes
  $self->output(@remapped);
 
}



sub make_genes {

  my ($self, $runnable, $reverse) = @_;
  my $genetype = $self->genetype;
  my $analysis_obj = $self->analysis;
  my $count = 0;
  my @genes;

  my $time  = time; chomp($time);
  my $contig = $self->vc;

  # are we working on the reverse strand?
  if(defined $reverse){
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
      $excount++;
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
  my ($self, $genes, $reverse) = @_;
  my $contig = $self->vc;
  if (defined $reverse){
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


1;

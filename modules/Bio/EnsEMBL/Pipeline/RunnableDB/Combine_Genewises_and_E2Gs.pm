#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs
#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs

=head1 SYNOPSIS

my $t_e2g = new Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs(
									       -db        => $db,   
									       -input_id  => $input_id
									      );

$t_e2g->fetch_input();
$t_e2g->run();
$t_e2g->output();
$t_e2g->write_output(); #writes to DB

=head1 DESCRIPTION

Combines predictions from Genewises with Est2genome predictions from cDNA alignments:


=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Pipeline::Tools::ExonUtils;
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise;
use Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator;
use Bio::SeqIO;

# all the parameters are read from GeneBuild config files
use Bio::EnsEMBL::Pipeline::Config::GeneBuild::General   qw (
							     GB_INPUTID_REGEX 	   
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							     GB_GW_DBNAME
							     GB_cDNA_DBHOST
							     GB_cDNA_DBUSER
							     GB_cDNA_DBNAME
							     GB_cDNA_DBPASS
							     GB_COMB_DBHOST
							     GB_COMB_DBUSER
							     GB_COMB_DBNAME
							     GB_COMB_DBPASS
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted qw (
							     GB_TARGETTED_GW_GENETYPE
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Similarity qw (
							      GB_SIMILARITY_GENETYPE
							     );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Combined qw (
							    GB_cDNA_GENETYPE
							    GB_COMBINED_MAX_INTRON
							    GB_COMBINED_GENETYPE
							   );

#use diagnostics;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  # IMPORTANT:
  # SUPER creates db, which is created on run_GeneBuild_runnable
  # this db is a reference to GB_DBHOST@GB_DBNAME containing
  # the features and the dna, so here it is used as refdb only

  # db with the genewises
  my $genewise_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $GB_GW_DBHOST,
						   '-user'   => $GB_GW_DBUSER,
						   '-pass'   => $GB_GW_DBPASS,
						   '-dbname' => $GB_GW_DBNAME,
						  );
  
  
  $genewise_db->dnadb($self->db);
  $self->genewise_db($genewise_db);

  # db with the cdnas
  my $cdna_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_cDNA_DBHOST,
     '-user'   => $GB_cDNA_DBUSER,
     '-dbname' => $GB_cDNA_DBNAME,
     '-pass' => $GB_cDNA_DBPASS,
    ); 
  
  $cdna_db->dnadb($self->db);
  $self->cdna_db($cdna_db);


  # db where we will write the output (different from any of the above):
  my $comb_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    '-host'   => $GB_COMB_DBHOST,
						    '-user'   => $GB_COMB_DBUSER,
						    '-dbname' => $GB_COMB_DBNAME,
						    '-pass'   => $GB_COMB_DBPASS,
                                                    ); 
  
  $comb_db->dnadb($self->db);
  $self->output_db($comb_db);

  
  #$self->db($genedb);
 
  return $self;
}

############################################################

sub fetch_input{
  my ($self,@args) = @_;
  
  my $entry = $self->input_id;
  my $chr_name;
  my $start;
  my $end;
  
  # input format is given by GeneConf::$GB_INPUTID_REGEX and is usually of the form '(^\S+\.\S+)\.(\d+)-(\d+)',
  
  if(!($entry =~ /$GB_INPUTID_REGEX/ ) ){
    $self->throw("Not a valid input id... $entry");
  }
  $chr_name = $1;
  $start    = $2;
  $end      = $3;
  print STDERR "input_id id : $chr_name .  $start - $end\n";  
  
  # genewises db
  my $slice_adaptor = $self->genewise_db->get_SliceAdaptor();
  my $slice         = $slice_adaptor->fetch_by_chr_start_end($chr_name,$start,$end);
  $self->query($slice);

  # get genewise genes
  my @similarity_genes = @{$self->query->get_all_Genes_by_type($GB_SIMILARITY_GENETYPE)};
  my @targetted_genes  = @{$self->query->get_all_Genes_by_type($GB_TARGETTED_GW_GENETYPE)};
  print STDERR "got " . scalar(@similarity_genes) . " similarity genewise genes\n";
  print STDERR "got " . scalar(@targetted_genes) . " targetted genewise genes\n";
  $self->gw_genes( @similarity_genes, @targetted_genes );
  
  # cdnas db
  my $slice_adaptor2 = $self->cdna_db->get_SliceAdaptor();
  my $cdna_vc        = $slice_adaptor2->fetch_by_chr_start_end($chr_name,$start,$end);
  $self->cdna_vc($cdna_vc);
  
  # get cdnas 
  my @e2g = @{$self->cdna_vc->get_all_Genes_by_type($GB_cDNA_GENETYPE)};
  print STDERR "got " . scalar(@e2g) . " cdnas ($GB_cDNA_GENETYPE)\n";
  
  # filter cdnas
  my @newe2g = $self->_filter_cdnas(@e2g);
  $self->e2g_genes(@newe2g);
  print STDERR "got " . scalar($self->e2g_genes) . " cdnas after filtering\n";
}


############################################################

sub run {
  my ($self,@args) = @_;
  
  # merge exons with frameshifts into a big exon
  my @merged_gw_genes = $self->_merge_gw_genes;
  print STDERR "got " . scalar(@merged_gw_genes) . " merged genewise genes\n";  
  
  # first of all, sort genewises by exonic length and genomic length  
  @merged_gw_genes = sort { my $result = ( $self->_transcript_length_in_gene($b) <=>
					   $self->_transcript_length_in_gene($a) );
			    unless ($result){
			      return ( $self->_transcript_exonic_length_in_gene($b) <=>
				       $self->_transcript_exonic_length_in_gene($a) );
			    }
			    return $result;
			    
			  } @merged_gw_genes;
  
  print STDERR "now have " . scalar(@merged_gw_genes) . " merged genes\n";
		   
  # keep track of the e2g genes that have been used already
  my %used_e2g;
  
  # find one-2-one matching between proteins and cdnas
 GENEWISE:
  foreach my $gw (@merged_gw_genes){
 
    # should be only 1 transcript
    my @gw_tran  = @{$gw->get_all_Transcripts};
    my @gw_exons = @{$gw_tran[0]->get_all_Exons}; # ordered array of exons
    my $strand   = $gw_exons[0]->strand;
    
    if($gw_exons[$#gw_exons]->strand != $strand){
      $self->warn("first and last gw exons have different strands - can't make a sensible combined gene\n");
      $self->unmatched_genewise_genes($gw);
      next GENEWISE;
    }
    
    # find the matching cdna
    my ($matching_cdnas,$utr_length_hash) = $self->match_protein_to_cdna($gw);
    my @matching_e2gs                     = @{$matching_cdnas};
    my %utr_length_hash                   = %$utr_length_hash;

    ##############################
    # now we can sort by length of the added UTR
    # using %utr_length_hash
    # not yet in use
    ##############################
    
    if(!@matching_e2gs){
      # store non matched genewise
      print STDERR "no matching cDNA for " . $gw->dbID ."\n";
      $self->unmatched_genewise_genes($gw);      
      next GENEWISE;
  }
    
    # from the matching ones, take the best fit
    
    my @list;
    foreach my $e2g ( @matching_e2gs ){
      # we check exon_overlap and fraction of overlap in gw and e2g:
      my @gw_trans  = @{$gw->get_all_Transcripts};
      my @e2g_trans = @{$e2g->get_all_Transcripts};
      my ($exon_overlap, $extent_overlap) = 
	Bio::EnsEMBL::Pipeline::GeneComparison::TranscriptComparator->_compare_Transcripts( $gw_trans[0], $e2g_trans[0] );
      push (@list, [$exon_overlap, $extent_overlap, $e2g]);
    }

    # sort them
    @list = sort{ 
      # by number of exon overlap
      my $result = ( $$b[0] <=> $$a[0] );
      unless ($result){
	# else, by extent of the overlap
	$result =  ( $$b[1] <=> $$a[1] );
      }
      unless ($result){
	# else, by genomic length of the cdna
	$result = (  $self->_transcript_length_in_gene($$b[2]) <=>
		     $self->_transcript_length_in_gene($$a[2]) );
      }
      unless ($result){
	# else, by exonic extent of the cdna
	$result = (  $self->_transcript_exonic_length_in_gene($$b[2]) <=>
		     $self->_transcript_exonic_length_in_gene($$a[2]) );
      }
    } @list;
    
    #test:
    print STDERR "matching cdnas:\n";
    foreach my $overlap ( @list ){
      print STDERR "cdna: ".$$overlap[2]->dbID.
	", exon_overlap: ".$$overlap[0].
	  ", extent_overlap: ".$$overlap[1].
	    ", extent_UTR: ".$utr_length_hash{$$overlap[2]}."\n";
    }
    
    my $count = 0;
    my $howmany = scalar(@list);
    my $e2g_match;
    do{
      my $this = shift @list;
      $count++;
      $e2g_match = $$this[2];
      
    } while( $used_e2g{$e2g_match} ); 	      
    
    unless ( $e2g_match){
      print STDERR "No cdna found for gw_gene".$gw->dbID." ";
      if ( $howmany == 0 ){
	print STDERR "(no cdna matched this gw gene)\n";
      }
      if ( ($count - 1) == $howmany && $howmany != 0 ){
	print STDERR "(all matching cdnas were already used)\n";
      }

      $self->unmatched_genewise_genes($gw);
      next GENEWISE;
    }
    
    $used_e2g{$e2g_match} = 1;
    
    
    print STDERR "combining gw gene : " . $gw->dbID.":\n";
#    Bio::EnsEMBL::Pipeline::Tools::GeneUtils->_print_Gene($gw);
    print STDERR "with e2g gene " . $e2g_match->dbID . ":\n";
#    Bio::EnsEMBL::Pipeline::Tools::GeneUtils->_print_Gene($e2g_match);
    
    my $combined_transcript = $self->combine_genes($gw, $e2g_match);
    if ( $combined_transcript ){
      $combined_transcript = $self->_transfer_evidence($combined_transcript, $e2g_match);
      $self->make_gene($combined_transcript);
    }
    else{
      print STDERR "no combined gene built from " . $gw->dbID ."\n";
      $self->unmatched_genewise_genes($gw);
      next GENEWISE;
    }
  }
  
  # remap to raw contig coords
  my @remapped = $self->remap_genes();
  $self->output(@remapped);
}

############################################################

=head2

  Description: This method checks the cdnas.             
               ( see Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils::_check_transcript() and _check_introns for more details )
ImportantNote: translation is not checked as these cdnas come without translation
   Arg       : a Bio::EnsEMBL::Gene array
   Return    : a Bio::EnsEMBL::Gene array

=cut

sub _filter_cdnas{
  my ($self,@e2g) = @_;
  my @newe2g;
  
  print STDERR "filtering ".scalar(@e2g)." cdnas\n";
 cDNA_GENE:
  foreach my $e2g (@e2g) {
    
  cDNA_TRANSCRIPT:
    foreach my $tran (@{$e2g->get_all_Transcripts}) {
      
      # rejecting on basis of intron length may not be valid here - it may not be that simple in the same way as it isn;t that simple in Targetted & Similarity builds

      next cDNA_TRANSCRIPT unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($tran,$self->query) && Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_introns($tran,$self->query));

      print STDERR "keeping trans_dbID:" . $tran->dbID . "\n";
      push(@newe2g,$e2g);
    }
  }
  return @newe2g;
}

############################################################
    
# method to calculate the exonic length of a transcript which is inside a gene
# this assumes that these genewise genes contain one transcript each
    
sub _transcript_exonic_length_in_gene{
    my ($self,$gene) = @_;
    my @trans = @{$gene->get_all_Transcripts};
    my $tran = $trans[0];
    my $exonic_length = 0;
    foreach my $exon ( @{$tran->get_all_Exons} ){
	$exonic_length += ($exon->end - $exon->start + 1);
    }
    return $exonic_length;
}

############################################################

# method to calculate the length of a transcript in genomic extent, 
# this assumes that these genewise genes contain one transcript each

sub _transcript_length_in_gene{
    my ($self,$gene) = @_;
    my @trans = @{$gene->get_all_Transcripts};
    my @exons= @{$trans[0]->get_all_Exons};
    my $genomic_extent = 0;
    if ( $exons[0]->strand == -1 ){
      @exons = sort{ $b->start <=> $a->start } @exons;
      $genomic_extent = $exons[0]->end - $exons[$#exons]->start + 1;
    }
    elsif( $exons[0]->strand == 1 ){
      @exons = sort{ $a->start <=> $b->start } @exons;
      $genomic_extent = $exons[$#exons]->end - $exons[0]->start + 1;
    }
    return $genomic_extent;
}

############################################################

sub write_output {
  my($self) = @_;
  
  # write genes in the database: GB_COMB_DBNAME@GB_COMB_DBHOST
  my $gene_adaptor = $self->output_db->get_GeneAdaptor;
  
 GENE: 
  foreach my $gene ($self->output) {	
    if(!$gene->analysis || 
       $gene->analysis->logic_name ne $self->analysis->logic_name){
      $gene->analysis($self->analysis);
    }
    eval {
      $gene_adaptor->store($gene);
      print STDERR "wrote gene dbID " . $gene->dbID . "\n";
#      foreach my $t ( @{$gene->get_all_Transcripts} ){
#	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($t);
#      }
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
#      foreach my $t ( @{$gene->get_all_Transcripts} ){
#	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($t);
#      }
    }
  }
}

############################################################

# make some lovely genes
sub make_gene{
  my ($self,@transcripts) = @_;
  
  # the genetype should be given in Bio::EnsEMBL::Pipeline::GeneConf
  my $genetype = $GB_COMBINED_GENETYPE;
  unless ( $genetype ){
    $self->throw("You must define GB_COMBINED_GENETYPE in Bio::EnsEMBL::Pipeline::GeneConf");
  }
  
  # an analysis should be passed in via the RunnableDB.m parent class:
  my $analysis = $self->analysis;
  unless ($analysis){
    $self->throw("You have to pass an analysis to this RunnableDB through new()");
  }

  my @genes;
  my $count=0;
  
  foreach my $trans(@transcripts){
    $trans->sort;
    
    unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript( $trans ) ){
      print STDERR "rejecting transcript\n";
      return;
    }

    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->add_Transcript($trans);
    $gene->analysis($analysis);

    # do not modify the analysis of the supporting features
    # they should be the original ones: cdna, targetted_genewise or similarity_genewise
    
    if($self->validate_gene($gene)){
      push (@genes,$gene);
      $count++;
    }
  }
  
  print STDERR "Produced genes:\n";
#  foreach my $gene (@genes){
#    Bio::EnsEMBL::Pipeline::Tools::GeneUtils->_print_Gene($gene);
#  }
  
  $self->combined_genes(@genes);
  
}

############################################################

=head2 match_protein_to_cdna

Description: this method tried to find the cdnas that can be merged with the genewise genes.
             Basically it checks for exact internal splice matching in the 5' and 3' exons of the genewise gene.
             In order to match a starting genewise exon with an cDNA exon, we need to have 
             a. exactly coinciding exon boundaries
             b. either cdna exon has start <= genewise exon start, 
             OR cdna exon start lies within $exon_slop bp of genewise exon start AND 
             the cdna transcript will add extra UTR exons. 
             Substitute "end" for "start" for 3prime ends of transcripts  
             BUT do not allow through any cDNA that will result just in a 
             shortened peptide without additional UTR exons.
             
=cut

sub match_protein_to_cdna{
  my ($self, $gw) = @_;

  my %UTR_hash;
  
  print STDERR "\nSearching cDNA for gw gene dbID: ".$gw->dbID."\n";
  my @matching_e2g;
  my @gw_tran = @{$gw->get_all_Transcripts};
  
  my @gw_exons = @{$gw_tran[0]->get_all_Exons};
  my $strand   = $gw_exons[0]->strand;
  if($gw_exons[$#gw_exons]->strand != $strand){
    $self->warn("first and last gw exons have different strands - can't make a sensible combined gene\n with ".$gw_tran[0]->dbId );
      return;
  }
  if ( @gw_exons ){
    if ($strand == 1 ){
      @gw_exons = sort { $a->start <=> $b->start } @gw_exons;
    }
    else{
      @gw_exons = sort { $b->start <=> $a->start } @gw_exons;
    }
  }
  else{
    $self->warn("gw gene without exons: ".$gw->dbID.", skipping it");
    return undef;
  }

  my $exon_slop = 20;
  
 cDNA:
  foreach my $e2g($self->e2g_genes){
    my @egtran  = @{$e2g->get_all_Transcripts};
    my @eg_exons = @{$egtran[0]->get_all_Exons};

    $strand   = $eg_exons[0]->strand; 
    if($eg_exons[$#eg_exons]->strand != $strand){
	$self->warn("first and last e2g exons have different strands - skipping transcript ".$egtran[0]->dbID);
	next cDNA;
    }
    if ($strand == 1 ){
      @eg_exons = sort { $a->start <=> $b->start } @eg_exons;
    }
    else{
      @eg_exons = sort { $b->start <=> $a->start } @eg_exons;
    }
    
    my $fiveprime_match  = 0;
    my $threeprime_match = 0;
    my $left_exon;
    my $right_exon;
    my $left_diff  = 0;
    my $right_diff = 0;
    my $UTR_length;

    
    # Lets deal with single exon genes first
    if ($#gw_exons == 0) {
      foreach my $current_exon (@eg_exons) {
	
	if($current_exon->strand != $gw_exons[0]->strand){
	  next cDNA;
	}
	
	# don't yet deal with genewise leakage for single exon genes
	if ($gw_exons[0]->end   <= $current_exon->end && 
	    $gw_exons[0]->start >= $current_exon->start){
	  $fiveprime_match  = 1;
	  $threeprime_match = 1;
	  
	  $left_exon   = $current_exon;
	  $right_exon  = $current_exon; 
	  $left_diff   = $gw_exons[0]->start - $current_exon->start;
	  $right_diff  = $current_exon->end  - $gw_exons[0]->end;   
	}	
      }
      # can match either end, or both
      if($fiveprime_match && $threeprime_match){
	$UTR_length = $self->_compute_UTRlength($egtran[0],$left_exon,$left_diff,$right_exon,$right_diff);
	$UTR_hash{$e2g} = $UTR_length;
	push(@matching_e2g, $e2g);
      }
      
      # Now the multi exon genewises
    } 
    else {
      
    cDNA_EXONS:
      foreach my $current_exon (@eg_exons) {
	
	if($current_exon->strand != $gw_exons[0]->strand){
	  next cDNA;
	}
	
	if($gw_exons[0]->strand == 1){
	  
	  #FORWARD:

	  # 5prime
	  if ($gw_exons[0]->end == $current_exon->end && 
	      # either e2g exon starts before genewise exon
	      ($current_exon->start <= $gw_exons[0]->start || 
	       # or e2g exon is a bit shorter but there are spliced UTR exons as well
	       (abs($current_exon->start - $gw_exons[0]->start) <= $exon_slop && $current_exon != $eg_exons[0]))){
	    
	    $fiveprime_match = 1;
	    $left_exon = $current_exon;
	    $left_diff = $gw_exons[0]->start - $current_exon->start;
	  }
	  # 3prime
	  elsif($gw_exons[$#gw_exons]->start == $current_exon->start &&
		# either e2g exon ends after genewise exon
		($current_exon->end >= $gw_exons[$#gw_exons]->end ||
		 # or there are UTR exons to be added
		 (abs($current_exon->end - $gw_exons[$#gw_exons]->end) <= $exon_slop && 
		  $current_exon != $eg_exons[$#eg_exons]))){
	    
	    $threeprime_match = 1;
	    $right_exon = $current_exon;
	    $right_diff  = $current_exon->end  - $gw_exons[0]->end;   
	  }
	}
	elsif($gw_exons[0]->strand == -1){
	
	  #REVERSE:

	  # 5prime
	  if ($gw_exons[0]->start == $current_exon->start &&
	      # either e2g exon ends after gw exon
	      ($current_exon->end >= $gw_exons[0]->end ||
	       # or there are UTR exons to be added
	       (abs($current_exon->end - $gw_exons[0]->end) <= $exon_slop &&
		$current_exon != $eg_exons[0]))){
	    print STDERR "fiveprime reverse match\n";
	    
	    $fiveprime_match = 1;
	    $right_exon  = $current_exon;
	    $right_diff  = $current_exon->end  - $gw_exons[0]->end;   
	  }
	  #3prime
	  elsif ($gw_exons[$#gw_exons]->end == $current_exon->end &&
		 # either e2g exon starts before gw exon
		 ($current_exon->start <= $gw_exons[$#gw_exons]->start ||
		  # or there are UTR exons to be added
		  (abs($current_exon->start - $gw_exons[$#gw_exons]->start) <= $exon_slop &&
		   $current_exon != $eg_exons[$#eg_exons]))){
	    print STDERR "threeprime reverse match\n";

	    $threeprime_match = 1;
	    $left_exon = $current_exon;
	    $left_diff = $gw_exons[0]->start - $current_exon->start;
	  }
	}
      }
      if($fiveprime_match || $threeprime_match){
	$UTR_length = $self->_compute_UTRlength($egtran[0],$left_exon,$left_diff,$right_exon,$right_diff);
	$UTR_hash{$e2g} = $UTR_length;
	push(@matching_e2g, $e2g);
	
	# test
	foreach my $egtran ( @{$e2g->get_all_Transcripts} ){
	  print STDERR "Found cDNA match trans_dbID:".$egtran->dbID."\n";
	}
      }
    }
  } 
  return (\@matching_e2g,\%UTR_hash);
}

############################################################

sub _compute_UTRlength{
 my ($self,$transcript,$left_exon, $left_diff,$right_exon,$right_diff) = @_;
 my $strand = $transcript->start_Exon->strand;
 my @exons = sort { $a->start <=> $b->start } @{ $transcript->get_all_Exons };

 my $UTRlength = 0;
 my $in_UTR = 1;

 foreach my $exon ( @exons ){
   if ( defined $left_exon && $exon == $left_exon ){
     $UTRlength += $left_diff;
     $in_UTR     = 0;
   }
   elsif( defiend $right_exons && $exon == $right_exon ){
     $UTRlength += $right_diff;
     $in_UTR     = 1;
   }
   elsif( $in_UTR == 1 ){
     $UTRlength += $exon->length;
   }
 }
 return $UTRlength;

}

############################################################

=head2 _merge_gw_genes

 Function: merges adjacent exons if they are frameshifted; stores component exons as subSeqFeatures of the merged exon

=cut

sub _merge_gw_genes {
  my ($self) = @_;

  my @merged;
  my $count = 1;
  
 GW_GENE:
  foreach my $gwg ($self->gw_genes){
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type('combined');
    $gene->dbID($gwg->dbID);
    my @pred_exons;
    my $ecount = 0;
    
    # order is crucial
    my @trans = @{$gwg->get_all_Transcripts};
    if(scalar(@trans) != 1) { $self->throw("expected one transcript for $gwg\n"); }
    
    # check the sanity of the transcript
    next GW_GENE unless ( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($trans[0],$self->query));
    
    ### we follow here 5' -> 3' orientation ###
    $trans[0]->sort;
   
    #print STDERR "checking for merge: ".$trans[0]->dbID."\n";
    #print STDERR "translation:\n";
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($trans[0]);

    my $cloned_translation = new Bio::EnsEMBL::Translation;
    
    my @gw_exons = @{$trans[0]->get_all_Exons};
    my $strand   = $gw_exons[0]->strand; 
    my $previous_exon;

  EXON:      
    foreach my $exon(@gw_exons){
	
	## genewise frameshift? we merge here two exons separated by max 10 bases into a single exon
	#if ($ecount && $pred_exons[$ecount-1]){
	#  $previous_exon = $pred_exons[$ecount-1];
	#}
	
	$ecount++;
	
	my $separation = 0;
	my $merge_it   = 0;
	
	## we put every exon following a frameshift into the first exon before the frameshift
	## following the ordering 5' -> 3'
	if (defined($previous_exon)){
	    
	    #print STDERR "previous exon: ".$previous_exon->start."-".$previous_exon->end."\n";
	    #print STDERR "current exon : ".$exon->start."-".$exon->end."\n";
	    
	    if ($strand == 1){
		$separation = $exon->start - $previous_exon->end - 1;
	    }
	    elsif( $strand == -1 ){
		$separation = $previous_exon->end - $exon->start - 1;
	    }
	    if ($separation <=10){
		$merge_it = 1;
	    }	
	}
	
	if ( defined($previous_exon) && $merge_it == 1){
	    
	    # combine the two
	    
	    # the first exon (5'->3' orientation always) is the containing exon,
	    # which gets expanded and the other exons are added into it
	    print STDERR "merging $exon into $previous_exon\n";
	    print STDERR $exon->start."-".$exon->end." into ".$previous_exon->start."-".$previous_exon->end."\n";
	    
	    $previous_exon->end($exon->end);
	    $previous_exon->add_sub_SeqFeature($exon,'');
	    
	    # if this is end of translation, keep that info:
	    if ( $exon == $trans[0]->translation->end_Exon ){
		$cloned_translation->end_Exon( $previous_exon );
		$cloned_translation->end($trans[0]->translation->end);
	    }
	    
	    my %evidence_hash;
	    foreach my $sf( @{$exon->get_all_supporting_features}){
		if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
		    next;
		}
		#print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
		$evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
		$previous_exon->add_supporting_features($sf);
	    }
	    next EXON;
	} 
	else{
	    # make a new Exon - clone $exon
	    my $cloned_exon = Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_clone_Exon($exon);
	    $cloned_exon->add_sub_SeqFeature($exon,'');
	    
	    # if this is start/end of translation, keep that info:
	    if ( $exon == $trans[0]->translation->start_Exon ){
		$cloned_translation->start_Exon( $cloned_exon );
		$cloned_translation->start($trans[0]->translation->start);
	    }
	    if ( $exon == $trans[0]->translation->end_Exon ){
		$cloned_translation->end_Exon( $cloned_exon );
		$cloned_translation->end($trans[0]->translation->end);
	    }
	    
	    $previous_exon = $cloned_exon;
	    push(@pred_exons, $cloned_exon);
	}
    }
    
    # transcript
    my $merged_transcript   = new Bio::EnsEMBL::Transcript;
    $merged_transcript->dbID($trans[0]->dbID);
    foreach my $pe(@pred_exons){
	$merged_transcript->add_Exon($pe);
    }
    
    $merged_transcript->sort;
    $merged_transcript->translation($cloned_translation);
    
    #print STDERR "merged_transcript:\n";
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($merged_transcript);

    # and gene
    $gene->add_Transcript($merged_transcript);
    push(@merged, $gene);
    $count++;

    # store match between merged and original gene so we can easily retrieve the latter if we need to
    $self->merged_unmerged_pairs($gene,$gwg);
  } # end GW_GENE
  
  return @merged;
}

############################################################

sub combine_genes{
  my ($self, $gw, $e2g) = @_;

  my $modified_peptide = 0;
  my @combined_transcripts  = ();
  
  # should be only 1 transcript
  my @gw_tran  = @{$gw->get_all_Transcripts};
  #  $gw_tran[0]->sort;
  my @gw_exons = @{$gw_tran[0]->get_all_Exons}; # ordered array of exons
  my @egtran = @{$e2g->get_all_Transcripts};
  #  $egtran[0]->sort;
  my @e2g_exons  = @{$egtran[0]->get_all_Exons}; # ordered array of exons
  
  # OK, let's see if we need a new gene
  # base it on the existing genewise one
  my $newtranscript = new Bio::EnsEMBL::Transcript;
  foreach my $exon(@gw_exons){
    $newtranscript->add_Exon($exon);
  }
  
  my $translation   = new Bio::EnsEMBL::Translation;
  $translation->start($gw_tran[0]->translation->start);
  $translation->end($gw_tran[0]->translation->end);
  $translation->start_Exon($gw_tran[0]->translation->start_Exon);
  $translation->end_Exon($gw_tran[0]->translation->end_Exon);
  
  $newtranscript->translation($translation);
  $newtranscript->translation->start_Exon($newtranscript->start_Exon);
  $newtranscript->translation->end_Exon($newtranscript->end_Exon);
  
  my $eecount = 0;
  my $modified_peptide_flag;
 
 EACH_E2G_EXON:
  foreach my $ee (@e2g_exons){
      
      # check strands are consistent
      if ($ee->strand != $gw_exons[0]->strand){
	  $self->warn("gw and e2g exons have different strands - can't combine genes\n") ;
	  next GENEWISE;
      }
      
      # single exon genewise prediction?
      if(scalar(@gw_exons) == 1) {
	  
	  ($newtranscript,$modified_peptide_flag) = $self->transcript_from_single_exon_genewise( $ee, 
												 $gw_exons[0], 
												 $newtranscript, 
												 $translation, 
												 $eecount, 
												 @e2g_exons);
      }
      
      else {
	  ($newtranscript,$modified_peptide_flag) = $self->transcript_from_multi_exon_genewise($ee, 
											       $newtranscript, 
											       $translation, 
											       $eecount,
											       $gw, 
											       $e2g)
	  }
      if ( $modified_peptide_flag ){
	  $modified_peptide = 1;
      }
      
      # increment the exon
      $eecount++;
      
  } # end of EACH_E2G_EXON
  
  
  ##############################
  # expand merged exons
  ##############################
  
  # the new transcript is made from a merged genewise gene
  # check the transcript and expand frameshifts in all but original 3' gw_exon
  # (the sub_SeqFeatures have been flushed for this exon)
  if (defined($newtranscript)){
      
      # test
      #print STDERR "before expanding exons, newtranscript: $newtranscript\n"; 
      #$self->_print_Transcript($newtranscript);
    
      foreach my $ex (@{$newtranscript->get_all_Exons}){
      
      if(scalar($ex->sub_SeqFeature) > 1 ){
	my @sf    = $ex->sub_SeqFeature;
	my $first = shift(@sf);
	
	$ex->end($first->end);
	
	# add back the remaining component exons
	foreach my $s(@sf){
	    $newtranscript->add_Exon($s);
	    $newtranscript->sort;
	}
	# flush the sub_SeqFeatures
	$ex->flush_sub_SeqFeature;
    }
  }
    
    #unless($self->compare_translations($gw_tran[0], $newtranscript) ){
    #  print STDERR "translation has been modified\n";
    #}


    # check that the result is fine
    unless( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Transcript($newtranscript,$self->query) ){
	print STDERR "problems with this combined transcript, return undef";
	return undef;
    }
    unless( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($newtranscript) ){
	print STDERR "problems with this combined translation, return undef";
	return undef;
    }
      
    # check translation is the same as for the genewise gene we built from
    my $foundtrans = 0;
    
    # the genewise translation can be modified due to a disagreement in a
    # splice site with cdnas. This can happen as neither blast nor genewise can
    # always find very tiny exons.
    # we then recalculate the translation using genomewise:
    
    my $newtrans;
    if ( $modified_peptide ){
      my $strand = $newtranscript->start_Exon->strand;
      
      print STDERR "before genomewise:\n";
#      Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($newtranscript);
      
      $newtrans = $self->_recalculate_translation($newtranscript,$strand); 
      
      print STDERR "after genomewise:\n";
#    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($newtrans);
      
#      unless($self->compare_translations($gw_tran[0], $newtrans) ){
#      	print STDERR "translation has been modified\n";
#      }
      
      # if the genomewise results gets stop codons, return the original transcript:
      unless( Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($newtrans) ){
	print STDERR "Arrgh, stop codons, returning the original transcript\n";
	$newtrans = $newtranscript;
      }
    }
    else{
      $newtrans = $newtranscript;
    }
    return $newtrans;
  }
  else{
      $self->warn("No combination could be built\n");
      return undef;
  }
}

############################################################

sub transcript_from_single_exon_genewise {
    my ($self, $eg_exon, $gw_exon, $transcript, $translation, $exoncount, @e2g_exons) = @_;
    
    # save out current translation end - we will need this if we have to unmerge frameshifted exons later
    my $orig_tend = $translation->end;
    
    # stay with being strict about gw vs e2g coords - may change this later ...
    # the overlapping e2g exon must at least cover the entire gw_exon
    if ($gw_exon->start >= $eg_exon->start && $gw_exon->end <= $eg_exon->end){
	
	my $egstart = $eg_exon->start;
	my $egend   = $eg_exon->end;
	my $gwstart = $gw_exon->start;
	my $gwend   = $gw_exon->end;
	
	# print STDERR "single exon gene, " . $gw_exon->strand  .  " strand\n";	    
	
	# modify the coordinates of the first exon in $newtranscript
	my $ex = $transcript->start_Exon;
	
	$ex->start($eg_exon->start);
	$ex->end($eg_exon->end);
	
	# need to explicitly set the translation start & end exons here.
	$translation->start_Exon($ex);
	
	# end_exon may be adjusted by 3' coding exon frameshift expansion. Ouch.

	$translation->end_Exon($ex);

	# need to deal with translation start and end this time - varies depending on strand
	
	############################################################
	#FORWARD:
	if($gw_exon->strand == 1){
	    my $diff = $gwstart - $egstart;
	    my $tstart = $translation->start;
	    my $tend = $translation->end;
	    
	    #print STDERR "diff: ".$diff." translation start : ".$tstart." end: ".$tend."\n";
	    #print STDERR "setting new translation to start: ".($tstart+$diff)." end: ".($tend+$diff)."\n";
	    $translation->start($tstart + $diff);
	    $translation->end($tend + $diff);
	    
	    $self->throw("Forward strand: setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
	    $self->throw("Forward strand: setting dodgy translation end: " . $translation->end . " exon_length: " . $translation->end_Exon->length . "\n") unless $translation->end <= $translation->end_Exon->length;
	  }
	############################################################
	#REVERSE:
	elsif($gw_exon->strand == -1){
	  my $diff   = $egend - $gwend;
	  my $tstart = $translation->start;
	  my $tend   = $translation->end;
	  $translation->start($tstart+$diff);
	  $translation->end($tend + $diff);
	  
	  $self->throw("Reverse strand: setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
	  $self->throw("Reverse strand: setting dodgy translation end: " . $translation->end . " exon_length: " . $translation->end_Exon->length . "\n") unless $translation->end <= $translation->end_Exon->length;
	}
	
	# expand frameshifted single exon genewises back from one exon to multiple exons
	if(scalar($ex->sub_SeqFeature) > 1){
	  print STDERR "frameshift in a single exon genewise\n";
	  my @sf = $ex->sub_SeqFeature;
	  
	  # save current start and end of modified exon
	  my $cstart = $ex->start;
	  my $cend   = $ex->end;
	  my $exlength = $ex->length;
	  
	  # get first exon - this has same id as $ex
	    my $first = shift(@sf);
	  $ex->end($first->end); # NB end has changed!!!
	  # don't touch start. 
	  # no need to modify translation start
	  
	    # get last exon
	  my $last = pop(@sf);
	  $last->end($cend);
	  $transcript->add_Exon($last);
	  # and adjust translation end - the end is still relative to the merged gw exon
	  $translation->end_Exon($last);
	  # put back the original end translation
	  $translation->end($orig_tend);
	  
	  # get any remaining exons
	  foreach my $s(@sf){
	    $transcript->add_Exon($s);
	    $transcript->sort;
	  }
	  # flush the sub_SeqFeatures
	  $ex->flush_sub_SeqFeature;
	}
	
	# need to add back exons, both 5' and 3'
	$self->add_5prime_exons($transcript, $exoncount, @e2g_exons);
	$self->add_3prime_exons($transcript, $exoncount, @e2g_exons);
	
      }
    return ($transcript,0);
  }

# this method will actually do the combination of both cdna and genewise gene.
# Note that if there is a match on one end but not on the other, the
# code will extend one end, but will leave the other as it is in the
# genewise genes. This will explit cdna matches that look fine on one end
# and we disregard the mismatching part.

############################################################

sub transcript_from_multi_exon_genewise {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;
  
  # $current_exon is the exon one the e2g_transcript we are in at the moment
  # $exoncount is the position of the e2g exon in the array

  # save out current translation->end - we'll need it if we have to expand 3prime exon later
 # my $orig_tend = $translation->end;

  my @gwtran  = @{$gw_gene->get_all_Transcripts};
  $gwtran[0]->sort;
  my @gwexons = @{$gwtran[0]->get_all_Exons};
  
  my @egtran  = @{$eg_gene->get_all_Transcripts};
  $egtran[0]->sort;
  my @egexons = @{$egtran[0]->get_all_Exons};
  
  # in order to match a starting genewise exon with an e2g exon, we need to have 
  # a. exactly coinciding exon ends
  # b. exon starts lying within $exon_slop bp of each other. 
  # previously we had required e2g start to be strictly <= gw start, but this will lose us some valid UTRs
  # substitute "end" for "start" for 3' ends of transcripts
  
  # compare to the first genewise exon
  
  if($gwexons[0]->strand == 1){
    return $self->transcript_from_multi_exon_genewise_forward($current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene);
  }
  elsif( $gwexons[0]->strand == -1 ){
    return $self->transcript_from_multi_exon_genewise_reverse($current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene);
  }
}

############################################################

sub transcript_from_multi_exon_genewise_forward{
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;

  my $modified_peptide = 0;

  my @gwtran  = @{$gw_gene->get_all_Transcripts};
  $gwtran[0]->sort;
  my @gwexons = @{$gwtran[0]->get_all_Exons};
  
  my @egtran  = @{$eg_gene->get_all_Transcripts};
  $egtran[0]->sort;
  my @egexons = @{$egtran[0]->get_all_Exons};
  
  # save out current translation->end - we'll need it if we have to expand 3prime exon later
  my $orig_tend = $translation->end;

  ###################
  my $exon_slop = 20;
  ###################
   
  ############### 5_PRIME:
  if (#they have a coincident end
      $gwexons[0]->end == $current_exon->end && 
      
      # either e2g exon starts before genewise exon
      ($current_exon->start <= $gwexons[0]->start || 
       
       # or e2g exon is a bit shorter but there are spliced UTR exons as well
       (abs($current_exon->start - $gwexons[0]->start) <= $exon_slop && 
	$current_exon != $egexons[0]))){
    

    
    my $current_start = $current_exon->start;
    my $gwstart       = $gwexons[0]->start;
    
    # this exon will be the start of translation, convention: phase = -1
    my $ex = $transcript->start_Exon;
    $ex->phase(-1);
    
    # modify the coordinates of the first exon in $newtranscript if
    # e2g is larger on this end than gw.
    if ( $current_exon->start < $gwexons[0]->start ){
      $ex->start($current_exon->start);
    }
    elsif( $current_exon->start == $gwexons[0]->start ){
      $ex->start($gwstart);
      $ex->phase($gwexons[0]->phase);
    }
    # if the e2g exon starts after the gw exon, 
    # modify the start only if this e2g exon is not the first of the transcript
    elsif(  $current_start > $gwstart && $exoncount != 0 ) {
      $ex->start($current_exon->start);
    }
    
    # add all the exons from the est2genome transcript, previous to this one
    Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence($current_exon, $ex);
    $self->add_5prime_exons($transcript, $exoncount, @egexons);
    
    # fix translation start 
    if($gwstart >= $current_start){
      # take what it was for the gw gene, and add on the extra
      my $tstart = $translation->start;
      #print STDERR "Forward 5': original translation start: $tstart ";
      $tstart += ($gwstart - $current_start);
      $translation->start($tstart);
      #print STDERR "re-setting translation start to: $tstart\n";
    }
    ############################################################
    # only trust a smaller cdna exon if it is not the first of the transcript
    # (it could be a truncated cdna)
    elsif($gwstart < $current_start && $exoncount != 0){
      
      $modified_peptide = 1;
      print STDERR "SHORTENING GENEWISE TRANSLATION\n";
      # genewise has leaked over the start. Tougher call - we need to take into account the 
      # frame here as well
      print STDERR "gw exon starts: $gwstart < new start: $current_start\n";
      print STDERR "modifying exon, as cdna exon is not the first of transcript-> exoncount = $exoncount\n";
      
      # $diff is the number of bases we chop from the genewise exon
      my $diff   = $current_start - $gwstart;
      my $tstart = $translation->start;
      $self->warn("this is a case where gw translation starts at $tstart > 1") if ($tstart>1);
      print STDERR "gw translation start: ".$tstart."\n";
      print STDERR "start_exon: ".$translation->start_Exon->start.
	"-".$translation->start_Exon->end.
	  " length: ".($translation->start_Exon->end - $translation->start_Exon->start + 1).
	    " phase: ".$translation->start_Exon->phase.
	      " end_phase: ".$translation->start_Exon->end_phase."\n";
      
      
      if($diff % 3 == 0) { 
	# we chop exactily N codons from the beginning of translation
	$translation->start(1); 
	}
      elsif ($diff % 3 == 1) { 
	# we chop N codons plus one base 
	$translation->start(3); 
      }
      elsif ($diff % 3 == 2) { 
	# we chop N codons plus 2 bases
	$translation->start(2); 
      }
      else {
	$translation->start(1);
	$self->warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");
      }
    }
    
    ############################################################
    
    else{
      print STDERR "gw exon starts: $gwstart > new start: $current_start";
      print STDERR "but cdna exon is the first of transcript-> exoncount = $exoncount, so we don't modify it\n";
    }
    $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
    
  } # end 5' exon
  
  ############### 3_PRIME:
  elsif (# they have coincident start
	 $gwexons[$#gwexons]->start == $current_exon->start && 
	 
	 # either e2g exon ends after genewise exon
	 ($current_exon->end >= $gwexons[$#gwexons]->end ||
	  
	    # or we allow to end before if there are UTR exons to be added
	    (abs($current_exon->end - $gwexons[$#gwexons]->end) <= $exon_slop && 
	     $current_exon != $egexons[$#egexons]))){   
      
      my $end_translation_shift = 0;
      
      # modify the coordinates of the last exon in $newtranscript
      # e2g is larger on this end than gw.
      my $ex = $transcript->end_Exon;
      
      # this exon is the end of translation, convention: end_phase = -1
      $ex->end_phase(-1);
      
      if ( $current_exon->end > $gwexons[$#gwexons]->end ){
	$ex->end($current_exon->end);
      }
      elsif( $current_exon->end == $gwexons[$#gwexons]->end ){
	$ex->end($gwexons[$#gwexons]->end);
	$ex->end_phase($gwexons[$#gwexons]->end_phase);
      }
      # if the e2g exon ends before the gw exon, 
      # modify the end only if this e2g exon is not the last of the transcript
      elsif ( $current_exon->end < $gwexons[$#gwexons]->end && $exoncount != $#egexons ){
	
	$modified_peptide = 1;
	print STDERR "SHORTENING GENEWISE TRANSLATION\n";
	  ## fix translation end iff genewise has leaked over - will need truncating
	  my $diff   = $gwexons[$#gwexons]->end - $current_exon->end;
	  print STDERR "diff: $diff\n";
	  my $tend   = $translation->end;
	  
	  my $gw_exon_length   = $gwexons[$#gwexons]->end - $gwexons[$#gwexons]->start + 1;
	  my $cdna_exon_length = $current_exon->end - $current_exon->start + 1;
	  print STDERR "gw exon length  : $gw_exon_length\n";
	  print STDERR "cdna exon length: $cdna_exon_length\n";
	  
	  my $length_diff = $gw_exon_length - $cdna_exon_length;
	  print STDERR "length diff: ".$length_diff."\n"; # should be == diff
	  
	  $ex->end($current_exon->end);
	  
	  if($diff % 3 == 0) { 
	    # we chop exactily N codons from the end of the translation
	    # so it can end where the cdna exon ends
	    $translation->end($cdna_exon_length);
	    $end_translation_shift = $length_diff;
	    
	  }
	  elsif ($diff % 3 == 1) { 
	    # we chop N codons plus one base 
	    # it should end on a full codon, so we need to end translation 2 bases earlier:
	    $translation->end($cdna_exon_length - 2); 
	    $end_translation_shift = $length_diff + 2;
	  }
	  elsif ($diff % 3 == 2) { 
	    # we chop N codons plus 2 bases
	    # it should end on a full codon, so we need to end translation 1 bases earlier:
	    $translation->end($cdna_exon_length - 1); 
	    $end_translation_shift = $length_diff + 1;
	  }
	  else {
	    # absolute genebuild paranoia 8-)
	    $translation->end($cdna_exon_length);
	    $self->warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");
	  }
	  print STDERR "Forward: translation end set to : ".$translation->end."\n";
      
      }
      # need to explicitly set the translation end exon for translation to work out
      my $end_ex = $transcript->end_Exon;
      $translation->end_Exon($end_ex);
      
      # strand = 1
      my $expanded = $self->expand_3prime_exon($ex, $transcript, 1);
      
      if($expanded){
	# set translation end to what it originally was in the unmerged genewise gene
	# taking into account the diff
	print STDERR "Forward: expanded 3' exon, re-setting end of translation from ".$translation->end." to orig_end ($orig_tend)- ( length_diff + shift_due_to_phases ) ($end_translation_shift)".($orig_tend - $end_translation_shift)."\n";
	$translation->end($orig_tend - $end_translation_shift);
      }
      
      
      # finally add any 3 prime e2g exons
      Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence($current_exon, $ex);
      $self->add_3prime_exons($transcript, $exoncount, @egexons);
      
    } # end 3' exon    
    
  return ($transcript,$modified_peptide);
}

##################################################

sub transcript_from_multi_exon_genewise_reverse{
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;
  
  my $modified_peptide = 0;

  my @gwtran  = @{$gw_gene->get_all_Transcripts};
  $gwtran[0]->sort;
  my @gwexons = @{$gwtran[0]->get_all_Exons};
  
  my @egtran  = @{$eg_gene->get_all_Transcripts};
  $egtran[0]->sort;
  my @egexons = @{$egtran[0]->get_all_Exons};
  
  # save out current translation->end - we'll need it if we have to expand 3prime exon later
  my $orig_tend = $translation->end;


  ###################
  my $exon_slop = 20;
  ###################
  
  ####################### 5_PRIME:
  if ($gwexons[0]->start == $current_exon->start && 
      # either e2g exon ends after gw exon
      ($current_exon->end >= $gwexons[0]->end ||
       # or there are UTR exons to be added
       (abs($current_exon->end - $gwexons[0]->end) <= $exon_slop &&
	$current_exon != $egexons[0]))){
    
    # sort out translation start
    my $tstart = $translation->start;
    if($current_exon->end >= $gwexons[0]->end){
      # take what it was for the gw gene, and add on the extra
      $tstart += $current_exon->end - $gwexons[0]->end;
      $translation->start($tstart);
    }
    elsif( $current_exon->end < $gwexons[0]->end && $current_exon != $egexons[0] ){
      # genewise has leaked over the start. Tougher call - we need to take into account the 
      # frame here as well
      $modified_peptide = 1;
      print STDERR "SHORTENING GENEWISE TRANSLATION\n";
      print STDERR "In Reverse strand. gw exon ends: ".$gwexons[0]->end." > cdna exon end: ".$current_exon->end."\n";
      print STDERR "modifying exon, as cdna exon is not the first of transcript-> exoncount = $exoncount\n";
      
      
      my $diff = $gwexons[0]->end - $current_exon->end;
      my $gwstart = $gwexons[0]->end;
      my $current_start = $current_exon->end;
      my $tstart = $translation->start;
      
      #print STDERR "before the modification:\n";
      #print STDERR "start_exon: ".$translation->start_Exon."\n";
      #print STDERR "gw translation start: ".$tstart."\n";
      #print STDERR "start_exon: ".$translation->start_Exon->start.
      #  "-".$translation->start_Exon->end.
      #    " length: ".($translation->start_Exon->end - $translation->start_Exon->start + 1).
      #      " phase: ".$translation->start_Exon->phase.
      #	" end_phase: ".$translation->start_Exon->end_phase."\n";
      
      
      if    ($diff % 3 == 0) { $translation->start(1); }
      elsif ($diff % 3 == 1) { $translation->start(3); }
      elsif ($diff % 3 == 2) { $translation->start(2); }
      else {
	$translation->start(1);
	$self->warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");}
    }
    
    
    $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
    
    
    # this exon is the start of translation, convention: phase = -1
    my $ex = $transcript->start_Exon;
    $ex->phase(-1);
    
    # modify the coordinates of the first exon in $newtranscript
    if ( $current_exon->end > $gwexons[0]->end){ 
      $ex->end($current_exon->end);
      $ex->phase(-1);
    }
    elsif (  $current_exon->end == $gwexons[0]->end){
      $ex->end($gwexons[0]->end);
      $ex->phase($gwexons[0]->phase);
    }
    elsif (  $current_exon->end < $gwexons[0]->end && $current_exon != $egexons[0] ){
      $ex->end($current_exon->end);
    }
    
    # need to explicitly set the translation start exon for translation to work out
    $translation->start_Exon($ex);
    
    #print STDERR "after the modification:\n";
    #print STDERR "start_Exon: ".$translation->start_Exon."\n";
    #print STDERR "gw translation start: ".$translation->start."\n";
    #print STDERR "start_exon: ".$translation->start_Exon->start.
    #	"-".$translation->start_Exon->end.
    #  " length: ".($translation->start_Exon->end - $translation->start_Exon->start + 1).
    #  " phase: ".$translation->start_Exon->phase.
    #      " end_phase: ".$translation->start_Exon->end_phase."\n";
    
    Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence($current_exon, $ex);
    $self->add_5prime_exons($transcript, $exoncount, @egexons);
    
  }
  # end 5' exon
  
  ###################### 3_PRIME:
  elsif ($gwexons[$#gwexons]->end == $current_exon->end && 
	 # either e2g exon starts before gw exon
	 ($current_exon->start <= $gwexons[$#gwexons]->start ||
	  # or there are UTR exons to be added
	  (abs($current_exon->start - $gwexons[$#gwexons]->start) <= $exon_slop &&
	   $current_exon != $egexons[$#egexons]))){
    
    
    my $end_translation_shift = 0;
    
    # this exon is the end of translation, convention: end_phase = -1
    my $ex = $transcript->end_Exon;
      $ex->end_phase(-1);
    
    # modify the coordinates of the last exon in $newtranscript
    if ( $current_exon->start < $gwexons[$#gwexons]->start ){
      # no need to modify translation->end as the 'end' of this exon has not changed
      $ex->start($current_exon->start);
      $ex->end_phase(-1);
    }
    elsif( $current_exon->start == $gwexons[$#gwexons]->start){
      $ex->start($gwexons[$#gwexons]->start);
      $ex->end_phase($gwexons[$#gwexons]->end_phase);
    }
    # if the e2g exon starts after the gw exon, 
    # modify the end only if this e2g exon is not the last of the transcript
    elsif ( $current_exon->start > $gwexons[$#gwexons]->start && $exoncount != $#egexons ){

      $modified_peptide = 1;
      print STDERR "SHORTENING GENEWISE TRANSLATION\n";
      print STDERR "In Reverse strand: gw exon start: ".$gwexons[$#gwexons]->start." < cdna exon start: ".$current_exon->start."\n";
      print STDERR "modifying exon, as cdna exon is not the last of transcript-> exoncount = $exoncount, and #egexons = $#egexons\n";
      
	## adjust translation
	my $diff   = $current_exon->start - $gwexons[$#gwexons]->start;
	print STDERR "diff: $diff\n";
	my $tend   = $translation->end; 
	
	my $gw_exon_length   = $gwexons[$#gwexons]->end - $gwexons[$#gwexons]->start + 1;
	my $cdna_exon_length = $current_exon->end - $current_exon->start + 1;
	print STDERR "gw exon length  : $gw_exon_length\n";
	print STDERR "cdna exon length: $cdna_exon_length\n";
	
	my $length_diff = $gw_exon_length - $cdna_exon_length;

	# modify the combined exon coordinate to be that of the cdna
	$ex->start($current_exon->start);

	if($diff % 3 == 0) { 
	  # we chop exactily N codons from the end of the translation
	  # so it can end where the cdna exon ends
	  $translation->end($cdna_exon_length); 
	  $end_translation_shift = $length_diff;
	}
	elsif ($diff % 3 == 1) { 
	  # we chop N codons plus one base 
	  # it should end on a full codon, so we need to end translation 2 bases earlier:
	  $translation->end($cdna_exon_length - 2);
	  $end_translation_shift = $length_diff + 2;
	}
	elsif ($diff % 3 == 2) { 
	  # we chop N codons plus 2 bases
	  # it should end on a full codon, so we need to end translation 1 bases earlier:
	  $translation->end($cdna_exon_length - 1);
	  $end_translation_shift = $length_diff + 1;
	}
	else {
	  # absolute genebuild paranoia 8-)
	  $translation->end($cdna_exon_length);
	  $self->warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");
	}
	#print STDERR "Reverse: translation end set to : ".$translation->end."\n";
      
      }	
      # strand = -1
      my $expanded = $self->expand_3prime_exon($ex, $transcript,-1);
      
      # need to explicitly set the translation end exon for translation to work out
      my $end_ex = $transcript->end_Exon;
      $translation->end_Exon($end_ex);
      
      if($expanded){
	# set translation end to what it originally was in the unmerged genewise gene
	print STDERR "Reverse: expanded 3' exon, re-setting translation exon ".$translation->end." to original end( $orig_tend ) - shifts_due_to_phases_etc ( $end_translation_shift ) :".($orig_tend - $end_translation_shift)."\n";
	$translation->end($orig_tend - $end_translation_shift);
      }
      Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_transfer_supporting_evidence($current_exon, $ex);
      $self->add_3prime_exons($transcript, $exoncount, @egexons);
      
    } # end 3' exon
  
  return ($transcript,$modified_peptide); 
}

############################################################

sub add_5prime_exons{
my ($self, $transcript, $exoncount, @e2g_exons) = @_;

      # add all the exons from the est2genome transcript, previous to this one
      # db handle will be screwed up, need to mak new exons from these
      my $c = 0;
      my $modified = 0;
      while($c < $exoncount){
	my $newexon = new Bio::EnsEMBL::Exon;
	my $oldexon = $e2g_exons[$c];
	$newexon->start($oldexon->start);
	$newexon->end($oldexon->end);
	$newexon->strand($oldexon->strand);

	# these are all 5prime UTR exons
	$newexon->phase(-1);
	$newexon->end_phase(-1);
	$newexon->contig($oldexon->contig);
	my %evidence_hash;
	#print STDERR "adding evidence at 5':\n";
	foreach my $sf( @{$oldexon->get_all_supporting_features} ){
	  if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	    next;
	  }
	  $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	  #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	  $newexon->add_supporting_features($sf);
	}
#	print STDERR "Adding 5prime exon " . $newexon->start . " " . $newexon->end . "\n";
	$transcript->add_Exon($newexon);
	$modified = 1;
	$transcript->sort;
	$c++;
      }
      if ($modified == 1){
       $transcript->translation->start_Exon->phase(-1);
      }
      
}

# $exon is the terminal exon in the genewise transcript, $transcript. We need
# to expand any frameshifts we merged in the terminal genewise exon. 
# The expansion is made by putting $exon to be the last (3' end) component, so we modify its
# start but not its end. The rest of the components are added. The translation end will have to be modified,
# this happens in the method _transcript_from_multi_exon....

############################################################

sub expand_3prime_exon{
  my ($self, $exon, $transcript, $strand) = @_;
  
  if(scalar($exon->sub_SeqFeature) > 1){
    print STDERR "expanding 3'prime frameshifted exon $exon in strand $strand: ".
      $exon->start."-".$exon->end." phase: ".$exon->phase." end_phase: ".$exon->end_phase."\n";
    my @sf = $exon->sub_SeqFeature;
    
    my $last = pop(@sf);
     #print STDERR "last component: ".$last->start."-".$last->end." phase ".$last->phase." end_phase ".$last->end_phase."\n";
    
    #print STDERR "setting exon $exon start: ".$last->start." phase: ".$last->phase."\n";  
    $exon->start($last->start); # but don't you dare touch the end!
    $exon->dbID($last->dbID);
    $exon->phase($last->phase);
    
    # add back the remaining component exons
    foreach my $s(@sf){
      #print STDERR "adding exon: ".$s->start."-".$s->end."\n";
      $transcript->add_Exon($s);
      $transcript->sort;
    }
    # flush the sub_SeqFeatures so we don't try to re-expand later
    $exon->flush_sub_SeqFeature;
    return 1;
  }
  
  # else, no expansion
  return 0;
}


############################################################

# $exoncount tells us which position in the array 
# of e2g exons corresponds to the end of the genewise transcript so we add back 
# exons 3' to that position.
# $exon and $transcript are references to Exon and Transcript objects.

sub add_3prime_exons {
my ($self, $transcript, $exoncount, @e2g_exons) = @_;
# need to deal with frameshifts - 3' exon is a special case as its end might have changed

      # add all the exons from the est2genome transcript, subsequent to this one
      my $c = $#e2g_exons;
      my $modified = 0;
      while($c > $exoncount){
	my $newexon = new Bio::EnsEMBL::Exon;
	my $oldexon = $e2g_exons[$c];
	$newexon->start($oldexon->start);
	$newexon->end($oldexon->end);
	$newexon->strand($oldexon->strand);
	
	# these are all exons with UTR:
	$newexon->phase(-1);
	$newexon->end_phase(-1);
	$newexon->contig($oldexon->contig);
	#print STDERR "adding evidence in 3':\n";
	my %evidence_hash;
	foreach my $sf( @{$oldexon->get_all_supporting_features }){
	  if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	    next;
	  }
	  $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	  #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	  $newexon->add_supporting_features($sf);
	}
	#	print STDERR "Adding 3prime exon " . $newexon->start . " " . $newexon->end . "\n";
	$transcript->add_Exon($newexon);
	$modified = 1;
        $transcript->sort;
	$c--;
      }
      if ($modified == 1){
	$transcript->translation->end_Exon->end_phase(-1);
      }

}
############################################################

=head2 compare_translations

  Description: it compares the peptides from two transcripts.
               It returns 1 if both translations are the same or when
               one is a truncated version of the other. It returns 0 otherwise.
               Also returns 0 if either of the transcripts has stops.
  ReturnType : BOOLEAN

=cut

sub compare_translations{
  my ($self, $genewise_transcript, $combined_transcript) = @_;
    
  my $seqout = new Bio::SeqIO->new(-fh => \*STDERR);
  
  my $genewise_translation;
  my $combined_translation;
  
  eval {
      $genewise_translation = $genewise_transcript->translate;
  };
  if ($@) {
    print STDERR "Couldn't translate genewise gene\n";
  }
  else{
    print STDERR "genewise: \n";             
    $seqout->write_seq($genewise_translation);
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Translation($genewise_transcripts[0]);
  }
  
  $@ = '';
  
  eval{
    $combined_translation = $combined_transcript->translate;
  };
  
  if ($@) {
    print STDERR "Couldn't translate combined gene:[$@]\n";
    return 0;
  }
  else{
    print STDERR "combined: \n";             
    $seqout->write_seq($combined_translation);
  }	 
  
  if ( $genewise_translation && $combined_translation ){
      my $gwseq  = $genewise_translation->seq;
      my $comseq = $combined_translation->seq;
      
      if($gwseq eq $comseq) {
	  print STDERR "combined translation is identical to genewise translation\n";
	  return 1;
      }
      elsif($gwseq =~ /$comseq/){
	  print STDERR "combined translation is a truncated version of genewise translation\n";
	  return 1;
      }
      elsif($comseq =~ /$gwseq/){
	  print STDERR "interesting: genewise translation is a truncated version of the combined-gene translation\n";
	  return 1;
      }
  }
  
  return 0;
}

############################################################

=head2 remap_genes

=cut

sub remap_genes {
  my ($self) = @_;
  my @newf;  
  my $contig = $self->query;
  
  my @genes = $self->combined_genes;
  push(@genes, $self->unmatched_genewise_genes);

  my $genecount = 0;
 GENE:  
  foreach my $gene (@genes) {
    $genecount++;
    my @t = @{$gene->get_all_Transcripts};
    my $tran = $t[0];
    
    # check that it translates
    unless(Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($tran)){
      print STDERR "rejecting gene\n";
      next GENE;
    }
    
    #print STDERR "**************about to remap:**********************\n";
    my $transcount = 0;
    foreach my $transcript ( @{$gene->get_all_Transcripts} ){
      $transcount++;
      $transcript->type( $genecount."_".$transcount );

      # set start and stop codons
      $transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_start_codon($transcript);
      $transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon($transcript);
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
    }

    eval {
      $gene->transform;	
    };
    
    #print STDERR "****************after remapping**************\n";
    foreach my $t ( @{$gene->get_all_Transcripts} ){
      #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($t);
    }
    if ($gene){
      push(@newf,$gene);
    }
    else{
      print STDERR "transform didn't give anything back on the gene:\n";
#      foreach my $t ( @{$gene->get_all_Transcripts} ){
#	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($t);
#      }
    }
    # did we throw exceptions?
    if ($@) {
      print STDERR "Couldn't reverse map gene:  [$@]\n";
#      foreach my $t ( @{$gene->get_all_Transcripts} ){
#	Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($t);
#      }
    }
  }
  
  return @newf;
}

############################################################

=head2 validate_gene

 Title   : validate_gene
 Usage   : $self->validate_gene($gene)
 Function: checks start and end coordinates of each exon of each transcript are sane
 Example : 
 Returns : 1 if gene is valid, otherwise zero
 Args    : $gene: Bio::EnsEMBL::Gene


=cut

sub validate_gene{
  my ($self, $gene) = @_;

  # should be only a single transcript
  my @transcripts = @{$gene->get_all_Transcripts};
  if(scalar(@transcripts) != 1) {
    my $msg = "Rejecting gene - should have one transcript, not " . scalar(@transcripts) . "\n";
    $self->warn($msg);
    return 0;
  }
  
  foreach my $transcript(@transcripts){
    foreach my $exon(@{$transcript->get_all_Exons}){
      unless ( Bio::EnsEMBL::Pipeline::Tools::ExonUtils->_validate_Exon($exon)){
	my $msg = "Rejecting gene because of invalid exon\n";
	$self->warn($msg);
	return 0;
      }
    }
  }
  
  return 1;
}

############################################################

=head2 _recalculate_translation

  
 Arg[1]: a transcript object
 Arg[2]: the strand where the transcript sits
 Return: a brand new transcript object
 Description: a transcript is used as evidence for genomewise
              to recalculate the ORF. The idea is to use this when
              the peptide has been shortened, due to a genewise model
              being incompatible with the cdna splicing. This can happen when genewise cannot find very short exons               
              and attaches them to one of the flanking exons.
              We tell genomewise to keep the splice boundaries pretty much
              static, so that we preserve the original splicing structure.
           
=cut

sub _recalculate_translation{
  my ($self,$mytranscript,$strand) = @_;

  my $this_is_my_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_clone_Transcript($mytranscript);
  my $transcript;

  my $slice = $self->query;
  my $inverted_slice = $slice->invert;
  
  # the genomic sequence to be used with Genomewise (a PrimarySeq)

  # genomewise doesn't know about strands, we need to put everything in forward strand:
  my $genomic_sequence;

  if ( $strand == -1 ){

      print STDERR "In reverse strand: inverting gene\n";
    $genomic_sequence = $inverted_slice;
    
    my $mygene               = Bio::EnsEMBL::Gene->new();
    $mygene->add_Transcript($mytranscript);
    my $gene                 = $mygene->transform($inverted_slice);
    my @inverted_transcripts = @{$gene->get_all_Transcripts};
    $transcript              = $inverted_transcripts[0];
    
    #print STDERR "transcript inverted in recalculate():\n";
    #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($transcript);
  }
  else{
      print STDERR "In forward strand:\n";
    $transcript       = $mytranscript;
    $genomic_sequence = $slice;
  }
  
  my @transcripts;
  push(@transcripts,$transcript);
  
  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::MiniGenomewise->new(
								       -genomic     => $genomic_sequence,
								       -transcripts => \@transcripts,
								       -smell       => 0,
								      );
  eval{
      $runnable->run;
  };
  if ($@){
      print STDERR $@;
  }
  my @trans = $runnable->output;
  
  unless ( scalar(@trans) == 1 ){
      $self->warn("Something went wrong running Genomewise. Got ".scalar(@trans).
		  " transcripts. returning without modifying the translation\n");
      return $mytranscript;
  }
  
  my $newtranscript;
  #print STDERR "transcript from Genomewise:\n";
  #Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Transcript($trans[0]);
    
  # if in the reverse strand, put it back in the original slice
  if ( $strand == -1 ){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->add_Transcript($trans[0]);   
    my $newgene        = $gene->transform($slice);
    my @newtranscripts = @{$newgene->get_all_Transcripts};
    $newtranscript     = $newtranscripts[0];
  }
  else{
    $newtranscript = $trans[0];
  }
  
  # check that everything is sane:
  unless (Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($newtranscript)){
      print STDERR "problem with the translation. Returning the original transcript\n";
      return $this_is_my_transcript;
  }
  return $newtranscript;
}

############################################################

=head2 _transfer_evidence

  
 Arg[1]: reference to Bio::EnsEMBL::Tanscript $combined_transcript
 Arg[2]: reference to Bio::EnsEMBL::Transcript $cdna_transcript
 Return: Bio::EnsEMBL::Transcript
 Description: transfers cdna evidence to combined transcript
           
=cut 

sub _transfer_evidence {
  my ($self, $combined_transcript, $cdna_transcript) = @_;
  foreach my $combined_exon(@{$combined_transcript->get_all_Exons}){
    foreach my $cdna_exon(@{$cdna_transcript->get_all_Exons}){
      # exact match or overlap? 

# exact match
#
#      if($combined_exon->start == $cdna_exon->start &&
#         $combined_exon->end == $cdna_exon->end &&
#         $combined_exon->strand == $cdna_exon->strand){
#         print STDERR "exact match " . $combined_exon->dbID . " with " . $cdna_exon->dbID . "; transferring evidence\n";
#         Bio::EnsEMBL::Pipeline::Tools::ExonUtils-> _transfer_supporting_evidence($cdna_exon, $combined_exon);
#      }

       # overlap - feature boundaries may well be wonky
       if($combined_exon->overlaps($cdna_exon)){
         print STDERR "overlap: " . $combined_exon->dbID . " with " . $cdna_exon->dbID . "; transferring evidence\n";
         Bio::EnsEMBL::Pipeline::Tools::ExonUtils-> _transfer_supporting_evidence($cdna_exon, $combined_exon);
      }
    }
  }
  return $combined_transcript;
}


############################################################
#
# GET/SET METHODS
#
############################################################

sub cdna_vc {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_cdna_vc'} = $value;
    }
    return $obj->{'_cdna_vc'};

}

############################################################

# Function: get/set for e2g gene array

sub e2g_genes {
  my ($self, @genes) = @_;
  
  if (!defined($self->{'_e2g_genes'})) {
    $self->{'_e2g_genes'} = [];
  }
  
  if (@genes) {
    push(@{$self->{'_e2g_genes'}},@genes);
  }
  
  return @{$self->{'_e2g_genes'}};
}

############################################################

#  Function: get/set for genewise gene array

sub gw_genes {
  my ($self, @genes) = @_;
  if (!defined($self->{'_gw_genes'})) {
    $self->{'_gw_genes'} = [];
  }
  
  if (scalar(@genes)) {
    push(@{$self->{'_gw_genes'}},@genes);
  }
  
  return @{$self->{'_gw_genes'}};
}

############################################################

# Function: get/set for combined gene array

sub combined_genes {
  my ($self, @genes) = @_;
  
  if (!defined($self->{'_combined_genes'})) {
    $self->{'_combined_genes'} = [];
  }
  
  if (@genes) {
    push(@{$self->{'_combined_genes'}},@genes);
  }
  
  return @{$self->{'_combined_genes'}};
}

############################################################

# Function: get/set for unmatched genewise gene array

sub unmatched_genewise_genes {
  my ($self, @genes) = @_;
  
  if (!defined($self->{'_unmatched_genewise_genes'})) {
    $self->{'_unmatched_genewise_genes'} = [];
  }

  # we need to store the unmerged version of the genewise gene  
  my %pairs = %{$self->merged_unmerged_pairs()};
  foreach my $merged(@genes){
    my $unmerged = $pairs{$merged};    
    push(@{$self->{'_unmatched_genewise_genes'}},$unmerged);
  }
  
  return @{$self->{'_unmatched_genewise_genes'}};
}

############################################################

# Function: get/set for pairs of merged and unmerged genewise genes
# key is merged gene, value is unmerged
sub merged_unmerged_pairs {
  my ($self, $merged_gene, $unmerged_gene) = @_;
  
  if (!defined($self->{'_merged_unmerged_pairs'})) {
    $self->{'_merged_unmerged_pairs'} = {};
  }
  
  if ($unmerged_gene && $merged_gene) {
    $self->{'_merged_unmerged_pairs'}{$merged_gene}= $unmerged_gene;
  }
  
  # hash ref
  return $self->{'_merged_unmerged_pairs'};
}

############################################################

sub output{
  my ($self,@genes) = @_;
  
  if (!defined($self->{'_output'})) {
    $self->{'_output'} = [];
  }
  
  if(@genes){
    push(@{$self->{'_output'}},@genes);
  }
  
  return @{$self->{'_output'}};
}

############################################################

sub cdna_db {
    my( $self, $cdna_db ) = @_;
    
    if ($cdna_db) 
    {
	$cdna_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
	    || $self->throw("Input [$cdna_db] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
	$self->{_cdna_db} = $cdna_db;
    }
    return $self->{_cdna_db};
}
############################################################

sub genewise_db {
    my( $self, $genewise_db ) = @_;
    
    if ($genewise_db) 
    {
	$genewise_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
	    || $self->throw("Input [$genewise_db] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
	$self->{_genewise_db} = $genewise_db;
    }
    return $self->{_genewise_db};
}

############################################################

sub output_db {
    my( $self, $output_db ) = @_;
    
    if ($output_db) 
    {
	$output_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
	    || $self->throw("Input [$output_db] isn't a Bio::EnsEMBL::DBSQL::DBAdaptor");
	$self->{_output_db} = $output_db;
    }
    return $self->{_output_db};
}

############################################################


1;

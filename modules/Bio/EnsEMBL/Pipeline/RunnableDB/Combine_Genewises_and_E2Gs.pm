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
                                                                      '-db_obj'      => $dbobj,
                                                                      '-golden_path' => $gp,
                                                                      '-input_id'    => $input_id
                                                                    );

$t_e2g->fetch_input();
$t_e2g->run();
$t_e2g->output();
$t_e2g->write_output(); #writes to DB

=head1 DESCRIPTION

Combines predictions from Genewises with Est2genome predictions from cDNA alignments - ripped straight out of TargettedGeneE2G.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::Combine_Genewises_and_E2Gs;

use vars qw(@ISA);
use strict;

# no idea what this is:
#use Storable qw(dclone);

# Object preamble - inheriets from Bio::Root::RootI


use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Gene;
use Bio::SeqIO;

# all the parameters are read from GeneConf.pm
use Bio::EnsEMBL::Pipeline::GeneConf;


@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  # IMPORTANT:
  # SUPER creates dbobj, which is created on run_GeneBuild_runnable
  # this dbobj is a reference to GB_DBHOST@GB_DBNAME containing
  # the features and the dna, so here it is used as refdb only

  # db with the genewises
  my $genewise_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						   '-host'   => $GB_GW_DBHOST,
						   '-user'   => $GB_GW_DBUSER,
						   '-pass'   => $GB_GW_DBPASS,
						   '-dbname' => $GB_GW_DBNAME,
						  );
  
  
  $genewise_db->dnadb($self->dbobj);
  $self->genewise_db($genewise_db);

  # db with the cdnas
  my $cdna_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    '-host'   => $GB_cDNA_DBHOST,
						    '-user'   => $GB_cDNA_DBUSER,
						    '-dbname' => $GB_cDNA_DBNAME,
						    ); 
  
  $cdna_db->dnadb($self->dbobj);
  $self->cdna_db($cdna_db);


  # db where we will write the output (different from any of the above):
  my $comb_db =  new Bio::EnsEMBL::DBSQL::DBAdaptor(
						    '-host'   => $GB_COMB_DBHOST,
						    '-user'   => $GB_COMB_DBUSER,
						    '-dbname' => $GB_COMB_DBNAME,
						    '-pass'   => $GB_COMB_DBPASS,
                                                    ); 
  
  $comb_db->dnadb($self->dbobj);
  $self->output_db($comb_db);

  



  #$self->dbobj($genedb);
 
  return $self;
}

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
  my $chrname;
  my $start;
  my $end;
  
  # input format: chrname.start-end
  if(!($entry =~ /(.*)\.(.*)\-(.*)/)) {
    $self->throw("Not a valid input id... $entry");
  }
  $chrname    = $1;
  $start   = $2 - 10000;
  $end     = $3 + 10000;
  print STDERR "To compensante for the last run of TargettedGeneWise we add 10000b on either side!\n";
  print STDERR "Chromosome id : $chrname\n";
  print STDERR "Range         : $start - $end\n";  
  
  # genewises
  my $sgpa = $self->genewise_db->get_StaticGoldenPathAdaptor();
  my $vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$start,$end);
  $self->vc($vc);
  
  # cdna
  $sgpa = $self->cdna_db->get_StaticGoldenPathAdaptor();
  my $cdna_vc = $sgpa->fetch_VirtualContig_by_chr_start_end($chrname,$start,$end);
  $self->cdna_vc($cdna_vc);
  
}


=head2 run

 Title   : run
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub run {
  my ($self,@args) = @_;
  
  # get genewise genes
  my @similarity_genes = $self->vc->get_Genes_by_Type($GB_SIMILARITY_GENETYPE,'evidence');
  my @targetted_genes  = $self->vc->get_Genes_by_Type($GB_TARGETTED_GW_GENETYPE,'evidence');
  print STDERR "got " . scalar(@similarity_genes) . " similarity genewise genes\n";
  print STDERR "got " . scalar(@targetted_genes) . " targetted genewise genes\n";
  $self->gw_genes( @similarity_genes, @targetted_genes );

  # get e2g genes  
  my @e2g = $self->cdna_vc->get_Genes_by_Type('exonerate_e2g','evidence');

  print STDERR "got " . scalar(@e2g) . " exonerate_e2g genes\n";
  my @newe2g;

 cDNA_GENE:
  foreach my $e2g (@e2g) {
  
  cDNA_TRANSCRIPT:
    foreach my $tran ($e2g->each_Transcript) {
      my $found = 0;
      my @exons = $tran->get_all_Exons;
      @exons = sort {$a->start <=> $b->start} @exons;
      my $seqname;

    cDNA_EXON:
      for (my $i = 1; $i <= $#exons; $i++) {

	# reject trnascripts with long introns
	my $intron = $exons[$i]->start - $exons[$i-1]->end - 1;
	if ($intron > $GB_COMBINED_MAX_INTRON) {
	  print STDERR "rejecting trans_dbID: ".$tran->dbID." for long intron: ". $intron.">".$GB_COMBINED_MAX_INTRON."\n";
	  next cDNA_TRANSCRIPT;
	}
	# check contig consistency
	if ( !( $exons[$i-1]->seqname eq $exons[$i]->seqname ) ){
	  print STDERR "transcript ".$tran->dbID." is partly outside the contig, skipping it...\n";
	  next cDNA_TRANSCRIPT;
	}
      }
      print STDERR "keeping trans_dbID:" . $tran->dbID . "\n";
      push(@newe2g,$e2g);
    }
  }
  
  $self->e2g_genes(@newe2g);
  
  print STDERR "got " . scalar($self->e2g_genes) . " sensible e2g genes\n";
  
  # find which gw matches which e2gs
  my @merged_gw_genes = $self->_merge_gw_genes;
  print STDERR "got " . scalar(@merged_gw_genes) . " merged gw genes\n";  
  
  # check:
  #print STDERR "after merging:\n";
  #foreach my $gw (@merged_gw_genes){
  #  my @trans = $gw->each_Transcript;
  #  print STDERR "genewise ".$gw->dbID." with ".scalar( $trans[0]->get_all_Exons )." exons\n";
  #}

  # first of all, sort gw genes by exonic length and genomic length  
  @merged_gw_genes = sort { my $result = ( $self->_transcript_exonic_length_in_gene($b) <=>
					   $self->_transcript_exonic_length_in_gene($a) );
			    unless ($result){
			      return ( $self->_transcript_length_in_gene($b) <=>
				       $self->_transcript_length_in_gene($a) );
			    }
			    return $result;
			  
			  } @merged_gw_genes;
  
  # keep track of the e2g genes that have been used already
  # we aim for a one-2-one matching between proteins and cdnas
  my %used_e2g;
		   
 GENEWISE:
  foreach my $gw (@merged_gw_genes){
      # should be only 1 transcript
      my @gw_tran  = $gw->each_Transcript;
      my @gw_exons = $gw_tran[0]->get_all_Exons; # ordered array of exons
      my $strand   = $gw_exons[0]->strand;
      
      if($gw_exons[$#gw_exons]->strand != $strand){
	  $self->warn("first and last gw exons have different strands - can't make a sensible combined gene\n");
	  next GENEWISE;
      }
      
      my @matching_e2gs = $self->match_gw_to_e2g($gw);
      
      next GENEWISE unless scalar(@matching_e2gs);
      
      # from the matching ones, take the best fit
      my @list;
      foreach my $e2g ( @matching_e2gs ){
	# we check exon_overlap and fraction of overlap in gw and e2g:
	my ($exon_overlap, $extent_overlap) = $self->_check_overlap( $gw, $e2g );
	push (@list, [$exon_overlap, $extent_overlap,$e2g]);
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
	  # else, by exon extent of the cdna
	  $result = (  $self->_transcript_exonic_length_in_gene($$b[2]) <=>
		       $self->_transcript_exonic_length_in_gene($$a[2]) );
	}
	unless ($result){
	  # else, by genomic length of the cdna
	  $result = (  $self->_transcript_length_in_gene($$b[2]) <=>
		       $self->_transcript_length_in_gene($$a[2]) );
	}
      } @list;
      
      #test:
      print STDERR "matching cdnas:\n";
      foreach my $overlap ( @list ){
	print STDERR "cdna: ".$$overlap[2]->dbID.
	  ", exon_overlap: ".$$overlap[0].
	    ", extent_overlap: ".$$overlap[1]."\n";
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
	  if ( ($count-1) == $howmany && $howmany != 0 ){
	      print STDERR "(all matching cdnas were already used)\n";
	  }
	  next GENEWISE;
      }
      
      $used_e2g{$e2g_match} = 1;
      
      ## pick longest gene match for each gw (exon length - though this may not be the best way)
      #my $chosen_e2g;
      #my $longest = 0;
      #foreach my $e2g(@matching_e2gs){
      #  my $length = 0;
      #  foreach my $exon ($e2g->get_all_Exons){
#	$length += $exon->end - $exon->start + 1; 
#      }
#      if ($length > $longest){
#	$chosen_e2g = $e2g;
#	$longest = $length;
#      }
#    }
#    # there is one transcript per gene
#    my @e2g_tran = $chosen_e2g->each_Transcript; 
#    if ( @gw_tran == 1 &&  @e2g_tran == 1){
#      print STDERR "combining : " . $gw_tran[0]->dbID . " with " . $e2g_tran[0]->dbID . "\n";
#    }
      #    else{
      #      $self->warn("genes with more thatn one transcript -> gw: ".scalar(@gw_tran)." e2g: ".scalar(@e2g_tran));
      #    }
      
      print STDERR "combining gw gene : " . $gw->dbID.":\n";
      foreach my $tran ( $gw->each_Transcript){
	  $tran->sort;
	  foreach my $exon ($tran->get_all_Exons){
	      print STDERR $exon->start."-".$exon->end." ";
	  }
	  print STDERR "\n";
      }
      print STDERR "with e2g gene " . $e2g_match->dbID . ":\n";
      foreach my $tran ( $e2g_match->each_Transcript){
	  $tran->sort;
	  foreach my $exon ($tran->get_all_Exons){
	      print STDERR $exon->start."-".$exon->end." ";
	  }
	  print STDERR "\n";
      }
      
      # build combined genes
      $self->combine_genes($gw, $e2g_match);
  }
  
  
  # remap to raw contig coords
  my @remapped = $self->remap_genes();
  
  $self->output(@remapped);
}

# method to calculate the exonic length of a transcript which is inside a gene
# this assumes that these genewise genes contain one transcript each

sub _transcript_exonic_length_in_gene{
    my ($self,$gene) = @_;
    my @trans = $gene->each_Transcript;
    my $tran = $trans[0];
    my $exonic_length = 0;
    foreach my $exon ($tran->get_all_Exons){
      $exonic_length += ($exon->end - $exon->start + 1);
    }
    return $exonic_length;
}

# method to calculate the length of a transcript in genomic extent, 
# this assumes that these genewise genes contain one transcript each

sub _transcript_length_in_gene{
    my ($self,$gene) = @_;
    my @trans = $gene->each_Transcript;
    my @exons= $trans[0]->get_all_Exons;
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

# method to calculate the amount of overlap between two transcripts

sub _check_overlap{
    my ($self, $gw_gene, $e2g_gene) = @_;
    my @gw_trans  = $gw_gene->each_Transcript;
    my @e2g_trans = $e2g_gene->each_Transcript;
    
    my $exon_overlap = 0;
    my $extent_overlap = 0;
    foreach my $gw_exon ( $gw_trans[0]->get_all_Exons ){
	foreach my $e2g_exon ( $e2g_trans[0]->get_all_Exons ){
	    if ( $gw_exon->overlaps( $e2g_exon ) ){
		$exon_overlap++;
		if ( $gw_exon->start >= $e2g_exon->start ){
		    if ( $gw_exon->end < $e2g_exon->end ){
			$extent_overlap += ( $gw_exon->end - $gw_exon->start + 1 );
		    }
		    elsif ( $gw_exon->end >= $e2g_exon->end ){
			$extent_overlap += ( $e2g_exon->end - $gw_exon->start + 1 );
		    }
		}
		elsif( $gw_exon->start < $e2g_exon->start ){
		    if ( $gw_exon->end < $e2g_exon->end ){
			$extent_overlap += ( $gw_exon->end - $e2g_exon->start + 1 );
		    }
		    elsif ( $gw_exon->end >= $e2g_exon->end ){
			$extent_overlap += ( $e2g_exon->end - $e2g_exon->start + 1 );
		    }
		}
	    }
	}
    }
    #my $gw_length = $self->_transcript_length_in_gene( $gw_gene );
    #my $e2g_length = $self->_transcript_length_in_gene( $e2g_gene );
    
    return( $exon_overlap, $extent_overlap );
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


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut

  
sub write_output {
    my($self) = @_;
    
    # write genes in the database: GB_COMB_DBNAME@GB_COMB_DBHOST
    my $gene_adaptor = $self->output_db->get_GeneAdaptor;
    
  GENE: 
    foreach my $gene ($self->output) {	
	
	eval {
	    $gene_adaptor->store($gene);
	    print STDERR "wrote gene dbID " . $gene->dbID . "\n";
	}; 
	if( $@ ) {
	    print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
	}
	
    }
    
}

=head2 e2g_genes

 Title   : e2g_genes
 Usage   :
 Function: get/set for e2g gene array
 Example :
 Returns : 
 Args    :


=cut

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

=head2 gw_genes

 Title   : gw_genes
 Usage   :
 Function: get/set for genewise gene array
 Example :
 Returns : 
 Args    :


=cut

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

=head2 combined_genes

 Title   : combined_genes
 Usage   :
 Function: get/set for combined gene array
 Example :
 Returns : 
 Args    :


=cut

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

=head2 combine_genes

 Title   : combine_genes
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub combine_genes{
  my ($self, $gw, $e2g) = @_;
  
  my $genetype = $GB_COMBINED_GENETYPE;
  if(!defined ($genetype) || $genetype eq ''){
    $genetype = 'combined_gw_e2g';
    $self->warn("Setting genetype to $genetype\n");
  }
  # get the appropriate analysis from the AnalysisAdaptor
  my $anaAdaptor = $self->output_db->get_AnalysisAdaptor;
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
       -module          => 'Combine_Genewises_and_E2Gs',
      );
  }
  
  # do the actual combining of predictions
  my @newtrans = $self->_make_newtranscript($genetype, $analysis_obj, $gw, $e2g);
  
  $analysis_obj->gff_feature('gene');
  
  # make some lovely genes
  my @genes;
  my $count=0;
  
  foreach my $trans(@newtrans){
    $trans->sort;
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->type($genetype);
    $gene->add_Transcript($trans);
    $gene->analysis($analysis_obj);
    
    if($self->validate_gene($gene)){
      push (@genes,$gene);
      $count++;
    }
  }

  print STDERR "Produced genes:\n";
  foreach my $gene (@genes){
      foreach my $tran ( $gene->each_Transcript ){
	  foreach my $exon ( $tran->get_all_Exons ){
	      print STDERR "exon: ".$exon->start."-".$exon->end." phase: ".$exon->phase." end_phase ".$exon->end_phase."\n";
	      #foreach my $sf ($exon->each_Supporting_Feature){
	      #  print STDERR "evidence: ".$sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	      #}
	  }
	  print STDERR "\n";
      }
  }
  
  $self->combined_genes(@genes);
  
}

=head2 match_gw_to_e2g

 Title   : match_gw_to_e2g
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub match_gw_to_e2g{
  my ($self, $gw) = @_;
  
  print STDERR "\nSearching cDNA for gw gene dbID: ".$gw->dbID."\n";
  my @matching_e2g;
  my @gw_tran = $gw->each_Transcript;
  
  my @gw_exons = $gw_tran[0]->get_all_Exons;
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
  
  # in order to match a starting genewise exon with an e2g exon, we need to have 
  # a. exactly coinciding exon boundaries
  # b. either e2g exon has start <= gw exon start, 
  # OR e2g exon start lies within $exon_slop bp of gw exon start AND the e2g transcript will add extra UTR exons. 
  # previously we had required e2g start to be strictly <= gw start, but this will lose us some valid UTRs
  # substitute "end" for "start" for 3' ends of transcripts
  # BUT don't allow through any e2gs that will result just in a shortened prediction without additional UTR exons.
  my $exon_slop = 20;
  
 E2G:
  foreach my $e2g($self->e2g_genes){
    my @egtran  = $e2g->each_Transcript;
    my @eg_exons = $egtran[0]->get_all_Exons;

    $strand   = $eg_exons[0]->strand; 
    if($eg_exons[$#eg_exons]->strand != $strand){
	$self->warn("first and last e2g exons have different strands - skipping transcript ".$egtran[0]->dbID);
	next E2G;
    }
    if ($strand == 1 ){
      @eg_exons = sort { $a->start <=> $b->start } @eg_exons;
    }
    else{
      @eg_exons = sort { $b->start <=> $a->start } @eg_exons;
    }
    
    my $fiveprime_match = 0;
    my $threeprime_match = 0;
    
    # Lets deal with single exon genes first
    if ($#gw_exons == 0) {
      foreach my $current_exon (@eg_exons) {
	
	if($current_exon->strand != $gw_exons[0]->strand){
	  next E2G;
	}
	
	# don't yet deal with genewise leakage for single exon genes
	if ($gw_exons[0]->end   <= $current_exon->end && 
	    $gw_exons[0]->start >= $current_exon->start){
	  $fiveprime_match = 1;
	  $threeprime_match = 1;
	}
	
      }
      #      if($fiveprime_match && $threeprime_match){
      # can match either end, or both
      if($fiveprime_match || $threeprime_match){
	push(@matching_e2g, $e2g);
      }
      
      # Now the multi exon genewises
    } else {
      
    E2G_EXONS:
	foreach my $current_exon (@eg_exons) {
	    if($current_exon->strand != $gw_exons[0]->strand){
		next E2G;
	    }
	    
	    if($gw_exons[0]->strand == 1){
		
		#FORWARD:
		if ($gw_exons[0]->end == $current_exon->end && 
		    # either e2g exon starts before genewise exon
		    ($current_exon->start <= $gw_exons[0]->start || 
		     # or e2g exon is a bit shorter but there are spliced UTR exons as well
		     (abs($current_exon->start - $gw_exons[0]->start) <= $exon_slop && 
		      $current_exon != $eg_exons[0]))){
		    
		    $fiveprime_match = 1;
		}
		elsif($gw_exons[$#gw_exons]->start == $current_exon->start &&
		      # either e2g exon ends after genewise exon
		      ($current_exon->end >= $gw_exons[$#gw_exons]->end ||
		       # or there are UTR exons to be added
		       (abs($current_exon->end - $gw_exons[$#gw_exons]->end) <= $exon_slop && 
			$current_exon != $eg_exons[$#eg_exons]))){
		    
		    $threeprime_match = 1;
		}
	    }
	    
	    elsif($gw_exons[0]->strand == -1){
	      #REVERSE:
		if ($gw_exons[0]->start == $current_exon->start &&
		    # either e2g exon ends after gw exon
		    ($current_exon->end >= $gw_exons[0]->end ||
		     # or there are UTR exons to be added
		     (abs($current_exon->end - $gw_exons[0]->end) <= $exon_slop &&
		      $current_exon != $eg_exons[0]))){
		    print STDERR "fiveprime reverse match\n";
		    
		    $fiveprime_match = 1;
		}
		elsif ($gw_exons[$#gw_exons]->end == $current_exon->end &&
		       # either e2g exon starts before gw exon
		       ($current_exon->start <= $gw_exons[$#gw_exons]->start ||
			# or there are UTR exons to be added
			(abs($current_exon->start - $gw_exons[$#gw_exons]->start) <= $exon_slop &&
			 $current_exon != $eg_exons[$#eg_exons]))){
		    print STDERR "threeprime reverse match\n";
		    $threeprime_match = 1;
		}
	    }
	}
#      if($fiveprime_match && $threeprime_match){
	if($fiveprime_match || $threeprime_match){
	    push(@matching_e2g, $e2g);
	    
	    # test
	    foreach my $egtran ( $e2g->each_Transcript ){
	      print STDERR "Found cDNA match trans_dbID:".$egtran->dbID."\n";
	      #foreach my $exon ($egtran->get_all_Exons){
	      #  print STDERR $exon->start."-".$exon->end."  ";
	      #}
		#print STDERR "\n";
	    }
	}
    }
} 
  
  
  return @matching_e2g;
  
}
  

=head2 _merge_gw_genes

 Title   : _merge_gw_genes
 Usage   :
 Function: merges adjacent exons if they are frameshifted; stores component exons as subSeqFeatures of the merged exon
 Example :
 Returns : 
 Args    :


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
    my @trans = $gwg->each_Transcript;
    if(scalar(@trans) != 1) { $self->throw("expected one transcript for $gwg\n"); }
    
    ### we follow here 5' -> 3' orientation ###
    $trans[0]->sort;
   
    my @gw_exons = $trans[0]->get_all_Exons;
    
    # check contig consistency:
    for(my $i=1; $i<=$#gw_exons;$i++){
      if ( !( $gw_exons[$i-1]->seqname eq $gw_exons[$i]->seqname ) ){
	print STDERR "transcript (gene_id; ".$gwg->dbID.") is partly outside the contig, skipping it...\n";
	next GW_GENE;
      }
    }
    
    my $strand = $gw_exons[0]->strand; 
    my $previous_exon;

  EXON:      
    foreach my $exon(@gw_exons){
            
      ## genewise frameshift? we merge here two exons separated by max 10 bases into a single exon
      if ($ecount && $pred_exons[$ecount-1]){
	$previous_exon = $pred_exons[$ecount-1];
      }
      
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
	$previous_exon->end($exon->end);
	$previous_exon->add_sub_SeqFeature($exon,'');
	
	my %evidence_hash;
	foreach my $sf($exon->each_Supporting_Feature){
	  if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	      next;
	    }
	  #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	  $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	  $previous_exon->add_Supporting_Feature($sf);
	}
	  next EXON;
      } 
      else{
	# make a new Exon - clone $exon
	my $cloned_exon = new Bio::EnsEMBL::Exon;
	$cloned_exon->start($exon->start);
	$cloned_exon->end($exon->end);
	$cloned_exon->strand($exon->strand);
	$cloned_exon->phase($exon->phase);
	$cloned_exon->end_phase($exon->end_phase);
	$cloned_exon->contig_id($exon->contig_id);
	
	$cloned_exon->attach_seq($self->vc->primary_seq);
	$cloned_exon->add_sub_SeqFeature($exon,'');
	
	#print STDERR "in merged gw_gene, adding evidence in cloned exon:\n";
	my %evidence_hash;
	foreach my $sf($exon->each_Supporting_Feature){
	  if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	    next;
	  }
	  #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	  $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	  $cloned_exon->add_Supporting_Feature($sf);
	}
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

    my $cloned_translation = new Bio::EnsEMBL::Translation;
    $cloned_translation->start($trans[0]->translation->start);
    $cloned_translation->start_exon($trans[0]->translation->start_exon);
    $cloned_translation->end($trans[0]->translation->end);
    $cloned_translation->end_exon($trans[0]->translation->end_exon);

    $merged_transcript->translation($cloned_translation);
    
    # and gene
    $gene->add_Transcript($merged_transcript);
    push(@merged, $gene);
    $count++;
  }
  
  return @merged;
}

=head2 _make_newtranscript

 Title   : _make_newtranscript
 Usage   :
 Function: makes new transcript by combining the genewise and est2genome predictions. Its a monster.
 Example :
 Returns : 
 Args    :


=cut

sub _make_newtranscript {
  my ($self, $genetype, $analysis_obj, $gw, $e2g) = @_;
  my @combined_transcripts  = ();
  
  # should be only 1 transcript
  my @gw_tran  = $gw->each_Transcript;
  #  $gw_tran[0]->sort;
  my @gw_exons = $gw_tran[0]->get_all_Exons; # ordered array of exons
  my @egtran = $e2g->each_Transcript;
  #  $egtran[0]->sort;
  my @e2g_exons  = $egtran[0]->get_all_Exons; # ordered array of exons
  
  # OK, let's see if we need a new gene
  # base it on the existing genewise one
  my $newtranscript = new Bio::EnsEMBL::Transcript;
  foreach my $exon(@gw_exons){
    $newtranscript->add_Exon($exon);
  }

  my $translation   = new Bio::EnsEMBL::Translation;
  $translation->start($gw_tran[0]->translation->start);
  $translation->end($gw_tran[0]->translation->end);
  $translation->start_exon($gw_tran[0]->translation->start_exon);
  $translation->end_exon($gw_tran[0]->translation->end_exon);
  
  #print STDERR "in _make_newtranscript() translation end: ".$translation->end."\n";
  #print STDERR "translation: ".$translation."\n";
  #print STDERR "translation end exon: ".
  #  $translation->end_exon->start."-".$translation->end_exon->end."\n";



  $newtranscript->translation($translation);
  
  $newtranscript->translation->start_exon($newtranscript->start_exon);
  $newtranscript->translation->end_exon($newtranscript->end_exon);
  
  my $eecount = 0;
 EACH_E2G_EXON:
  foreach my $ee (@e2g_exons){
    
    # check strands are consistent
    if ($ee->strand != $gw_exons[0]->strand){
      $self->warn("gw and e2g exons have different strands - can't combine genes\n") ;
      next GENEWISE;
    }
    
    # single exon genewise prediction?
    if(scalar(@gw_exons) == 1) {
      print STDERR "Making new transcript from single exon genewise\n";

      $newtranscript = $self->transcript_from_single_exon_genewise( $ee, 
								    $gw_exons[0], 
								    $newtranscript, 
								    $translation, 
								    $eecount, 
								    @e2g_exons);
    }
    
    else {
      $newtranscript = $self->transcript_from_multi_exon_genewise($ee, 
								  $newtranscript, 
								  $translation, 
								  $eecount,
								  $gw, 
								  $e2g)
    }
    
    # increment the exon
    $eecount++;
    
  } # end foreach my $ee
  
  # the new transcript is made from a merged genewise gene
  # check the transcript and expand frameshifts in all but original 3' gw_exon
  # (the sub_SeqFeatures have been flushed for this exon)
  if (defined($newtranscript)){
    
    #print STDERR "before expanding exons, newtranscript: $newtranscript\n"; 
    #print STDERR "start_exon       : ".$newtranscript->translation->start_exon."\n";
    #print STDERR "start_exon start : ".$newtranscript->translation->start_exon->start."\n";
    #print STDERR "start_exon end   : ".$newtranscript->translation->start_exon->end."\n";
    #print STDERR "start translation: ".$newtranscript->translation->start."\n";
    #print STDERR "end_exon         : ".$newtranscript->translation->end_exon."\n";
    #print STDERR "end_exon start   : ".$newtranscript->translation->end_exon->start."\n";
    #print STDERR "end_exon end     : ".$newtranscript->translation->end_exon->end."\n";
    #print STDERR "end translation  : ".$newtranscript->translation->end."\n";
    foreach my $ex($newtranscript->get_all_Exons){
      
      # test
      #print STDERR $ex->start."-".$ex->end." phase: ".$ex->phase." end_phase: ".$ex->end_phase."\n";
      
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
    # test
    #print STDERR "after expanding exons, newtranscript: $newtranscript\n"; 
    #print STDERR "start_exon       : ".$newtranscript->translation->start_exon."\n";
    #print STDERR "start_exon start : ".$newtranscript->translation->start_exon->start."\n";
    #print STDERR "start_exon end   : ".$newtranscript->translation->start_exon->end."\n";
    #print STDERR "start translation: ".$newtranscript->translation->start."\n";
    #print STDERR "end_exon         : ".$newtranscript->translation->end_exon."\n";
    #print STDERR "end_exon start   : ".$newtranscript->translation->end_exon->start."\n";
    #print STDERR "end_exon end     : ".$newtranscript->translation->end_exon->end."\n";
    #print STDERR "end translation  : ".$newtranscript->translation->end."\n";
    my $count = 0;
    my $previous_ex;
    foreach my $ex($newtranscript->get_all_Exons){
      #print STDERR $ex->start."-".$ex->end." phase: ".$ex->phase." end_phase: ".$ex->end_phase."\n";
      
      # check phases
      if ( $previous_ex ){
	if ( $previous_ex->end_phase != $ex->phase && $previous_ex->end_phase != -1 && $ex->phase != -1 ){
	  $self->warn("inconsistent phases");
	}
      }
      $previous_ex = $ex;
    }

    # dclone messes up database handles
    foreach my $ex($newtranscript->get_all_Exons){
      $ex->attach_seq($self->vc);
      $ex->contig_id($self->vc->id);
      # add new analysis object to the supporting features
      foreach my $sf($ex->each_Supporting_Feature){
	$sf->analysis($analysis_obj);
	$sf->source_tag($genetype);
      }
    }
    
    # double check sane translation start & end
    if($newtranscript->translation->start < 1){
      print STDERR "dodgy translation start - defaulting to 1\n";
      $newtranscript->translation->start(1);
    }
    
    if($newtranscript->translation->end < 1 || 
       $newtranscript->translation->end > $newtranscript->translation->end_exon->length ){
      print STDERR "dodgy translation end " . $newtranscript->translation->end . " - defaulting to end_exon length " . $newtranscript->translation->end_exon->length. "\n";
      $newtranscript->translation->end($newtranscript->translation->end_exon->length);
    }
    


    # check translation is the same as for the genewise gene we're built from
    my $foundtrans = 0;
  GWG:
    foreach my $gwg($self->gw_genes) {
      if($gw->dbID == $gwg->dbID){
	print STDERR "comparing transcripts\n";
	my @tran = $gwg->each_Transcript;
	$foundtrans = $self->compare_transcripts($gwg, $newtranscript);
	
	if ($foundtrans == 1){
	  push (@combined_transcripts, $newtranscript); 	
	  last GWG;
	}
      }
    }
    
 
    if(!$foundtrans){
      $self->warn("prediction with UTRs is not the same as genewise prediction - discarding it\n");
    }
    
  }
  
  
  return @combined_transcripts;
  
}

sub compare_transcripts{
  my ($self, $genewise_gene, $combined_transcript) = @_;
  my @genewise_transcripts = $genewise_gene->each_Transcript;
  
  my $seqout = new Bio::SeqIO->new(-fh => \*STDERR);
  
  if(scalar(@genewise_transcripts != 1)) {
    $self->warn("Panic! Got " . scalar(@genewise_transcripts) . " transcripts, expecting only 1!\n");
    return 0;
  }
  
  my $genewise_translation;
  my $combined_translation;
  
  eval {
      $genewise_translation = $genewise_transcripts[0]->translate;
  };
  
  if ($@) {
      print STDERR "Couldn't translate genewise gene\n";
  }
  else{
    print STDERR "genewise: \n";             
    #$seqout->write_seq($genewise_translation);
  }
  
  $@ = '';
  
  eval{
      $combined_translation = $combined_transcript->translate;
  };
  
  
  if ($@) {
    print STDERR "Couldn't translate combined gene:[$@]\n";
    return 0;
  }
  elsif($combined_translation->seq =~ /\*/){
    print STDERR "combined translation has stops\n";
    return 0;
  }
  else{
    print STDERR "combined: \n";             
    #$seqout->write_seq($combined_translation);
  }	 

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

  return 0;
}


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
	my $ex = $transcript->start_exon;
	
	$ex->start($eg_exon->start);
	$ex->end($eg_exon->end);
	
	# need to explicitly set the translation start & end exons here.
	$translation->start_exon($ex);
	
	# end_exon may be adjusted by 3' coding exon frameshift expansion. Ouch.
	$translation->end_exon($ex);
	
	#checks
	#print STDERR "Single exon genewise:\n";
	#print STDERR "gw    exon: ".$gwstart."-".$gwend." length: ".($gwend-$gwstart+1)."\n";
	#print STDERR "e2g   exon: ".$egstart."-".$egend." length: ".($egend-$egstart+1)."\n";
	#print STDERR "start exon: ".$ex->start."-".$ex->end." length: ".($ex->end-$ex->start+1)." ($ex)\n";
	#print STDERR "translation start_exon: ".$translation->start_exon." end_exon ".$translation->end_exon."\n";

	# need to deal with translation start and end this time - varies depending on strand
	
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
	    $self->throw("Forward strand: setting dodgy translation end: " . $translation->end . " exon_length: " . $translation->end_exon->length . "\n") unless $translation->end <= $translation->end_exon->length;
	}
	#REVERSE:
	elsif($gw_exon->strand == -1){
	    my $diff   = $egend - $gwend;
	    my $tstart = $translation->start;
	    my $tend   = $translation->end;
	    $translation->start($tstart+$diff);
	    $translation->end($tend + $diff);
	    
	    $self->throw("Reverse strand: setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
	    $self->throw("Reverse strand: setting dodgy translation end: " . $translation->end . " exon_length: " . $translation->end_exon->length . "\n") unless $translation->end <= $translation->end_exon->length;
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
	  $translation->end_exon($last);
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
	$self->add_5prime_exons(\$transcript, $exoncount, @e2g_exons);
	$self->add_3prime_exons(\$transcript, $exoncount, @e2g_exons);
	
      }
    return $transcript;
  }

# this method will actually do the combination of both cdna and genewise gene.
# Note that if there is a match on one end but not on the other, the
# code will extend one end, but will leave the other as it is in the
# genewise genes. This will explit cdna matches that look fine on one end
# and we disregard the mismatching part.

sub transcript_from_multi_exon_genewise {
  my ($self, $current_exon, $transcript, $translation, $exoncount, $gw_gene, $eg_gene) = @_;
  
  # $current_exon is the exon one the e2g_transcript we are in at the moment
  # $exoncount is the position of the e2g exon in the array

  # save out current translation->end - we'll need it if we have to expand 3prime exon later
  my $orig_tend = $translation->end;

  my @gwtran  = $gw_gene->each_Transcript;
  $gwtran[0]->sort;
  my @gwexons = $gwtran[0]->get_all_Exons;
  
  my @egtran  = $eg_gene->each_Transcript;
  $egtran[0]->sort;
  my @egexons = $egtran[0]->get_all_Exons;

  # in order to match a starting genewise exon with an e2g exon, we need to have 
  # a. exactly coinciding exon ends
  # b. exon starts lying within $exon_slop bp of each other. 
  # previously we had required e2g start to be strictly <= gw start, but this will lose us some valid UTRs
  # substitute "end" for "start" for 3' ends of transcripts

  ###################
  my $exon_slop = 20;
  ###################

  # compare to the first genewise exon
  

  ################### FORWARD:
  if($gwexons[0]->strand == 1){
      
      ############### 5_PRIME:
      if ($gwexons[0]->end == $current_exon->end && 
	  # either e2g exon starts before genewise exon
	  ($current_exon->start <= $gwexons[0]->start || 
	   # or e2g exon is a bit shorter but there are spliced UTR exons as well
	   (abs($current_exon->start - $gwexons[0]->start) <= $exon_slop && 
	    $current_exon != $egexons[0]))){
	  
	  my $current_start = $current_exon->start;
	  my $gwstart = $gwexons[0]->start;
	  
	  # this exon will be the start of translation, convention: phase = -1
	  my $ex = $transcript->start_exon;
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
	  $self->add_5prime_exons(\$transcript, $exoncount, @egexons);
	  
	  # fix translation start 
	  if($gwstart >= $current_start){
	    # take what it was for the gw gene, and add on the extra
	    my $tstart = $translation->start;
	    #print STDERR "Forward 5': original translation start: $tstart ";
	    $tstart += ($gwstart - $current_start);
	    $translation->start($tstart);
	    #print STDERR "re-setting translation start to: $tstart\n";
	  }
	  # only trust a smaller cdna exon if it is not the first of the transcript
	  # (it could be a truncated cdna)
	  elsif($gwstart < $current_start && $exoncount != 0){
	      
	      # genewise has leaked over the start. Tougher call - we need to take into account the 
	      # frame here as well
	      print STDERR "gw exon starts: $gwstart < new start: $current_start\n";
	      print STDERR "modifying exon, as cdna exon is not the first of transcript-> exoncount = $exoncount\n";
	      
	      # $diff is the number of bases we chop from the genewise exon
	      my $diff   = $current_start - $gwstart;
	      my $tstart = $translation->start;
	      $self->warn("this is a case where gw translation starts at $tstart > 1") if ($tstart>1);
	      
	      print STDERR "gw translation start: ".$tstart."\n";
	      print STDERR "start_exon: ".$translation->start_exon->start.
		"-".$translation->start_exon->end.
		  " length: ".($translation->start_exon->end - $translation->start_exon->start + 1).
		    " phase: ".$translation->start_exon->phase.
		      " end_phase: ".$translation->start_exon->end_phase."\n";
	      
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
	  else{
	      print STDERR "gw exon starts: $gwstart > new start: $current_start";
	      print STDERR "but cdna exon is the first of transcript-> exoncount = $exoncount, so we don't modify it\n";
	  }
	  $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
	  
      } # end 5' exon
      
      ############### 3_PRIME:
      elsif ($gwexons[$#gwexons]->start == $current_exon->start && 
	     # either e2g exon ends after genewise exon
	     ($current_exon->end >= $gwexons[$#gwexons]->end ||
	      # or we allow to end before if there are UTR exons to be added
	      (abs($current_exon->end - $gwexons[$#gwexons]->end) <= $exon_slop && 
	       $current_exon != $egexons[$#egexons]))){   
	  
	my $end_translation_shift = 0;

	# modify the coordinates of the last exon in $newtranscript
	# e2g is larger on this end than gw.
	my $ex = $transcript->end_exon;
	
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
	my $end_ex = $transcript->end_exon;
	$translation->end_exon($end_ex);

	# strand = 1
	my $expanded = $self->expand_3prime_exon($ex, $transcript, 1);
	
	if($expanded){
	  # set translation end to what it originally was in the unmerged genewise gene
	  # taking into account the diff
	  print STDERR "Forward: expanded 3' exon, re-setting end of translation from ".$translation->end." to orig_end ($orig_tend)- ( length_diff + shift_due_to_phases ) ($end_translation_shift)".($orig_tend - $end_translation_shift)."\n";
	  $translation->end($orig_tend - $end_translation_shift);
	}
	
	# finally add any 3 prime e2g exons
	$self->add_3prime_exons(\$transcript, $exoncount, @egexons);
	
      } # end 3' exon    
      
  }
  ########################### REVERSE:
  elsif($gwexons[0]->strand == -1){
      
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
	print STDERR "In Reverse strand. gw exon ends: ".$gwexons[0]->end." > cdna exon end: ".$current_exon->end."\n";
	print STDERR "modifying exon, as cdna exon is not the first of transcript-> exoncount = $exoncount\n";
	
	
	my $diff = $gwexons[0]->end - $current_exon->end;
	my $gwstart = $gwexons[0]->end;
	my $current_start = $current_exon->end;
	my $tstart = $translation->start;
	
	#print STDERR "before the modification:\n";
	#print STDERR "start_exon: ".$translation->start_exon."\n";
	#print STDERR "gw translation start: ".$tstart."\n";
	#print STDERR "start_exon: ".$translation->start_exon->start.
	#  "-".$translation->start_exon->end.
	#    " length: ".($translation->start_exon->end - $translation->start_exon->start + 1).
	#      " phase: ".$translation->start_exon->phase.
	#	" end_phase: ".$translation->start_exon->end_phase."\n";
	
	
	if    ($diff % 3 == 0) { $translation->start(1); }
	elsif ($diff % 3 == 1) { $translation->start(3); }
	elsif ($diff % 3 == 2) { $translation->start(2); }
	else {
	  $translation->start(1);
	  $self->warn("very odd - $diff mod 3 = " . $diff % 3 . "\n");}
      }


      $self->throw("setting very dodgy translation start: " . $translation->start.  "\n") unless $translation->start > 0;
      
      
      # this exon is the start of translation, convention: phase = -1
      my $ex = $transcript->start_exon;
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
      $translation->start_exon($ex);
      
      #print STDERR "after the modification:\n";
      #print STDERR "start_exon: ".$translation->start_exon."\n";
      #print STDERR "gw translation start: ".$translation->start."\n";
      #print STDERR "start_exon: ".$translation->start_exon->start.
      #	"-".$translation->start_exon->end.
      #  " length: ".($translation->start_exon->end - $translation->start_exon->start + 1).
      #  " phase: ".$translation->start_exon->phase.
      #      " end_phase: ".$translation->start_exon->end_phase."\n";
      
            
      $self->add_5prime_exons(\$transcript, $exoncount, @egexons);
      
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
      my $ex = $transcript->end_exon;
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
      my $end_ex = $transcript->end_exon;
      $translation->end_exon($end_ex);
      
      if($expanded){
	# set translation end to what it originally was in the unmerged genewise gene
	print STDERR "Reverse: expanded 3' exon, re-setting translation exon ".$translation->end." to original end( $orig_tend ) - shifts_due_to_phases_etc ( $end_translation_shift ) :".($orig_tend - $end_translation_shift)."\n";
	$translation->end($orig_tend - $end_translation_shift);
      }
      
      $self->add_3prime_exons(\$transcript, $exoncount, @egexons);
      
    } # end 3' exon
  }  
  
  return $transcript; 
}

sub add_5prime_exons{
my ($self, $transcript, $exoncount, @e2g_exons) = @_;

      # add all the exons from the est2genome transcript, previous to this one
      # db handle will be screwed up, need to make new exons from these
      my $c = 0;
      while($c < $exoncount){
	my $newexon = new Bio::EnsEMBL::Exon;
	my $oldexon = $e2g_exons[$c];
	$newexon->start($oldexon->start);
	$newexon->end($oldexon->end);
	$newexon->strand($oldexon->strand);

	# these are all 5prime UTR exons
	$newexon->phase(-1);
	$newexon->end_phase(-1);
	$newexon->contig_id($oldexon->contig_id);
	$newexon->attach_seq($self->vc);
	my %evidence_hash;
	#print STDERR "adding evidence at 5':\n";
	foreach my $sf($oldexon->each_Supporting_Feature){
	  if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	    next;
	  }
	  $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	  #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	  $newexon->add_Supporting_Feature($sf);
	}
#	print STDERR "Adding 5prime exon " . $newexon->start . " " . $newexon->end . "\n";
	$$transcript->add_Exon($newexon);
	$$transcript->sort;
	$c++;
      }
      
}

# $exon is the terminal exon in the genewise transcript, $transcript. We need
# to expand any frameshifts we merged in the terminal genewise exon. 
# The expansion is made by putting $exon to be the last (3' end) component, so we modify its
# start but not its end. The rest of the components are added. The translation end will have to be modified,
# this happens in the method _transcript_from_multi_exon....

sub expand_3prime_exon{
  my ($self, $exon, $transcript, $strand) = @_;
  
  if(scalar($exon->sub_SeqFeature) > 1){
    print STDERR "expanding 3'prime frameshifted exon $exon in strand $strand: ".
      $exon->start."-".$exon->end." phase: ".$exon->phase." end_phase: ".$exon->end_phase."\n";
    my @sf = $exon->sub_SeqFeature;
    
    #if ($strand == 1){
#      my $first = shift(@sf);
#      print STDERR "first component: ".$first->start."-".$first->end." phase ".$first->phase." end_phase ".$first->end_phase."\n";
#      print STDERR "setting exon $exon end: ".$first->end." end_phase: ".$first->end_phase."\n";        
#      $exon->end($first->end); 
#      $exon->dbID($first->dbID);
#      $exon->end_phase($first->end_phase);
    
#      # add back the remaining component exons
#      foreach my $s(@sf){
#	print STDERR "adding exon: ".$s->start."-".$s->end."\n";
#	$transcript->add_Exon($s);
#	$transcript->sort;
#      }

#      # sort out the translation end:

#    }

    
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


# $exoncount tells us which position in the array 
# of e2g exons corresponds to the end of the genewise transcript so we add back 
# exons 3' to that position.
# $exon and $transcript are references to Exon and Transcript objects.

sub add_3prime_exons {
my ($self, $transcript, $exoncount, @e2g_exons) = @_;
# need to deal with frameshifts - 3' exon is a special case as its end might have changed

      # add all the exons from the est2genome transcript, subsequent to this one
      my $c = $#e2g_exons;
      while($c > $exoncount){
	my $newexon = new Bio::EnsEMBL::Exon;
	my $oldexon = $e2g_exons[$c];
	$newexon->start($oldexon->start);
	$newexon->end($oldexon->end);
	$newexon->strand($oldexon->strand);
	
	# these are all exons with UTR:
	$newexon->phase(-1);
	$newexon->end_phase(-1);
	$newexon->contig_id($oldexon->contig_id);
	$newexon->attach_seq($self->vc);
	#print STDERR "adding evidence in 3':\n";
	my %evidence_hash;
	foreach my $sf($oldexon->each_Supporting_Feature){
	  if ( $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} ){
	    next;
	  }
	  $evidence_hash{$sf->hseqname}{$sf->hstart}{$sf->hend}{$sf->start}{$sf->end} = 1;
	  #print STDERR $sf->start."-".$sf->end."  ".$sf->hstart."-".$sf->hend."  ".$sf->hseqname."\n";
	  $newexon->add_Supporting_Feature($sf);
	}
	#	print STDERR "Adding 3prime exon " . $newexon->start . " " . $newexon->end . "\n";
	$$transcript->add_Exon($newexon);
	$$transcript->sort;
	$c--;
      }

}


=head2 remap_genes

 Title   : remap_genes
 Usage   :
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub remap_genes {
  my ($self) = @_;
  my @newf;  
  my $contig = $self->vc;
  
  my @genes = $self->combined_genes;

GENE:  
  foreach my $gene (@genes) {
      my @t = $gene->each_Transcript;
      my $tran = $t[0];
      
      # check that it translates - not the est2genome genes
      if($gene->type eq 'TGE_gw' || $gene->type eq 'combined_gw_e2g'){
	  
	  my $translates = $self->check_translation($tran);
	  next GENE unless $translates;
      }
      
      eval {
	  my $genetype = $gene->type;
	  my $newgene = $contig->convert_Gene_to_raw_contig($gene);
	  # need to explicitly add back genetype and analysis.
	  $newgene->type($genetype);
	  $newgene->analysis($gene->analysis);
	  
	  # sort out supporting feature coordinates
	  foreach my $tran ($newgene->each_Transcript) {
	      foreach my $exon($tran->get_all_Exons) {
		  foreach my $sf($exon->each_Supporting_Feature) {
		      # this should be sorted out by the remapping to rawcontig ... strand is fine
		      if ($sf->start > $sf->end) {
			  my $tmp = $sf->start;
			  $sf->start($sf->end);
			  $sf->end($tmp);
		      }
		  }
	      }
	  }
	  
	  # is this a special case single coding exon gene with UTRS?
	  if($tran->translation->start_exon() eq $tran->translation->end_exon() 
	     && $gene->type eq 'combined_gw_e2g'){
	    #print STDERR "single coding exon, with UTRs\n";
	    
	    # problems come about when we switch from + strand on the vc to - strand on raw contig.
	    my $vc_strand;
	    
	    foreach my $exon($tran->get_all_Exons) {
	      if ($exon eq $tran->translation->start_exon()) {
		$vc_strand = $exon->strand;
		last;
	      }
	    }
	      
	      foreach my $tran ($newgene->each_Transcript) {
		  foreach my $exon ($tran->get_all_Exons) {
		      
		      # oh dear oh dear oh dear
		      # this is still giving some problems
		      if ($exon eq $tran->translation->start_exon()) {
			  if($vc_strand == 1 && $exon->strand == -1){
			      print STDERR "vc strand 1, raw strand -1 - flipping translation start/end\n";
			      $self->warn("something very strange is about to happen to the  exon coordinates for transcript ". $tran->dbID);
			      print STDERR "exon start: ".$exon->start.
				  " changing to ".$exon->end - ($tran->translation->end -1)."\n";
			      print STDERR "exon end  : ".$exon->end.
				  " changing to ".$exon->end - ($tran->translation->start -1)."\n";
			      $exon->end($exon->end - ($tran->translation->start -1));
			      $exon->start($exon->end - ($tran->translation->end -1));
			  }
		      }
		  }
	      }
	      
	  } # end special case single coding exon
	  
	  # final exon coord sanity check
      foreach my $exon($newgene->get_all_Exons){
	# make sure we deal with stickies!
	if($exon->isa("Bio::EnsEMBL::StickyExon")){
	  foreach my $ce($exon->each_component_Exon){
	    # exon start and end must both be within the raw contig!!!
	    if($ce->start < 1){
	      $self->throw("can't set exon->start < 1 (" . $ce->start . ") - discarding gene\n");
	    }
	    
	    if($ce->end > $ce->contig->primary_seq->length){
	      $self->throw("exon extends beyond end of contig - discarding gene\n");
	    }
	  }
	}
	else{
	  # regular exon
	  # exon start and end must both be within the raw contig!!!
	  if($exon->start < 1){
	    $self->throw("can't set exon->start < 1 (" . $exon->start . ") - discarding gene\n");
	  }
	  
	  if($exon->end > $exon->contig->primary_seq->length){
	    $self->throw("exon extends beyond end of contig - discarding gene\n");
	  }
	}
      }
      # if we get to here, the gene is fine, so push it onto the array to be returned
      push(@newf,$newgene);

    };
    
    # did we throw exceptions?
    if ($@) {
      print STDERR "Couldn't reverse map gene:  [$@]\n";
    }

  }
  
  return @newf;
}

=head2 check_translation

 Title   : check_translation
 Usage   :
 Function: checks that transcript translates
 Example :
 Returns : 1 if transcript translates with no stops; otherwise 0
 Args    :


=cut

sub check_translation {
  my ($self, $transcript) = @_;
  my $tseq;

  eval{
    $tseq = $transcript->translate;
  };

  if((!defined $tseq) || ($@)){
    my $msg = "problem translating :\n$@\n";
    $self->warn($msg);
    return 0;
  }

  if ($tseq->seq =~ /\*/ ) {
    $self->warn("discarding gene because of stops\n");
    return 0;
  }

  # if we get to here, all is well
  return 1;

}

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype, $count, $analysis_obj)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::DB::ContigI,
           $genetype: string, $count: integer
           $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript{
  my ($self, $gene, $contig, $genetype, $count, $analysis_obj)=@_;
  $genetype = 'combined_gw_e2g' unless defined ($genetype);
  $count = 1 unless defined ($count);

  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";}
  unless ($contig->isa ("Bio::EnsEMBL::DB::ContigI"))
    {print "$contig must be Bio::EnsEMBL::DB::ContigI\n";}

  my $transcript   = new Bio::EnsEMBL::Transcript;
  my $translation  = new Bio::EnsEMBL::Translation;    
  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
    
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    # make an exon
    my $exon = new Bio::EnsEMBL::Exon;
    
    $exon->contig_id($contig->id);
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase($exon_pred->phase);
    $exon->attach_seq($contig);
    
    # sort out supporting evidence for this exon prediction
    foreach my $subf($exon_pred->sub_SeqFeature){
      $subf->feature1->source_tag($genetype);
      $subf->feature1->primary_tag('similarity');
      $subf->feature1->score(100);
      $subf->feature1->analysis($analysis_obj);
	
      $subf->feature2->source_tag($genetype);
      $subf->feature2->primary_tag('similarity');
      $subf->feature2->score(100);
      $subf->feature2->analysis($analysis_obj);
      
      $exon->add_Supporting_Feature($subf);
    }
    
    push(@exons,$exon);
    
    $excount++;
  }
  
  if ($#exons < 0) {
    print STDERR "Odd.  No exons found\n";
  } 
  else {
    
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
  my @transcripts = $gene->each_Transcript;
  if(scalar(@transcripts) != 1) {
    my $msg = "Rejecting gene - should have one transcript, not " . scalar(@transcripts) . "\n";
    $self->warn($msg);
    return 0;
  }

  foreach my $transcript(@transcripts){
    foreach my $exon($transcript->get_all_Exons){
      if(!$self->validate_exon($exon)){
	my $msg = "Rejecting gene because of invalid exon\n";
	$self->warn($msg);
	return 0;
      }
    }
  }
  
  return 1;
}

=head2 validate_exon

 Title   : validate_exon
 Usage   : $self->validate_exon($exon)
 Function: checks start and end coordinates of exon are sane
 Example : 
 Returns : 1 if exon is valid, otherwise zero
 Args    : $exon: Bio::EnsEMBL::Exon


=cut

sub validate_exon{
  my ($self, $exon) = @_;

  if($exon->start < 0){
    my $msg = "rejecting exon, start < 0 : " . $exon->start . "\n";
    $self->warn($msg);
    return 0;
  }

  elsif($exon->start > $exon->end){
    my $msg = "rejecting exon, start > end : " . $exon->start . " > " . $exon->end . "\n";
    $self->warn($msg);
    return 0;
  }

  elsif($exon->start == $exon->end){
    my $msg = "naughty exon, start == end : " . $exon->start . " == " . $exon->end . " - letting it through\n";
    $self->warn($msg);
    return 1;
  }
  
  return 1;
}

=head2 cdna_vc

 Title   : cdna_vc
 Usage   : $obj->cdna_vc($newval)
 Function: 
 Returns : value of cdna_vc
 Args    : newvalue (optional)


=cut

sub cdna_vc {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_cdna_vc'} = $value;
    }
    return $obj->{'_cdna_vc'};

}

1;

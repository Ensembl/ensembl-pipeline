#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGenewise.pm
#
# Cared for by Ensembl <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGenewise.pm - Targetted genewise Runnable DB

=head1 SYNOPSIS

my $tgw = new Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGenewise
    (  -db => $db,
       -input_id => $input_id);

  $tgw->fetch_input;
  $tgw->run();
  $tgw->output();
  $tgw->write_output(); # write to db

=head1 DESCRIPTION

This object manages the data fetching, running, output parsing, and data storing of Targetted Genewise in the Ensembl pipeline.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGenewise;

use vars qw(@ISA);
use strict;
# Object preamble
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::SeqIO;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Sequences qw (
							     GB_PROTEIN_INDEX
							     GB_PROTEIN_SEQFETCHER
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted qw (
							     GB_TARGETTED_SINGLE_EXON_COVERAGE
							     GB_TARGETTED_MULTI_EXON_COVERAGE
							     GB_TARGETTED_MAX_INTRON
							     GB_TARGETTED_MIN_SPLIT_COVERAGE
							     GB_TARGETTED_GW_GENETYPE
							     GB_TARGETTED_MASKING
							     GB_TARGETTED_SOFTMASK
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases qw (
							     GB_GW_DBHOST
							     GB_GW_DBUSER
							     GB_GW_DBPASS
							     GB_GW_DBNAME
                   GB_GW_DBPORT                                          
							    );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($output_db) = $self->_rearrange([qw(OUTPUT_DB)], @args);

  # makes it easier to run standalone if required
  
  # protein sequence fetcher
  if(!defined $self->seqfetcher) {
    my $seqfetcher = $self->make_seqfetcher($GB_PROTEIN_INDEX, $GB_PROTEIN_SEQFETCHER);
    $self->seqfetcher($seqfetcher);
  }
  throw("no output database defined can't store results $!") unless($output_db);
  $self->output_db($output_db);
  # IMPORTANT
  # SUPER creates db, which is a reference to GB_DBHOST@GB_DBNAME containing
  # features and dna
  # Here it is used as refdb only and we need to make a connection to GB_GW_DBNAME@GB_GW_DBHOST
 

  return $self;
}

=head2 make_seqfetcher

 Title   : make_seqfetcher
 Usage   :
 Function: get/set for sequence fetcher
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    : $indexname - string, $seqfetcher_class - string


=cut

sub make_seqfetcher{
  my ( $self, $index, $seqfetcher_class ) = @_;
  my $seqfetcher;

  (my $class = $seqfetcher_class) =~ s/::/\//g;
  require "$class.pm";

  if(defined $index && $index ne ''){
    my @db = ( $index );
    # make sure that your class is compatible with the index type
    $seqfetcher = "$seqfetcher_class"->new('-db' => \@db, );
  }
  else{
    throw("Can't make seqfetcher\n");
  }

  return $seqfetcher;
}

=head2 protein_id

 Title   : protein_id
 Usage   :
 Function: get/set
 Example :
 Returns : 
 Args    :


=cut

sub protein_id {
    my( $self, $value ) = @_;    
    if ($value) {
        $self->{'_protein_id'} = $value;
    }
    return $self->{'_protein_id'};
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
  my $name;
  my $protein_id; 
  print STDERR "\n\nInput id = ".$entry."\n";
  # chr12:10602496,10603128:Q9UGV6:
  #  print STDERR $entry."\n";
  ($name, $protein_id) = split /\|/, $self->input_id;
  my @array = split(/:/,$name);

  if(@array != 6) {
    throw("Malformed slice name [$name].  Format is " .
          "coord_system:version:start:end:strand");
  }
  
  my ($cs_name, $cs_version, $seq_region, $start, $end, $strand) = @array;
  # we want to give genewise a bit more genomic than the one found by pmatch, 
  if($start > $end){
    my $tmp_start = $end;
    $end = $start;
    $start = $tmp_start;
  }
  #print STDERR "Have pmatch results ".$start." ".$end." ".protein_id."\n";
  my $new_start  = $start - 10000;
  my $new_end    = $end   + 10000;
  
  #print STDERR "Fetching ".$seq_region." ".$start." ".$end."\n";
  if($new_start < 1){
    $new_start = 1;
  }
  my $sliceadp = $self->db->get_SliceAdaptor();
  my $slice = $sliceadp->fetch_by_region($cs_name,$seq_region,
                                         $new_start,$new_end,
                                         $strand, $cs_version);
  
  if($slice->end() > $slice->seq_region_length) {
    #refetch slice, since we went past end, second call is fast
   $new_end = $slice->seq_region_length();
    $slice = $sliceadp->fetch_by_region($cs_name, $seq_region,
                                        $new_start, $new_end,
                                        $strand, $cs_version);
  }


  #print STDERR "Have ".$slice->name." sequence to run\n";
  $self->query($slice);
  my $seq;
  if(@$GB_TARGETTED_MASKING){
    $seq = $slice->get_repeatmasked_seq($GB_TARGETTED_MASKING, $GB_TARGETTED_SOFTMASK);
  }else{
    $seq = $slice;
  }
  $self->protein_id($protein_id);
  #print STDERR $protein_id."\n";
  #print STDERR "running on targetted ".$protein_id." and ".$slice->name."length ".$slice->length."\n";

  # genewise runnable
  # repmasking?

  #print STDERR "Have slice ".$new_start." ".$new_end." ".$seq->length."\n";
  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( '-genomic'        => $seq,
								    '-ids'            => [ $protein_id ] ,
								    '-seqfetcher'     => $self->seqfetcher,
								    '-check_repeated' => 1);
 
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

sub run {
   my ($self,@args) = @_;

   #print STDERR "run runnable\n";
   eval{
     $self->runnable->run();
   };
   if($@){
     $self->throw("Error in BlastMiniGenewise run: \n[$@]\n");
   }
   
   $self->convert_gw_output;
   #print STDERR "converted output\n";
   # clean up tmpfile
   my $tmpfile = $self->{'_tmpfile'};
   unlink $tmpfile;
   #print STDERR "deleted temp files\n";
   # remap genes to raw contig coords
   my @remapped = $self->remap_genes();
   #print STDERR "remapped output\n";
   #print STDERR "have ".@remapped." remapped gene\n";
   $self->output(@remapped);
   #print STDERR "defined output\n";
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
    
   if(@genes){
     push(@{$self->{'_output'}},@genes);
   }

   return @{$self->{'_output'}};
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
      $obj->{'_runnable'} = $value;
    }
    return $obj->{'_runnable'};
}

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   
    Args    :   none

=cut


sub write_output {
  my($self) = @_;
  #print STDERR "writing genes\n";
  my $gene_adaptor = $self->output_db->get_GeneAdaptor;
  my @genes = $self->output;
  #print STDERR "have ".@genes." genes\n";
 GENE: foreach my $gene ($self->output) {	
    # do a per gene eval...
    eval {
      $gene_adaptor->store($gene);
      #print STDERR "wrote gene dbID " . $gene->dbID . "\n";
    }; 
    if( $@ ) {
      print STDERR "UNABLE TO WRITE GENE\n\n$@\n\nSkipping this gene\n";
    }
  }
  return @genes;
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

=head2 convert_gw_output

 Title   : convert_gw_output
 Usage   :
 Function: converts output from Genewise into genes
 Example :
 Returns : 
 Args    :


=cut

sub convert_gw_output {
  my ($self) = @_;
  my $count = 1;
  my $genetype = $GB_TARGETTED_GW_GENETYPE;
  if(!$genetype){
    $genetype = 'TGE_gw';
    $self->warn("Setting genetype to $genetype\n");
  }
  my @results  = $self->runnable->output;
  #print STDERR "BlastMiniGenewise produced ".@results." results\n";

  # Throw here if zero results? Suggests something v. bad has happened 
  # - usually corrupt sequence file means sequences not fetched. We should 
  # never fail to fetch sequences ina a targetted run!
  if(!@results){
    $self->warn("BMG didn't produce any results for ".$self->input_id." ".
               $self->protein_id);
    return;
  }
  # get the appropriate analysis from the AnalysisAdaptor
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;

  my $analysis_obj = $self->analysis;
  if(!$analysis_obj){
    $analysis_obj = $anaAdaptor->fetch_by_logic_name($genetype);
    #print STDERR "have adaptor and analysis objects\n";
  }
  if (!$analysis_obj) {
    # make a new analysis object
    $analysis_obj = new Bio::EnsEMBL::Analysis
      (-db              => 'NULL',
       -db_version      => 1,
       -program         => $genetype,
       -program_version => 1,
       -gff_source      => $genetype,
       -gff_feature     => 'gene',
       -logic_name      => $genetype,
       -module          => 'TargettedGenewise',
      );
  }

  #print STDERR "about to make genes\n";
  my @genes = $self->make_genes($count, $genetype, $analysis_obj, \@results);
  
  # check for stops?
  #print STDERR "have made ".@genes." genes\n";
  #print STDERR "RUNNABLEDB code produced ".@genes." genes\n\n";
  my $exon_count;
  foreach my $g(@genes){
    foreach my $t(@{$g->get_all_Transcripts}){
      $exon_count += @{$t->get_all_Exons};
      # foreach my $e(@{$t->get_all_Exons}){
      # print STDERR "exon ".$e->start." ".$e->end." ".$e->strand."\n";
      # }
    }
  }
  print STDERR "Have ".$exon_count." exons from ".$self->input_id.
    "'s analysis\n";
  $self->gw_genes(@genes);
  
}

=head2 make_genes

 Title   : make_genes
 Usage   : $self->make_genes($count, $genetype, $analysis_obj, $runnable)
 Function: converts the output from $runnable into Bio::EnsEMBL::Genes in
           $contig(VirtualContig) coordinates. The genes have type $genetype, 
           and have $analysis_obj attached. Each Gene has a single Transcript, 
           which in turn has Exons(with supporting features) and a Translation
 Example : 
 Returns : array of Bio::EnsEMBL::Gene
 Args    : $count: integer, $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis, 
           $runnable: Bio::EnsEMBL::Pipeline::RunnableI


=cut

sub make_genes {
  my ($self, $count, $genetype, $analysis_obj, $results) = @_;
  my $contig = $self->query;
  my @genes;
  ##print STDERR "making genes\n";
  throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") 
    unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");
  #print STDERR "have ".@$results." transcript\n";
 MAKE_GENE:  foreach my $tmpf (@$results) {
    my $transcript = $self->make_transcript($tmpf,$self->query,$genetype,$count, $analysis_obj);
    #print STDERR "have validated transcript\n";	
    # validate transcript - validate_transcript returns an array ref
    my $valid_transcripts = $self->validate_transcript($transcript);
    next MAKE_GENE unless defined $valid_transcripts;

    my $gene;
    # make one gene per valid transcript
    foreach my $valid (@$valid_transcripts){

      # add a start codon if appropriate
      $valid = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_start_codon($valid);

      # add a stop codon if appropriate
      $valid = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon($valid);

      $gene   = new Bio::EnsEMBL::Gene;
      $gene->type($genetype);
      $gene->analysis($analysis_obj);
      $gene->add_Transcript($valid);
      push(@genes,$gene);
    }

  }
  #print STDERR "have made ".@genes." genes\n";
  return @genes;
}

=head2 validate_transcript

 Title   : validate_transcript 
 Usage   : my @valid = $self->validate_transcript($transcript)
 Function: Validates a transcript - rejects if mixed strands, 
                                    rejects if low coverage, 
                                    splits if long introns and insufficient coverage of parental protein
                                    rejects unless exon coordinates are sane
 Returns : Ref to @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut

sub validate_transcript {
  my ($self, $transcript) = @_;
  
  my @valid_transcripts;
  #print STDERR "validting transcripts\n";
  my $valid = 1;
  my $split = 0;
  
  
  # check exon phases:
  my @exons = @{$transcript->get_all_Exons};
  #print "there are ".@exons." exons\n";
  #$transcript->sort;
  for (my $i=0;$i<(scalar(@exons-1));$i++){
    my $end_phase = $exons[$i]->end_phase;
    my $phase    = $exons[$i+1]->phase;
    if ( $phase != $end_phase ){
      $self->warn("rejecting transcript with inconsistent phases( $phase - $end_phase) ");
      return undef;
    }
  }

  # check coverage of parent protein
  my $threshold = $GB_TARGETTED_SINGLE_EXON_COVERAGE;
       if(scalar(@exons) > 1){
	 $threshold = $GB_TARGETTED_MULTI_EXON_COVERAGE;
       }

  if(!defined $threshold){
    print STDERR "You must define GB_TARGETTED_SINGLE_EXON_COVERAGE and GB_TARGETTED_MULTI_EXON_COVERAGE in Config::GeneBuild::Targetted.pm\n";
    return undef;
  }

  my $coverage  = $self->check_coverage($transcript);
  if ($coverage < $threshold){
    $self->warn ("Coverage of ". $self->protein_id . " is only $coverage - will be rejected\n");
    return undef;
  }
  
  #print STDERR "Coverage of ". $self->protein_id . " is $coverage%\n";

  my $previous_exon;
  foreach my $exon (@{$transcript->get_all_Exons}){
    if(!$self->validate_exon($exon)){
      print STDERR "Rejecting gene because of invalid exon\n";
      return undef;
    }
       
    # check intron size
    if (defined($previous_exon)) {
      my $intron;
      
      if ($exon->strand == 1) {
				$intron = abs($exon->start - $previous_exon->end - 1);
      } else {
				$intron = abs($previous_exon->start - $exon->end - 1);
      }
      
#      if ($intron > 250000 && $coverage < 95) {

      if ($intron > $GB_TARGETTED_MAX_INTRON && $coverage < $GB_TARGETTED_MIN_SPLIT_COVERAGE ) {
	print STDERR "Intron too long $intron  for transcript " . $transcript->{'temporary_id'} . " with coverage $coverage\n";
	$split = 1;
	$valid = 0;
      }
      
      # check sensible strands
      if ($exon->strand != $previous_exon->strand) {
	print STDERR "Mixed strands for gene " . $transcript->{'temporary_id'} . "\n";
	return undef;
      }
    }
    $previous_exon = $exon;
  }
  
  if ($valid) {
    # make a new transcript that's a copy of all the important parts of the old one
    # but without all the db specific gubbins
    my $newtranscript  = Bio::EnsEMBL::Transcript->new;
    my $newtranslation = Bio::EnsEMBL::Translation->new;

    $newtranscript->translation($newtranslation);
    $newtranscript->translation->start_Exon($transcript->translation->start_Exon);
    $newtranscript->translation->end_Exon  ($transcript->translation->end_Exon);
    $newtranscript->translation->start     ($transcript->translation->start);
		$newtranscript->translation->end       ($transcript->translation->end);
		
		foreach my $exon (@{$transcript->get_all_Exons}){
      $newtranscript->add_Exon($exon);
      foreach my $sf (@{$exon->get_all_supporting_features}){
				$sf->seqname($exon->dbID);
      }
    }
		
    push(@valid_transcripts,$newtranscript);
  }
  elsif ($split){
    # split the transcript up.
    my $split_transcripts = $self->split_transcript($transcript);
    push(@valid_transcripts, @$split_transcripts);
  }

  if(scalar(@valid_transcripts)){
    return \@valid_transcripts;
  }
  else { 
    return undef;
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
  my $contig = $self->query;

  my @genes = $self->gw_genes;
 # print STDERR "REMAPPING GENES\n";
GENE:  foreach my $gene (@genes) {

    my @t = @{$gene->get_all_Transcripts};
    my $tran = $t[0];
    #print STDERR "Have transcript ".$tran."\n";
    #foreach my $exon(@{$tran->get_all_Exons}){
      #print STDERR "Exon ".$exon->start." ".$exon->end." ".
    #    $exon->strand."\n";
    #}
    #print STDERR "gene has ".@t." transcripts\n";
    # check that it translates
    if($gene->type eq $GB_TARGETTED_GW_GENETYPE){
      
      my $translates = $self->check_translation($tran);
      if(!$translates){
        my $msg = "discarding gene - translation has stop codons\n";
        $self->warn($msg);
        next GENE;
      }
    }
    eval {
      my $genetype = $gene->type;
      #$gene->transform;
      # need to explicitly add back genetype and analysis.
      #$newgene->type($genetype);
      #$newgene->analysis($gene->analysis);
      push(@newf,$gene);

      # sort out supporting feature coordinates
      foreach my $tran (@{$gene->get_all_Transcripts}) {
	#print STDERR "transcript has ".$tran->get_all_Exons." exons\n";
	foreach my $exon (@{$tran->get_all_Exons}) {
#	  print STDERR "exon has ".$exon->get_all_supporting_features." supporting features\n";
	  foreach my $sf (@{$exon->get_all_supporting_features}) {
#	   print STDERR "have ".$sf."\n";
	    # this should be sorted out by the remapping to rawcontig ... strand is fine
	    if ($sf->start > $sf->end) {
	      my $tmp = $sf->start;
	      $sf->start($sf->end);
	      $sf->end($tmp);
	    }
	  }
	}
      }
    };

    # did we throw exceptions?
    if ($@) {
      print STDERR "Couldn't reverse map gene:  [$@]\n";
      #$self->throw("couldn't reverse map gene $@");
    }
  }

  return @newf;
}

=head2 check_translation

 Title   : check_translation
 Usage   :
 Function: 
 Example :
 Returns : 1 if transcript translates with no stops (excluding terminal stops), otherwise 0
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
  #print "translation ".$tseq->seq."\n";
  if ($tseq->seq =~ /\*/ ) {
    return 0;
  }
  else{
    return 1;
  }
}



=head2 check_coverage

 Title   : check_coverage
 Usage   :
 Function: checks how much of the parent protein is covered by the genewise prediction
 Example :
 Returns : %coverage of parent protein
 Args    :


=cut

sub check_coverage{
  my ($self, $transcript) = @_;
  my $pstart = 0;
  my $pend = 0;
  my $protname = $self->protein_id;
  my $plength;
  my $fetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher;
  #print STDERR "checking coverage\n";
  my $matches = 0;

  foreach my $exon (@{$transcript->get_all_Exons}) {
    $pstart = 0;
    $pend   = 0;
    #my $exonadp = $exon->adaptor;
    #if(!$exonadp){
    #  die "no exon adaptor defined : $!";
    #}else{
    #  print "exon adaptor is ".$exonadp."\n";
    #}
    my @sfs = @{$exon->get_all_supporting_features};
    #print STDERR "have ".@sfs." supporting features\n";
    foreach my $f(@sfs){
     #print STDERR "have ".$f." from get_all_supporting_features\n";
      if (!defined($protname)){
	$protname = $f->hseqname;
      }
      if($protname ne $f->hseqname){
	warn("$protname ne " . $f->hseqname . "\n");
      }
      
      if((!$pstart) || $pstart > $f->hstart){
	$pstart = $f->hstart;
      }
      
      if((!$pend) || $pend < $f->hend){
	$pend= $f->hend;
      }
    }
    $matches += ($pend - $pstart + 1);
  }
  
  my $seq; 
  eval{
    $seq = $self->seqfetcher->get_Seq_by_acc($protname);
  };
  if ($@) {
    throw("Error fetching sequence for [$protname]: [$@]\n");
  }
  
  throw("No sequence fetched for [$protname]\n") unless defined $seq;
  
  $plength = $seq->length;

  if(!defined($plength) || $plength == 0){
    warn("no sensible length for $protname - can't get coverage\n");
    return 0;
  }

  my $coverage = $matches/$plength;
  $coverage *= 100;
  return $coverage;
}

=head2 make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $contig, $genetype, $count, $analysis_obj)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene, 
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature 
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::RawContig,
           $genetype: string, $count: integer
           $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub make_transcript{
  my ($self, $gene, $contig, $genetype, $count, $analysis_obj)=@_;
  $genetype = 'TGE_gw' unless defined ($genetype);
  $count = 1 unless defined ($count);
  #print STDERR "making transcript\n";
  unless ($gene->isa ("Bio::EnsEMBL::SeqFeatureI"))
    { print "$gene must be Bio::EnsEMBL::SeqFeatureI\n"; }
  unless ($contig->isa ("Bio::EnsEMBL::Slice"))
    { print "$contig must be Bio::EnsEMBL::Slice\n"; }

  my $transcript   = Bio::EnsEMBL::Transcript->new;
  my $translation  = Bio::EnsEMBL::Translation->new;    
  $transcript->translation($translation);

  my $excount = 1;
  my @exons;
  #print "have ".scalar($gene->sub_SeqFeature)." exons\n";
  foreach my $exon_pred ($gene->sub_SeqFeature) {
    my $exon = Bio::EnsEMBL::Exon->new;
   
    $exon->start($exon_pred->start);
    $exon->end  ($exon_pred->end);
    $exon->strand($exon_pred->strand);
    
    $exon->phase($exon_pred->phase);
    $exon->end_phase($exon_pred->end_phase);
    
    $exon->slice($contig);
    #$exon->adaptor($self->db->get_ExonAdaptor);
    
    # sort out supporting evidence for this exon prediction
    my @sf = $exon_pred->sub_SeqFeature;
    #print STDERR "Making Supporting Features in TargettedGenewise\n";
    #foreach my $f(@sf){
    #  print STDERR "supporting feature ".$f->gffstring."\n";
    #}
    #print STDERR "\n\n";
    if(@sf){
      my $align = new Bio::EnsEMBL::DnaPepAlignFeature(-features => \@sf); 
    
      $align->seqname($contig->seq_region_name);
      $align->slice($contig);
      $align->score(100);
      $align->analysis($analysis_obj);
      #print STDERR "adding ".$align." to exon\n";
      $exon->add_supporting_features($align);
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
     #for forward strand:
    
    #start_translation: position on the translation->start_exon coordinate system where
    #the translation starts (counting from the left)
    
    #end_translation  : position on the translation->end_Exon coordinate system where
    #the translation ends (counting from the left)
    
    #for reverse strand:
    
    #start_translation: position on the translation->start_exon coordinate system where
    #the translation starts (counting from the right, which is the direction of translation now)

    #end_translation  : position on the translation->end_Exon coordinate system where
    #the translation ends (counting from the right)
  
    $translation->start_Exon($exons[0]);
    $translation->end_Exon  ($exons[$#exons]);
    
     # phase is relative to the 5' end of the transcript (start translation)
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

=head2 split_transcript

 Title   : split_transcript 
 Usage   : my @splits = $self->split_transcript($transcript)
 Function: splits a transcript into multiple transcripts at long introns. Rejects single exon 
           transcripts that result. 
 Returns : Ref to @Bio::EnsEMBL::Transcript
 Args    : Bio::EnsEMBL::Transcript

=cut


sub split_transcript{
  my ($self, $transcript) = @_;
  #$transcript->sort;
  my @split_transcripts   = ();

  if(!($transcript->isa("Bio::EnsEMBL::Transcript"))){
    $self->warn("[$transcript] is not a Bio::EnsEMBL::Transcript - cannot split");
    return (); # empty array
  }
  
  my $prev_exon;
  my $exon_added = 0;
  my $curr_transcript = Bio::EnsEMBL::Transcript->new;
  my $translation     = Bio::EnsEMBL::Translation->new;
  $curr_transcript->translation($translation);

EXON:   foreach my $exon (@{$transcript->get_all_Exons}){


    $exon_added = 0;
      # is this the very first exon?
    if($exon == $transcript->start_Exon){

      $prev_exon = $exon;
      
      # set $curr_transcript->translation start and start_exon
      $curr_transcript->add_Exon($exon);
      $exon_added = 1;
      $curr_transcript->translation->start_Exon($exon);
      $curr_transcript->translation->start($transcript->translation->start);
      push(@split_transcripts, $curr_transcript);
      next EXON;
    }
    
    if ($exon->strand != $prev_exon->strand){
      return (); # empty array
    }

    # We need to start a new transcript if the intron size between $exon and $prev_exon is too large
    my $intron = 0;
    if ($exon->strand == 1) {
      $intron = abs($exon->start - $prev_exon->end + 1);
    } else {
      $intron = abs($exon->end   - $prev_exon->start + 1);
    }
    
    if ($intron > $GB_TARGETTED_MAX_INTRON) {
      $curr_transcript->translation->end_Exon($prev_exon);

      # need to account for end_phase of $prev_exon when setting translation->end
      $curr_transcript->translation->end($prev_exon->end - $prev_exon->start + 1 - $prev_exon->end_phase);
      
      # start a new transcript 
      my $t  = Bio::EnsEMBL::Transcript->new;
      my $tr = Bio::EnsEMBL::Translation->new;
      $t->translation($tr);

      # add exon unless already added, and set translation start and start_exon
      $t->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      $t->translation->start_Exon($exon);
     
      if ($exon->phase == 0) {
	$t->translation->start(1);
      } elsif ($exon->phase == 1) {
	$t->translation->start(3);
      } elsif ($exon->phase == 2) {
	$t->translation->start(2);
      }

      # start exon always has phase 0
      $exon->phase(0);

      # this new transcript becomes the current transcript
      $curr_transcript = $t;

      push(@split_transcripts, $curr_transcript);
    }

    if($exon == $transcript->end_Exon){
      $curr_transcript->add_Exon($exon) unless $exon_added;
      $exon_added = 1;

      $curr_transcript->translation->end_Exon($exon);
      $curr_transcript->translation->end($transcript->translation->end);
      
    }

    else{
      # just add the exon
      $curr_transcript->add_Exon($exon) unless $exon_added;
    }
    foreach my $sf (@{$exon->get_all_supporting_features}){
	  $sf->feature1->seqname($exon->slice->seq_region_name);

      }
    # this exon becomes $prev_exon for the next one
    $prev_exon = $exon;

  }

  # discard any single exon transcripts
  my @final_transcripts = ();
  my $count = 1;
  
  foreach my $st(@split_transcripts){
    #$st->sort;
    my @ex = @{$st->get_all_Exons};
    if(scalar(@ex) > 1){
      $st->{'temporary_id'} = $transcript->dbID . "." . $count;
      $count++;
      push(@final_transcripts, $st);
    }
  }

  return \@final_transcripts;

}

=head2 output_db

 Title   : output_db
 Usage   : needs to be moved to a genebuild base class
 Function: 
           
 Returns : 
 Args    : 

=cut

sub output_db {
    my( $self, $output_db ) = @_;
    
    if ($output_db) {
      $output_db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")
        || throw("Input [$output_db] isn't a Bio::EnsEMBL::".
                 "DBSQL::DBAdaptor");
      $self->{_output_db} = $output_db;
    }
    if(!$self->{_output_db}){
      $self->{_output_db}= new Bio::EnsEMBL::DBSQL::DBAdaptor
        (
         '-host'   => $GB_GW_DBHOST,
         '-user'   => $GB_GW_DBUSER,
         '-pass'   => $GB_GW_DBPASS,
         '-dbname' => $GB_GW_DBNAME,
         '-port' => $GB_GW_DBPORT,
        );
    }
    return $self->{_output_db};
}



1;

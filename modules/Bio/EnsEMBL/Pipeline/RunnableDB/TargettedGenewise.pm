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
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Sequences qw (
							     GB_PROTEIN_INDEX
							     GB_PROTEIN_SEQFETCHER
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted;

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
  $output_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     '-host'   => $GB_GW_DBHOST,
     '-user'   => $GB_GW_DBUSER,
     '-pass'   => $GB_GW_DBPASS,
     '-dbname' => $GB_GW_DBNAME,
     '-port' => $GB_GW_DBPORT,
    ) if(!$output_db);
  
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
  ($name, $protein_id) = split /\,/, $self->input_id;
  my @array = split(/:/,$name);

  if(@array != 6) {
    throw("Malformed slice name [$name].  Format is " .
          "coord_system:version:seq_region:start:end:strand");
  }
  
  my ($cs_name, $cs_version, $seq_region, $start, $end, $strand) = @array;
  # we want to give genewise a bit more genomic than the one found by pmatch, 
  if($start > $end){
    my $tmp_start = $end;
    $end = $start;
    $start = $tmp_start;
  }
  #print STDERR "Have pmatch results ".$start." ".$end." ".protein_id."\n";
  my $new_start  = $start - $GB_TARGETTED_TERMINAL_PADDING;
  my $new_end    = $end   + $GB_TARGETTED_TERMINAL_PADDING;
  
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
  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( '-genomic'          => $seq,
								    '-ids'              => [ $protein_id ] ,
								    '-seqfetcher'       => $self->seqfetcher,
								    '-endbias'          => 1,
								    '-gap'              => $GB_TARGETTED_GENEWISE_GAP,
								    '-extension'        => $GB_TARGETTED_GENEWISE_EXTENSION,
								    '-matrix'           => $GB_TARGETTED_GENEWISE_MATRIX,
								    '-terminal_padding' => $GB_TARGETTED_TERMINAL_PADDING,
								    '-exon_padding'     => $GB_TARGETTED_EXON_PADDING,
								    '-minimum_intron'   => $GB_TARGETTED_MINIMUM_INTRON,
								    '-check_repeated'   => 1,
                                                                    '-fullseq'           => $GB_TARGETTED_FULLSEQ, 
                                                                     );
 
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

   eval{
     $self->runnable->run();
   };
   if($@){
     throw("Error in BlastMiniGenewise run: \n[$@]\n");
   }

   $self->convert_gw_output;

   my $tmpfile = $self->{'_tmpfile'};
   unlink $tmpfile if ($tmpfile);

   $self->output($self->gw_genes);
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
  print STDERR "writing genes\n";
  my $gene_adaptor = $self->output_db->get_GeneAdaptor;
  my @genes = $self->output;
  #print STDERR "have ".@genes." genes\n";
 GENE: foreach my $gene ($self->output) {	
    # do a per gene eval...
    eval {
      $gene_adaptor->store($gene);
      print STDERR "wrote gene dbID " . $gene->dbID . "\n";
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
    warning("Setting genetype to $genetype\n");
  }

  my @genes  = $self->runnable->output;
  print STDERR "BlastMiniGenewise produced ".@genes." genes\n";

  # Throw here if zero results? Suggests something v. bad has happened
  # - usually corrupt sequence file means sequences not fetched. We should
  # never fail to fetch sequences in a targetted run!
  if(!@genes){
    warning("BMG didn't produce any results for ".$self->input_id." ".
	     $self->protein_id);
    return;
  }
  # get the appropriate analysis from the AnalysisAdaptor
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;

  my $analysis_obj = $self->analysis;
  if(!$analysis_obj){
    $analysis_obj = $anaAdaptor->fetch_by_logic_name($genetype);
  }

  if (!$analysis_obj) {
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

  my @processed_genes = $self->process_genes($count, $genetype, $analysis_obj, \@genes);

  $self->gw_genes(@processed_genes);
}

=head2 process_genes

 Title   : process_genes
 Usage   : $self->process_genes($count, $genetype, $analysis_obj, $runnable)
 Function: checks validity of genes/transcripts retruned from genewise run
 Example :
 Returns : array of Bio::EnsEMBL::Gene
 Args    : $count: integer, $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis, 
  $results: ref to array of genes


=cut

sub process_genes {
  my ($self, $count, $genetype, $analysis_obj, $results) = @_;
  my $contig = $self->query;
  my @genes;

  throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") 
    unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");
  my @seqfetchers;
  push (@seqfetchers, $self->seqfetcher);

 PROCESS_GENE:
  foreach my $unprocessed_gene (@$results) {
    foreach my $transcript(@{$unprocessed_gene->get_all_Transcripts}){

      my $valid_transcripts = 
	Bio::EnsEMBL::Pipeline::Tools::GeneUtils->validate_Transcript($transcript,
								      $self->query,
								      $GB_TARGETTED_MULTI_EXON_COVERAGE,
								      $GB_TARGETTED_SINGLE_EXON_COVERAGE,
								      $GB_TARGETTED_MAX_INTRON,
								      $GB_TARGETTED_MIN_SPLIT_COVERAGE,
								      \@seqfetchers
								     );

      next PROCESS_GENE unless defined $valid_transcripts;

      foreach my $checked_transcript (@$valid_transcripts){

	# add a start codon if appropriate
	$checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_start_codon($checked_transcript);
	
	# add a stop codon if appropriate
	$checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon($checked_transcript);
	
	# attach analysis to supporting features
	foreach my $exon(@{$checked_transcript->get_all_Exons}){
	  foreach my $sf(@{$exon->get_all_supporting_features}){
	    $sf->analysis($analysis_obj);
	  }
	}
	foreach my $tsf (@{$checked_transcript->get_all_supporting_features}){
	  $tsf->analysis($analysis_obj);
	}

	my $gene   = new Bio::EnsEMBL::Gene;
	$gene->type($genetype);
	$gene->analysis($analysis_obj);
	$gene->add_Transcript($checked_transcript);
	push(@genes,$gene);
      }
    }

  }

  return @genes;
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
    my $msg = 
    warning("problem translating :\n$@\n");
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

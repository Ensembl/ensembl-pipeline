#
# Ensembl module for Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise.pm
#
# Cared for by Ensembl <ensembl-dev@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise.pm - Targetted genewise Runnable DB

=head1 SYNOPSIS

my $tgw = new Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise
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


package Bio::EnsEMBL::Pipeline::RunnableDB::TargettedGeneWise;

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
use Bio::EnsEMBL::Pipeline::Tools::GeneUtils;



use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Sequences qw (
							     GB_PROTEIN_INDEX
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Targetted qw (
							     GB_TARGETTED_SINGLE_EXON_COVERAGE
							     GB_TARGETTED_MULTI_EXON_COVERAGE
							     GB_TARGETTED_MAX_INTRON
							     GB_TARGETTED_MIN_SPLIT_COVERAGE
							     GB_TARGETTED_GW_GENETYPE
							    );

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::Scripts   qw (
							     GB_SIZE
							    );

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($output_db) = $self->_rearrange([qw(OUTPUT_DB)], @args);
  
  # protein sequence fetcher
  if(!defined $self->seqfetcher) {
    my $seqfetcher = $self->make_seqfetcher($GB_PROTEIN_INDEX);
    $self->seqfetcher($seqfetcher);
  }
  $self->throw("no output database defined can't store results $!") unless($output_db);
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
 Function: get/set
 Example :
 Returns : Bio::DB::RandomAccessI
 Args    :


=cut

sub make_seqfetcher{
  my ( $self, $index ) = @_;
  my $seqfetcher;

  if(defined $index && $index ne ''){
    my @db = ( $index );
    $seqfetcher = new Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs(
								  '-db' => \@db,
								 );
  }
  else{
    $self->throw("Can't make seqfetcher\n");
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
  my $chr_name;
  my $start;
  my $end;
  my $protein_id; 

  # chr12:10602496,10603128:Q9UGV6:
#  print STDERR $entry."\n";
  if( !($entry =~ /([^\:]+):(\d+),(\d+):([^\:]+):/)) {
      $self->throw("Not a valid input id... $entry");
  }
  
  $chr_name    = $1;
  $protein_id = $4;
  $start   = $2;
  $end     = $3;
  if ($2 > $3) { # let blast sort it out
      $start  = $3;
      $end    = $2;
  }

  # we want to give genewise a bit more genomic than the one found by pmatch, 
  # but we don't want to exceed the multiple of $GB_SIZE,
  # if transcripts cross this boundary, they will get mangled afterwards when we
  # read them in GeneBuilder and store them again using $GB_SIZE's
  
  my $chunk_size = $GB_SIZE;
  my $new_start  = $start - 10000;
  my $new_end    = $end   + 10000;
  
  my $chunk_start = (int($start/$GB_SIZE)) * $GB_SIZE + 1;
  my $chunk_end   = $chunk_start + $GB_SIZE - 1;
  
  # check that end passed in is not larger than end of chunk
  if ( $end > $chunk_end ){
    $self->throw("the end of the protein match passed in ($end) is larger than the genomic chunk we are in ($chunk_end)");
  }
  $new_start = (($start - 10000) < $chunk_start) ? $chunk_start : ($start - 10000);
  $new_end   = (($end + 10000)   > $chunk_end)   ? $chunk_end   : ($end + 10000);
  #print STDERR "fetching slice ".$chr_name." ".$new_start." ".$new_end." \n";
  my $sliceadp = $self->db->get_SliceAdaptor();
  my $slice = $sliceadp->fetch_by_chr_start_end($chr_name,$new_start,$new_end);
  
  $self->query($slice);
  $self->protein_id($protein_id);
  print STDERR $protein_id."\n";

  # check slice seq
#  print ">slice\n" . $slice->seq . "\n";
  print "trying to fetch $chr_name.$new_start-$new_end\n";

  # genewise runnable
  # repmasking?

  my $r = Bio::EnsEMBL::Pipeline::Runnable::BlastMiniGenewise->new( '-genomic'    => $slice,
								    '-ids'        => [ $protein_id ] ,
								    '-seqfetcher' => $self->seqfetcher);
 
#  $self->runnable($r);
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
   ($self->runnable)[0]->run();
   
   $self->convert_gw_output;
   #print STDERR "converted output\n";
   # clean up tmpfile
   my $tmpfile = $self->{'_tmpfile'};
   unlink $tmpfile;
   #print STDERR "deleted temp files\n";
   # remap genes to raw contig coords
   my @remapped = $self->remap_genes();
   #print STDERR "remapped output\n";
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

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   
    Args    :   none

=cut


sub write_output {
  my($self) = @_;
  
  my $gene_adaptor = $self->output_db->get_GeneAdaptor;
  my @genes = $self->output;
  print STDERR "have ".@genes." genes\n";
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
  if(!defined $genetype || $genetype eq ''){
    $genetype = 'TGE_gw';
    $self->warn("Setting genetype to $genetype\n");
  }
  my @results  = ($self->runnable)[0]->output;
  print STDERR "have ".@results." from blastmini genewise\n";
  # get the appropriate analysis from the AnalysisAdaptor
  my $anaAdaptor = $self->db->get_AnalysisAdaptor;

  my $analysis_obj = $anaAdaptor->fetch_by_logic_name($genetype);
  #print STDERR "have adaptor and analysis objects\n";

  if ( !defined $analysis_obj ) {
    # make a new analysis object
    $analysis_obj = new Bio::EnsEMBL::Analysis
      (-db              => 'NULL',
       -db_version      => 1,
       -program         => $genetype,
       -program_version => 1,
       -gff_source      => $genetype,
       -gff_feature     => 'gene',
       -logic_name      => $genetype,
       -module          => 'TargettedGeneWise',
      );
  }

  #print STDERR "about to make genes\n";
  my @genes = $self->make_genes($count, $genetype, $analysis_obj, \@results);
  
  # check for stops?
  #print STDERR "have made ".@genes." genes\n";
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

  $self->throw("[$analysis_obj] is not a Bio::EnsEMBL::Analysis\n") 
    unless defined($analysis_obj) && $analysis_obj->isa("Bio::EnsEMBL::Analysis");

 MAKE_GENE:  foreach my $tmpf (@$results) {
    my $transcript = Bio::EnsEMBL::Pipeline::Tools::GeneUtils->SeqFeature_to_Transcript($tmpf,$self->query, $analysis_obj, $self->output_db, 0);
    my $valid_transcripts = $self->validate_transcript($transcript);
    next MAKE_GENE unless defined $valid_transcripts;

    my $gene;
    # make one gene per valid transcript
    foreach my $valid (@$valid_transcripts){
      $gene   = new Bio::EnsEMBL::Gene;
      $gene->type($genetype);
      $gene->analysis($analysis_obj);
      $gene->add_Transcript($valid);
      push(@genes,$gene);
    }

  }

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
  
  # check coverage of parent protein
  my $threshold = $GB_TARGETTED_SINGLE_EXON_COVERAGE;
  if(scalar(@{$transcript->get_all_Exons}) > 1){
    $threshold = $GB_TARGETTED_MULTI_EXON_COVERAGE;
  }
  
  if(!defined $threshold){
    print STDERR "You must define GB_TARGETTED_SINGLE_EXON_COVERAGE and GB_TARGETTED_MULTI_EXON_COVERAGE in Config::GeneBuild::Targetted.pm\n";
    return undef;
  }
  
  my @seqfetchers = ($self->seqfetcher);

  my $valid_transcripts = 
    Bio::EnsEMBL::Pipeline::Tools::GeneUtils->validate_Transcript($transcript,
								  $self->query,
								  $threshold,
								  $GB_TARGETTED_MAX_INTRON,
								  $GB_TARGETTED_MIN_SPLIT_COVERAGE,
								  \@seqfetchers,
								 );

  return $valid_transcripts;
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

GENE:  foreach my $gene ($self->gw_genes) {

    my @t = @{$gene->get_all_Transcripts};
    my $tran = $t[0];

    # check that it translates
    if($gene->type eq $GB_TARGETTED_GW_GENETYPE){
      
      my $translates = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_Translation($tran);
      if(!$translates){
	my $msg = "discarding gene - translation has stop codons\n";
	$self->warn($msg);
	next GENE;
      }
  }
    eval {
      my $genetype = $gene->type;
      $gene->transform;
      push(@newf,$gene);

      # sort out supporting feature coordinates
      foreach my $tran (@{$gene->get_all_Transcripts}) {
	foreach my $exon (@{$tran->get_all_Exons}) {
	  foreach my $sf (@{$exon->get_all_supporting_features}) {
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


=head2 output_db

 Title   : output_db
 Usage   : needs to be moved to a genebuild base class
 Function: 
           
 Returns : 
 Args    : 

=cut

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


1;

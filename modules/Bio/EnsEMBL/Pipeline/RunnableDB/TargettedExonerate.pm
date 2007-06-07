package Bio::EnsEMBL::Pipeline::RunnableDB::TargettedExonerate;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
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

use Bio::EnsEMBL::Pipeline::Config::GeneBuild::TargettedExonerate;

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

=head2 protein

 Title   : protein
 Usage   :
 Function: get/set for protein Bio::Seq
 Example :
 Returns : 
 Args    :


=cut


sub protein {
  my( $self, $id ) = @_;    

  if (defined($self->{'_protein'})) {
    return $self->{'_protein'};
  }   

  my $seqfetcher = $self->seqfetcher;    
  
  my $seq;
  print STDERR "Fetching ".$id." sequence\n";
  eval {
    $seq = $seqfetcher->get_Seq_by_acc($id);
  };
  
  if ($@) {
    $self->throw("Problem fetching sequence for [$id]: [$@]\n");
  }
  if(!$seq){
    print STDERR "have had problems fetching sequence for ".$id."\n";
  }
  $self->{'_protein'} = $seq;

  return $seq;
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
  print STDERR "\n\nInput id = ".$entry."\n";
  # chromosome:NCBI36:10:103520079:103524963:1,P55075

  my ($name, $protein_id) = split /\,/, $self->input_id;
  my @array = split(/:/,$name);

  if(@array != 6) {
    throw("Malformed slice name [$name].  Format is " .
          "coord_system:version:seq_region:start:end:strand");
  }
  
  my ($cs_name, $cs_version, $seq_region, $start, $end, $strand) = @array;
  # we want to give exonerate a bit more genomic than the one found by pmatch,
  if($start > $end){
    my $tmp_start = $end;
    $end = $start;
    $start = $tmp_start;
  }
  print STDERR "Have pmatch results ".$start." ".$end." ".$protein_id."\n";
  my $new_start  = $start - $GB_TARGETTEDXRATE_TERMINAL_PADDING;
  my $new_end    = $end   + $GB_TARGETTEDXRATE_TERMINAL_PADDING;

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


  print STDERR "Have ".$slice->name." sequence to run\n";
  $self->genomic($slice);
  my $seq;
  if(@$GB_TARGETTEDXRATE_MASKING){
    $seq = $slice->get_repeatmasked_seq($GB_TARGETTEDXRATE_MASKING, $GB_TARGETTEDXRATE_SOFTMASK);
  }else{
    $seq = $slice;
  }

  $self->protein($protein_id);

  # make target(genomic) and query(protein) tmp sequencefiles
  my $genfile = $self->write_sequence_to_file($self->genomic);
  my $pepfile = $self->write_sequence_to_file($self->protein);

  $self->genfile($genfile);
  $self->pepfile($pepfile);


  print STDERR "running on targetted ".$protein_id." and ".$slice->name." length ".$slice->length."\n";

  my $r = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
                                                                     -analysis => $self->analysis,
                                                                     -target_file => $genfile,
                                                                     -query_type => 'protein',
                                                                     -query_file => $pepfile,
                                                                     -options => "--model protein2genome --bestn $GB_TARGETTEDXRATE_BESTN ".
                                                                                 "--maxintron $GB_TARGETTEDXRATE_MAX_INTRON ",
                                                                    );
 
  $self->runnable($r);

}

sub run{
  my ($self) = @_;
  my @results;
  
  throw("Can't run - no runnable objects") unless ($self->runnable);
  
  $self->runnable->run;     
  push ( @results, @{$self->runnable->output} );

  unlink $self->genfile;
  unlink $self->pepfile;

  my $genes = $self->make_genes(\@results);
 
  

  # check genes
  print STDERR "TargettedExonerate: starting with " . scalar(@$genes) . " unchecked transcripts\n";
  my $checked_genes = $self->check_genes($genes);
  print STDERR "TargettedExonerate: now have " . scalar(@$checked_genes) . " checked transcripts\n";
  $self->output($checked_genes);
#  $self->write_output;
}



sub runnable{
  my ($self, $runnable) = @_;

  if(!$self->{'_runnables'}){
    $self->{'_runnables'} = undef;
  }
  if($runnable){
    $self->{'_runnables'} = $runnable;
  }

  return $self->{'_runnables'};
}

sub make_genes{
  my ($self,$transcripts) = @_;

  my (@genes);

  my $slice_adaptor = $self->db->get_SliceAdaptor;
  
  my %genome_slices;

  foreach my $tran ( @$transcripts ){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->type($self->analysis->logic_name);
    
    ############################################################
    # put a slice on the transcript
    
    my $slice_id = $tran->start_Exon->seqname;      
    if (not exists $genome_slices{$slice_id}) {
      # assumes genome seqs were named in the Ensembl API Slice naming 
      # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
      $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
    }
    my $slice = $genome_slices{$slice_id};
    
    foreach my $exon (@{$tran->get_all_Exons}){
      $exon->slice($slice);
      foreach my $evi (@{$exon->get_all_supporting_features}){
        $evi->slice($slice);
        $evi->analysis($self->analysis);
      }
    }
    foreach my $evi (@{$tran->get_all_supporting_features}) {
      $evi->slice($slice);
      $evi->analysis($self->analysis);
    }

    $tran->slice($slice);
    $gene->add_Transcript($tran);
    push( @genes, $gene);
  }


  # NEED TO MAKE TRANSLATIONS - look at EST_GENEBUILDER/TranscriptCoalescer

  return \@genes;
}

sub check_genes{
  my ($self, $unchecked_genes) = @_;
  my $contig = $self->genomic;
  my @checked_genes;

  my @seqfetchers;
  push (@seqfetchers, $self->seqfetcher);

 PROCESS_GENE:
  foreach my $unchecked_gene (@$unchecked_genes) {
    foreach my $transcript(@{$unchecked_gene->get_all_Transcripts}){

      my $valid_transcripts = 
	Bio::EnsEMBL::Pipeline::Tools::GeneUtils->validate_Transcript($transcript,
								      $self->genomic,
								      $GB_TARGETTEDXRATE_MULTI_EXON_COVERAGE,
								      $GB_TARGETTEDXRATE_SINGLE_EXON_COVERAGE,
								      $GB_TARGETTEDXRATE_MAX_INTRON,
								      $GB_TARGETTEDXRATE_MIN_SPLIT_COVERAGE,
								      \@seqfetchers,
								     );

      next PROCESS_GENE unless defined $valid_transcripts;

    CHECKED_TRANSCRIPT:
      foreach my $checked_transcript (@$valid_transcripts){

        # balance complexity with coverage - this is not the right place for this, and possibly not the right way to do it ...
        # Note that this check is not done in _validate_Transcript above as the method is not passed a 'low complexity'
        # threshold value
        my $coverage     =  Bio::EnsEMBL::Pipeline::Tools::GeneUtils->_check_coverage($transcript, \@seqfetchers);
        if(defined $coverage) {
          next CHECKED_TRANSCRIPT unless Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_check_low_complexity($transcript, $coverage);
        }

	# add a start codon if appropriate
	$checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_start_codon($checked_transcript);
	
	# add a stop codon if appropriate
	$checked_transcript = Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->set_stop_codon($checked_transcript);
	
	# attach analysis to supporting features
	foreach my $exon(@{$checked_transcript->get_all_Exons}){
	  foreach my $sf(@{$exon->get_all_supporting_features}){
	    $sf->analysis($unchecked_gene->analysis);
	  }
	}
	foreach my $tsf (@{$checked_transcript->get_all_supporting_features}){
	  $tsf->analysis($unchecked_gene->analysis);
	}

	my $gene   = new Bio::EnsEMBL::Gene;
	$gene->type($unchecked_gene->type);
	$gene->analysis($unchecked_gene->analysis);
	$gene->add_Transcript($checked_transcript);
	push(@checked_genes,$gene);
      }
    }

  }

  return \@checked_genes;
}

sub output{
  my ($self, $output) = @_;
  
  if(!$self->{'_output'}){
    $self->{'_output'} = [];
  }
  if($output && @$output){
    push(@{$self->{'_output'}}, @$output);
  }
  
  return @{$self->{'_output'}};
}

sub genfile{
  my ($self, $genfile) = @_;
  
  if($genfile){
    $self->{'_genfile'} = $genfile;
  }
  
  return $self->{'_genfile'};
}

sub pepfile{
  my ($self, $pepfile) = @_;
  
  if($pepfile){
    $self->{'_pepfile'} = $pepfile;
  }
  
  return $self->{'_pepfile'};
}

sub write_output{
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

sub get_tmp_file {
    my ($self,$dir,$stub,$ext) = @_;
   
    if ($dir !~ /\/$/) {
        $dir = $dir . "/";
    }

    # This is not good

    my $num = int(rand(100000));
    my $file = $dir . $stub . "." . $num . "." . $ext;

    while (-e $file) {
        $num = int(rand(100000));
        $file = $dir.$stub . "." . $num . "." . $ext;
    }

    return $file;
}
   

sub write_sequence_to_file {
    my ($self, $seqobj) = @_;
  
    if (!defined($seqobj)) {
	throw("Must enter a Bio::Seq or a Bio::PrimarySeq object to the write_sequence_to_file");
    }
    if (!$seqobj->isa("Bio::Seq") && !$seqobj->isa("Bio::PrimarySeqI")) {
        throw("Must enter a Bio::Seq or a Bio::PrimarySeqI object to the write_sequence_to_file. Currently [$seqobj]");
    }

    my $file      = $self->get_tmp_file("/tmp","seq","fa");
    my $clone_out = Bio::SeqIO->new(-file => ">$file" , '-format' => 'Fasta');
      
    $clone_out->write_seq($seqobj);

    return $file;
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

=head2 genomic

    Title   :   genomic
    Usage   :   $self->genomic($genomic);
    Function:   Get/set genomic
    Returns :   
    Args    :   

=cut

sub genomic {
    my ($self, $genomic) = @_;

    if (defined($genomic)){ 
	$self->{'_genomic'} = $genomic; 
    }

    return $self->{'_genomic'}
}


1;

#
# Written by Eduardo Eyras
#
# Copyright GRL/EBI 2002
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::Blat

=head1 SYNOPSIS
$database  = a full path location for the directory containing the target (genomic usually) sequence,
@sequences = a list of Bio::Seq objects,
$exonerate = a location for the binary,
$options   = a string with options ,

  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::NewExonerate->new(
								 -database    => $database,
								 -query_seqs  => \@sequences,
								 -query_type => 'dna',
			                                         -target_type=> 'dna',
                                                                 -exonerate   => $exonerate,
								 -options     => $options,
								);

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
 where @results is an array of SeqFeatures, each one representing an aligment (e.g. a transcript), 
 and each feature contains a list of alignment blocks (e.g. exons) as sub_SeqFeatures, which are
 in fact feature pairs.
 
=head1 DESCRIPTION

NewExonerate takes a Bio::Seq (or Bio::PrimarySeq) object and runs Exonerate
against a set of sequences.  The resulting output file is parsed
to produce a set of features.

 here are a few examples of what it can do at this stage:

1. Aligning cdnas to genomic sequence:
   exonerate --exhaustive no --model est2genome cdna.fasta genomic.masked.fasta
   ( this is the default )

2. Behaving like est2genome:
   exonerate --exhaustive yes --model est2genome cdna.fasta genomic.masked.fasta

3. Behaving like blastn:
   exonerate --model affine:local dna.fasta genomic.masked.fasta

4. Smith-Waterman:
   exonerate --exhaustive --model affine:local query.fasta target.fasta

5. Needleman-Wunsch:
   exonerate --exhaustive --model affine:global query.fasta target.fasta

6. Generate ungapped Protein <---> DNA alignments:
   exonerate --gapped no --showhsp yes protein.fasta genome.fasta


=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Pipeline::Runnable::NewExonerate;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($database, $query_seqs, $query_type, $target_type, $exonerate, $options) = $self->_rearrange([qw(
												       DATABASE
												       QUERY_SEQS
												       QUERY_TYPE
												       TARGET_TYPE
												       EXONERATE
												       OPTIONS
												      )
												   ], @args);
  

  $self->{_output} = [];
  # must have a target and a query sequences
  unless( $query_seqs ){
    $self->throw("Exonerate needs a query_seqs: $query_seqs");
  }
  $self->query_seqs(@{$query_seqs});
  
  # you can pass a sequence object for the target or a database (multiple fasta file);
  if( $database ){
    $self->database( $database );
  }
  else{
    $self->throw("Exonerate needs a target - database: $database");
  }

  # Target type: dna  - DNA sequence
  #              prot - protein sequence
  #              dnax - DNA sequence translated in six frames to protein
  #              The default is dna
  if ($target_type){
    $self->target_type($target_type);
  }
  else{
    print STDERR "Defaulting target type to dna\n";
    $self->target_type('dna');
  }

  # Query type: dna  - DNA sequence
  #             rna  - RNA sequence
  #             prot - protein sequence
  #             dnax - DNA sequence translated in six frames to protein
  #             rnax - DNA sequence translated in three frames to protein
  #             The default is dna
  if ($query_type){
    $self->query_type($query_type);
  }
  else{
    print STDERR "Defaulting query type to dna\n";
    $self->query_type('dna');
  }

  # can choose which exonerate to use
  
  #$self->exonerate('exonerate-0.6.7') unless $exonerate;
  $self->exonerate($self->find_executable($exonerate));
  #$self->exonerate('exonerate-0.6.7');

  # can add extra options as a string
  $self->options("   --saturatethreshold 800 --exhaustive no --model est2genome --ryo \"RESULT: %S %p %g %V\\n\" "); 
  
  #if ($options){
  #  $self->options($options);
  #}
  return $self;
}

############################################################
#
# Analysis methods
#
############################################################

=head2 run

Usage   :   $obj->run($workdir, $args)
Function:   Runs exonerate script and puts the results into the file $self->results
            It calls $self->parse_results, and results are stored in $self->output
=cut
  
sub run {
  my ($self) = @_;
    
  my $dir         = $self->workdir();
  my $exonerate   = $self->exonerate;
  my @query_seqs  = $self->query_seqs;
  my $query_type  = $self->query_type;
  my $target_type = $self->target_type;

  # set the working directory (usually /tmp)
  $self->workdir('/tmp') unless ($self->workdir());
  print STDERR" working directory ".$self->workdir()."\n";

  # results go here:
  $self->results($self->workdir()."/results.$$");
  
  # target sequence
  my $target;
  if ( $self->database ){
    $target = $self->database;
  }
  #elsif( $self->target_seq ){
  #  
  #  # write the target sequence into a temporary file then
  #  $target = "$dir/target_seq.$$";
  #  open( TARGET_SEQ,">$target") || $self->throw("Could not open $target $!");
  #  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
  #		 '-fh'     => \*TARGET_SEQ);
  #  $seqout->write_seq($self->target_seq);
  #  close( TARGET_SEQ );
  #}
  
  # write the query sequence into a temporary file then
  my $query = "$dir/query_seqs.$$";
  
  open( QUERY_SEQ,">$query") || $self->throw("Could not open $query $!");
  
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*QUERY_SEQ);
  
  # we write each Bio::Seq sequence in the fasta file $query
  my %length;
  foreach my $query_seq ( @query_seqs ){
      #print STDERR "length( ".$query_seq->display_id.") = ".$query_seq->length."\n";
      $length{ $query_seq->display_id} = $query_seq->length;
      $seqout->write_seq($query_seq);
  }
  close( QUERY_SEQ );
  
  my $command ="exonerate-0.6.7 ".
      $self->options.
	  " --querytype $query_type --targettype $target_type --query $query --target $target > ".$self->results; 
  print STDERR "running exonerate: $command\n";
  
  # system calls return 0 (true in Unix) if they succeed
  $self->throw("Error running exonerate\n") if (system ($command));
  
  $self->output( $self->parse_results( \%length ));
  
  # remove interim files (but do not remove the database if you are using one)
  unlink $query;
  if ( $self->genomic){
    unlink $target;
  }
}

############################################################

=head2 parse_results

 Usage   :   $obj->parse_results
 Function:   reads the Blat output (in PSL format or psLayout ) which has been written to
             a local file $self->results. can accept filenames, filehandles or pipes (\*STDIN)
 Returns :   a list of Features (each alignment), each feature with a list of sub_Seqfeatures (each block of the alignment)
 Args    :   optional filename

=cut

sub parse_results {
  my ($self, $length) = @_;
  
  my $filehandle;
  my %length = %$length;
  my @features_within_features;
  
  my $resfile = $self->results();
  
  if (-e $resfile) {
    if (-z $self->results) {
      print STDERR "Exonerate didn't find any matches\n";
      return; 
    } 
    else {
      open (OUT, "<$resfile") or $self->throw("Error opening ", $resfile, " \n");
      $filehandle = \*OUT;
    }
  }
  else { #it'a a filehandle
    $filehandle = $resfile;
  }
  
  ############################################################
  # store each alignment as a transcript with supporting features
  my @transcripts;

  ############################################################
  # extract values
  while (<$filehandle>){
      
      print $_;
      
      ############################################################
      # the output is of the format:
      #
      # --ryo "RESULT: %S %p %V\n"
      # 
      # It shows the alignments in "sugar" + percent_id + "vulgar blocks" format. 
      #
      # Sugar contains 9 fields
      # ( <qy_id> <qy_start> <qy_len> <qy_strand> <tg_id> <tg_start> <tg_len> <tg_strand> <score> ), 
      # 
      # The vulgar (Verbose Useful Labelled Gapped Alignment Report) blocks are a series 
      # of <label, query_length, target_length> triplets. The label may be one of the following: 
      #
      # M Match 
      # G Gap 
      # C Codon gap 
      # N Non-equivalenced region 
      # 5 5' splice site 
      # 3 3' splice site 
      # I Intron 
      #
      # example:
      # RESULT: AW793782.1 25 250 + 10_NT_035057.1-141075 104447 126318 + 652 82.88 M 30 30 5 0 2 I 0 21645 3 0 2 M 35 35 G 1 0 M 16 16 G 1 0 M 134 134 G 1 0 M 7 7
      # 
      # This gives rise to:
      # M 30 30  ---> exon
      # 5 0 2 
      # I 0 21645 --> intron
      # 3 0 2 
      # M 35 35    \
      # G 1 0      |
      # M 16 16    |
      # G 1 0      |-> exon
      # M 134 134  |
      # G 1 0      |
      # M 7 7     /
      # 
      chomp;
      my ( $tag, $q_id, $q_start, $q_length, $q_strand, $t_id, $t_start, $t_length, $t_strand, $score, $gene_orientation, $perc_id, @blocks) = split;
      
      # the SUGAR 'start' coordinates are 1 less than the actual position on the sequence
      $q_start++;
      $t_start++;

      next unless ( $tag && $tag eq 'RESULT:' );
      
      ############################################################
      # initialize the feature
      my (%query, %target);
      ( $query{score}  , $target{score} )  = ($score,$score);
      ( $query{percent}, $target{percent}) = ( $perc_id, $perc_id );
      ( $query{source} , $target{source} ) = ('exonerate','exonerate');
      ( $query{start}  , $target{start})   = ( $q_start, $t_start );
      ( $query{name}   , $target{name})    = ( $q_id, $t_id);
      if ( $q_strand eq '+' ){ $q_strand = 1; }
      else{ $q_strand = -1; }
      if ( $t_strand eq '+' ){ $t_strand = 1; }
      else{ $t_strand = -1; }

      if ( $gene_orientation eq '+' ){ $gene_orientation = 1 }
      elsif( $gene_orientation eq '-' ){ $gene_orientation = -1 }
      else{ $gene_orientation = 0 }

      ( $query{strand} , $target{strand})  = ( $q_strand, $t_strand );
      
      

      my $transcript = Bio::EnsEMBL::Transcript->new();      
      my $exon = Bio::EnsEMBL::Exon->new();
      my @features;
      my $in_exon = 0;
      my $target_gap_length = 0;
      
      ############################################################
      # exons are delimited by M - I
      # supporting features are delimited by M - G  and  M - I
    TRIAD:
      for( my $i=0; $i<=$#blocks; $i+=3){
	  
	  # do not look at splice sites now
	  if ( $blocks[$i] eq '5' || $blocks[$i] eq '3' ){
	      next TRIAD;
	  }

	  ############################################################
	  # match
	  if ( $blocks[$i] eq 'M' ){
	      if ( $in_exon == 0 ){
		  $in_exon = 1;
		  
		  # start a new exon
		  $exon = Bio::EnsEMBL::Exon->new();
		  $exon->seqname( $t_id );
		  $exon->start( $target{start} );
		  $exon->end( $exon->start +  $blocks[$i+2] - 1 );
		  $exon->strand( $t_strand );
		  $exon->phase( 0 );
		  $exon->end_phase( 0 );
	      }
	      
	      # start a new feature pair
	      $query{end}  = $query{start}  + $blocks[$i+1] - 1;
	      $target{end} = $target{start} + $blocks[$i+2] - 1;
	      print STDERR "query : $query{start} -  $query{end}\n";
	      print STDERR "target: $target{start} -  $target{end}\n";
	      
	      my $feature_pair = $self->create_FeaturePair(\%target, \%query);
	      push( @features, $feature_pair);
	      
	      # re-set the start:
	      $query{start}  = $query{end}  + 1;
	      $target{start} = $target{end} + 1;
	  }
	  ############################################################
	  # gap 
	  if ( $blocks[$i] eq 'G' ){
	      if ( $in_exon ){
		  # keep the same exon
		  
		  # if the gap is in the query, we move the target
		  if ( $blocks[$i+2] ){
		      $target{start} += $blocks[$i+2];
		      # also move the exon:
		      $exon->end( $exon->end + $blocks[$i+2]); 
		  }
		  # if the gap is in the target, we move the query
		  if ( $blocks[$i+1] ){
		      $query{start} += $blocks[$i+1];
		      $target_gap_length += $blocks[$i+1]
		  }
	      }
	  }
	  ############################################################
	  # intron
	  if( $blocks[$i] eq 'I' ){
	      if ( $in_exon ){
		  # emit the current exon
		  $transcript->add_Exon($exon);
		  
		  # add the supporting features
		  my $supp_feature;
		  eval{
		      $supp_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@features);
		  };
		  print STDERR "intron: adding evidence : ".$supp_feature->cigar_string."\n";
		  $exon->add_supporting_features( $supp_feature );
		  
		  $in_exon = 0;
		  
		  @features = ();
		  		  
		  # reset the start in the target only
		  my $intron_length = $blocks[$i+2];
		  $target{start} += $intron_length + 1;
	      }
	  }
      } # end of TRIAD
      
      ############################################################
      # emit the last exon and the last set of supporting features
      # and add the supporting features
      #print STDERR "created features:\n";
      #foreach my $f (@features){
      #print STDERR $f->gffstring."\n";
      #}

      my $supp_feature;
      eval{
	  $supp_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@features);
      };
      print STDERR "outside: adding evidence : ".$supp_feature->cigar_string."\n";
      $exon->add_supporting_features( $supp_feature );
      $transcript->add_Exon($exon);
      
      my $coverage = sprintf "%.2f", ( $q_length - $target_gap_length ) / $length{$q_id};
      print STDERR "coverage: $coverage\n";
      
      foreach my $exon ( @{$transcript->get_all_Exons} ){
	  foreach my $evidence ( @{$exon->get_all_supporting_features} ){
	      $evidence->score($coverage);
	  }
      }

      # test
      print STDERR "produced transcript:\n";
    Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils->_print_Evidence($transcript);
      
      ############################################################
      # reject single-exon alignments which cover less than half the cdna/est:
      if ( scalar( @{$transcript->get_all_Exons} ) ){
	  my $mapped_length = $transcript->get_all_Exons->[0]->end - $transcript->get_all_Exons->[0]->start + 1;
	  if ( $mapped_length <= $length{$q_id}/2 ){
	      next;
	  }
      }
      
      push( @transcripts, $transcript );
      
  } # end of while loop
  
  close $filehandle;
  
  #get rid of the results file
  unlink $self->results;
  
  return @transcripts;
}

############################################################
#
# get/set methods
#
############################################################

sub query_seqs {
  my ($self, @seqs) = @_;
  if (@seqs){
    unless ($seqs[0]->isa("Bio::PrimarySeqI") || $seqs[0]->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    push(@{$self->{_query_seqs}}, @seqs) ;
  }
  return @{$self->{_query_seqs}};
}

############################################################

sub genomic {
  my ($self, $seq) = @_;
  if ($seq){
    unless ($seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_genomic} = $seq ;
  }
  return $self->{_genomic};
}

############################################################

sub exonerate {
  my ($self, $location) = @_;
  if ($location) {
    $self->throw("Exonerate not found at $location: $!\n") unless (-e $location);
    $self->{_exonerate} = $location ;
  }
  return $self->{_exonerate};
}

############################################################

sub options {
  my ($self, $options) = @_;
  if ($options) {
    $self->{_options} = $options ;
  }
  return $self->{_options};
}

############################################################

sub output {
    my ($self, @output) = @_;
  if (@output) {
      unless( $self->{_output} ){
	  $self->{_output} = [];
      }
      push( @{$self->{_output}}, @output );
  }
  return @{$self->{_output}};
}

############################################################

sub database {
  my ($self, $database) = @_;
  if ($database) {
    $self->{_database} = $database;
  }
  return $self->{_database};
}
############################################################

sub query_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna' || $type eq 'protein' ){
      $self->throw("not the right query type: $type");
    }
    $self->{_query_type} = $type;
  }
  return $self->{_query_type};
}

############################################################

sub target_type {
  my ($self, $mytype) = @_;
  if (defined($mytype) ){
    my $type = lc($mytype);
    unless( $type eq 'dna'|| $type eq 'protein' ){
      $self->throw("not the right target type: $type");
    }
    $self->{_target_type} = $type ;
  }
  return $self->{_target_type};
}

############################################################


1;


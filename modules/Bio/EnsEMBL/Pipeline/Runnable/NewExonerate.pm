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
  $self->options("   --exhaustive no --model est2genome --ryo \"RESULT: %S %p %V\\n\" "); 
  #$self->options("   --exhaustive no --model est2genome  --showvulgar");
  
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
  foreach my $query_seq ( @query_seqs ){
    $seqout->write_seq($query_seq);
  }
  close( QUERY_SEQ );
    
  my $command ="exonerate-0.6.7 ".
      $self->options.
	  " --querytype $query_type --targettype $target_type --query $query --target $target > ".$self->results; 
  print STDERR "running exonerate: $command\n";
  
  # system calls return 0 (true in Unix) if they succeed
  $self->throw("Error running exonerate\n") if (system ($command));
  
  $self->output( $self->parse_results );
  
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
  my ($self,$filehandle) = @_;
  
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
  
  #extract values
  while (<$filehandle>){
    
      print $_;
      
      # --ryo "RESULT: %S %p %V\n"
      # 
      # Shows the alignments in "sugar" + perc_it + "vulgar blocks" format. 
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
      my ( $tag, $q_id, $q_start, $q_length, $q_strand, $t_id, $t_start, $t_length, $t_strand, $score, $perc_id, @blocks) = split;
      
      next unless ( $tag eq 'RESULT:' );
      
      my $seq_end   = $seq_start + $seq_length - 1;
      my $slice_end = $slice_start + $slice_length  - 1;
      
      my $superfeature = Bio::EnsEMBL::SeqFeature->new();
      

      my (%query, %target);
      ( $query{score}  , $target{score} )  = ($score,$score);
      ( $query{percent}, $target{percent}) =  ( $perc_id, $perc_id );
      ( $query{source} , $target{source} ) = ('exonerate','exonerate');
      ( $query{start}  , $target{start})   = ( $q_start, $t_start );
      
      my ( $exon, $intron ) = (0,0);
      for( my $i=0; $i<=$#blocks; $i+=3){
	  
	  # do not look at splice sites now
	  if ( $blocks[$i] eq '5' || $blocks[$i] eq '3' ){
	      next;
	  }
	  
	  # intron
	  if ( $blocks[$i] eq 'I' ){
	      # emit the current exon
	      my $feature_pair = $self->create_FeaturePair(\%feat1, \%feat2);
	      $superfeature->add_sub_SeqFeature( $feature_pair,'EXPAND');


	      # and reset the start



	  # gap
	  if ( $blocks[$i] eq 'G' ){
	      # take the value of the gap
	      my $gap = $blocks[$i+1] || $blocks[$i+2];
	      
	      # add it to the end
	      if ($exon ){
		  $query{end}  += $gap;
		  $target{end} += $gap;
	      }
	  }
	      
	  # match
	  if ( $blocks[$i] eq 'M' ){
	      if ( $exon == 0 ){
		  $query{end}  = $query{start}  + $blocks[$i+1] - 1;
		  $target{end} = $target{start} + $blocks[$i+2] - 1;
	      }
		  $exon = 1;
		  $intron = 0;
	      
	      
	  }
	  elsif( $blocks[$i] eq 'I' ){
	      $exon = 0;
	      $intron = 1;
	  }
	  
	  
      }


      # create as many features as blocks there are in each output line
#    my (%feat1, %feat2);
#    $feat1{name} = $t_name;
#    $feat2{name} = $q_name;
    
#    $feat2{strand} = 1;
#    $feat1{strand} = $strand;
    
#    # all the block sizes add up to $matches + $mismatches + $rep_matches
    
#    # percentage identity =  ( matches not in repeats + matches in repeats ) / ( alignment length )
#    #print STDERR "calculating percent_id and score:\n";
#    #print STDERR "matches: $matches, rep_matches: $rep_matches, mismatches: $mismatches, q_length: $q_length\n";
#    #print STDERR "percent_id = 100x".($matches + $rep_matches)."/".( $matches + $mismatches + $rep_matches )."\n";
#    my $percent_id = sprintf "%.2f", ( 100 * ($matches + $rep_matches)/( $matches + $mismatches + $rep_matches ) );
    
#    # or is it ...?
#    ## percentage identity =  ( matches not in repeats + matches in repeats ) / query length
#    #my $percent_id = sprintf "%.2d", (100 * ($matches + $rep_matches)/$q_length );
    
#    # we put basically score = coverage = ( $matches + $mismatches + $rep_matches ) / $q_length
#    #print STDERR "score = 100x".($matches + $mismatches + $rep_matches)."/".( $q_length )."\n";
    
#    my $score   = sprintf "%.2f", ( 100 * ( $matches + $mismatches + $rep_matches ) / $q_length );
    
#    # size of each block of alignment (inclusive)
#    my @block_sizes     = split ",",$block_sizes;
    
#    # start position of each block (you must add 1 as psl output is off by one in the start coordinate)
#    my @q_start_positions = split ",",$q_starts;
#    my @t_start_positions = split ",",$t_starts;
    
#    $superfeature->seqname($q_name);
#    $superfeature->score( $score );
#    $superfeature->percent_id( $percent_id );

#    # each line of output represents one possible entire aligment of the query (feat1) and the target(feat2)
#    for (my $i=0; $i<$block_count; $i++ ){
      
#      $feat2 {start} = $q_start_positions[$i] + 1;
#      $feat2 {end}   = $feat2{start} + $block_sizes[$i] - 1;
      
#      $feat1 {start} = $t_start_positions[$i] + 1;
#      $feat1 {end}   = $feat1{start} + $block_sizes[$i] - 1;
      
#      # we put all the features with the same score and percent_id
#      $feat2 {score}   = $score;
#      $feat1 {score}   = $feat2 {score};
#      $feat2 {percent} = $percent_id;
#      $feat1 {percent} = $feat2 {percent};
      
#      # other stuff:
#      $feat1 {db}         = undef;
#      $feat1 {db_version} = undef;
#      $feat1 {program}    = 'blat';
#      $feat1 {p_version}  = '1';
#      $feat1 {source}     = 'blat';
#      $feat1 {primary}    = 'similarity';
#      $feat2 {source}     = 'blat';
#      $feat2 {primary}    = 'similarity';
      
#      my $feature_pair = $self->create_FeaturePair(\%feat1, \%feat2);
#      $superfeature->add_sub_SeqFeature( $feature_pair,'EXPAND');
#    }
#    push(@features_within_features, $superfeature);
  
  } # CLOSE WHILE LOOP
  close $filehandle;
  
  #get rid of the results file
  unlink $self->results;
  
  return @features_within_features;
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


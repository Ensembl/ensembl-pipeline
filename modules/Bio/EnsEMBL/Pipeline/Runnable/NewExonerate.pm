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

Bio::EnsEMBL::Pipeline::Runnable::NewExonerate

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

  my ($database, $query_seqs, $query_type, $target_type, $exonerate, $options) = 
      $self->_rearrange([qw(
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

  ############################################################
  # Target type: The default is dna
  if ($target_type){
    $self->target_type($target_type);
  }
  else{
    print STDERR "Defaulting target type to dna\n";
    $self->target_type('dna');
  }

  ############################################################
  # Query type: The default is dna
  if ($query_type){
    $self->query_type($query_type);
  }
  else{
    print STDERR "Defaulting query type to dna\n";
    $self->query_type('dna');
  }

  ############################################################
  # we are using the latest version: exonerate-0.6.7
  $self->exonerate('/usr/local/ensembl/bin/exonerate-0.6.7');

  ############################################################
  # options
  my $basic_options = " --exhaustive FALSE --model est2genome --ryo \"RESULT: %S %p %g %V\\n\" "; 
  #my $basic_options = "  --ryo \"RESULT: %S %p %g %V\\n\" "; 
  
  # can add extra options as a string
  if ($options){
    $basic_options .= $options;
  }
  $self->options( $basic_options);
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

  my $verbose = 0;
    
  my $dir         = $self->workdir();
  my $exonerate   = $self->exonerate;
  my @query_seqs  = $self->query_seqs;
  my $query_type  = $self->query_type;
  my $target_type = $self->target_type;

  # set the working directory (usually /tmp)
  $self->workdir('/tmp') unless ($self->workdir());
  
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
  
  ############################################################
  # write the query sequence into a temporary file then
  my $query = "$dir/query_seqs.$$";
  open( QUERY_SEQ,">$query") || $self->throw("Could not open $query $!");
  my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
			       '-fh'     => \*QUERY_SEQ);
  
  # we write each Bio::Seq sequence in the fasta file $query
  my %length;
  foreach my $query_seq ( @query_seqs ){
      ## calculate the length
      #print STDERR "length( ".$query_seq->display_id.") = ".$query_seq->length."\n";
      $length{ $query_seq->display_id} = $query_seq->length;
      $seqout->write_seq($query_seq);
  }
  close( QUERY_SEQ );
  
  my $command ="exonerate-0.6.7 ".
      $self->options.
	  " --querytype $query_type --targettype $target_type --query $query --target $target |";

  print STDERR "running exonerate: $command\n" if $verbose;

  open( EXO, $command ) || $self->throw("Error running exonerate $!");
  
  # system calls return 0 (true in Unix) if they succeed
  #$self->throw("Error running exonerate\n") if (system ($command));
  
    
  ############################################################
  # store each alignment as a transcript with supporting features
  my @transcripts;

  ############################################################
  # parse results - avoid writing to disk the output
  while (<EXO>){
      
    print STDERR $_ if $verbose;
      
      ############################################################
      # the output is of the format:
      #
      # --ryo "RESULT: %S %p %g %V\n"
      # 
      # It shows the alignments in "sugar" + percent_id + gene orientation + "vulgar blocks" format. 
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
    my ( $tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score, $perc_id, $gene_orientation, @blocks) = split;
    
    # the SUGAR 'start' coordinates are 1 less than the actual position on the sequence
    $q_start++;
    $t_start++;
    
    next unless ( $tag && $tag eq 'RESULT:' );
    
    ############################################################
    # initialize the feature
    my (%query, %target);
    my (%rev_query, %rev_target);

    ( $query{score}  , $target{score} )  = ($score,$score);
    ( $query{percent}, $target{percent}) = ( $perc_id, $perc_id );
    ( $query{source} , $target{source} ) = ('exonerate','exonerate');
    ( $query{name}   , $target{name})    = ( $q_id, $t_id);
 
    if ( $q_strand eq '+' ){ $q_strand = 1; }
    else{ $q_strand = -1; }
    if ( $t_strand eq '+' ){ $t_strand = 1; }
    else{ $t_strand = -1; }

    ( $query{strand} , $target{strand})  = ( $q_strand, $t_strand );

    ############################################################
    # coordinates are corrected later according to strand
    $query{start}  = $q_start;
    $target{start} = $t_start;
        
    ############################################################
    # orientation is not used at the moment
    if ( $gene_orientation eq '+' ){   $gene_orientation = 1 }
    elsif( $gene_orientation eq '-' ){ $gene_orientation = -1 }
    else{ $gene_orientation = 0 }
    
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
	  $exon->strand( $t_strand );
	  $exon->phase( 0 );
	  $exon->end_phase( 0 );
	}
	# for every match we increase the end coordinate of the current exon
	
	#print STDERR "setting exon end = ". ($exon->start +  $blocks[$i+2] - 1 )."\n";
	$exon->end( $exon->start +  $blocks[$i+2] - 1 );
	
	# start a new feature pair
	$query{end}  = $query{start}  + $blocks[$i+1] - 1;
	$target{end} = $target{start} + $blocks[$i+2] - 1;
	
	if ($verbose){
	  print STDERR "FEATURE:\n";
	  print STDERR "query length: ".$length{$q_id}."\n";
	  print STDERR "query : $query{start} -  $query{end}\n";
	  print STDERR "target: $target{start} -  $target{end}\n";
	}
	
	############################################################
	# if the query has been inverted we coun from the end and
	# inver start and end to inforce start < end
	%rev_query  = %query;
	%rev_target = %target;
	if ( $q_strand == -1 ){
	  $rev_query{end}   = $length{$q_id} - $query{start} + 1;
	  $rev_query{start} = $length{$q_id} - $query{end} + 1;
	  
	  print STDERR "rev_query{end}   = ".($length{$q_id})." - ".$query{start}." + 1\n" if $verbose;
	  print STDERR "rev_query{start} = ".($length{$q_id})." - ".$query{start}." + 1\n" if $verbose;
	  my $feature_pair = $self->create_FeaturePair(\%rev_target, \%rev_query);
	  print STDERR "adding feature: ".$feature_pair->gffstring."\n" if $verbose;
	  push( @features, $feature_pair);
	}
	else{	      	    
	  my $feature_pair = $self->create_FeaturePair(\%target, \%query);
	  print STDERR "adding feature: ".$feature_pair->gffstring."\n" if $verbose;
	  push( @features, $feature_pair);
	}
	
	# re-set the start:
	$query{start}  = $query{end}  + 1;
	$target{start} = $target{end} + 1;
	#if ($verbose){
	#  print STDERR "Re-SET POINTERS:\n";
	#  print STDERR "query : $query{start} -  $query{end}\n";
	#  print STDERR "target: $target{start} -  $target{end}\n";
	#}
	
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
	  if ($verbose){
	    print STDERR "GAP:\n";
	    print STDERR "query : $query{start} -  $query{end}\n";
		print STDERR "target: $target{start} -  $target{end}\n";
	  }
	}
      }
      ############################################################
      # intron
      if( $blocks[$i] eq 'I' ){
	if ( $in_exon ){
	  # emit the current exon
	  if ($verbose){
	    print STDERR "EXON: ".$exon->start."-".$exon->end."\n";
	  }
	  
	  $transcript->add_Exon($exon);
	  
	  # add the supporting features
	  my $supp_feature;
	  eval{
	    $supp_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@features);
	  };
	  print STDERR "intron: adding evidence : ".$supp_feature->gffstring."\n" if $verbose;
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
    
    if ($verbose){
      print STDERR "EXON: ".$exon->start."-".$exon->end."\n";
    }
    if ( scalar(@features) ){
      my $supp_feature;
      eval{
	$supp_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@features);
      };
      if ($@){
	print STDERR $@."\n";
      }
      print STDERR "outside: adding evidence : ".$supp_feature->gffstring."\n" if $verbose;
      $exon->add_supporting_features( $supp_feature );
    }
    else{
      print STDERR "No more features to add" if $verbose;
    }
    $transcript->add_Exon($exon);


    ############################################################
    # compute coverage
    # q_start reported by the sugar/cigar lines is one less 
    # as exonerate counts between lines
    my $aligned_length = $q_end - ($q_start + 1) + 1;
    my $coverage = sprintf "%.2f", 100 * ( $aligned_length - $target_gap_length ) / $length{$q_id};
    
    print STDERR "coverage = ( $aligned_length - $target_gap_length ) / ".$length{$q_id}." =  $coverage\n" if $verbose;
    
    foreach my $exon ( @{$transcript->get_all_Exons} ){
      foreach my $evidence ( @{$exon->get_all_supporting_features} ){
	$evidence->score($coverage);
      }
    }
        
    push( @transcripts, $transcript );
    
  } # end of while loop
  
  $self->output( @transcripts );

  close(EXO);  
    
  ############################################################
  # remove interim files (but do not remove the database if you are using one)
  unlink $query;
  if ( $self->genomic){
    unlink $target;
  }
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


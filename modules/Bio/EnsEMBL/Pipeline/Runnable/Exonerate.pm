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

Bio::EnsEMBL::Pipeline::Runnable::Exonerate

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Pipeline::Runnable::NewExonerate->new(
								 -query_seqs     => \@q_seqs,
                                                             [or -query_file     => $q_file]   
								 -query_type     => 'dna',
								 -target_seqs    => \@t_seqs,
                                                             [or -target_file    => $t_file]   
			                                         -target_type    => 'dna',
                                                                 -exonerate      => $exonerate,
								 -options        => $options,
								);

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
 where @results is an array of SeqFeatures, each one representing 
 an alignment (e.g. a transcript), and each feature contains a 
 list of alignment blocks (e.g. exons) as sub_SeqFeatures, which are
 in fact feature pairs.
 
=head1 DESCRIPTION

Exonerate takes a Bio::Seq (or Bio::PrimarySeq) object and runs 
Exonerate against a set of sequences.  The resulting output file 
is parsed to produce a set of features.

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

package Bio::EnsEMBL::Pipeline::Runnable::Exonerate;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Pipeline::RunnableI;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Pipeline::RunnableI);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($query_file, $query_seqs, $query_type,
      $target_file, $target_seqs, $target_type,
      $exonerate, $options, $verbose) = 
          $self->_rearrange([qw(
                                QUERY_FILE
                                QUERY_SEQS
                                QUERY_TYPE
                                TARGET_FILE
                                TARGET_SEQS
                                TARGET_TYPE
                                EXONERATE
                                OPTIONS
                                VERBOSE
			   )
			 ], @args);
  $self->_verbose($verbose) if $verbose;
  $self->{_output} = [];

  if (defined($query_seqs)) {     
    $self->throw("You must supply an array reference with -query_seqs") 
        if ref($query_seqs) ne "ARRAY";
    $self->query_seqs($query_seqs);
  }
  elsif (defined $query_file) {
    $self->throw("The given query file does not exist") if ! -e $query_file;
    $self->query_file($query_file);
  }
  
  if (defined($target_seqs)) {     
    $self->throw("You must supply an array reference with -target_seqs") 
        if ref($target_seqs) ne "ARRAY";
    $self->target_seqs($target_seqs);
  }
  elsif (defined $target_file) {
    $self->throw("The given database does not exist") if ! -e $target_file;
    $self->target_file($target_file);
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
  # We default exonerate-0.6.7
  if ($exonerate) {
    $self->exonerate($exonerate);
  } else {
    $self->exonerate('/usr/local/ensembl/bin/exonerate-0.6.7');
  }

  ############################################################
  # options
  my $basic_options = "--ryo \"RESULT: %S %pi %ql %tl %g %V\\n\" "; 

  # can add extra options as a string
  if ($options){
    $basic_options .= " ".$options;
  }
  $self->options($basic_options);

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

  # Set name of results file
  $self->results($self->workdir . "/results.$$");

  # Write query sequences to file if necessary
  $self->_write_sequences;

  # Build exonerate command

  my $command =$self->exonerate . " " .$self->options .
      " --querytype " . $self->query_type .
      " --targettype " . $self->target_type .
      " --query " . $self->query_file .
      " --target " . $self->target_file. " --forwardcoordinates FALSE";
  

  # Execute command and parse results

  print STDERR "Exonerate command : $command\n"
    if $self->_verbose;

  open( EXO, "$command |" ) || $self->throw("Error running exonerate $!");

  $self->parse_results(\*EXO);

  close(EXO);

  return 1
}

sub parse_results {
  my ($self, $fh) = @_;

  my %strand_lookup = ('+' => 1, '-' => -1, '.' => 0);

  # Each alignment will be stored as a transcript with 
  # exons and supporting features.  Initialise our
  # transcript.

  my @transcripts;

  # Parse output looking for lines beginning with 'RESULT:'.
  # Each line represents a distinct match to one sequence
  # containing multiple 'exons'.

 TRANSCRIPT:
  while (<$fh>){

    print STDERR $_ if $self->_verbose;

    next unless /^RESULT:/;

    chomp;

    my ($tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end,
	$t_strand, $score, $perc_id, $q_length, $t_length, $gene_orientation,
	@align_components) = split;
    
    # Increment all our start coordinates.  Exonerate has a 
    # coordinate scheme that counts _between_ nucleotides.
    
    $q_start++; $t_start++;

    $t_strand = $strand_lookup{$t_strand};
    $q_strand = $strand_lookup{$q_strand};
    $gene_orientation = $strand_lookup{$gene_orientation};

    # Read vulgar information and extract exon regions.
    my $proto_exons = $self->_extract_distinct_exons($t_start,
                                                     $t_strand,
                                                     $t_length,
						     $q_start, 
						     $q_strand,
                                                     $q_length,
						     \@align_components);

    # now we have extracted the exons and the coordinates are with 
    # reference to the forward strand of the query and target, we can 
    # use the gene_orienation to flip the strands if necessary
    if ($gene_orientation == -1 and $t_strand == 1) {
      $t_strand *= -1;
      $q_strand *= -1;
    }

    # Start building our pair of features.
    my %target = (percent => $perc_id,
		  source  => 'exonerate',
		  name    => $t_id,
		  strand  => $t_strand);

    my $transcript = Bio::EnsEMBL::Transcript->new();

    # Build FeaturePairs for each region of query aligned to a single
    # Exon.  Create a DnaDnaAlignFeature from these FeaturePairs and then
    # attach this to our Exon.

    my $aligned_query_residues = 0;

    foreach my $proto_exon (@$proto_exons){

      if(not defined $proto_exon->[0]){
        next;
      }

      my @feature_pairs;
      my $exon_end    = undef; # Where the overall exon end is tallied.
      my $exon_start  = undef; # Where the overall exon start is tallied.
      my $exon_phase      = $proto_exon->[0]->{phase}; 
      my $exon_end_phase  = $proto_exon->[0]->{end_phase}; 

      foreach my $exon_frag (@$proto_exon){
	$target{start} = $exon_frag->{target_start};
	$target{end}   = $exon_frag->{target_end};

        $aligned_query_residues += $exon_frag->{query_end} - $exon_frag->{query_start} + 1;

	# Watch for the ends of our Exon.
	$exon_end = $target{end}
	  if $target{end} > $exon_end;
	$exon_start = $target{start}
	  if (($target{start} < $exon_start)||(! $exon_start));

	# Build the query feature for our feature pair
	my %query = (name    => $q_id,
		     start   => $exon_frag->{query_start},
		     end     => $exon_frag->{query_end},
		     strand  => $q_strand,
		     percent => $perc_id,
		     source  => 'exonerate',);

	$self->throw("Hit coordinate < 1.")
	  if (($query{start} < 1) || ($query{end} < 1));

	my $feature_pair = 
	  $self->create_FeaturePair($target{start}, 
				    $target{end}, 
				    $target{strand},
				    $query{start},
				    $query{end},
				    $query{strand},
				    $query{name},
				    $query{score},
				    $query{percent},
				    0,
				    $target{name});


	push @feature_pairs, $feature_pair;
      }

      unless (@feature_pairs){
	$self->warn("Found an exonerate match that didn't " . 
		    "have any exonic regions.");
	next TRANSCRIPT
      }

      # Build our exon and set its key values.
      my $exon = Bio::EnsEMBL::Exon->new();

      $exon->seqname($t_id);
      $exon->start($exon_start);
      $exon->end($exon_end);
      $exon->phase($exon_phase);
      $exon->end_phase($exon_end_phase);
      $exon->strand($target{strand});

      # Use our feature pairs for this exon to create a single 
      # supporting feature (with cigar line).
      my $supp_feature;

      eval{
	$supp_feature =
	  Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@feature_pairs);
      };

      if ($@){
	$self->warn($@);
	next TRANSCRIPT;
      }

      $exon->add_supporting_features($supp_feature);

      # Finally, add Exon to growing transcript.
      $transcript->add_Exon($exon);
    }

    # Calculate query coverage, now that we have seen all matches.
    my $coverage = sprintf "%.2f",
	100 * ($aligned_query_residues) / $q_length;

    # Retrospectively set the coverage for each item of 
    # supporting evidence.
    foreach my $exon ( @{$transcript->get_all_Exons} ){
      foreach my $evidence ( @{$exon->get_all_supporting_features} ){
	$evidence->score($coverage);
      }
    }

    push @transcripts, $transcript 
      if scalar @{$transcript->get_all_Exons};
  }

  $self->output(@transcripts);

  return 1;
}

sub _extract_distinct_exons {
  my ($self, 
      $target_coord, $target_strand, $target_length,
      $query_coord,  $query_strand, $query_length,
      $vulgar_components) = @_;

  # This method works along the length of a vulgar line 
  # exon-by-exon.  Matches that comprise an exon are 
  # grouped and an array of 'proto-exons' is returned.
  # Coordinates from the vulgar line are extrapolated 
  # to actual genomic/query coordinates.

  my @exons;
  my $exon_number = 0;
  my $cumulative_query_coord  = $query_coord;
  my $cumulative_target_coord = $target_coord;

  my @last_vulgar_component;
  while (@$vulgar_components){
    $self->throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @$vulgar_components . "] items left to process.")
      unless scalar @$vulgar_components >= 3;

    my $type          = shift @$vulgar_components;
    my $query_match_length  = shift @$vulgar_components;
    my $target_match_length = shift @$vulgar_components;

    $self->throw("Vulgar string does not start with a match.  Was not " . 
		 "expecting this.")
      if ((scalar @exons == 0) && ($type ne 'M'));

    next if $type eq "S"; 
    # partial codon; these will be dealt with by "M" case

    if ($type eq 'M'){
      my %hash;

      $hash{type}         = $type;
      $hash{phase}        = 0;
      $hash{end_phase}    = 0;
      
      # adjust for split codons
      if (@last_vulgar_component and $last_vulgar_component[0] eq "S") {
        $query_match_length  += $last_vulgar_component[1];
        $target_match_length += $last_vulgar_component[2];
        # the "3 -" is necessary to convert to ensembl phase convention
        $hash{phase} = 3 - $last_vulgar_component[1];
      }
      if (@$vulgar_components and $vulgar_components->[0] eq "S") {
        $query_match_length  += $vulgar_components->[1];
        $target_match_length += $vulgar_components->[2];        
        # luckily, this one already fits into ensembl phase convention
        $hash{end_phase} = $vulgar_components->[1];
      }

      if ($target_strand != -1) {
        $hash{target_start} = $cumulative_target_coord;
        $hash{target_end}   = $cumulative_target_coord + $target_match_length - 1;
      } else {
        $hash{target_end}   = $target_length - $cumulative_target_coord + 1;
        $hash{target_start} = $hash{target_end} - $target_match_length + 1;
      }

      if ($query_strand != -1) {
	$hash{query_start} = $cumulative_query_coord;
	$hash{query_end}   = $cumulative_query_coord + $query_match_length - 1;
      } else {
	$hash{query_end}   = $query_length - $cumulative_query_coord + 1;
	$hash{query_start} = $hash{query_end} - $query_match_length + 1;
      }

      push @{$exons[$exon_number]}, \%hash;
    }

    $cumulative_target_coord += $target_match_length;
    $cumulative_query_coord += $query_match_length;

    # in protein mode, any insertion on the genomic side should be treated as 
    # an intron
    if ($type eq "I" or
        $type eq "F" or
        ($type eq "G" and 
         $target_match_length > 0 and
         $self->target_type eq "dna" and 
         $self->query_type eq "protein")) {
      $exon_number++;
    }

    @last_vulgar_component = ($type, $query_match_length, $target_match_length);
  }

  return \@exons;
}


sub _write_sequences {
  my ($self) = @_;

  if ($self->query_seqs) {
    my $query_file = $self->workdir . "/exonerate_q.$$";
    my $seqout = Bio::SeqIO->new('-format' => 'fasta',
                                 '-file'     => ">$query_file");   
    foreach my $seq ( @{$self->query_seqs} ) {
      $seqout->write_seq($seq);
    }
    # register the file for deletion
    $self->file($query_file);
    $self->query_file($query_file);
  }

  if ($self->target_seqs) {
    my $target_file = $self->workdir . "/exonerate_t.$$";
    my $seqout = Bio::SeqIO->new('-format' => 'fasta',
                                 '-file'   => ">$target_file");   
    foreach my $seq ( @{$self->target_seqs} ) {
      $seqout->write_seq($seq);
    }
    # register the file for later deletion
    $self->file($target_file);
    $self->target_file($target_file);
  } 

}

############################################################
#
# get/set methods
#
############################################################

sub query_seqs {
  my ($self, $seqs) = @_;
  if ($seqs){
    unless ($seqs->[0]->isa("Bio::PrimarySeqI") || $seqs->[0]->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_query_seqs} = $seqs;
  }
  return $self->{_query_seqs};
}

############################################################

sub target_seqs {
  my ($self, $seqs) = @_;
  if ($seqs){
    unless ($seqs->[0]->isa("Bio::PrimarySeqI") || $seqs->[0]->isa("Bio::SeqI")){
      $self->throw("query seq must be a Bio::SeqI or Bio::PrimarySeqI");
    }
    $self->{_target_seqs} = $seqs;
  }
  return $self->{_target_seqs};
}

############################################################

sub query_file {
    my ($self, $file) = @_;

    if ($file) {
	$self->{_query_file} = $file;
    }
    return $self->{_query_file};
}

############################################################

sub target_file {
    my ($self, $file) = @_;

    if ($file) {
	$self->{_target_file} = $file;
    }
    return $self->{_target_file};
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
      $self->throw("Not the right target type: $type");
    }
    $self->{_target_type} = $type ;
  }
  return $self->{_target_type};
}

############################################################


sub _verbose {
  my $self = shift;
  
  if (@_){
    $self->{_verbose} = shift;
  }
  
  return $self->{_verbose}
}


1;


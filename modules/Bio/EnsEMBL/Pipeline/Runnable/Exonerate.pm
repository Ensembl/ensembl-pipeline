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
$database  = a full path location for the directory containing 
             the target (genomic usually) sequence,
@sequences = a list of Bio::Seq objects,
$exonerate = a location for the binary,
$options   = a string with options ,

  my $runnable = 
    Bio::EnsEMBL::Pipeline::Runnable::Exonerate->new(
		     -database      => $database,
		     -query_seqs    => \@sequences,
		     -query_type    => 'dna',
		     -target_type   => 'dna',
		     -exonerate     => $exonerate,
		     -options       => $options);

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

  my ($database,
      $query_seqs,
      $query_type,
      $target_type,
      $exonerate,
      $options,
      $verbose) = $self->_rearrange([qw(
					DATABASE
					QUERY_SEQS
					QUERY_TYPE
					TARGET_TYPE
					EXONERATE
					OPTIONS
					VERBOSE
				       )
				    ], @args);

  $self->_verbose($verbose) if $verbose;

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
  # We default exonerate-0.6.7
  if ($exonerate) {
      $self->exonerate($exonerate);
  } else {
      $self->exonerate('/usr/local/ensembl/bin/exonerate-0.6.7');
  }

  ############################################################
  # options
  my $basic_options = "--ryo \"RESULT: %S %pi %ql %g %V\\n\" "; 

  # can add extra options as a string
  if ($options){
    $basic_options .= " ".$options;
  }
  $self->options($basic_options);

  return $self;
}

sub DESTROY {
  my $self = shift;

  unlink $self->_query_file;
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

  # Write query sequences to file

  $self->_write_sequences;

  # Build exonerate command

  my $command =$self->exonerate . " " .$self->options .
      " --querytype " . $self->query_type .
	" --targettype " . $self->target_type .
	  " --query " . $self->_query_file .
	    " --target " . $self->database;


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
	$t_strand, $score, $perc_id, $q_length, $gene_orientation,
	@align_components) = split;

    # Increment all our start coordinates.  Exonerate has a 
    # coordinate scheme that counts _between_ nucleotides.

    $q_start++; $t_start++;

    unless ($t_strand eq '+'){
      $self->warn("Target strand is not positive [$t_strand].  Was " . 
		  " not expecting this.");
      next TRANSCRIPT
    }

    # Start building our pair of features.
    my %target = (percent => $perc_id,
		  source  => 'exonerate',
		  name    => $t_id,
		  strand  => $strand_lookup{$t_strand},);

    # Hmm, what to do with the gene orientation?  Nothing for now.
    $gene_orientation = $strand_lookup{$gene_orientation};

    # Initialise the transcript feature that we are about to 
    # populate with exons.
    my $transcript = Bio::EnsEMBL::Transcript->new();

    # Read vulgar information and extract exon regions.

    my $q_coord = $q_strand eq '+' ? $q_start : $q_end;

    my $proto_exons = $self->_extract_distinct_exons($t_start, 
						     $q_coord, 
						     $q_strand,
						     \@align_components);

    # Build FeaturePairs for each region of query aligned to a single
    # Exon.  Create a DnaDnaAlignFeature from these FeaturePairs and then
    # attach this to our Exon.

    my $target_gaps = 0;

    foreach my $proto_exon (@$proto_exons){
      my @feature_pairs;
      my $exon_end    = undef; # Where the overall exon end is tallied.
      my $exon_start  = undef; # Where the overall exon start is tallied.
      my $prev_end; # This is used to count gaps in the target/genomic sequence.

      foreach my $exon_frag (@$proto_exon){
	$target{start} = $exon_frag->{target_start};
	$target{end}   = $exon_frag->{target_end};

	# Account for genomic sequence gaps.
	$target_gaps += $exon_frag->{target_start} - $prev_end - 1
	  if defined $prev_end;
	$prev_end = $target{end};

	# Watch for the ends of our Exon.
	$exon_end = $target{end}
	  if $target{end} > $exon_end;
	$exon_start = $target{start}
	  if (($target{start} < $exon_start)||(! $exon_start));

	# Build the query feature for our feature pair
	my %query = (name    => $q_id,
		     start   => $exon_frag->{query_start},
		     end     => $exon_frag->{query_end},
		     strand  => $strand_lookup{$q_strand},
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
      $exon->strand($strand_lookup{$t_strand});
      $exon->phase(0);
      $exon->end_phase(0);


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
	100 * ($q_end - $q_start + 1 - $target_gaps) / $q_length;

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

  return 1
}

sub _extract_distinct_exons {
  my ($self, $target_start, $query_coord, 
      $query_strand, $vulgar_components) = @_;

  # This method works along the length of a vulgar line 
  # exon-by-exon.  Matches that comprise an exon are 
  # grouped and an array of 'proto-exons' is returned.
  # Coordinates from the vulgar line are extrapolated 
  # to actual genomic/query coordinates.

  my @exons;
  my $exon_number = 0;
  my $cumulative_query_coord  = $query_coord;
  my $cumulative_target_coord = $target_start;

  while (@$vulgar_components){
    $self->throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @$vulgar_components . "] items left to process.")
      unless scalar @$vulgar_components >= 3;

    my $type          = shift @$vulgar_components;
    my $query_length  = shift @$vulgar_components;
    my $target_length = shift @$vulgar_components;

    $self->throw("Vulgar string does not start with a match.  Was not " . 
		 "expecting this.")
      if ((scalar @exons == 0) && ($type ne 'M'));

    if ($type eq 'M'){
      my %hash;

      $hash{type}         = $type;
      $hash{target_start} = $cumulative_target_coord;
      $hash{target_end}   = $cumulative_target_coord + $target_length - 1;

      if ($query_strand ne '-') {
	$hash{query_start} = $cumulative_query_coord;
	$hash{query_end}   = $cumulative_query_coord  + $query_length - 1;
      } else {
	$hash{query_start} = $cumulative_query_coord - $query_length + 1;
	$hash{query_end}   = $cumulative_query_coord;
      }

      push @{$exons[$exon_number]}, \%hash;
    }

    $cumulative_target_coord += $target_length;

    if ($query_strand ne '-') {
      $cumulative_query_coord  += $query_length;
    } else {
      $cumulative_query_coord  -= $query_length;
    }

    $exon_number++ 
      if $type eq 'I';
  }

  return \@exons;
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
    $self->throw("Genomic sequence database file [" . 
		 $self->{_database} . "] does not exist.")
      unless -e $self->{_database};
  }

  $self->throw("No genomic sequence database provided for Exonerate.")
    unless $self->{_database};

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
      $self->throw("Not the right target type: $type");
    }
    $self->{_target_type} = $type ;
  }
  return $self->{_target_type};
}

############################################################


sub _write_sequences {
  my $self = shift;

  my $seqout = Bio::SeqIO->new('-file'   => ">" . $self->_query_file,
			       '-format' => 'Fasta',
			      );

  foreach my $query_seq ( $self->query_seqs ){
    $seqout->write_seq($query_seq);
  }
}

sub _query_file {
  my $self = shift;

  if (@_) {
    $self->{_query_file} = shift;
  }

  unless ($self->{_query_file}){
    $self->{_query_file} = $self->workdir . "/query_seqs.$$"
  }

  return $self->{_query_file};
}

sub _verbose {
  my $self = shift;

  if (@_){
    $self->{_verbose} = shift;
  }

  return $self->{_verbose}
}


1;


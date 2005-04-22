=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::RunnableDB::DuplicationFinder;

=head1 SYNOPSIS

  # Almost all config is set in Bio::EnsEMBL::Pipeline:Config::GeneDupl

  my $df = 
    Bio::EnsEMBL::Pipeline::RunnableDB::DuplicationFinder->new(
      -input_id   => 'ENSG00000183362',
      );

  $df->fetch_input();
  $df->run();
  $df->write_output();

=head1 DESCRIPTION

RunnableDB used to drive the identification of putative gene paralogues.  The
matching Runnable module is Bio::EnsEMBL::Pipeline::GeneDuplication::Finder.

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::RunnableDB::DuplicationFinder;

use strict;
use Bio::EnsEMBL::Pipeline::RunnableDB;
use Bio::EnsEMBL::Pipeline::GeneDuplication::Finder;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Pipeline::Config::GeneDupl qw(GD_OPTIONS
						GD_BLAST_EXE
						GD_BLAST_VARIANT
						GD_CODEML_EXE
						GD_WORK_DIR
						GD_DISTANCE_METHOD
						GD_OUTPUT_METHOD
						GD_OUTPUT_DBNAME
						GD_OUTPUT_DBHOST
						GD_OUTPUT_DBUSER
						GD_OUTPUT_DBPASS
						GD_OUTPUT_DBPORT
						GD_OUTPUT_TYPE);

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Pipeline::RunnableDB);

sub new {
  my ($class, @args) = @_;

  # Object creation is not done via the parent class to 
  # avoid the need for a database adaptor and analysis
  # object to be passed at construction.  These are not
  # needed if the output format is text.

  my $self = {};
  bless $self, $class;

  my ($db,
      $input_id, 
      $seqfetcher, 
      $analysis) = &rearrange([qw(DB
				  INPUT_ID
				  SEQFETCHER
				  ANALYSIS )], @args);

  $self->input_id($input_id) if defined $input_id;
  $self->db($db)             if defined $db;
  $self->seqfetcher          if defined $seqfetcher;
  $self->analysis            if defined $analysis;

  return $self;
}

sub fetch_input {
  my $self = shift;

  # Check that an input_id exists.

  unless (defined $self->input_id) {
    throw("Input_id is unset.  Unable to proceed without it.")
  }

  # Clean up after any previous runs.

  $self->output(undef);

  # Set options, if this appears necessary.
  # If the input_id seems to come from a different species to that 
  # already configured, set config options relevant to new 
  # input_id/species.

  unless (defined $self->input_id && 
	  defined $self->_ingroup_regex && 
	  $self->input_id =~ /^$self->_ingroup_regex/) {
    $self->_set_options;
    $self->_verify_options;
  }

  # Build runnable

  $self->_build_runnable;

  # Retrieve input sequence from blast database.  Throw an
  # error if the sequence is not present.

  my $runnable = ($self->runnable)[0];

  unless ($self->input_id){
    throw("Unable to run without an input id.")
  }

  my $input_seq = 
    $runnable->_seq_fetcher->fetch($self->input_id);

  unless ($input_seq){
    throw("Failed to fetch input sequence [". $self->input_id ."].");
  }

  $self->_input_seq($input_seq);

  return 1
}

sub run{
  my $self = shift;

  my $runnable = ($self->runnable)[0];

  my $result;

### Annoying loop added to retry paml runs if they fail. Remove once paml is replaced.
 RETRY:
  for (my $i = 0; $i < 3; $i++) {
###

  eval {
    $result = $runnable->run($self->_input_seq);
  };

  if ($@){
    next RETRY
  }

  last
}

  if ($@) {
    # Trying to work around a bloody PAML bug:
    if ($self->_identical_seqs_bug($runnable->alignment)){
      warning("Identical sequences causing mayhem with PAML run.")
    } else {
      throw ("Problem running PAML :\n$@")
    }
  }

  if ($result){
    $self->output($result);
  }

  return 1
}

sub write_output{
  my $self = shift;

  if ($GD_OUTPUT_METHOD eq 'files') {
    $self->_write_output_as_text;
  } elsif ($GD_OUTPUT_METHOD eq 'db') {
    $self->_write_output_to_database;
  } else {
    throw("Unrecognised output option in Config : " . 
      "GD_OUTPUT_METHOD = [$GD_OUTPUT_METHOD]")
  }

  return 1
}

sub output {
  my $self = shift;

  if (@_) {
    $self->{_output} = shift;
  }

  return $self->{_output}
}

sub _set_options {
  my $self = shift;

  my $id_prefix;

  foreach my $prefix (keys %$GD_OPTIONS){
print STDERR "Input id : " .$self->input_id . " Prefix : " . $prefix . "\n";
    if ($self->input_id =~ /^$prefix/) {
      $id_prefix = $prefix;
      last
    }
  }

  unless (defined $id_prefix) {
    throw("Input sequence id [".$self->input_id."] does not match any " . 
	  "regex/id-prefix in config file.")
  }

  $self->{_blastdb_file}     = $GD_OPTIONS->{$id_prefix}->{GD_BLASTDB_FILE};
  $self->{_hit_coverage}     = $GD_OPTIONS->{$id_prefix}->{GD_HIT_COVERAGE};
  $self->{_hit_identity}     = $GD_OPTIONS->{$id_prefix}->{GD_HIT_IDENTITY};
  $self->{_distance_cutoff}  = $GD_OPTIONS->{$id_prefix}->{GD_DISTANCE_CUTOFF};
  $self->{_ingroup_regex}    = $GD_OPTIONS->{$id_prefix}->{GD_INGROUP_REGEX};
  $self->{_outgroup_regexes} = $GD_OPTIONS->{$id_prefix}->{GD_OUTGROUP_REGEXES};
  $self->{_transl_regex}     = $GD_OPTIONS->{$id_prefix}->{GD_TRANSL_REGEX};
  $self->{_genetic_code}     = $GD_OPTIONS->{$id_prefix}->{GD_GENETIC_CODE};
  $self->{_output_dir}       = $GD_OPTIONS->{$id_prefix}->{GD_OUTPUT_DIR};

  return 1
}

sub _verify_options {
  my $self = shift;

  # Verify options set in config file.

    # Check options for output type are compatible

  if ($GD_OUTPUT_METHOD eq 'files') {
    unless (-d $self->_output_dir) {
      throw("Output directory [" . $self->_output_dir . "] appears not " .
	    "to exist or is inaccessible.")
    }	
  } elsif ($GD_OUTPUT_METHOD eq 'db') {
    unless (defined $GD_OUTPUT_DBNAME && 
	    defined $GD_OUTPUT_DBHOST && 
	    defined $GD_OUTPUT_DBUSER) {
      throw("Details of output database need to be provided " .
	    "in order for output to be written to it.")
    }

    my $outdb = 
      Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
        -host   => $GD_OUTPUT_DBHOST,
	-dbname => $GD_OUTPUT_DBNAME,
	-port   => $GD_OUTPUT_DBPORT,
	-user   => $GD_OUTPUT_DBUSER,
	-pass   => $GD_OUTPUT_DBPASS);

    unless (defined $outdb && 
	    $outdb->isa("Bio::EnsEMBL::Compara::DBSQL::DBAdaptor")){
      throw("Unable to connect to output database using options " .
	    "specified in config ")
    }
  } else {
    throw("Unrecognized output method [$GD_OUTPUT_METHOD].  " .
	  "Accepted options are \'files\' or \'db\'.")
  }

    # Check BLAST-related options.

#  unless (defined $self->_blastdb_file && -e $self->_blastdb_file) {
#    throw("Input BLAST database is either not specified or does " . 
#	  "not exist [" . $self->_blastdb_file . "].")
#  }

  if (defined $GD_BLAST_EXE && $GD_BLAST_EXE ne '' && 
      ! -x $GD_BLAST_EXE) {
    throw("Blast binary specified is not executable [$GD_BLAST_EXE]")
  }

  unless (defined $GD_BLAST_VARIANT &&
     ($GD_BLAST_VARIANT eq 'ncbi' || 
      $GD_BLAST_VARIANT eq 'wu_new' || 
      $GD_BLAST_VARIANT eq 'wu_old' )) {
    throw("Unrecognised BLAST index type [$GD_BLAST_VARIANT].  Should " .
	  "be on of either \'ncbi\', \'wu_new\' or \'wu_old\'")
  }

    # Check codeml, if specified

  if (defined $GD_CODEML_EXE && $GD_CODEML_EXE ne '' 
      && ! -x $GD_CODEML_EXE) {
    throw("Codeml binary specified is not executable [$GD_CODEML_EXE]")
  }

    # Check user regexes.

  unless (defined $self->_ingroup_regex && $self->_ingroup_regex ne '') {
    throw("No ingroup regular expression specified in config")
  }

  if (defined $self->_outgroup_regexes && ! scalar @{$self->_outgroup_regexes}) {
    throw("Output regular expressions should be specified as " . 
	  "an anonymous array of strings.")
  }

  unless (defined $self->_transl_regex && $self->_transl_regex ne '') {
    throw("No regular expression for matching translation stable " . 
	  "ids specified in config")
  }

    # Warnings if coverage, identity and genetic distance cutoffs not
    # set.

  unless (defined $self->_hit_coverage) {
    warning("Hit coverage cutoff has not been set, using default " .
	    "value (which might be anything).")
  }

  unless (defined $self->_hit_identity) {
    warning("Hit identity cutoff has not been set, using default " .
	    "value (which might be anything).")
  }

  unless (defined $self->_distance_cutoff) {
    warning("Genetic distance cutoff has not been set, using default " .
	    "value (which might be anything).")
  }

    # Check genetic code.

  unless (defined $self->_genetic_code) {
    throw("Genetic code has not been set in config file.  Game over.")
  }

    # Check whether to return output alignments in nucleotide or amino 
    # acid coordinates.

  unless (defined $GD_OUTPUT_TYPE && 
	  ($GD_OUTPUT_TYPE eq 'nucleotide' ||
	   $GD_OUTPUT_TYPE eq 'aminoacid')) {
    throw("Must specify mode for result presentation.  Options are " . 
	  "\'nt\' and \'aa\', but you said [$GD_OUTPUT_TYPE]")
  }

  return 1
}

sub _blastdb_file {
  my $self = shift;

  return $self->{_blastdb_file}
}

sub _hit_coverage {
  my $self = shift;

  return $self->{_hit_coverage}
}

sub _hit_identity {
  my $self = shift;

  return $self->{_hit_identity}
}

sub _distance_cutoff {
  my $self = shift;

  return $self->{_distance_cutoff}
}

sub _ingroup_regex {
  my $self = shift;

  return $self->{_ingroup_regex}
}

sub _outgroup_regexes {
  my $self = shift;

  return $self->{_outgroup_regexes}
}

sub _transl_regex {
  my $self = shift;

  return $self->{_transl_regex}
}

sub _genetic_code {
  my $self = shift;

  return $self->{_genetic_code}
}

sub _output_dir {
  my $self = shift;

  return $self->{_output_dir}
}

sub _build_runnable {
  my $self = shift;

  my $gene_dupl 
    = Bio::EnsEMBL::Pipeline::GeneDuplication::Finder->new(
	 '-blastdb'                => $self->_blastdb,
         '-blast_program'          => $GD_BLAST_EXE,
         '-blast_index_type'       => $GD_BLAST_VARIANT,
         '-codeml'                 => $GD_CODEML_EXE,
	 '-hit_coverage'           => $self->_hit_coverage,
	 '-hit_identity'           => $self->_hit_identity,
	 '-distance_cutoff'        => $self->_distance_cutoff,
	 '-regex_query_species'    => $self->_ingroup_regex,
	 '-regex_outgroup_species' => $self->_outgroup_regexes,
	 '-genetic_code'           => $self->_genetic_code,
	 '-work_dir'               => $GD_WORK_DIR,
	 '-distance_method'        => $GD_DISTANCE_METHOD,
	);

  $self->runnable($gene_dupl);
}

sub _write_output_as_text {
  my $self = shift;

  my $runnable = ($self->runnable)[0];

  my $result = $self->output;

  # Open output file

  my $output_filename = $self->_make_output_file;

  open (OUT, ">$output_filename") 
    or throw("Cannot write to file [$output_filename]");


  # Return early if there were no accepted matches.
  unless (defined $result) {
    print OUT $self->input_id . " : No homologous matches were " .
      "found that satisfied the match criteria.";
###
print STDOUT $self->input_id . " : No homologous matches were " . "found that satisfied the match criteria.";
###
    return 1
  }

  # Work though our results, generating vital stats where necessary.

  foreach my $match (@{$result}) {

    if ($GD_OUTPUT_TYPE eq 'nucleotide') {

### Must put pair align stuff back in here!!!!
# Need similarity, cigars, starts, ends

      print OUT join("\t",
		     $match->{ka},
		     $match->{ks},
		     $self->_distance_cutoff,
		     $match->{alignment}->[0]->display_id,
		     $self->_build_cigar($match->{alignment}->[0]),
		     $match->{alignment}->[1]->display_id,
		     $self->_build_cigar($match->{alignment}->[1]),
		     $match->{identity},
		     $match->{coverage}) . "\n";
    } else {

      my $transl_regex = $self->_transl_regex;

      $match->{alignment}->[0]->desc =~ /($transl_regex\w*)/;
      my $query_translation_id = $1;
      $match->{alignment}->[1]->desc =~ /($transl_regex\w*)/;
      my $match_translation_id = $1;

      unless ($query_translation_id ne '' && $match_translation_id ne '') {
	warning("Unable to determine translation stable ids.\n" . 
		"Translation stable id regex is [" . $transl_regex . "]\n" . 
		"Query seq desc line is [" . $match->{alignment}->[0]->desc . "]\n" .
		"Match seq desc line is [" . $match->{alignment}->[1]->desc . "]")
      }

      my $pair_align = $self->_alignment_vitals_aa(@{$match->{alignment}});

      print OUT join ("\t",
		      $match->{ka},
		      $match->{ks},
		      0, #$match->{N},
		      0, #$match->{S},
		      0, #$match->{lnL},
		      $self->_distance_cutoff,
		      $match->{alignment}->[0]->display_id,
		      $query_translation_id,
		      $pair_align->{_query_cigar},
		      $pair_align->{_query_start},
		      $pair_align->{_query_end},
		      $pair_align->{_query_coverage},
		      $pair_align->{_query_identity},
		      $pair_align->{_query_similarity},
		      $match->{alignment}->[1]->display_id,
		      $match_translation_id,
		      $pair_align->{_match_cigar},
		      $pair_align->{_match_start},
		      $pair_align->{_match_end},
		      $pair_align->{_match_coverage},
		      $pair_align->{_match_identity},
		      $pair_align->{_match_similarity})
	. "\n";

    }

  }

  close(OUT);

  return 1
}

sub _write_output_to_database {
  my $self = shift;

  throw("Output to database mode is not yet implemented.");
}

sub _alignment_vitals_aa {
  my ($self, $query, $match) = @_;

  my %similarity_lookup 
    = ('A' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'V' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'F' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'P' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'M' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'I' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'L' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'W' => ['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'],
       'D' => ['D', 'E'],
       'E' => ['D', 'E'],
       'R' => ['R', 'H', 'K'],
       'K' => ['R', 'H', 'K'],
       'S' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'],
       'T' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'],
       'Y' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'],
       'H' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q', 'R', 'H', 'K'],
       'C' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'],
       'N' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'],
       'G' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'],
       'Q' => ['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'],
       'B' => ['B'],
       'Z' => ['Z']);

  $query->alphabet('DNA');
  $match->alphabet('DNA');

  my $query_protein 
    = $query->translate(undef, undef, undef, $self->_genetic_code);

  my $match_protein 
    = $match->translate(undef, undef, undef, $self->_genetic_code);

  my %aligned_pair;

  $aligned_pair{_query_cigar} = $self->_build_cigar($query_protein);
  $aligned_pair{_match_cigar} = $self->_build_cigar($match_protein);

  $aligned_pair{_query_start} = 1;
  $aligned_pair{_match_start} = 1;

  my @query = split //, $query_protein->seq;
  my @match = split //, $match_protein->seq;

  if (scalar @query != scalar @match) {
    throw("Aligned sequences are not the same length")
  }

  my $query_length = 0;
  my $match_length = 0;
  my $identical_aa = 0;
  my $covered_aa   = 0;
  my $similar_aa   = 0;

  for (my $i = 0; $i < scalar @query; $i++){
    $query_length++ if $query[$i] ne 'X';
    $match_length++ if $match[$i] ne 'X';

    $covered_aa++ if (($query[$i] ne 'X') and ($match[$i] ne 'X'));

    if (($query[$i] eq $match[$i]) and
	($query[$i] ne 'X')){
      $identical_aa++;
    } elsif (($query[$i] ne 'X') and
	     (grep /$match[$i]/, @{$similarity_lookup{$query[$i]}})) {
      $similar_aa++;
    }

  }

  $aligned_pair{_query_end}        = $query_length;
  $aligned_pair{_match_end}        = $match_length;
  $aligned_pair{_query_identity}   = sprintf "%3.2f", ($identical_aa/$covered_aa)*100;
  $aligned_pair{_match_identity}   = sprintf "%3.2f", ($identical_aa/$covered_aa)*100;
  $aligned_pair{_query_coverage}   = sprintf "%3.2f", ($covered_aa/$query_length)*100;
  $aligned_pair{_match_coverage}   = sprintf "%3.2f", ($covered_aa/$match_length)*100;
  $aligned_pair{_query_similarity} = sprintf "%3.2f", (($similar_aa+$identical_aa)/$covered_aa)*100;
  $aligned_pair{_match_similarity} = sprintf "%3.2f", (($similar_aa+$identical_aa)/$covered_aa)*100;

  return \%aligned_pair;
}

sub _build_cigar {
  my ($self, $gapped_seq) = @_;

  my @seq = split //, $gapped_seq->seq;

  my @cigar;
  my $mode = 'M';
  $mode = 'D' 
    if ($seq[0] eq 'X' || $seq[0] eq '-');
  my $prev_mode = $mode;
  my $pointer_pusher = 0;

  foreach my $unit (@seq){
    if ($unit eq 'X' || $unit eq '-'){
      $mode = 'D'
    } else {
      $mode = 'M'
    }

    if ($mode ne $prev_mode){
      $pointer_pusher++;
      $prev_mode = $mode;
    }

    $cigar[$pointer_pusher]->{mode} = $mode
      unless $cigar[$pointer_pusher]->{mode};

    $cigar[$pointer_pusher]->{length}++;
  }

  my $cigar = '';
  foreach my $cigarette (@cigar){
    if ($cigarette->{length} > 1){
      $cigar .=  $cigarette->{length} . $cigarette->{mode};
    } else {
      $cigar .= $cigarette->{mode};
    }
  }

  return $cigar;
}

sub _identical_seqs_bug {
  my ($self, $aligned_seqs) = @_;

  return 0 unless $aligned_seqs;

  my $verdict = 0;

  for (my $i = 0; $i < scalar @$aligned_seqs; $i++){
    for (my $j = $i + 1; $j < scalar @$aligned_seqs; $j++){

      next unless $aligned_seqs->[$i]->seq eq $aligned_seqs->[$j]->seq;

      $verdict = 1;
      last
    }
  }

  return $verdict
}

sub _make_output_file {
  my $self = shift;

  unless (-d $self->_output_dir) {
    throw("Output directory [" . $self->_output_dir . "] does not exist.")
  }

  # Parse off the last two digits (or at a pinch, the last two alphanumeric
  # characters) from our input id.  This will be used as an output directory
  # name (so that the output directories are finely grained and managable).
  $self->input_id =~ /(\w)[\W\_]?(\w)$/;

  my $sub_dir = $self->_output_dir . "/" . $1 . $2;
  $sub_dir =~ s/\/\//\//g;  # !

  unless (-d $sub_dir) {
    my $success = mkdir $sub_dir;

    unless ($success) {
      throw("Can not make sub-directory [$sub_dir]")
    }
  }

  return $sub_dir . '/' .$self->input_id . ".txt";
}

sub _blastdb {
  my $self = shift;

  unless (defined $self->{_blastdb}){
    $self->{_blastdb} = 
      Bio::EnsEMBL::Pipeline::Runnable::BlastDB->new(
        -dbfile        => $self->_blastdb_file,
        -molecule_type => 'dna',
        -index_type    => 'wu_new');

    $self->{_blastdb}->make_fetchable_index(1);
    $self->{_blastdb}->db_formatted(1);
  }

  return $self->{_blastdb}
}

sub _input_seq {
  my $self = shift;

  if (@_) {
    $self->{_input_seq} = shift;

    unless ($self->{_input_seq} && 
	    $self->{_input_seq}->isa("Bio::Seq")) {
      throw("Input sequence is not a Bio::Seq.  It is a [" . 
	    $self->{_input_seq} . "]")
    }
  }

  return $self->{_input_seq}
}

1;

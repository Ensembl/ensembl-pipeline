package Bio::EnsEMBL::Pipeline::GeneDuplication::Finder;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB;
use Bio::EnsEMBL::Pipeline::GeneDuplication::PAML;
use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;
use Bio::EnsEMBL::Pipeline::GeneDuplication::Result;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);


my $DEFAULT_DISTANCE_METHOD  = 'NeiGojobori';
my $DEFAULT_BLAST_PROGRAM    = 'wublastn';
my $DEFAULT_DISTANCE_CUTOFF  = 1.000;

@ISA = qw(Bio::EnsEMBL::Root Bio::EnsEMBL::Pipeline::RunnableI);

### Constructor ###

=head2 new

  Args       : -dbfile       - (string) Full path to a fasta formatted file of 
                               nucleotide sequences.
               -blastdb      - A Bio::EnsEMBL::Pipeline::Runnable::BlastDB 
                               object for which the run method has been invoked.
               -query_seq    - A Bio::Seq which is comprised of nucleotides.
               -blast_program           - Manually set the blast program (and
                               full path) that should be used.  Make sure the
                               program matches the index type chosen.  Defaults
                               to wublastn.
               -blast_index_type        - The distribution of blast to use.  See 
                               Bio::EnsEMBL::Pipeline::Runnable::BlastDB->index_type
                               documentation.  Defaults to 'wu_new'.
               -hit_identity - (optional) Hit identity.  Defaults to 0.80
               -hit_coverage - (optional) Hit coverage.  Defults to 0.80
               -work_dir     - (optional) Dir where working files are 
                                          placed.  Defaults to /tmp.
               -codeml       - Full path to codeml executable.
               -genetic_code - 0 for Universal, 1 for Mitochondrial.  See
                               the Bio::Seq->translate method for a full 
                               list of options.
               -regex_query_species     - A regular expression that will parse
                               some portion of the ids of the query species (and
			       not the ids of outgroup species).  Eg. 'ENSG'.
               -regex_outgroup_species  - A ref to an array of regular expressions
                               that will parse the ids of the various outgroup 
                               species (but not the ids of the query species).
                               E.g. ['ENSMUSG', 'ENSRNOG']
               -distance_cutoff         - A genetic distance cutoff to apply for
                               occasions where an outgroup species is not set, or
                               an outgroup species match was not found.
               -distance_method         - The method used to calculate the genetic
                               distance and Ka/Ks ratios.  Options are 'NeiGojobori'
                               and 'ML', ML being the maximum likelihood method
                               of Yang.
  Example    : none
  Description: Constructs new object
  Returntype : Bio::EnsEMBL::Pipeline::GeneDuplication::Finder
  Exceptions : Throws if database file not specified or does not exist.
  Caller     : General

=cut

sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;

  my ($query,
      $blastdb,
      $blast_program,
      $blast_index_type,
      $work_dir,
      $codeml,
      $genetic_code,
      $regex_query_species,
      $regex_outgroup_species,
      $identity_cutoff,
      $coverage_cutoff,
      $distance_cutoff,
      $distance_method) = rearrange([qw(QUERY
					BLASTDB
					BLAST_PROGRAM
					BLAST_INDEX_TYPE
					WORK_DIR
					CODEML_EXECUTABLE
					GENETIC_CODE
					REGEX_QUERY_SPECIES
					REGEX_OUTGROUP_SPECIES
					HIT_IDENTITY
					HIT_COVERAGE
					DISTANCE_CUTOFF
					DISTANCE_METHOD
				       )],@args);

  $self->_work_dir($work_dir) if $work_dir;

  if ($blastdb && $blastdb->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB")){
    $self->_blastdb($blastdb);
  } else {
    throw ("Need a Bio::EnsEMBL::Pipeline::Runnable::BlastDB object.");
  }

  $self->_query_seq($query)                               if $query;
  $self->_codeml($codeml)                                 if $codeml;
  $self->_genetic_code($genetic_code);
  $self->_regex_query_species($regex_query_species)       if $regex_query_species;
  $self->_regex_outgroup_species($regex_outgroup_species) if $regex_outgroup_species;
  $self->_identity_cutoff($identity_cutoff)               if $identity_cutoff;
  $self->_coverage_cutoff($coverage_cutoff)               if $coverage_cutoff;
  $self->_blast_program($blast_program)                   if $blast_program;

  $self->_identity_cutoff(80) unless $self->_identity_cutoff;
  $self->_coverage_cutoff(80) unless $self->_coverage_cutoff;

  $self->_distance_method($distance_method) if $distance_method;
  $self->_distance_cutoff($distance_cutoff) if $distance_cutoff;

  return $self
}


### Public methods ###


=head2 run

  Args[1]    : [optional] Input sequence (Bio::Seq)
  Example    : none
  Description: Top level method that executes the gene duplicate 
               finding algorithm.
  Returntype : Bio::EnsEMBL::Pipeline::GeneDuplication::Result
  Exceptions : Warns if PAML run fails.
               Throws if PAML returns multiple result sets (unlikely).
  Caller     : General

=cut

sub run {
  my ($self, $input_seq) = @_;

  $self->_query_seq($input_seq)
    if $input_seq;

  # Derive a list of sequences that look like duplications of 
  # the query sequence.

  my $accepted_hits 
    = $self->_find_recent_duplications($self->_query_seq);

  unless (scalar @$accepted_hits > 1) {
    print "No homologous matches were found that satisfied the match criteria.\n";
    return 0
  }

  return $accepted_hits;
}


### Hit Chooser Methods ###

=head2 _find_recent_duplications

  Args[1]    : Bio::Seq query sequence
  Example    : none
  Description: This is the implementation of the main algorithm of 
               this module.
  Returntype :
  Exceptions :
  Caller     : General

=cut

sub _find_recent_duplications {
  my ($self, $query) = @_;

  # Perform blast search with this query sequence.
  my $bplite_report = $self->_run_blast($query);

  # Take blast report and filter for hits with good coverage
  # and sufficiently high identity.
  my $good_hits = $self->_preliminary_filter($bplite_report);

  # Check that good hits have been found.
  unless (scalar(keys %$good_hits)){
    print "Did not find hits to query sequence.\n";
    return []
  }

  # From the set of promising hits, which will possibly include
  # hits to other species, filter using outgroup sequences (if
  # any) and the genetic distance between each hit and the query
  # sequence.
  my $recent_duplications = $self->_phylogenetic_filter($good_hits);

  # Return a complex data object that contains the pairwise
  # match details for each final accepted hit.
  return $recent_duplications;
}


=head2 _run_blast

  Args[1]    : Bio::Seq query sequence
  Example    : none
  Description:
  Returntype :
  Exceptions :
  Caller     : General

=cut

sub _run_blast {
  my ($self, $query) = @_;

  # We are looking for duplicates of the query gene
  # in our blast database.
  $self->_query_seq($query) if $query;

  my $bplite_report;

  eval{
    $bplite_report = $self->_blast_obj->run;
  };

  throw ("Blast did not run successfully.  Blast program was [".
	       $self->_blast_program."].  Index type was [".
	       $self->_blast_index_type."].")
    if $@;

  throw ("Blast process did not return a report.")
    unless ($bplite_report->isa("Bio::Tools::BPlite"));

  return $bplite_report
}


=head2 _preliminary_filter

  Args[1]    :
  Example    : none
  Description:
  Returntype :
  Exceptions :
  Caller     : General

=cut

sub _preliminary_filter {
  my ($self, $bplite_report) = @_;

  # Process our blast report.  
  #
  # For each blast hit:
  #   * throw away self matches (if any)
  #   * filter hits by identity and coverage
  #   * align WHOLE genes of promising hits, re-filter by 
  #       identity and coverage
  #   * calculate genetic distance between each subject and the query

  my %good_hits;

 DISTINCT_HIT:
  while (my $sbjct = $bplite_report->nextSbjct){

    # Mangle the BPLite::Sbjct object for its own good.  Quite 
    # often the hit ids parsed by BPLite include the whole 
    # Fasta header description line.  This is problematic if 
    # sequence ids need to be compared or a hash is keyed on 
    # this sequence id.  Here we simply lop the description
    # from each Sbjct->name, if there is one.
    $sbjct = $self->_fix_sbjct($sbjct);

    # It appears that the BPLite::Sbjct object only allows 
    # HSPs to be accessed once (as this process is closely 
    # tied to the parsing of the Blast report).  Hence, here
    # we loop through them all here and store them in an 
    # array.

    my @hsps;

    while (my $hsp = $sbjct->nextHSP) {
      push (@hsps, $hsp);
    }

    # Skip hit if it is a match to self.

    next DISTINCT_HIT
      if ($self->_query_seq->display_id eq $sbjct->name);

    # First, filter hits by identity

    my $hit_identity = $self->_hit_identity(\@hsps);

    next DISTINCT_HIT
      unless ($hit_identity >= $self->_identity_cutoff);

    # Second, filter hits by coverage

    next DISTINCT_HIT
      unless ($self->_appraise_hit_coverage($sbjct, \@hsps));

    # Third, this has become a serious hit.  Make a pairwise 
    # alignment of the query and hit.

    my @seqs = ($self->_fetch_seq($self->_query_seq->display_id),
		$self->_fetch_seq($sbjct->name));

    throw ("Didnt correctly obtain two sequences for alignment.")
      unless scalar @seqs == 2;

    my $align = $self->_pairwise_align(\@seqs);

      # and filter by coverage and identity for the whole 
      # aligned sequences.

    my ($min_coverage, $identity) 
      = $self->_calculate_coverage_and_identity($align);

    unless ($min_coverage >= $self->_coverage_cutoff &&
	    $identity     >= $self->_identity_cutoff) {

      next DISTINCT_HIT;
    }

    # Forth, calculate the genetic distance while we are here.

    my ($ka, $ks, $n, $s) 
      = $self->_run_pairwise_paml($align);

    # Finally, store all this useful information.  Perhaps a proper
    # object would make this easier to access later?

die "undefined N and S" unless defined $n and defined $s;

    $good_hits{$sbjct->name} = {'Ka'       => $ka,
				'Ks'       => $ks,
				'N'        => $n,
				'S'        => $s,
				'coverage' => $min_coverage,
				'identity' => $identity,
				'alignment'=> $align};

  }

    return \%good_hits
}


=head2 _phylogenetic_filter

  Args[1]    : 
  Example    : none
  Description:
  Returntype :
  Exceptions :
  Caller     : General

=cut

sub _phylogenetic_filter {
  my ($self, $good_hits) = @_;

  my %hits_by_species;

 HIT:
  foreach my $hit_id (keys %$good_hits) {

    # Now, partition hits according to their species.  The species
    # from which the subject is derived is determined by a
    # regular expression match to the sequence id.

    my $query_regex = $self->_regex_query_species;

    if ($hit_id =~ /$query_regex/) {

      push(@{$hits_by_species{$self->_regex_query_species}}, $hit_id);

      next HIT;

    } else {

      foreach my $regex (@{$self->{_regex_outgroup_species}}) {
	if ($hit_id =~ /$regex/){
	  push (@{$hits_by_species{$regex}}, $hit_id);
	  next HIT;
	}
      }

    }

    warning("Didnt match hit id to any regex [". $hit_id ."].");
  }

  # Return now if there are no intra-species hits.

  unless (defined @{$hits_by_species{$self->_regex_query_species}} &&
	 scalar @{$hits_by_species{$self->_regex_query_species}}) {
    return []
  }

  # Sort our hits by their distance to the query sequence.

  my %sorted_hits_by_species;

  foreach my $species (keys %hits_by_species) {

    my @sorted_hits 
      = sort {$good_hits->{$a}->{Ks} <=> $good_hits->{$b}->{Ks}} 
	@{$hits_by_species{$species}};
    $sorted_hits_by_species{$species} = \@sorted_hits;
  }

  # Accept all query species hits with a distance less than
  # the distance to the most related outgroup species.

  $self->outgroup_distance($self->_distance_cutoff);
  my $closest_outgroup_distance = $self->_distance_cutoff;

  foreach my $regex (@{$self->_regex_outgroup_species}){
    next
      unless $sorted_hits_by_species{$regex};

    my $closest_hit_id       = $sorted_hits_by_species{$regex}->[0];

    throw("Missing distance value")
      unless defined $good_hits->{$closest_hit_id}->{Ks} ||
	$good_hits->{$closest_hit_id}->{Ks} == 0;

    my $closest_hit_distance = $good_hits->{$closest_hit_id}->{Ks};

    if (($closest_hit_distance < $closest_outgroup_distance) && 
	($closest_hit_distance > 0)) {

      $closest_outgroup_distance = $closest_hit_distance;
    }
  }
  $self->outgroup_distance($closest_outgroup_distance);

  my @accepted_hits;

  my @query_species_hits = @{$sorted_hits_by_species{$self->_regex_query_species}};

  foreach my $hit_id (@query_species_hits) {
    if ($good_hits->{$hit_id}->{Ks} <= $closest_outgroup_distance &&
	$good_hits->{$hit_id}->{Ks} >= 0) {

      push (@accepted_hits, $good_hits->{$hit_id});
    } elsif ($good_hits->{$hit_id}->{Ks} >= 0) {
      last;
    }
  }

  return \@accepted_hits
}



### Utility Methods ###

=head2 _run_pairwise_paml

  Args       : An arrayref to a set of aligned Bio::Seq objects.
  Example    : none
  Description: Uses an array of aligned Bio::Seq objects to execute
               PAML in pairwise mode.
  Returntype : Bio::Tools::Phylo::PAML
  Exceptions : none
  Caller     : $self->run

=cut

sub _run_pairwise_paml {
  my ($self, $aligned_seqs) = @_;

  # Cope with a bloody PAML quirk.  Identical sequences
  # cause a fatal error.

  if ($aligned_seqs->[0]->seq eq $aligned_seqs->[1]->seq) {
    return (0, 0, 0, 0)
  }

  # The business of running codeml.

  my $retries = 0;
  my $result;
  my $paml; # Object must stay in scope until parsing is complete.

  while (! defined $result && $retries < 3) {

    $paml = Bio::EnsEMBL::Pipeline::GeneDuplication::PAML->new(
			     '-work_dir'     => $self->_work_dir,
			     '-executable'   => $self->_codeml,
			     '-aligned_seqs' => $aligned_seqs,
			     '-runmode'      => '-2',
			     '-seqtype'      => '1',
			     '-model'        => '0',
			     '-nssites'      => '0',
			     '-icode'        => ($self->_genetic_code) - 1
			    );

    my $parser = $paml->run_codeml;

    # Often a codeml run might look successful, but the Bioperl parser
    # will not be able to read the output created.  This is usually
    # due to a silent codeml failure and a lack of output.  Hence,
    # here an attempt is made to check the success of the parser.

    eval {
      $result = $parser->next_result();
    };

    if ($@ || ! defined $result) {
      warning("Problem deriving PAML result object.\n$@")
    }

    $retries++
  }

  unless (defined $result) {
    throw ("Result has become undefined.  This is probably a problem " .
	   "with Bio::Tools::Phylo::PAML.")
  }

  my $matrix;

  eval {
    if ($self->_distance_method eq 'NeiGojobori') {
      $matrix = $result->get_NGmatrix()
    } else {
      $matrix = $result->get_MLmatrix()
    }
  };

  if ($@){
    throw ("PAML failed to give a file that could be parsed.\n" .
	     "This is a PAML/codeml error.\n" .
	   "The offending aligned sequences are :\n>" . 
	   $aligned_seqs->[0]->display_id . "\n" . $aligned_seqs->[0]->seq . "\n>" . 
	   $aligned_seqs->[1]->display_id . "\n" . $aligned_seqs->[1]->seq . "\n$@");
  }

  $matrix->[0]->[1]->{N} = 0 unless defined $matrix->[0]->[1]->{N};
  $matrix->[0]->[1]->{S} = 0 unless defined $matrix->[0]->[1]->{S};

  return ($matrix->[0]->[1]->{dN},
	  $matrix->[0]->[1]->{dS},
	  $matrix->[0]->[1]->{N},
	  $matrix->[0]->[1]->{S})
}


=head2 _pairwise_align

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _pairwise_align {
  my ($self, $seqs) = @_;

  throw ("Pairwise alignment was only expecting two sequences.")
    unless ((scalar @$seqs) == 2);

  my $cba 
    = Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment->new(
	-genetic_code => 1);

  $cba->sequences($seqs);

  my $aligned_seqs = $cba->run_alignment;

  return $aligned_seqs
}

=head2 _calculate_coverage_and_identity

  Args[1]    :
  Example    :
  Description:
  Returntype : Array (min_coverage, identity)
  Exceptions :
  Caller     :

=cut


sub _calculate_coverage_and_identity {
  my ($self, $align) = @_;

  my $seq1 = $align->[0]->seq;
  my $seq2 = $align->[1]->seq;

  my $align_len = length($seq1);
  my $seq_length1 = 0;
  my $seq_length2 = 0;
  my $matches     = 0;
  my $identical   = 0;

  for (my $i = 0; $i < $align_len; $i++) {
    my $base1 = substr($seq1, $i, 1);
    my $base2 = substr($seq2, $i, 1);

    if ($base1 ne '-') {
      $seq_length1++;
    }

    if ($base2 ne '-') {
      $seq_length2++;
    }

    if ($base1 ne '-' && $base2 ne '-'){
      $matches++;
    }

    if ($base1 eq $base2 && $base1 ne '-'){
      $identical++;
    }
  }

  my $coverage1 = $matches/$seq_length1;
  my $coverage2 = $matches/$seq_length2;

  my ($min_coverage) = sort {$a <=> $b} ($coverage1, $coverage2);

  my $identity  = $identical/$matches;

  return ($min_coverage, $identity)
}


=head2 _appraise_hit_coverage

  Args[1]    :
  Example    :
  Description:
  Returntype : 0 or 1
  Exceptions :
  Caller     :

=cut

sub _appraise_hit_coverage{
  my ($self, $sbjct, $hsps) = @_;

  # Select our sequence length as being that of the 
  # longer sequence.

  my $sbjct_length 
    = $self->_fetch_seq($sbjct->name)->length;

  my ($longest_length) = sort {$a <=> $b} ($self->_query_seq->length, $sbjct_length);


  # Look at all the hits along the length of the 
  # query and tally the collective coverage of the hits.

  my @query_coverage;

  foreach my $hsp (@$hsps) {

    for (my $base_position = $hsp->query->start; 
	 $base_position <= $hsp->query->end;
	 $base_position++){
      $query_coverage[$base_position]++;
    }
  }

  my $covered_bases = 0;

  foreach my $base (@query_coverage){
    $covered_bases++ if $base;
  }

  # Return true if good coverage exists.

  return 1 if ($covered_bases >= $self->_coverage_cutoff * $longest_length);

  # Otherwise return false.
  return 0;
}


=head2 _fix_sbjct

  Args[1]    :
  Example    :
  Description: A work-around for a BPLite::Sbjct annoyance.  The 
               sbjct->name object returns the whole fasta description 
               line for a subject hit.  If the input fasta sequence 
               file includes more than an id on the description line, 
               this will be passed back every time the name method is 
               called.  This is a real pest is you are trying to 
               match ids via a regex or use the ids as hash keys.
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _fix_sbjct {
  my ($self, $sbjct) = @_;

  my $sbjct_name = $sbjct->name;

  $sbjct_name =~ s/\W*(\S+).*/$1/;

  # BAD!
  $sbjct->{NAME} = $sbjct_name;

  return $sbjct;
}


=head2 _hit_identity

  Args[1]    :
  Example    :
  Description:
  Returntype : 
  Exceptions :
  Caller     :

=cut

sub _hit_identity{
  my ($self, $hsps) = @_;

  my $tally_query_length = 0;
  my $tally_matched_bases = 0;

  foreach my $hsp (@$hsps) {
    $tally_query_length  += length($hsp->querySeq);
    $tally_matched_bases += $hsp->positive;
  }

  if ($tally_matched_bases && $tally_query_length){
    return $tally_matched_bases/$tally_query_length;
  }

  return 0
}


=head2 outgroup_distance

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub outgroup_distance {
  my $self = shift; 

  if (@_) {
    $self->{_outgroup_distance} = shift;
  }

  return $self->{_outgroup_distance};
}


### Sequence fetching ###

=head2 _fetch_seq

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _fetch_seq {
  my ($self, $id) = @_;

  throw ("Cant fetch sequence without an id.")
    unless $id;

  $self->{_cache} = {}
    unless $self->{_cache};

  if ($self->{_cache}->{$id}){
    return $self->{_cache}->{$id}
  }

  my $seq = $self->_seq_fetcher->fetch($id);

  throw ("Sequence fetch failed for id [$id].")
    unless ($seq && $seq->isa("Bio::Seq"));

  $self->{_cache}->{$seq->display_id} = $seq;

  return $seq
}


=head2 _fetch_seqs

  Args       : An arrayref of string sequence ids.
  Example    : none
  Description: An alias to $self->_fetch_seq, but handles 
               multiple sequences.
  Returntype : Arrayref of Bio::Seq
  Exceptions : none
  Caller     : $self->run

=cut

sub _fetch_seqs {
  my ($self, $seq_ids) = @_;

  my @seqs;

  foreach my $seq_id (@$seq_ids){
    push (@seqs, $self->_fetch_seq($seq_id));
  }

  return \@seqs;
}


=head2 _force_cache

  Args[1]    : Bio::Seq
  Example    : $self->_force_cache($seq);
  Description: Allows a sequence to be manually added to the seqfetcher 
               cache.  This is useful for coping with user supplied 
               sequences (eg. passed as a query sequence) that dont 
               exist in any database.
  Returntype : 1
  Exceptions : Warns if sequence already exists in cache.  Throws if
               a defined sequence isnt supplied.
  Caller     : 

=cut

sub _force_cache {
  my ($self, $seq) = @_;

  throw ("Trying to add something odd to the sequence cache [$seq].")
    unless (defined $seq);

  if ($self->{_cache}->{$seq->display_id}){
    warning('Sequence [' . $seq->display_id . 
	    '] already exists in cache, but will replace.');
  }

  $self->{_cache}->{$seq->display_id} = $seq;

  return 1
}


=head2 _seq_fetcher

  Args       : (optional) A seqfetcher of any variety, as long as 
               it has a 'fetch' method.
  Example    : none
  Description: Holds SeqFetcher object.
  Returntype : Bio::EnsEMBL::Pipeline::SeqFetcher::xxx
  Exceptions : none
  Caller     : $self->_candidate_hits, $self->_fetch_seqs

=cut

sub _seq_fetcher {
  my $self = shift;

  $self->{_seq_fetcher} = shift if @_;

  if (! $self->{_seq_fetcher}){
    $self->{_seq_fetcher} = 
      Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB->new(
				     -db => $self->_blastdb);
  }

  return $self->{_seq_fetcher};
}


### Blast-related Getter/Setters ###

=head2 _dbfile

  Args       : (optional) String
  Example    : none
  Description: Holds the filename of the blast database.
  Returntype : String.
  Exceptions : none
  Caller     : $self->new, $self->_hit_chooser.

=cut

sub _dbfile {
  my $self = shift;

  return $self->_blastdb->dbfile
}


=head2 _blast_obj

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _blast_obj {
  my $self = shift;

  if (@_){
    $self->{_blast_obj} = shift;
    return
  }

  unless ($self->{_blast_obj}){

    # Create a new blast object with our mixed species input database.

    $self->{_blast_obj} 
      = Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast->new(
		 -program         => $self->_blast_program,
		 -blastdb         => $self->_blastdb,
		 -queryseq        => $self->_query_seq,
		 -options         => '-hspmax 200',
		 -workdir         => $self->_work_dir,
		 -identity_cutoff => $self->_identity_cutoff);
  }


  return $self->{_blast_obj};
}


=head2 _blastdb

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _blastdb {
  my $self = shift;

  if (@_) {
    $self->{_blastdb} = shift;

    throw ("Blast database must be a Bio::EnsEMBL::Pipeline::Runnable::BlastDB.")
      unless ($self->{_blastdb}->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));

    throw ("Blast database has not been formatted.")
      unless $self->{_blastdb}->db_formatted;

    throw ("Blast database has been built without the " . 
	   "make_fetchable_index flag set (and this is " .
	   "a problem because the database can not be " . 
	   "used for sequence fetching).")
      unless $self->{_blastdb}->make_fetchable_index
  }

  throw ("Blast database object not set.")
    unless ($self->{_blastdb});

  return $self->{_blastdb};
}


=head2 _blast_program

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _blast_program {
  my $self = shift;

  if (@_){
    $self->{_blast_program} = shift;
    return
  }

  $self->{_blast_program} = $DEFAULT_BLAST_PROGRAM
    unless $self->{_blast_program};

  return $self->{_blast_program};
}


=head2 _blast_index_type

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _blast_index_type {
  my $self = shift;

  return $self->_blastdb->index_type;
}


### Getter/Setters ###

=head2 _query_seq

  Args       : (optional) Bio::Seq
  Example    : none
  Description: Holds the query sequence for which we are searching for duplicates.
  Returntype : Bio::Seq
  Exceptions : none
  Caller     : $self->run

=cut

sub _query_seq {
  my $self = shift;

  if (@_) {
    $self->{_query_seq} = shift;

    throw ("Query sequence is not a Bio::Seq object [" . 
	   $self->{_query_seq} . "]")
      unless $self->{_query_seq}->isa("Bio::Seq");

    throw ("Cant add query sequence to sequence cache manually.")
      unless $self->_force_cache($self->{_query_seq});

    return
  }

  throw ("Query sequence has not been set.")
    unless $self->{_query_seq};

  return $self->{_query_seq};
}


=head2 _work_dir

  Args       : (optional) String.
  Example    : none
  Description: Holds the path to the working directory.
  Returntype : String.
  Exceptions : none
  Caller     : $self->new, $self->_hit_chooser, $self->_run_pairwise_paml.

=cut

sub _work_dir {
  my $self = shift;

  if (@_) {
    $self->{_work_dir} = shift;
    return
  }

  throw ("Work directory not set.")
    unless $self->{_work_dir};

  return $self->{_work_dir};
}


=head2 _identity_cutoff

  Args       : (optional) an int or a float - a percentage value anyways.
  Example    : none
  Description: Holds identity cutoff percentage value.
  Returntype : A float value
  Exceptions : none
  Caller     : $self->new, $self->_hit_chooser.

=cut

sub _identity_cutoff {
  my $self = shift;

  if (@_) {
    $self->{_identity_cutoff} = shift;
    $self->{_identity_cutoff} /= 100 if $self->{_identity_cutoff} > 1;
    return
  }

  throw ("Blast match identity cutoff has not been set.")
    unless $self->{_identity_cutoff};

  return $self->{_identity_cutoff};
}


=head2 _coverage_cutoff

  Args       : (optional) an int or a float - a percentage value anyways.
  Example    : none
  Description: Holds coverage cutoff percentage value.
  Returntype : A float value.
  Exceptions : none
  Caller     : $self->new, $self->_hit_chooser.

=cut

sub _coverage_cutoff {
  my $self = shift; 

  if (@_) {
    $self->{_coverage_cutoff} = shift;
    $self->{_coverage_cutoff} /= 100 if $self->{_coverage_cutoff} > 1;
    return
  }

  throw ("The coverage cutoff has not been set.")
    unless $self->{_coverage_cutoff};

  return $self->{_coverage_cutoff};
}


=head2 _distance_cutoff

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _distance_cutoff {
  my $self = shift; 

  if (@_) {
    $self->{_distance_cutoff} = shift;
    return
  }

  $self->{_distance_cutoff} = $DEFAULT_DISTANCE_CUTOFF
    unless $self->{_distance_cutoff};

  return $self->{_distance_cutoff};
}


=head2 _regex_query_species

  Args       : String
  Example    : $self->_regex_query_species('ENSG')
  Description: Holds a regex string that will match the id of
               any sequence from the query species.
  Returntype : String or 1
  Exceptions : Throws when regex is not set.
  Caller     : $self->new, $self->_hit_chooser

=cut

sub _regex_query_species {
  my $self = shift;

  if (@_) {
    $self->{_regex_query_species} = shift;
    return
  }

  throw ("The regular expression used to match the sequence"
	       ." ids from the query species has not been set.")
    unless $self->{_regex_query_species};

  return $self->{_regex_query_species};
}


=head2 _regex_outgroup_species

  Args       : ref to an array of strings
  Example    : $self->_regex_outgroup_species(['ENSRNO', 'ENSMUS']);
  Description: Holds an array of regexs that will allow the sequence
               id of all non-query species sequences to be matched.
  Returntype : arrayref
  Exceptions : Warns if called while unset.
  Caller     : $self->new, $self->_hit_chooser

=cut

sub _regex_outgroup_species {
  my $self = shift; 

  if (@_) {
    $self->{_regex_outgroup_species} = shift;
    return
  }

  warning('No outgroup species regex provided.  ' .
	  'This may or may not be what you intend.')
    unless $self->{_regex_outgroup_species};

  return $self->{_regex_outgroup_species};
}


=head2 _genetic_code

  Args       : int
  Example    : $self->_genetic_code(1);
  Description: Holds an integer representing the genetic code.  To 
               choose the correct integer consult the documentation 
               used by the Bio::Seq->translate method.  1 is universal, 
               2 is vertebrate mitochondria.
  Returntype : int
  Exceptions : Warns if called while unset.
  Caller     : $self->new, $self->run, $self->_run_pairwise_paml

=cut

sub _genetic_code {
  my $self = shift; 

  if (@_) {
    $self->{_genetic_code} = shift;
    return
  }

  unless (defined $self->{_genetic_code}) {
    throw ('Genetic code unset.')
  }

  return $self->{_genetic_code};
}


=head2 _distance_method

  Args       : String
  Example    : 
  Description: 
  Returntype : 
  Exceptions : Throws if set to an unrecognised string.
  Caller     : 

=cut

sub _distance_method {
  my $self = shift; 

  if (@_) {
    $self->{_distance_method} = shift;

    unless ($self->{_distance_method} =~ /NeiGojobori/i |
	    $self->{_distance_method} =~ /ML/i){
      throw ("Distance method must be set to either " .
		   "NeiGojobori or ML, not [".
		   $self->{_distance_method}."]");
    }
  }

  $self->{_distance_method} = $DEFAULT_DISTANCE_METHOD
    unless $self->{_distance_method};

  return $self->{_distance_method};
}


=head2 _codeml

  Args       : [optional] String
  Example    : $self->_codeml('/path/to/codeml')
  Description: Holds the path to the codeml executable
  Returntype : String
  Exceptions : Throws if a full path is included, but the 
               file is not executable.
  Caller     : $self->new, $self->_run_pairwise_paml

=cut

sub _codeml {
  my $self = shift; 

  if (@_) {
    $self->{_codeml} = shift;
    return
  }

  $self->{_codeml} = 'codeml'
    unless $self->{_codeml};

  # If it looks like our executable comes with a full 
  # path, check that it will work.

  throw ("codeml executable not found or not " .
	       "executable. Trying to execute path " .
	       "[". $self->{_codeml} ."]")
    if ($self->{_codeml} =~ /^\//
	& !-x $self->{_codeml});

  return $self->{_codeml};
}

return 1

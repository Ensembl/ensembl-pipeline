package Bio::EnsEMBL::Pipeline::GeneDuplication;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB;
use Bio::EnsEMBL::Pipeline::GeneDuplication::Chooser;
use Bio::EnsEMBL::Pipeline::GeneDuplication::PAML;
use Bio::EnsEMBL::Pipeline::GeneDuplication::Result;
use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;
use Bio::Tools::Run::Alignment::Clustalw;

my $DEFAULT_DISTANCE_METHOD = 'NeiGojobori';

@ISA = qw(Bio::EnsEMBL::Root);

=head2 new

  Args       : -blastdb      - A Bio::EnsEMBL::Pipeline::Runnable::BlastDB 
                               object for which the run method has been invoked.
               -hit_identity - (optional) Hit identity.  Defaults to 0.80
               -hit_coverage - (optional) Hit coverage.  Defults to 0.80
               -work_dir     - (optional) Dir where working files are 
                                          placed.  Defaults to /tmp.
  Example    : none
  Description: Constructs new GeneDuplication object
  Returntype : Bio::EnsEMBL::Pipeline::GeneDuplication
  Exceptions : Throws if database file not specified or does not exist.
  Caller     : General

=cut

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($blastdb,
      $hit_identity,
      $hit_coverage,
      $regex_query_species,
      $regex_outgroup_species,
      $genetic_code,
      $distance_method,
      $codeml,
      $work_dir,)
    = $self->_rearrange([qw(BLASTDB
			    HIT_IDENTITY
			    HIT_COVERAGE
			    REGEX_QUERY_SPECIES
			    REGEX_OUTGROUP_SPECIES
			    GENETIC_CODE
			    DISTANCE_METHOD
			    CODEML_EXECUTABLE
			    WORK_DIR)],@args);

  $self->throw("Blast database must be a Bio::EnsEMBL::Pipeline::Runnable::BlastDB.") 
    unless ($blastdb && $blastdb->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));

  $self->throw("Blast database does not appear to be formatted.") 
    unless ($blastdb->db_formatted);

  $self->_blastdb($blastdb);

  $self->_working_dir($work_dir) if $work_dir;

  $self->_identity_cutoff($hit_identity) if $hit_identity;
  $self->_coverage_cutoff($hit_coverage) if $hit_coverage;

  $self->_identity_cutoff(80) unless $self->_identity_cutoff;
  $self->_coverage_cutoff(80) unless $self->_coverage_cutoff;

  $self->_regex_query_species($regex_query_species) 
    if $regex_query_species;
  $self->_regex_outgroup_species($regex_outgroup_species) 
    if $regex_outgroup_species;

  $self->_genetic_code($genetic_code) if $genetic_code;

  $self->_distance_method($distance_method)
    if $distance_method;

  $self->_codeml($codeml) if $codeml;

  return $self
}

=head2 run

  Args[1]    : Input sequence (Bio::Seq)
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

  die "Cant run without a Bio::Seq input sequence" 
    unless ($input_seq && $input_seq->isa("Bio::Seq"));

  $self->_query_sequence($input_seq);

  # Derive a list of sequences that look like duplications of 
  # the query sequence.

  my $seq_ids = $self->_candidate_hits($input_seq);

  # Add the id of our query seq, or it will be omitted.

  push (@$seq_ids, $input_seq->display_id);

  unless (scalar @$seq_ids > 1) {
    print "No homologous matches were found that satisfied the match criteria.\n";
    return 0
  }

  # Perform a codon based alignment of our sequences.

  my $seqs = $self->_fetch_seqs($seq_ids);

  my $cba = 
    Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment->new(
       -genetic_code => $self->_genetic_code);

  $cba->sequences($seqs);

  my $aligned_seqs = $cba->run_alignment;

  # Run PAML with these aligned sequences.

  my $parser = $self->_run_pairwise_paml($aligned_seqs);

  my @results;

  eval {
    # This is a bit stupid, but we dont know until here
    # whether our run has been successful.
    push (@results, $parser->next_result());
  };

  if ($@){
    $self->warn("PAML run was unsuccessful.\n$@");
    return 1;
  }

  while (my $result = $parser->next_result()) {
    push (@results, $result); # More stupidity.
  }

  unless (@results) {
    print "Duplications not found for this gene.\n";
    return 0
  }


  $self->throw("There are more than two sets of results returned from\n" .
	       "the PAML parser.  This was not expected.") 
    if scalar @results > 1;

  return $self->_extract_results($results[0])  
}


=head2 _extract_results

  Args[1]    : Bio::Tools::Phylo::PAML::Result
  Example    : none
  Description: Derive the PAML result object (and output matrix) into 
               the much simpler GeneDuplication::Result object.
  Returntype : Bio::EnsEMBL::Pipeline::GeneDuplication::Result
  Exceptions : none
  Caller     : $self->run

=cut


sub _extract_results {
  my ($self, $result) = @_;

  my $query_id = $self->_query_sequence->display_id;

  my $matrix;

  $matrix = $result->get_MLmatrix() 
    if $self->_distance_method =~ /ML/;
  $matrix = $result->get_NGmatrix() 
    if $self->_distance_method =~ /NeiGojobori/;

  $self->throw("Failed to retrieve a result matrix from ".
	       "the PAML result.")
    unless $matrix;

  my @otus = $result->get_seqs();

  my $result_obj = Bio::EnsEMBL::Pipeline::GeneDuplication::Result->new(
		       -id              => $query_id,
		       -distance_method => $self->_distance_method);


  for(my $i = 0; $i < scalar @otus; $i++){

    my $match_id;

    $match_id = $otus[$i]->display_id
      if ($otus[$i]->display_id ne $query_id);

    for (my $j = $i+1; $j < scalar @otus; $j++){

      next unless (($otus[$i]->display_id eq $query_id)||
		   ($otus[$j]->display_id eq $query_id));

      $match_id = $otus[$j]->display_id unless $match_id;

      $result_obj->add_match($match_id,
			    $matrix->[$i]->[$j]->{'dN'},
			    $matrix->[$i]->[$j]->{'dS'});
    }
  }

  return $result_obj
}



=head2 _candidate_hits

  Args[1]    : A query sequence (Bio::Seq)
  Example    : none
  Description: Method handles the use of the Chooser module
               and returns a list of likely duplicated sequences. 
  Returntype : reference to an array of id strings.
  Exceptions : none
  Caller     : $self->run

=cut

sub _candidate_hits {
  my ($self, $query) = @_;

  my $hit_chooser = $self->_hit_chooser;

  my $seq_ids = $hit_chooser->find_recent_duplications($query);

  return $seq_ids; #returns an arrayref of sequence ids from query species
}

=head2 _hit_chooser

  Args       : none
  Example    : none
  Description: A getter/setter that actually creates a GeneDuplication::Chooser
               object if one does not already exist. 
  Returntype : Bio::EnsEMBL::Pipeline::GeneDuplication::Chooser
  Exceptions : none
  Caller     : $self->_candidate_hits

=cut

sub _hit_chooser {
  my $self = shift;

  unless ($self->{_hit_chooser}) {

    $self->{_hit_chooser} = 
      Bio::EnsEMBL::Pipeline::GeneDuplication::Chooser->new(
		     -blastdb                => $self->_blastdb,
		     -work_dir               => $self->_working_dir,
		     -regex_query_species    => $self->_regex_query_species,
		     -regex_outgroup_species => $self->_regex_outgroup_species,
		     -identity_cutoff        => $self->_identity_cutoff,
		     -coverage_cutoff        => $self->_coverage_cutoff)
  }

  return $self->{_hit_chooser}
}


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

  my $paml = Bio::EnsEMBL::Pipeline::GeneDuplication::PAML->new(
			     '-work_dir'     => $self->_working_dir,
			     '-executable'   => $self->_codeml,
			     '-aligned_seqs' => $aligned_seqs,
			     '-runmode'      => '-2',
			     '-seqtype'      => '1',
			     '-model'        => '0',
			     '-nssites'      => '0',
			     '-icode'        => ($self->_genetic_code) - 1
			    );

  my $parser = $paml->run_codeml;

  return $parser;
}


### Sequence fetching ###

=head2 _fetch_seqs

  Args       : An arrayref of string sequence ids.
  Example    : none
  Description: Handles retrieval of sequences using the seqfetcher.
  Returntype : Arrayref of Bio::Seq
  Exceptions : none
  Caller     : $self->run

=cut


sub _fetch_seqs {
  my ($self, $seq_ids) = @_;

  my @seqs;

  foreach my $seq_id (@$seq_ids){
    # This is a bit of a hack here - we get the sequences from the
    # Chooser objects sequence fetcher.  This gets around problems
    # with the query sequence not being in the blast database.
    # Really, this module and Chooser.pm need to inherit from
    # a base class that would have the functionality for caching
    # sequences from the seq fetch.  Right now only the Chooser.pm
    # module has this and I'm loathe to add it here a second time.
#    push (@seqs, $self->_seq_fetcher->fetch($seq_id));
    push (@seqs, $self->_hit_chooser->_fetch_seq($seq_id));
  }

  return \@seqs;
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


### Getter/Setters ###

=head2 _query_sequence

  Args       : (optional) Bio::Seq
  Example    : none
  Description: Holds the query sequence for which we are searching for duplicates.
  Returntype : Bio::Seq
  Exceptions : none
  Caller     : $self->run, $self->_extract_results

=cut

sub _query_sequence {
  my $self = shift;

  if (@_) {
    $self->{_query_seq} = shift;
  }

  $self->throw("No query sequence attached.")
    unless $self->{_query_seq};

  return $self->{_query_seq};
}

=head2 _blastdb

  Args       : (optional) String
  Example    : none
  Description: Holds the filename of the blast database.
  Returntype : String.
  Exceptions : none
  Caller     : $self->new, $self->_hit_chooser.

=cut

sub _blastdb {
  my $self = shift;

  $self->{_blastdb} = shift if @_;

  return $self->{_blastdb}
}

=head2 _working_dir

  Args       : (optional) String.
  Example    : none
  Description: Holds the path to the working directory.
  Returntype : String.
  Exceptions : none
  Caller     : $self->new, $self->_hit_chooser, $self->_run_pairwise_paml.

=cut

sub _working_dir {
  my $self = shift;

  if (@_) {
    $self->{_working_dir} = shift;
  } 

  return $self->{_working_dir};
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
  }

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
  }

  return $self->{_coverage_cutoff};
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
  }

  $self->throw('No query species regex provided.')
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
  }

  $self->warn('No outgroup species regex provided.  ' .
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
  }

  $self->throw('Genetic code unset.')
    unless $self->{_genetic_code};

  return $self->{_genetic_code};
}


=head2 _distance_method

  Args       : String
  Example    : $self->_genetic_code(1);
  Description: Holds an integer representing the genetic code.  To 
               choose the correct integer consult the documentation 
               used by the Bio::Seq->translate method.  1 is universal, 
               2 is vertebrate mitochondria.
  Returntype : int
  Exceptions : Warns if called while unset.
  Caller     : $self->new, $self->run, $self->_run_pairwise_paml

=cut

sub _distance_method {
  my $self = shift; 

  if (@_) {
    $self->{_distance_method} = shift;

    unless ($self->{_distance_method} =~ /NeiGojobori/i |
	    $self->{_distance_method} =~ /ML/i){
      $self->throw("Distance method must be set to either " .
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
  }

  $self->{_codeml} = 'codeml'
    unless $self->{_codeml};

  # If it looks like our executable comes with a full 
  # path, check that it will work.

  $self->throw("codeml executable not found or not " .
	       "executable. Trying to execute path " .
	       "[". $self->{_codeml} ."]")
    if ($self->{_codeml} =~ /^\//
	& !-x $self->{_codeml});

  return $self->{_codeml};
}




return 1

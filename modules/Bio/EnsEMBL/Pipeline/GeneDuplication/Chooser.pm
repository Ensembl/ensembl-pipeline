package Bio::EnsEMBL::Pipeline::GeneDuplication::Chooser;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Pipeline::SeqFetcher::FetchFromBlastDB;
use Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment;
use Bio::EnsEMBL::Pipeline::GeneDuplication::PAML;
use Bio::EnsEMBL::Pipeline::Runnable::BlastDB;
use Bio::EnsEMBL::Pipeline::Runnable::MinimalBlast;


my $DEFAULT_BLAST_FLAV      = 'wublastn';
my $DEFAULT_DISTANCE_CUTOFF = 1;

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;
  
  my ($query_seq,
      $blastdb,
      $work_dir,
      $regex_query_species,
      $regex_outgroup_species,
      $identity_cutoff,
      $coverage_cutoff,
      $distance_cutoff,
      $genetic_code,
      $codeml) = $self->_rearrange([qw(QUERY_SEQ
						BLASTDB
						WORK_DIR
						REGEX_QUERY_SPECIES
						REGEX_OUTGROUP_SPECIES
						IDENTITY_CUTOFF
						COVERAGE_CUTOFF
						DISTANCE_CUTOFF
						GENETIC_CODE
						CODEML_EXECUTABLE)],@args);

  $self->throw("Blast database must be a Bio::EnsEMBL::Pipeline::Runnable::BlastDB.") 
    unless ($blastdb && $blastdb->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));

  $self->_blastdb($blastdb);

  $self->_work_dir($work_dir)                             if $work_dir;
  $self->_query_seq($query_seq)                           if $query_seq;
  $self->_regex_query_species($regex_query_species)       if $regex_query_species;
  $self->_regex_outgroup_species($regex_outgroup_species) if $regex_outgroup_species;
  $self->_identity_cutoff($identity_cutoff)               if $identity_cutoff;
  $self->_coverage_cutoff($coverage_cutoff)               if $coverage_cutoff;
  $self->_genetic_code($genetic_code)                     if $genetic_code;
  $self->_codeml($codeml)                                 if $codeml;

  return $self
}


=head2 DESTROY

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub DESTROY {
  my $self = shift;

  print "Cleaning up.\n";

  $self->_seq_fetcher->db->remove_index_files;  
}


=head2 find_recent_duplications

  Args[1]    : Bio::Seq query sequence
  Example    : none
  Description:
  Returntype :
  Exceptions :
  Caller     : General

=cut

sub find_recent_duplications {
  my ($self, $query) = @_;

  # We are looking for duplicates of the query gene
  # in our blast database.
  $self->_query_seq($query) if $query;

  my $bplite_report = $self->_blast_obj->run;

  $self->throw("Blast process did not return a report.")
    unless ($bplite_report->isa("Bio::Tools::BPlite"));

  return $self->_process_for_same_species_duplicates($bplite_report);
}


=head2 _process_for_same_species_duplicates 

  Args[1]    :
  Example    :
  Description: This is the main algorithmic implementation.  See the docs above for an explanation of what is going on here.
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _process_for_same_species_duplicates {
  my ($self, $bplite_report) = @_;

  # Process our blast report.  For each blast hit to our query sequence:
  #   * throw away self matches (if any)
  #   * filter by coverage
  #   * calculate genetic distance between each subject and the query
  #   * add subject sequence to correct species hash

  my %species_hash;
  my %hit_distance;
  my $have_an_outgroup = 0;
  my $report_empty = 1;


 PARTITION_HITS:
  while (my $sbjct = $bplite_report->nextSbjct){

    # Mangle the BPLite::Sbjct object for its own good.
    $sbjct = $self->_fix_sbjct($sbjct);

    $report_empty = 0;

    # Skip hit if it is a match to self.

    my $sbjct_id = $sbjct->name;
    $sbjct_id =~ s/\W*(\w+).*/$1/;

    next PARTITION_HITS 
      if ($self->_query_seq->display_id eq $sbjct_id);

    # First, filter by coverage

    next PARTITION_HITS 
      unless ($self->_appraise_hit_coverage($sbjct));

    # Second, filter by genetic distance.

    $hit_distance{$sbjct_id} 
      = $self->_calculate_pairwise_distance(
            $self->_query_seq->display_id, 
	    $sbjct_id,
	    'synonymous');

    # Third, partition hits according to their species.  The species
    # from which the subject is derived is determined by a
    # regular expression match to the sequence id.

    my $query_regex = $self->_regex_query_species;

    if ($sbjct_id =~ /$query_regex/) {

      push(@{$species_hash{$self->_regex_query_species}}, $sbjct);

      next PARTITION_HITS;
    } else {

      foreach my $regex (@{$self->{_regex_outgroup_species}}) {
	if ($sbjct_id =~ /$regex/){
	  $have_an_outgroup = 1;
	  push (@{$species_hash{$regex}}, $sbjct);
	  next PARTITION_HITS;
	}
      }

    }

    $self->throw("Didnt match hit id to any regex [$sbjct_id].");
  }

  $self->throw("Did not find any hits to query sequence.")
    if $report_empty;

  # Sort our hits by their distance to the query sequence.

  my %sorted_species_hits;

  foreach my $species (keys %species_hash) {

    my @sorted_hits 
      = sort {$hit_distance{$a->name} <=> $hit_distance{$b->name}} 
	@{$species_hash{$species}};

    $sorted_species_hits{$species} = \@sorted_hits;
  }

  # Accept all query species hits with a distance less than
  # the distance to the most related outgroup species.

  my $closest_outgroup_distance = $self->_distance_cutoff;

  foreach my $regex (@{$self->_regex_outgroup_species}){
    if (defined $sorted_species_hits{$regex}->[0]->name && 
	defined $hit_distance{$sorted_species_hits{$regex}->[0]->name} &&
	$hit_distance{$sorted_species_hits{$regex}->[0]->name} < $closest_outgroup_distance){

      $closest_outgroup_distance = $hit_distance{$sorted_species_hits{$regex}->[0]->name};
    }
  }

  my @accepted_ids;

  foreach my $sbjct (@{$sorted_species_hits{$self->_regex_query_species}}) {

    if ($hit_distance{$sbjct->name} <= $closest_outgroup_distance) {
      push (@accepted_ids, $sbjct->name);

    } else {
      last;
    }
  }

  return \@accepted_ids;
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

  $sbjct_name =~ s/\W*(\w+).*/$1/;

  # BAD!
  $sbjct->{NAME} = $sbjct_name;

  return $sbjct;
}


=head2 _calculate_pairwise_distance

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _calculate_pairwise_distance {
  my ($self, $input_id_1, $input_id_2, $distance_measure) = @_;

  $distance_measure = 'synonymous' unless $distance_measure;

  my @seqs = ($self->_fetch_seq($input_id_1), 
	      $self->_fetch_seq($input_id_2));

  return 0 if ($seqs[0]->display_id eq $seqs[1]->display_id);

  $self->throw("Didnt correctly obtain two sequences for alignment.")
    unless scalar @seqs == 2;

  my $align = $self->_pairwise_align(\@seqs);

  my $paml_parser;

  eval {
    $paml_parser = $self->_run_pairwise_paml($align);
  };

  if ($@){
    $self->throw("Pairwise use of PAML FAILED!\n$@");
    my $fh = $self->{_filehandle};
    print $fh "Error encountered while analysing $input_id_1 versus $input_id_2\n";
    return 0
  }

  my $result;
  my $NGmatrix;

  eval {
    $result = $paml_parser->next_result();
    $NGmatrix = $result->get_NGmatrix();
  };

  if ($@){
    $self->warn("PAML failed to give a file that could be parsed.  No doubt PAML threw an error!\n$@");

    my $fh = $self->{_filehandle};
    print $fh "Couldnt derive matrix from PAML run for  $input_id_1 versus $input_id_2\n";
    return 0
  }
  
#print "Synonymous distance is : " . $NGmatrix->[0]->[1]->{dS} . "\n";

  return $NGmatrix->[0]->[1]->{dN} if $distance_measure eq 'nonsynonymous';
  return $NGmatrix->[0]->[1]->{dS} if $distance_measure eq 'synonymous';
  return 0;
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

  $self->throw("Pairwise alignment was only expecting two sequences.")
    unless ((scalar @$seqs) == 2);

  my $cba 
    = Bio::EnsEMBL::Pipeline::GeneDuplication::CodonBasedAlignment->new(
	-genetic_code => 1);

  $cba->sequences($seqs);

  my $aligned_seqs = $cba->run_alignment;

  return $aligned_seqs
}


=head2 _run_pairwise_paml

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _run_pairwise_paml {
  my ($self, $aligned_seqs) = @_;

  my $paml = Bio::EnsEMBL::Pipeline::GeneDuplication::PAML->new(
                             '-work_dir'   => $self->_work_dir,
			     '-executable' => $self->_codeml,
			     '-aligned_seqs' => $aligned_seqs,
			     '-runmode'    => '-2', 
			     '-seqtype'    => '1',
			     '-model'      => '0',
			     '-nssites'    => '0',
			     '-icode'      => ($self->_genetic_code - 1));

  my $parser = $paml->run_codeml;

  return $parser;
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
  my ($self, $sbjct) = @_;

  # First, throw out hits that are way longer than
  # the query.

  my $sbjct_length 
    = $self->_fetch_seq($sbjct->name)->length;

  return 0 
    if ($sbjct_length > 
	((2 - $self->_coverage_cutoff) * $self->_query_seq->length));

  # If still here, look at all the hits along the length of the 
  # query and tally the collective coverage of the hits.

  my @query_coverage;

  while (my $hsp = $sbjct->nextHSP) {

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

  return 1 if ($covered_bases >= $self->_coverage_cutoff * $self->_query_seq->length);

  # Otherwise return false.
  return 0;
}


=head2 _blast_obj

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

### Blast object ###

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
		 -program         => $self->_blast_flav,
		 -blastdb         => $self->_blastdb,
		 -queryseq        => $self->_query_seq,
		 -options         => '',
		 -workdir         => $self->_work_dir,
		 -identity_cutoff => $self->_identity_cutoff);
  }


  return $self->{_blast_obj};
}




=head2 _blast_flav

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _blast_flav {
  my $self = shift;

  if (@_){
    $self->{_blast_flav} = shift;
    return
  }

  $self->{_blast_flav} = $DEFAULT_BLAST_FLAV
    unless $self->{_blast_flav};

  return $self->{_blast_flav};
}


### Sequence Fetching ###

=head2 _force_cache

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _force_cache {
  my ($self, $seq) = @_;

  $self->throw("Trying to add something odd to the sequence cache [$seq].")
    unless (defined $seq);

  if ($self->{_cache}->{$seq->display_id}){
    $self->warn('Sequence [' . $seq->display_id . 
		'] already exists in cache, but will replace.');
  }

  $self->{_cache}->{$seq->display_id} = $seq;

  return 1
}


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

  $self->throw("Cant fetch sequence without an id.")
    unless $id;

  $self->{_cache} = {}
    unless $self->{_cache};

  if ($self->{_cache}->{$id}){
    return $self->{_cache}->{$id}
  }

  my $seq = $self->_seq_fetcher->fetch($id);

  $self->{_cache}->{$seq->display_id} = $seq;

  return $seq
}


=head2 _seq_fetcher

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

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


### Getters/Setters ###

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
    $self->throw("Blast database must be a Bio::EnsEMBL::Pipeline::Runnable::BlastDB.")
      unless ($self->{_blastdb}->isa("Bio::EnsEMBL::Pipeline::Runnable::BlastDB"));
    return
  }

  $self->throw("Unable to return unset blastdb object.")
    unless $self->{_blastdb};

  return $self->{_blastdb};
}


=head2 _work_dir

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _work_dir {
  my $self = shift;

  if (@_) {
    $self->{_work_dir} = shift;
    return
  }

  $self->throw("Work directory not set.")
    unless $self->{_work_dir};

  return $self->{_work_dir};
}


=head2 _query_seq

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _query_seq {
  my $self = shift;

  if (@_) {
    $self->{_query_seq} = shift;

    $self->throw("Cant add query sequence to sequence cache manually.")
      unless $self->_force_cache($self->{_query_seq});

    return
  }

  $self->throw("Query sequence has not been set.")
    unless $self->{_query_seq};

  return $self->{_query_seq};
}


=head2 _regex_query_species

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _regex_query_species {
  my $self = shift;

  if (@_) {
    $self->{_regex_query_species} = shift;
    return
  }

  $self->throw("The regular expression used to match the sequence"
	       ." ids from the query species has not been set.")
    unless $self->{_regex_query_species};

  return $self->{_regex_query_species};
}


=head2 _regex_outgroup_species

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _regex_outgroup_species {
  my $self = shift;

  if (@_) {
    $self->{_regex_outgroup_species} = shift;
    return
  }

  $self->throw("The regular expression used to match any sequence"
	       ." ids from outgroup species has not been set.")
    unless $self->{_regex_outgroup_species};

  return $self->{_regex_outgroup_species};
}


=head2 _identity_cutoff

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _identity_cutoff {
  my $self = shift;

  if (@_) {
    $self->{_identity_cutoff} = shift;
    $self->{_identity_cutoff} /= 100 if $self->{_identity_cutoff} > 1;
    return
  }

  $self->throw("Blast match identity cutoff has not been set.")
    unless $self->{_identity_cutoff};

  return $self->{_identity_cutoff};
}


=head2 _coverage_cutoff

  Args[1]    :
  Example    :
  Description:
  Returntype :
  Exceptions :
  Caller     :

=cut

sub _coverage_cutoff {
  my $self = shift; 

  if (@_) {
    $self->{_coverage_cutoff} = shift;
    $self->{_coverage_cutoff} /= 100 if $self->{_coverage_cutoff} > 1;
    return
  }

  $self->throw("Blast match coverage cutoff has not been set.")
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


=head2 _genetic_code

  Args       : int
  Example    : $self->_genetic_code(1);
  Description: Holds an integer representing the genetic code.  To 
               choose the correct integer consult the documentation 
               used by the Bio::Seq->translate method.  1 is universal, 
               2 is vertebrate mitochondria.
  Returntype : int
  Exceptions : Warns if called while unset.
  Caller     : $self->new, $self->_run_pairwise_paml

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

return 1;

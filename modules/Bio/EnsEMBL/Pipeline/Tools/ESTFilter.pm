
# Cared for by Dan Andrews <dta@sanger.ac.uk>
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME
  
Bio::EnsEMBL::Pipeline::Tools::ESTFilter
 
=head1 SYNOPSIS

Instantiate the object:

my $est_filter = = Bio::EnsEMBL::Pipeline::ESTFilter->new();

While it is perfectly fine to run this module with the default 
settings, a few options do exist to tailor the sensitivity of 
the EST trimming.  These are designed to be set when the object
is instantiated:

my $est_filter = = Bio::EnsEMBL::Pipeline::ESTFilter->new(
                          '-min_length'                => 100,
			  '-three_prime_twilight_zone' =>  80,
			  '-bleb_length'               =>   8,
			  '-allowed_perc_n'            =>  10);

The options set the following values:

min_length                - the minimum EST length in base pairs

three_prime_twilight_zone - the length of 3-prime bases that are 
initially checked for poor quality sequence.  80 is the default 
setting and seems to work well.

bleb_length               - using hyper-technical jargon, a bleb 
is a run of bases all of the same kind, as is commonly seen at 
the end of a sequencing chromatogram.  This module aims to trim 
these from the sequence if they exist.  The default value is 8 
and it is recommended that values between 6 and 9 be used.

allowed_perc_n            - maximum allowable percentage of the 
EST sequence that can be comprised of Ns.


So, to actually filter a set of EST sequence from a file the 
following rationale should serve you well:

use Bio:SeqIO;

my $seqio = Bio::SeqIO->new(-format => 'fasta', 
			    -file   => 'unfiltered_ests.fa');

while (my $est = $seqio->next_seq){

  if ($est_filter->appraise($est)){ 
  # The appraise method returns 1 or 0, depending on whether
  # the EST is accepted or rejected by the filter.

  my $filtered_seq = $est_filter->filtered_sequence;
  # As the EST sequence could have been trimmed by the filter
  # before it was deemed acceptable, the sequence should be
  # derived from the ESTFilter object.

  print STDOUT ">" . $filtered_seq->display_id . "\n" . 
     $filtered_seq->seq . "\n";

  }
}

It is possible to retrieve information about how may ESTs have been
rejected or accepted, and the reasons for this:

print STDERR $est_filter->status . "\n";
# Returns the number of accepted/rejected ESTs that have 
# been processed since the filter was instantiated.

while (my ($error, $count) = each %{$est_filter->error_types}) {
  print STDERR $error . " - " . $count . "\n";
}
# Returns more detailed information as to why ESTs were rejected
# and of how many sequences were trimmed.


=head1 DESCRIPTION

This module conducts very simple filtering and trimming of input EST 
sequences.  In particular, this module has been implemented to trim
low quality sequence from the 3-prime of EST sequence reads.

=head1 CONTACT
  
Post general queries to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Pipeline::Tools::ESTFilter;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;
use Bio::Seq;

# Object preamble - inherits from Bio::Root::Object;

@ISA = qw(Bio::EnsEMBL::Root);

##### 'Public' methods #####

sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ($length_cutoff,
      $three_prime_twilight_zone,
      $bleb_length,
      $n_percentage) = $self->_rearrange([qw(MIN_LENGTH
					     THREE_PRIME_TWILIGHT_ZONE
					     BLEB_LENGTH
					     ALLOWED_PERC_N)],@args);
  

  # Set minimum est length cut-off
  $length_cutoff = 100 unless $length_cutoff;
  $self->_min_length($length_cutoff);

  # Set length of 3-prime region that is a candidate for low
  # quality sequence.  Experience and a bit of checking shows
  # that this tends to be around 80 bases (at least, the greatest
  # density of 'blebs' are found in the last 80 bases of human
  # EST datasets).
  # Set this to 80 unless you have better information.
  $three_prime_twilight_zone = 80 unless $three_prime_twilight_zone;
  $self->_three_prime_twilight_zone($three_prime_twilight_zone);

  # Set the maximum allowable length of a run of a single bases 
  # (i.e. a bleb).  A value between 6 and 9 is most sensitive.  
  # Dont use a value outside this range unless you have a 
  # reason to do so.
  $bleb_length = 8 unless $bleb_length;
  $self->_bleb_length($bleb_length);

  # Set the maximum allowable percentage of Ns in the sequence
  $n_percentage = 5 unless $n_percentage;
  $self->_max_perc_Ns($n_percentage);

  return $self;
}


sub appraise {
  my ($self, $est_sequence) = @_;
  
  $self->_est_sequence($est_sequence);

  $self->_trimmed_status(0);

  # First, bail out if the untouched sequence is too short.
  if (! $self->_check_length) {
    $self->status('reject');

    return 0
  }

  # Check for sequence oddities
  if (($self->_check_Ns)&&
      ($self->_check_single_base_blebs)) {
    $self->status('pass');

    return $self->_est_sequence;
    
  } else {
    # If oddities exist, see if they can be
    # rectified by trimming the sequence.
    $self->_apply_bleb_remover;

    if (($self->_check_length)&&
	($self->_check_Ns)&&
	($self->_check_single_base_blebs)) {
      $self->status('pass');

      return $self->_est_sequence;
    } else {
      # Odd sequence, chuck it.
      $self->status('reject');

      return 0;
    }    
  }
}

sub filtered_seq {
  my $self = shift;

  return $self->_est_sequence;
}


sub error_messages {
  my ($self, $string, $error_shortname) = @_;

  if ($string eq 'flush'){
    $self->{'_error_string'} = '';
  }

  if (($string)&&($string ne 'flush')) {
    $self->{'_error_string'} .= $string;

    $self->{'_error_tally'}->{$error_shortname}++;
  }

  return $self->{'_error_string'}
}


sub status {
  my ($self, $result) = @_;

  if ($result eq 'pass'){
    $self->{'pass'}++;
  }

  if ($result eq 'reject'){
    $self->{'reject'}++;
  }
  
  if (!$result){
    my $status_string = "Passed : " . $self->{'pass'} . 
      "  Rejected : " . $self->{'reject'};
    return $status_string
  }
}


sub error_types {
  my $self = shift;

  return $self->{'_error_tally'};
}


##### Internal methods #####

sub _check_length {
  my ($self) = @_;

  if ($self->_est_sequence->length < $self->_min_length){
    $self->_trimmed_status(1);  # Dont care if it has been or not.
    $self->error_messages("Sequence shorter than length cutoff", 
			  "Short sequence (rejected)");

    return 0
  }

  return 1
}


sub _check_Ns {
  my ($self) = @_;

  my $seq = $self->_est_sequence->seq;

  my $n_count = $seq =~ tr/[Nn]//;

  if ($n_count >= $self->_max_perc_Ns * $self->_est_sequence->length) {

    $self->error_messages("More than 10 percent of the sequence are Ns", 
			  "Too many Ns (rejected)")
      if $self->_trimmed_status;

    return 0;
  }

  return 1;
}


sub _check_single_base_blebs {
  my ($self) = @_;

  my $bleb_length = $self->_bleb_length;
  my $seq = $self->_est_sequence->subseq(($self->_est_sequence->length - $self->_three_prime_twilight_zone + 1), 
					 $self->_est_sequence->length);

  if (($seq =~ /A{$bleb_length}/i)||
      ($seq =~ /T{$bleb_length}/i)||
      ($seq =~ /C{$bleb_length}/i)||
      ($seq =~ /G{$bleb_length}/i)) {

    return 0 if !$self->_trimmed_status;
    # If we still have blebs after trimming there is not
    # much more that can be done.  It is possible that the
    # blebs represent real sequence, hence we cant really
    # throw the EST away.  Hence, if the EST has been trimmed
    # we return 1 and pass it.

  }

  return 1
}


sub _apply_bleb_remover {
  my ($self) = @_;

  my $window_length = 20; # bp
  my $window_overlap = 10; # bp
  my $bleb_score_cutoff = 16;

  my $est_untrimmed_length = $self->_est_sequence->length;
  my $cutoff_coordinate = $est_untrimmed_length; # Setting up the first window
  my $bleb_score = 0;

  while ($bleb_score < $bleb_score_cutoff){
    last if ($cutoff_coordinate - $window_length <= 1); 

    $bleb_score = $self->_calculate_blebbiness($self->_est_sequence->subseq(
			 $cutoff_coordinate - $window_length, 
			 $cutoff_coordinate));

    my $look_ahead = $window_length;  #####!!!!!!

    unless ($cutoff_coordinate - 2*($window_length) <= 1) {

      $look_ahead = $self->_calculate_blebbiness($self->_est_sequence->subseq(
			 ($cutoff_coordinate - 2*($window_length)), 
			 ($cutoff_coordinate - $window_length)));
    }

    $bleb_score = $look_ahead if $look_ahead < $bleb_score;

    $cutoff_coordinate -= $window_overlap;
  }

  $cutoff_coordinate += $window_overlap;
 
  if ($cutoff_coordinate != $est_untrimmed_length) {

    $self->error_messages("Total sequences trimmed (but not discarded)", 
			  "Sequences Trimmed (not rejected)");

    my $new_seq = Bio::Seq->new('-display_id' => $self->_est_sequence->display_id,
				'-seq' => $self->_est_sequence->subseq(1, $cutoff_coordinate)
			      );

    $self->_est_sequence($new_seq);

    $self->_trimmed_status(1);

  }

}


sub _calculate_blebbiness {
  my ($self, $seq) = @_;

  my @nt_stack = split //, $seq;

  my $prev_base;
  my $score = 0;

  while (scalar @nt_stack) {
   
    my $current_base = pop @nt_stack;

    $score++ if (($current_base ne $prev_base)&&($current_base =~ /[ATGCatgc]/));

    $prev_base = $current_base;
  }

  return $score;
}


##### Storage/Retrieval methods #####

sub _est_sequence {
  my $self = shift;

  if (@_) {
    $self->{'_est_sequence'} = shift;

    $self->throw("Sequence passed to ESTFilter is not a Bio::Seq") 
      unless $self->{'_est_sequence'}->isa("Bio::Seq");
  }

  return $self->{'_est_sequence'}
}

sub _min_length {
  my $self = shift;

  if (@_){
    $self->{'_min_length'} = shift;
  }

  return $self->{'_min_length'}
}


sub _trimmed_status {
  my $self = shift;

  if (@_){
    $self->{'_trimmed_status'} = shift;
  }

  return $self->{'_trimmed_status'};
}

sub _three_prime_twilight_zone {
  my $self = shift;

  if (@_){
    $self->{'_three_prime_bases'} = shift;
  }

  return $self->{'_three_prime_bases'}
}

sub _bleb_length {
  my $self = shift;

  if (@_){
    $self->{'_bleb_length'} = shift;

    if ($self->{'_bleb_length'} > 9 || $self->{'_bleb_length'} < 6) {
      $self->warn("The \'bleb length\' should be set >= 6 and <= 9 " .
		  "for highest sensitivity.");
    }

  }

  return $self->{'_bleb_length'}
}

sub _max_perc_Ns {
  my $self = shift;

  if (@_){
    # Covert integer percentage to a proportion
    
    my $int_perc = shift;
    my $proportion = $int_perc/100;
    
    $self->{'_max_perc_Ns'} = $proportion;
  }

  return $self->{'_max_perc_Ns'}
}

sub _allowed_stops {
  my $self = shift;

  if (@_){
    $self->{'_allowed_stops'} = shift;
  }

  return $self->{'_allowed_stops'}
}

sub _max_composition {
  my $self = shift;

  if (@_){
    # Covert integer percentage to a proportion
    
    my $int_perc = shift;
    my $proportion = $int_perc/100;

    $self->{'_max_composition'} = $proportion;
  }

  return $self->{'_max_composition'}
}


return 1;
